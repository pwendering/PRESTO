function [ecModel,ecModel_batch] = enhanceGEM(model,toolbox, unmod, name,modelVer)
% enhanceGEM
%
%   Main function for running the GECKO pipeline. It returns an ecModel and
%   its constrained version with an upper limit in the total protein pool
%   usage pseudoreaction (ecModel_batch) with calibrated Kcat parameters that
%   allow it to grow at a specified experimental growth rate.
%
%   model       a GEM MATLAB structure compatible with the COBRA or RAVEN
%               toolbox.
%   toolbox     string with the name of the prefered toolbox for model SBML
%               export (COBRA or RAVEN)
%   unmod         logical indicating wether the model without manual
%                   curations should be used (opt, default false
%   name        Desired name for the ecModel (opt, default '')
%   modelVer    modelVer of the original GEM (opt, default '')
%   ecModel     an ecModel MATLAB structure suitable for incorporation of   
%                      proteomics data as individual enzyme usage constraints.
%   ecModel_batch  an ecModel MATLAB structure with a global constraint on 
%                  the total protein pool usage pseudoreaction,
%                  proportional to the measured total protein content (Ptot)
%
%   Usage: [ecModel,ecModel_batch] = enhanceGEM(model,toolbox,name,modelVer)
%
%   Ivan Domenzain. Last edited: 2020-10-05
%
if nargin < 3
    raw=false;
end
if nargin < 4
    name    = '';
end
if nargin < 5
    modelVer = '';
end

cd change_model
%remove wrongly assigned proteins
if contains(name, 'yeast', 'IgnoreCase', true)
    model=preprocessYeast(model);
end
cd ..
%Convert model to RAVEN for easier visualization later on:
format short e
if isfield(model,'rules')
    %initCobraToolbox
    model = ravenCobraWrapper(model);
end

fprintf('\n***************************************************************')
fprintf('\n   GECKO: Adding enzyme constraints to a genome-scale model')
fprintf('\n***************************************************************\n\n')

%Get model-specific parameters
parameters = getModelParameters;

%Remove blocked rxns + correct model.rev:
cd change_model

%Remove blocked rxns + correct model.rev:
[model,name,modelVer] = preprocessModel(model,name,modelVer);

fprintf('\n==================')
fprintf('\nGenerating ecModel:')
fprintf('\n==================\n')

%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model);
[kcats, origintab]      = matchKcats(model_data,parameters.org_name);

%Save a table with origins for each kcat reaction
writematrix(origintab, ['../../models/', name, '/', name, '_kcatOrigins.txt'])

%Integrate enzymes in the model:
cd ../change_model
ecModel = readKcatData(model_data,kcats);
cd ../../models
if unmod
ecModel = saveECmodel(ecModel,toolbox, name, ['raw' name], modelVer);
else
%save(fullfile(name, ['raw' name '.mat']), 'ecModel')
%add manual modifications to the model
cd ../geckomat/change_model
[manecModel,modifications] = manualModifications(ecModel);
cd ../../models
manecModel = saveECmodel(manecModel,toolbox, name,  ['manmod' name], modelVer);
end
%save(fullfile(name, ['manmod' name '.mat']), 'manecModel')


%Constrain model to batch conditions:
fprintf('\n==============================================================')
fprintf('\nGenerating ecModel with shared pool assumption (ecModel_batch):')
fprintf('\n==============================================================\n')
cd ../geckomat/limit_proteins
disp('Generate batch models')
if unmod
[ecModel_batch, OptSigma, ~] = getConstrainedModel(ecModel, {},name, true);
[adpecModel_batch, OptSigma, rawchanges] = getConstrainedModel(ecModel, {},name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])
%save table with changes
    if ~isempty(rawchanges)
        writetable(rawchanges, ['../../models/', name, '/raw', name, '_kcatModifications.txt'])
    end
else
[manecModel_batch, OptSigma, ~] = getConstrainedModel(manecModel,modifications,name, true);
[adpmanecModel_batch, OptSigma, manmodchanges] = getConstrainedModel(manecModel,modifications,name);
disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])
%save table with changes
    if ~isempty(manmodchanges)
        writetable(manmodchanges, ['../../models/', name, '/manmod', name, '_kcatModifications.txt'])
    end
end

%Save output models:
fprintf('\n=============')
fprintf('\nSaving models:')
fprintf('\n=============\n')
cd ../../models
if unmod
%save unfitted models
ecModel_batch = saveECmodel(ecModel_batch,toolbox,name, ['raw' name '_batch'],modelVer);
%safve adapted models
adpecModel_batch = saveECmodel(adpecModel_batch,toolbox, name, ['adpraw' name '_batch'],modelVer);
else
manecModel_batch=saveECmodel(manecModel_batch, toolbox, name, ['manmod' name '_batch'], modelVer);
adpmanecModel_batch=saveECmodel(adpmanecModel_batch, toolbox, name,  ['adpmanmod' name '_batch'], modelVer);
end
cd ../geckomat

end