function batch_models = addGKOconstrains(models, org_name)
%Function to add protein pool constrain from GECKO models to PRESTO to
%compare performance
%INPUT:
% - cell models: A cell array of the PRESTO model with condition specific
%               biomass in each cell
% - orgname: The model organism

switch org_name
    case 'saccharomyces cerevisiae'
        configuration_yeast
        condNames=readChenetal(topDir, [], true);
    case 'escherichia coli'
        configuration_ecoli
        condNames=readDavidi2016([],topDir);
end

[loginf, ~]= getloginf(gkologFile, org_name);
if length(loginf.condNames)~=length(condNames)
    error("GECKO model creation log does not contain model for each condition. Matchin of sigma values impossible")
elseif length(models)~=length(condNames)
    error("Number of models does not match number of conditions. Matching of Sigma ismpossible")
end
%order GECKO model creation logfile according to conditions
[match, idx]=ismember(condNames, loginf.condNames);
if any(~match, 'all')
    error("Conditions in GECKO model creation log does not match experimental conditions")
end
loginf=loginf(idx,:);


%convert modeld to batch models
cd([geckoDir, '/geckomat/limit_proteins'])

batch_models=cell(1,length(models));
%get the amount of protein mass accounted for 
%by enzymes in the model
[f,~] = measureAbundance(models{1}.enzymes);
disp(['Fraction of model enzymes in total protein content (f) = ', string(f)])
%build initial model
load('../../databases/parameters.mat')
%check if COBRA toolbox is setup correctly
if ~strcmp(parameters.org_name, org_name)
    error(['matlab toolbox is not configured for organism ' org_name])
end
parameters.sigma=loginf.Sigma(1);
models{1}.MWs=models{1}.protMW/1000; %GECKO uses kg per mol 
batch_models{1}= constrainEnzymes(models{1}, f, true, [], sumProtein(models{1}));
batch_models{1}.csense=repmat('E', length(batch_models{1}.mets),1);

% restore default sigma
load('../../databases/parameters.mat')
parameters.sigma=0.5;
save('../../databases/parameters.mat', 'parameters')


%Adapt Biomass and protein pool constrain for all other conditions
for i=2:length(models)
    batch_models{i}=scaleBioMass(batch_models{1}, sumProtein(models{i}), GAM(i));
    batch_models{i}.ub(end)=f*loginf.Sigma(i)*sumProtein(models{i});
end
cd(topDir)
end