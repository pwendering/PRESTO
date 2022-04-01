function models = adjBaseModel(model,P,nutrExch,mod_GAM)
%% models = adjBaseModel(model,P,nutrExch,mod_GAM)
% Re-scale a model to multiple experimental conditions using the total
% protein content and nutrient uptake fluxes.
% The function assumes that reverse counterparts of reactions have a suffix
% '_REV'.
% Input:
%   struct model:       metabolic model
%   double P:           array of total protein contents for one or multiple
%                       conditions
%   table nutrExch:     table containing nutrient uptake fluxes for reactions
%                       specified as row names across conditions (columns)
%   double mod_GAM:         growth-associated maintenance, array that contains 
%                       a GAM value for ech condition
% Output:
%   cell models:        cell array containing adjusted model(s)

if ~( (numel(P) == size(nutrExch,2)) && (size(nutrExch,2) == numel(mod_GAM)) )
    error('adjBaseModel: incorrect input dimensions')
end
%detect modle organism 
if contains(model.name, 'Escherichia coli', 'IgnoreCase', true)
    org_name='escherichia coli';
    %This is suboptimal since a lot of unesescary variables like GAM are
    %loaded into the namespace... could be improved
    configuration_ecoli
elseif contains(model.name, 'Yeast', 'IgnoreCase', true)
    org_name='saccharomyces cerevisiae';
    configuration_yeast
end

% initialize models
models = cell(1,numel(P));
for i=1:numel(models)
    % re-scale biomass reaction (GECKO toolbox function)
    cd(fullfile(geckoDir, 'geckomat'))
    params=getModelParameters();
    if ~strcmp(params.org_name, org_name)
        error(['matlab toolbox is not configured for organism ' org_name])
    end
    %check if cobra toolbox confiuration matches org_name
    cd('limit_proteins/')
    models{i} = scaleBioMass(model,P(i),mod_GAM(i));
    cd(topDir)
    % change exchange reaction bounds (COBRA toolbox function)
    excRxns = nutrExch.Properties.RowNames;
    % set upper and lower bounds for all reactions to zero
    models{i} = changeRxnBounds(models{i},...
        models{i}.rxns(ismember(strrep(models{i}.rxns,'_REV',''),excRxns)),0,'b');
    % set all bounds with available information
    excRxns(nutrExch.(i)<0) = strcat(excRxns(nutrExch.(i)<0), '_REV');
    models{i} = changeRxnBounds(models{i},excRxns,...
        abs(nutrExch.(i)),'u');
end
end