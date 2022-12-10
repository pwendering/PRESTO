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

% initialize models
models = cell(1,numel(P));
geckoDir = char(regexp(which('getModelParameters'), '.*(?=\\geckomat)', 'match'));
if isempty(geckoDir)
    error('Organism-specific GECKO directory could not be located, please run configuration script')
end
topDir = char(regexp(pwd, '.*PRESTO', 'match'));

for i=1:numel(models)
    % re-scale biomass reaction (GECKO toolbox function)
    cd(fullfile(geckoDir, 'geckomat'))
    params = getModelParameters();
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