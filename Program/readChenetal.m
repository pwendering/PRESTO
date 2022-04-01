function[condNames, E, expVal, nutrExch, P] =readChenetal(topDir, model, rem_cond)
% Function to read in chemostat and proteomics data from the dataset of
% Chen et al. 2021 from Data/Chen
%INPUT:
% - string topDir:  Path of the working directory in which Data/Chen is
%                           localized
% - struc model:    (optional)An ecModel generated by GECKO to which protein abundancies
%                           should be matched. If this is empty no protein
%                           Abundancies will be returned (E=[])
% - boolean rem_cond:  (optional)if true outlier samples are removed
%OUTPUT:
% - cell condNames: Conditions names of proteomics conditions
% - double E:     Matrix of enzyme abundancies (columns = conditions, rows=
%                       model enzymes (missing values 0)
% - double expVal:  Experimental growth rates for each conditions
% - double nutrExch: Matrix bearing experimentally measured nutrient exchange
%                               rates. (Missing values =0)
% - double P:   Total cellular protein content in g/gDW per condition
%
fprintf('Reading experimental data from file...\n')
%% Check wether enzyme field has different ordering than protein pseudometabolites in matrix
configuration_yeast
protAbundanceData = readtable([topDir '/Data/Chen/abs_proteomics_yeast.tsv'], 'FileType', 'text', 'ReadRowNames',1);
condNames = protAbundanceData.Properties.VariableNames;

if nargin>1 && ~isempty(model)
    enzRxnIdx = find(contains(model.rxns,enzMetPfx)&~ismember(model.rxns, {'prot_pool_exchange'}));
    enzMetIdx = find(any(model.S(:,enzRxnIdx),2)&~ismember(model.mets, {'prot_pool'}));
    new_enzymes= cellfun(@(x)regexp(x,'[A-Z0-9]{6}','match'),model.metNames(enzMetIdx));
    if ~isempty(setdiff(new_enzymes, model.enzymes))||~isempty(setdiff(model.enzymes, new_enzymes))
        error('number of enzyme pseudometabolites differs from model.enzyme field. Unable to match proteomics data')
    elseif ~isequal(new_enzymes, model.enzymes)
        warning(['Model.enzyme field does not have same ordering as enzyme pseudometabolites in stochiometric matrix\n' ...
            'Ordering E matrix according to pseudometabolites'])
        model.enzymes=new_enzymes;
    end
    %Only read in protein abundance data if a model is given as input
    protIDs = protAbundanceData.Properties.RowNames;
    % map protein abundance data to proteins in the model
    idxModel2Data = cell2mat(cellfun(@(x)find(ismember(protIDs,x)),...
        model.enzymes,'un',0));
    % set up enzyme abundance matrix
    E = nan(numel(model.enzymes),size(protAbundanceData,2));
    E(ismember(model.enzymes,protIDs),:) = table2array(protAbundanceData(idxModel2Data,:));
    E(isnan(E)) = 0;
else
    E=[];
end
clear protAbundanceData

if nargin <3
    rem_cond=false;
end
growthData = readtable([topDir '/Data/Chen/growth_rates_yeast.tsv'], 'FileType', 'text');
% re-order column as in protein abundance file
growthData = [growthData(:,[1 2]), growthData(:,2+cellfun(@(x)find(...
    ismember(growthData.Properties.VariableNames(3:end),x)),...
    condNames))];
% experimental growth rates
expVal = table2array(growthData(1,3:end));
% nutrient uptake bounds
nutrExch = growthData(2:end,3:size(growthData,2));
nutrExch.Properties.RowNames = growthData.(1)(2:end);
clear growthData


% read data for total protein content
pTab = readtable([topDir '/Data/Chen/total_protein_yeast.tsv'], 'FileType', 'text', 'ReadRowNames', 1);
matchIdx = cellfun(@(x)find(ismember(pTab.Properties.RowNames,x)),condNames);
P = nan(size(condNames));
P(matchIdx) = table2array(pTab);
if ~all(P>0)
    warning('Total protein not available for all conditions')
    fprintf('==> Using maximum protein content for these conditions\n\n')
    P(isnan(P)) = max(P);
    if rem_cond
        %% Exclude inappropriate datasets
        % exclude Lathvee2017_Temp33,Lathvee2017_Temp36,Lathvee2017_Temp38
        idxExclude = contains(condNames,{'Lahtvee2017_Temp33',...
            'Lahtvee2017_Temp36','Lahtvee2017_Temp38'});
        condNames=condNames(~idxExclude);
        expVal = expVal(~idxExclude);
        P = P(~idxExclude);
        if ~isempty(model)
            E = E(:,~idxExclude);
        end
        nutrExch = nutrExch(:,~idxExclude);
    end
end
end