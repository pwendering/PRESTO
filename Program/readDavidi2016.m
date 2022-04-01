function [condNames, E, expVal, nutrExch, P] = readDavidi2016(model,topDir)
% Read data from text files created by readDataDavidi2016 and adjust the
% format for use in PRESTO.
% Input:
%   struct model:           (optional)enzyme-constraint metabolic model generated using
%                           the GECKO toolbox
%   char topDir:            (optional) top-level directory; default: '.'
% Output:
%   cellstr condNames:      names of experimental conditions
%   double E:               matrix of enzyme abundances (dimension: model
%                           enzymes x number of conditions)
%   double expVal:          experimentally measured growth rates for each condition
%   table nutrExch:         table containing available carbon sources and,
%                           if presnt, associated uptake rates
%   double P:               total cellular protein content in g/gDW per condition

if nargin < 2
    topDir = '.';
end


fprintf('Reading experimental data from file...\n')
configuration_ecoli
%% enzyme abundances
protAbundanceData = readtable(fullfile(topDir,'Data','Davidi','abs_proteomics_ecoli.tsv'),...
    'FileType', 'text', 'ReadRowNames', true);
condNames = protAbundanceData.Properties.VariableNames;
if ~isempty(model)
    enzRxnIdx = find(contains(model.rxns,enzMetPfx) & ~ismember(model.rxns, {'prot_pool_exchange'}));
    enzMetIdx = find(any(model.S(:,enzRxnIdx),2) & ~ismember(model.mets, {'prot_pool'}));

    new_enzymes = cellfun(@(x)regexp(x,'[A-Z0-9]{6}','match'),model.metNames(enzMetIdx));
    if ~isempty(setdiff(new_enzymes, model.enzymes))||~isempty(setdiff(model.enzymes, new_enzymes))
        error('number of enzyme pseudometabolites differs from model.enzyme field. Unable to match proteomics data')
    elseif ~isequal(new_enzymes, model.enzymes)
        warning(['Model.enzyme field does not have same ordering as enzyme pseudometabolites in stochiometric matrix\n' ...
            'Ordering E matrix according to pseudometabolites'])
        model.enzymes = new_enzymes;
    end

    % match E to gene IDs in the model
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

% read experimental growth rates
growthRateTab = readtable(fullfile(topDir,'Data','Davidi','growth_rates_ecoli.tsv'),...
    'FileType','text','ReadRowNames',true);
growthRateTab.Properties.RowNames = regexprep(growthRateTab.Properties.RowNames,...
    '[@=.+]','_');
% if necessary, re-order according to condNames from protein abundance table
matchIdx = cellfun(@(x)find(ismember(condNames,x)),growthRateTab.Properties.RowNames);
expVal = table2array(growthRateTab(matchIdx,:));
clear growthRateTab

pTab = readtable(fullfile(topDir,'Data','Davidi','total_protein_ecoli.tsv'),...
    'FileType','text','ReadRowNames',true);
pTab.Properties.RowNames = regexprep(pTab.Properties.RowNames,...
    '[@=.+]','_');
% if necessary, re-order according to condNames from protein abundance table
matchIdx = cellfun(@(x)find(ismember(condNames,x)),pTab.Properties.RowNames);
P = table2array(pTab(matchIdx,:));
if ~all(P>0)
    warning('Total protein not available for all conditions')
    fprintf('==> Using maximum protein content for these conditions\n\n')
    P(isnan(P)) = max(P);
end

% nutrient uptake data
nutrExch = readtable(fullfile(topDir,'Data','Davidi','csource_ecoli.tsv'),...
    'FileType','text','ReadRowNames',true);
% if necessary, re-order according to condNames from protein abundance table
matchIdx = cellfun(@(x)find(ismember(condNames,x)),nutrExch.Properties.VariableNames);
nutrExch = nutrExch(:,matchIdx);

% add M9 medium
% from doi: 10.1016/j.jbiotec.2009.10.007 (medium reference in Valgepea et
% al.)
% nutrients from MetaCyc M9 definition (w/o glycerol) and missing elements found by
% biomassPrecursorCheck

% medium used in Schmidt et al. 2015 further contains thiamine solution

M9ExchRxns = {'EX_na1_e_REV','EX_pi_e_REV','EX_cl_e_REV',...
    'EX_k_e_REV','EX_nh4_e_REV','EX_mg2_e_REV','EX_so4_e_REV','EX_mobd_e_REV',...
    'EX_mn2_e_REV','EX_ni2_e_REV','EX_zn2_e_REV','EX_cu2_e_REV','EX_ca2_e_REV',...
    'EX_fe2_e_REV','EX_fe3_e_REV','EX_cd2_e_REV','EX_cobalt2_e_REV','EX_h2o_e_REV',...
    'EX_o2_e_REV'};

M9 = array2table(repelem(1000,numel(M9ExchRxns),numel(condNames)),...
    'RowNames',M9ExchRxns,'VariableNames',nutrExch.Properties.VariableNames);

nutrExch = [nutrExch;M9];

end