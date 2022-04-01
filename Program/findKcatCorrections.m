function [kcatOrig,kcatCorr,protIDs,reportTable] = findKcatCorrections(origModel,corrModel,enzMetPfx)
%% [kcatOrig,kcatCorr,protIDs] = findCorrectionsInModel(origModel,corrModel,enzRxnPfx)
% Compare turnover values between the original GECKO model and after the
% kcat correction.
% Input:
%   struct origModel:       enzyme-constraint metabolic model generated
%                           using the GECKO2.0 toolbox
%   struct corrModel:       same model as origModel with changed entries
%                           for kcat values in enzyme mass balance constraints
%   char enzMetPfx:         prefix for enzyme mass balance constraints
%                           (i.e. values in model.metNames); default: 'prot_'
% Output:
%   double kcatOrig:        original kcats that have been corrected
%   double kcatCorr:        corrected values for the kcats in kcatOrig
%   cellstr protIDs:        enzyme metabolite names (i.e., protein IDs)
%                           associated to corrected kcats
%   table reportTable:      Table with information on the corrected kcats
%                           including Protein IDm, Reaction ID, original
%                           kcat, updated kcat, fold increase and
%                           delta(absolute increase)

% find enzyme metabolite indices
enzMetIdx = find(contains(origModel.metNames,enzMetPfx));
% initialize arrays
kcatOrig = [];
kcatCorr = [];
protIDs = cell(0,1);
rxnIDs = [];

% initialize report table
% loop over enzyme mass balance constraints
for i=1:numel(enzMetIdx)
    tmpKcatIdx = find(origModel.S(enzMetIdx(i),:)<0);
    % find kcats in both models
    tmpOrigKcat = -1./origModel.S(enzMetIdx(i),tmpKcatIdx)/3600;
    tmpCorrKcat = -1./corrModel.S(enzMetIdx(i),tmpKcatIdx)/3600;
    % if kcat has been updated, store both in returned array
    changeIdx = tmpCorrKcat>tmpOrigKcat;
    % also store associated enzyme metabolite name
    kcatOrig = [kcatOrig; tmpOrigKcat(changeIdx)'];
    kcatCorr = [kcatCorr; tmpCorrKcat(changeIdx)'];
    protIDs = [protIDs; repmat(origModel.metNames(enzMetIdx(i)),sum(changeIdx),1)];
    rxnIDs=[rxnIDs; tmpKcatIdx(changeIdx)'];
    
    
end

% create report table
reportTable = unique(...
    sortrows(...
        table(...
            strrep(protIDs,enzMetPfx,''),...
            rxnIDs,...
            kcatOrig,...
            kcatCorr,...
            kcatCorr./kcatOrig,...
            kcatCorr-kcatOrig,...
            'VariableNames',{'PROTEIN ID', 'REACTION ID', 'KCAT ORIG [s^-1]',...
            'KCAT UPDATED [s^-1]','FOLD-INCREASE','DELTA'}),...
        'DELTA','descend'),...
    'stable');
end