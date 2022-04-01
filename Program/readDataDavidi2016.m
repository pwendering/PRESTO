%% Read data from Davidi et al. 2016 study (doi: 10.1073/pnas.1514240113)

% define outfile names
protOutFileName = fullfile('Data', 'Davidi', 'abs_proteomics_ecoli.tsv');
growthOutFileName =  fullfile('Data', 'Davidi', 'growth_rates_ecoli.tsv');
ptotOutFileName = fullfile('Data', 'Davidi', 'total_protein_ecoli.tsv');
nutrExchOutFileName = fullfile('Data', 'Davidi', 'csource_ecoli.tsv');

% load conditions used by Davidi et al. 2016
load(fullfile('Data','Davidi','conditions31.mat'))
condNames = conditions31.cond;
expVal = conditions31.growth;
cReact = conditions31.carbonreac;
clear conditions31

% get references for the conditions used in Davidi et al. 2016
davidi_si = readtable(fullfile('Data','Davidi','pnas.1514240113.sd01.xlsx'),...
    'Sheet','growth conditions');
matchIdx = cellfun(@(x)find(ismember(condNames,x)),davidi_si.growthCondition);
refNames = davidi_si.reference(matchIdx);
nutrUptRates = davidi_si.uptakeRate_mmol_gCDW_H_(matchIdx);
clear davidi_si matchIdx

% enzyme abundances and gene IDs
load(fullfile('Data','Davidi','abundance_file.mat'))
geneIDs = cellstr(abundance_file.genes);

% take the maximum over replicates if present
trimConds = regexprep(abundance_file.cond,'_replicate\d+','');
uniqConds = unique(trimConds,'stable');
E = abundance_file.abun;
E = cell2mat(arrayfun(@(i)max(E(:,contains(trimConds,uniqConds{i})),[],2),...
            1:numel(uniqConds),'un',0));

% take re-ordered subset of matched conditions
matchIdx = cell2mat(cellfun(@(x)find(ismember(uniqConds,x)),condNames,'un',0));
E = E(:,matchIdx);
clear abundance_file uniqConds trimConds matchIdx

% total protein content (only measured in Valgepea et al. 2013)
% Values from Valgepea et al. study:
% μ=0.11 h-1	μ=0.21 h-1	μ=0.31 h-1	μ=0.40 h-1	μ=0.49 h-1
% P=61.0%        P=58.2%	P=55.4%     P=52.9%     P=50.4%
P = repelem(0.61,numel(condNames),1);

% create nutrient exchange table
uniqUptRxns = unique(cReact,'stable');
nutrExch = array2table(zeros(numel(uniqUptRxns),numel(condNames)),...
    'VariableNames',condNames,'RowNames',strcat(uniqUptRxns, '_REV'));

for i=1:numel(uniqUptRxns)
    nutrExch{i,ismember(cReact,uniqUptRxns(i))} = nutrUptRates(ismember(cReact,uniqUptRxns(i)))';
    if ~any(nutrExch{i,:}~=0 & ~isnan(nutrExch{i,:})) 
        nutrExch{i,isnan(nutrExch{i,:})} = 1000;
    else
        nutrExch{i,isnan(nutrExch{i,:})} = max(nutrExch{i,:});
    end
end


% translate KEGG gene IDs to UniProt IDs in the model
gn2upFile = fullfile('Data','Davidi','gene_name_uniprot_acc.csv');
% set webread timeout to 30s
weboptions.Timeout = 30;
if ~exist(gn2upFile,'file')
    upId = repmat({''},size(geneIDs));
    baseURL = 'https://www.uniprot.org';
    tool = 'uploadlists';
    
    % divide gene names into chunks of 500
    startIdx = 1:500:numel(geneIDs);
    for i=1:numel(startIdx)
        if i<numel(startIdx)
            endIdx = startIdx(i)+499;
        else
            endIdx = numel(geneIDs);
        end
        fprintf('Processing gene names %d to %d ...\n',startIdx(i),endIdx)
        tmpGeneIdx = startIdx(i):endIdx;
        tmpGnQuery = geneIDs(tmpGeneIdx);
        url = [baseURL '/' tool '/',...
            '?query=',strjoin(tmpGnQuery,','),...
            '&format=tab',...
            '&from=GENENAME',...
            '&to=ACC'...
            '&columns=id,organism'...
            ];
        % send request
        data = webread(url);
        % process output from webread to cell array with an entry for each
        % row
        rows = strsplit(strtrim(data),'\n');
        rows = rows(contains(rows,'Escherichia coli (strain K12)'));
        rows = cellfun(@strsplit,rows,'un', 0);
        % extract gene names to match with input gene names for the case
        % that not all IDs were matched
        tmpGnResponse = cellfun(@(x)x(end),rows);
        matchIdx = cellfun(@(x)find(ismember(tmpGnQuery,x)),tmpGnResponse);
        if numel(matchIdx)<numel(tmpGeneIdx)
            fprintf('\tGene ID could not be translated: %s\n',...
                tmpGnQuery{setdiff(1:numel(matchIdx),matchIdx)})
        end
        upId(tmpGeneIdx(matchIdx)) = cellfun(@(x)x(1),rows);
    end
    gn2up = cell2table([geneIDs upId],'VariableNames', {'GENENAME', 'ID'});
    writetable(gn2up,gn2upFile)
else
    gn2up = readtable(gn2upFile);
end

% write abundances to file
rowNames = gn2up.ID;
idxEmpty = cellfun(@isempty,rowNames);
rowNames(idxEmpty) = geneIDs(idxEmpty);
writetable(array2table(E, 'VariableNames', condNames,...
    'RowNames', rowNames), protOutFileName,'WriteRowNames', true,...
    'Delimiter', '\t', 'FileType', 'text')

% write growth rates to file
writetable(array2table(expVal,'RowNames',condNames), growthOutFileName,...
    'WriteVariableNames', false, 'WriteRowNames', true, 'Delimiter', '\t',...
    'FileType', 'text');

% write protein content to file
writetable(array2table(P,'RowNames',condNames), ptotOutFileName,...
    'WriteVariableNames', false, 'WriteRowNames', true,...
    'Delimiter','\t', 'FileType', 'text')

% write carbon source info to file
writetable(nutrExch,nutrExchOutFileName,...
    'WriteVariableNames', true, 'WriteRowNames', true,...
    'Delimiter','\t', 'FileType', 'text')

clear;clc