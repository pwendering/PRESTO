%% Read data from Chen et al. 2021 study (doi: 10.1073/pnas.2108391118)

% define output file names
protOutFileName = fullfile('Data','Chen','abs_proteomics_yeast.tsv');
growthOutFileName =  fullfile('Data','Chen','growth_rates_yeast.tsv');
ptotOutFileName = fullfile('Data','Chen', 'total_protein_yeast.tsv');
% translation file from gene names to UniProt identifiers
gn2upFile = fullfile('Data','Chen','gene_name_uniprot_acc.csv');
% Avogadro constant [1/mol]
N_A = 6.02214076e23;
% define data file names from Chen et al. 2021 (PNAS) DOI: 10.1073/pnas.2108391118
protDataFile = fullfile('Data', 'Chen', 'ProteomicsFlux.xlsx');
uniProtFile = fullfile('Data', 'Chen', 'UniProt.xlsx');
% read UniProt file (molecular mass in g/mol)
uniprot = readtable(uniProtFile);
% get names of sheets within protein abundance Excel file
sheetNames = sheetnames(protDataFile);
sheetNames = sheetNames(1:end-1); % exclude last sheet
% initialize array for absolute protein abundances (second to last sheet)
protAbcs = nan(size(uniprot,1),100);
geneNames = uniprot.Gene;
condNames = {};
% loop over sheets in Excel file
colCount = 0;
for i=1:numel(sheetNames)
    if i==1
        % first sheet contains exchange fluxes and growth rates
        growthData = readtable(protDataFile, 'Sheet', sheetNames{i});
        tmpCondIds = growthData.Properties.VariableNames(3:end);
        growthDataNonNum = growthData(:,[1 2]);
        
        % calculate average over replicates if needed
        uniqCondNames = unique(regexprep(tmpCondIds,'R\d+$',''),'stable');
        growthData = table2array(growthData(:,3:end));
        growthData = cell2mat(arrayfun(@(colIdx)...
            mean(growthData(:,contains(tmpCondIds,uniqCondNames{colIdx})),2,'omitnan'),...
            1:numel(uniqCondNames),'un',0));
        growthData = [growthDataNonNum array2table(growthData,...
            'VariableNames', uniqCondNames)];
        clear growthDataNonNum tmpCondIds uniqCondNames
    elseif i==2
        ptotData = readtable(protDataFile, 'Sheet', sheetNames{i});
    else
        % read abundances for current study
        tab = readtable(protDataFile, 'Sheet', sheetNames{i}, 'ReadRowNames',1);
        % find current gene names
        tmpGN = tab.Properties.RowNames;
        % five row names in 'DiBartolomeo2020' and one in 'Yu2021' are
        % alternative gene names, separated by ';'
        tmpGN = strtok(tmpGN,';');
        % get condition IDs 
        tmpCondIds = tab.Properties.VariableNames;
        
        % the data still contain replicates, so the maximum is calculated
        % across all replicates for a condition
        tab = table2array(tab);
        uniqCondNames = unique(regexprep(tmpCondIds,'R\d+$',''),'stable');
        tab = cell2mat(arrayfun(@(colIdx)...
            max(tab(:,contains(tmpCondIds,uniqCondNames{colIdx})),[],2),...
            1:numel(uniqCondNames),'un',0));
        nCond = size(tab,2);
        condNames = [condNames uniqCondNames];
        
        % map current gene names to the whole geneNames array
        matchIdx = cellfun(@(x)find(ismember(geneNames,x)),tmpGN);
        
        clear tmpCondIds uniqCondNames tmpGN
        % correct abundance units and append table to protAbcs
        switch sheetNames{i}
            case 'Lahtvee2017'
                % original unit: copy/pgCDW
                tab = 1e15*tab/N_A;
            case 'Yu2020'
                % original unit: fmol/mgCDW
                tab = tab/1e9;
            case 'DiBartolomeo2020'
                % original unit: g/gCDW
                tab = 1e3*diag(1./uniprot.Mass(matchIdx))*tab;
            case 'Yu2021'
                % original unit: fmol/mgCDW
                tab = tab/1e9;
        end
        % add abundances to matrix
        protAbcs(matchIdx,colCount+1:colCount+nCond) = tab;
        colCount = colCount+nCond;
    end
end

% remove all-zero rows and columns
geneNames(~any(protAbcs,2)) = [];
protAbcs(~any(protAbcs,2),:) = [];
protAbcs(:,~any(protAbcs)) = [];

% translate gene names to UniProt IDs using ID mapping service
if ~exist(gn2upFile,'file')
    upId = repmat({''},size(geneNames));
    baseURL = 'https://www.uniprot.org';
    tool = 'uploadlists';
    
    % divide gene names into chunks of 500
    startIdx = 1:500:numel(geneNames);
    for i=1:numel(startIdx)
        if i<numel(startIdx)
            endIdx = startIdx(i)+499;
        else
            endIdx = numel(geneNames);
        end
        fprintf('Processing gene names %d to %d ...\n',startIdx(i),endIdx)
        tmpGeneIdx = startIdx(i):endIdx;
        tmpGnQuery = geneNames(tmpGeneIdx);
        url = [baseURL '/' tool '/',...
            '?query=', strjoin(tmpGnQuery,','),...
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
        rows = rows(contains(rows,'Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker''s yeast)'));
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
    gn2up = cell2table([geneNames upId],'VariableNames', {'GENENAME', 'ID'});
    writetable(gn2up,gn2upFile)
else
    gn2up = readtable(gn2upFile);
end

% create new table with uniprot IDs times abundances
rowNames = gn2up.ID;
idxEmpty = cellfun(@isempty,rowNames);
rowNames(idxEmpty) = geneNames(idxEmpty);
writetable(array2table(protAbcs, 'VariableNames', condNames,...
    'RowNames', rowNames), protOutFileName,'WriteRowNames', true,...
    'Delimiter', '\t', 'FileType', 'text')
% write growth rates to file
writetable(growthData,growthOutFileName,'Delimiter','\t', 'FileType', 'text')
% write protein content to file
writetable(ptotData,ptotOutFileName,'Delimiter','\t', 'FileType', 'text')
clear