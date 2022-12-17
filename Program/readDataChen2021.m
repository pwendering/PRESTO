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
    URL = 'https://rest.uniprot.org/idmapping/run';
    
    % divide gene names into chunks
    chunksize = 500;
    startIdx = 1:chunksize:numel(geneNames);
    for i=1:numel(startIdx)
        if i<numel(startIdx)
            endIdx = startIdx(i)+chunksize-1;
        else
            endIdx = numel(geneNames);
        end
        fprintf('Processing gene names %d to %d ...\n',startIdx(i),endIdx)
        tmpGeneIdx = startIdx(i):endIdx;
        tmpGnQuery = geneNames(tmpGeneIdx);
        
        curl_cmd_jobid = ['curl --silent --request POST ' URL ' ',...
            '--form ids="', strjoin(tmpGnQuery,',') '" ',...
            '--form from="Gene_Name" ',...
            '--form to="UniProtKB"'];
        
        % send requestto retrieve job ID
        status = 1;
        trials = 0;
        while status ~= 0 && trials < 3
            trials = trials + 1;
            [status, job_id_json] = system(curl_cmd_jobid);
        end
        
        if status ~= 0
            error('Request not successful after %i attempts', trials)
        else
            job_id = jsondecode(job_id_json).jobId;
        end
        
        % wait until job is finished
        job_status = '';
        curl_cmd_jobstatus = ['curl --silent https://rest.uniprot.org/idmapping/status/' ...
            job_id];
        tic
        t = 0;
        while ~isequal(job_status, 'FINISHED') && t < 120
            t = toc;
            [status, job_status_json] = system(curl_cmd_jobstatus);
            if status == 0
                job_status = jsondecode(job_status_json).jobStatus;
            end
        end

        % get job results
        curl_cmd_jobresult = ['curl --silent https://rest.uniprot.org/idmapping/uniprotkb/results/' ...
            job_id '?fields=organism_name'];
        [status, job_result_json] = system(curl_cmd_jobresult);
        if status == 0
            job_result_struct = jsondecode(job_result_json).results;
            job_resuls_from = arrayfun(@(i)job_result_struct(i).from,...
                1:numel(job_result_struct),'UniformOutput',false)';
            job_resuls_to = arrayfun(@(i)job_result_struct(i).to.primaryAccession,...
                1:numel(job_result_struct),'UniformOutput',false)';
            job_resuls_org = arrayfun(@(i)job_result_struct(i).to.organism.scientificName,...
                1:numel(job_result_struct),'UniformOutput',false)';
            job_result_table = cell2table([job_resuls_from job_resuls_to job_resuls_org],...
                'VariableNames', {'from', 'to', 'organism'});
        else
            error('Job results could not be fetched')
        end
        
        % filter out non-organism specific hits
        keep_idx = cellfun(@(x)startsWith(x,...
            'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)'),...
            job_result_table.organism);
        job_result_table = job_result_table(keep_idx, :);

        tmpGnResponse = job_result_table.from;
        matchIdx = cellfun(@(x)find(ismember(tmpGnQuery,x)),tmpGnResponse);
        if numel(matchIdx)<numel(tmpGeneIdx)
            fprintf('\tGene ID could not be translated: %s\n',...
                tmpGnQuery{setdiff(1:numel(matchIdx),matchIdx)})
        end
        upId(tmpGeneIdx(matchIdx)) = job_result_table.to;
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