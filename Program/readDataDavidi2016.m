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
    URL = 'https://rest.uniprot.org/idmapping/run';
    
    % divide gene names into chunks
    chunksize = 500;
    startIdx = 1:chunksize:numel(geneIDs);
    for i=1:numel(startIdx)
        if i<numel(startIdx)
            endIdx = startIdx(i)+chunksize-1;
        else
            endIdx = numel(geneIDs);
        end
        fprintf('Processing gene names %d to %d ...\n',startIdx(i),endIdx)
        tmpGeneIdx = startIdx(i):endIdx;
        tmpGnQuery = geneIDs(tmpGeneIdx);
        
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
            'Escherichia coli (strain K12)'),...
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