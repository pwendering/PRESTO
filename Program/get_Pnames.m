function[protName]=get_Pnames(up_ID)
    %Function to obtaine gene names from uniprot ID using the web API
    %INPUT:
    % cell up_ID:   A cell array of character vectors containing Uniprot IDs
    %               of Proteins for which IDs should be retrieved
    %OUTPUT: 
    % cell protName: A cell array of character vectors containing the
    %               respective protein names found for up_ID
    URL = 'https://rest.uniprot.org/idmapping/run';

% divide gene names into chunks
chunk_size = 100;
startIdx = 1:chunk_size:numel(up_ID);
for i=1:numel(startIdx)
    if i<numel(startIdx)
        endIdx = startIdx(i)+(chunk_size-1);
    else
        endIdx = numel(up_ID);
    end
    fprintf('Processing IDs %d to %d ...\n',startIdx(i),endIdx)
    tmpIdx = startIdx(i):endIdx;
    tmpGnQuery = up_ID(tmpIdx);
    
    curl_cmd_jobid = ['curl -k --silent --request POST ' URL ' ',...
        '--form ids="', strjoin(tmpGnQuery,',') '" ',...
        '--form from="UniProtKB_AC-ID" ',...
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
    curl_cmd_jobstatus = ['curl -k --silent https://rest.uniprot.org/idmapping/status/' ...
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
    curl_cmd_jobresult = ['curl -k --silent "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/', ...
        job_id, '?fields=accession%2Cid%2Cprotein_name&format=tsv"'];
    [status, job_result_tsv] = system(curl_cmd_jobresult);
    if status == 0
        rows = strsplit(strtrim(job_result_tsv),'\n');
    else
        error('Job results could not be fetched')
    end
    if length(rows)-1 ~= length(up_ID)
        error('not all gene names could be retrieved from Uniprot server')
    end
    rows = cellfun(@(x) strsplit(x, '\t') ,rows(2:end), 'un', 0);
    protName=cellfun(@(x) extractBefore(x{4}, '(EC'), rows, 'un', 0);
end