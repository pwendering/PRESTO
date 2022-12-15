function KegginfTable=get_KEGGid(pIds, org_abbrev)
%Function to retrieve KEGG ids and pathway and orthology information 
%matched to Proteins in the model 
%INPUT
% - char pIds:  A character cell array of uniprot protein IDs as used in
%               GECKO
% - char org_abbrev: A single character vector giving the KEGG organism ID
%                   to query eg.('eco' - E.coli, 'sce' - yeast) 
%OUTPUT:
% - table KegginfTable: A table giving 1- queried pId 2- corresponding
%                       KEGGID 3- KEGG Pathway Info 4 - KEGG orthology Inf


keggId = cell(0,1);
GnResponse = cell(0,1);
URL = 'https://rest.uniprot.org/idmapping/run';

% divide gene names into chunks
chunk_size = 100;
startIdx = 1:chunk_size:numel(pIds);
for i=1:numel(startIdx)
    if i<numel(startIdx)
        endIdx = startIdx(i)+(chunk_size-1);
    else
        endIdx = numel(pIds);
    end
    fprintf('Processing IDs %d to %d ...\n',startIdx(i),endIdx)
    tmpIdx = startIdx(i):endIdx;
    tmpGnQuery = pIds(tmpIdx);
    
    curl_cmd_jobid = ['curl --silent --request POST ' URL ' ',...
        '--form ids="', strjoin(tmpGnQuery,',') '" ',...
        '--form from="UniProtKB_AC-ID" ',...
        '--form to="KEGG"'];
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
    curl_cmd_jobresult = ['curl --silent https://rest.uniprot.org/idmapping/results/' ...
        job_id];
    [status, job_result_json] = system(curl_cmd_jobresult);
    if status == 0
        job_result_struct = jsondecode(job_result_json).results;
        job_result_table = struct2table(job_result_struct);
    else
        error('Job results could not be fetched')
    end
    
    % filter out non-organism specific hits
    keep_idx = cellfun(@(x)startsWith(x, [org_abbrev ':']), job_result_table.to);
    job_result_table = job_result_table(keep_idx, :);
    % extract gene names to match with input gene names for the case
    % that not all IDs were matched
    GnResponse = [GnResponse; regexprep(job_result_table.from,'^\w+:','')];
    keggId = [keggId; job_result_table.to];
end

if length(GnResponse)<length(pIds)
    warning(['Not all gene names could be matched: ', strjoin(pIds(~ismember(pIds, GnResponse)), ', '), ...
        ' were not found in KEGG database'])
end
% now obtain KEGG orthologies and pathways for each gene
KO = cell(size(keggId));
PW = cell(size(keggId));
for i=1:numel(keggId)
    if ~isempty(keggId(i))
        % get gene entry
        res = webread(['http://rest.kegg.jp/get/' keggId{i}]);
        % KEGG orthology
        tmpKO = regexp(res,'K\d{5}','match');
        if ~isempty(tmpKO)
            KO(i) = cellstr(strjoin(tmpKO,'|'));
        end
        % pathways
        tmpPW = regexp(res,[org_abbrev '\d{5}  [\w ,-]+'],'match');
        if ~isempty(tmpPW)
            PW{i} = tmpPW;
        end
    end
end

KegginfTable=table(GnResponse, keggId, PW, KO, 'VariableNames', {'UniprotID', 'KEGGID', 'PathwayInfo', 'KEGGOrthoInfo'})
end
