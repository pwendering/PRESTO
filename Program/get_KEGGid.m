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
baseURL = 'https://www.uniprot.org';
tool = 'uploadlists';

% divide gene names into chunks of 500
startIdx = 1:100:numel(pIds);
for i=1:numel(startIdx)
    if i<numel(startIdx)
        endIdx = startIdx(i)+99;
    else
        endIdx = numel(pIds);
    end
    fprintf('Processing IDs %d to %d ...\n',startIdx(i),endIdx)
    tmpIdx = startIdx(i):endIdx;
    tmpGnQuery = pIds(tmpIdx);
    url = [baseURL '/' tool '/',...
        '?query=', strjoin(tmpGnQuery,','),...
        '&format=tab',...
        '&from=ACC',...
        '&to=KEGG_ID'...
        '&columns=id,organism'...
        ];
    % send request
    data = webread(url);
    % process output from webread to cell array with an entry for each
    % row
    rows = strsplit(strtrim(data),'\n');
    rows = cellfun(@strsplit,rows(2:end),'un', 0);
    % filter out non-organism specific hits
    rows = rows(cellfun(@(x)contains(x(2),[org_abbrev ':']),rows));
    % extract gene names to match with input gene names for the case
    % that not all IDs were matched
    GnResponse = [GnResponse; regexprep(cellfun(@(x)x(1),rows),'^\w+:','')'];
    keggId = [keggId;cellfun(@(x)x(2),rows)'];
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
