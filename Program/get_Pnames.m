function[protName]=get_Pnames(up_ID)
    %Function to obtaine gene names from uniprot ID using the web API
    %INPUT:
    % cell up_ID:   A cell array of character vectors containing Uniprot IDs
    %               of Proteins for which IDs should be retrieved
    %OUTPUT: 
    % cell protName: A cell array of character vectors containing the
    %               respective protein names found for up_ID
    baseURL = 'https://www.uniprot.org';
    tool = 'uploadlists';
    url = [baseURL '/' tool '/',...
        '?query=', strjoin(up_ID,','),...
        '&format=tab',...
        '&from=ACC',...
        '&to=ACC'...
        '&columns=id,protein_names'...
        ];
    % send request
    disp('Obtaining gene names from uniprot server')
    data = webread(url);
    rows = strsplit(strtrim(data),'\n');
    if length(rows)-1 ~= length(up_ID)
        error('not all gene names could be retrieved from Uniprot server')
    end
    rows = cellfun(@(x) strsplit(x, '\t') ,rows(2:end), 'un', 0);
    protName=cellfun(@(x) extractBefore(x{2}, '(EC'), rows, 'un', 0);
end