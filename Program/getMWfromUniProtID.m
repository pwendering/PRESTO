function MW = getMWfromUniProtID(IDs)
%% MW = getMWfromUniProtID(IDs)
% Use the UniProt API to retrieve molecular weights for each of the
% supplied IDs
% Input:
%   cellstr IDs:        UniProt identifiers
% Output:
%   double MW:          molecular weights [g/mol]
% construct request
%
% 22.03.2022 Philipp Wendering, University of Potsdam, philipp.wendering@gmail.com
baseURL = 'https://www.uniprot.org';
tool = 'uploadlists';

startIdx = 1:100:numel(IDs);
for i=1:numel(startIdx)
    if i<numel(startIdx)
        endIdx = startIdx(i)+99;
    else
        endIdx = numel(IDs);
    end
    fprintf('Processing IDs %d to %d ...\n',startIdx(i),endIdx)
    
    url = [baseURL '/' tool '/',...
        '?query=', strjoin(IDs(startIdx(i):endIdx),','),...
        '&format=tab',...
        '&from=ACC',...
        '&to=ACC'...
        '&columns=mass'...
        ];
    % send request
    data = webread(url);
    % parse response
    rows = strsplit(strtrim(data),'\n');
    rows = cellfun(@strsplit,rows,'un', 0);
    MW(startIdx(i):endIdx) = cellfun(@(x)str2double(x(1)),rows(2:end))';
end
end