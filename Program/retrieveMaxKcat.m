function K = retrieveMaxKcat(filename,query,dlm)
%% K = retrieveMaxKcat(filename,query,dlm)
% Find the maximum in vitro turnover number measured for a given query
% (can be organism name, higher taxonomic rank, substrate, E.C. number).
% Input:
%   char filename:  path to file that has at least four columns:
%                   EC number | substrate | taxonomic classification | kcat
%                   (format is used for kcat matching in the GECKO toolbox)
%   char query:     (optional) organism name, higher taxonomic rank, substrate,
%                   or E.C. number for which the maximum kcat should be
%                   retrieved; default: '' ==> returns overall maximum
%   char dlm:       (optional) delimiter; default: TAB
% Output:
%   double K:       maximum turnover number [s^-1]
% 
% 22.03.2022 Philipp Wendering, University of Potsdam, philipp.wendering@gmail.com

% check input
if ~ischar(filename)
    error('Filename must be of type char')
elseif nargin < 2
    query = '';
end

if iscellstr(query)
    query = char(query);
elseif ~ischar(query)
    error('Query must be of type char')
end

% default delimiter
if nargin < 3
    dlm = '\t';
end
% read file
fid = fopen(filename,'r');
data = textscan(fid,'EC%s %s %s %f %*[^\n]','Delimiter',dlm,'EmptyValue', NaN);
fclose(fid);
% scan text columns for query
idx = cell2mat(arrayfun(@(i)count(data{i},query,'ignorecase',1),1:numel(data)-1,'un',0));
if sum(any(idx)) >= 1
    % if the query was found in multiple columns, choose the one where most
    % hits were found
    [~,I] = max(sum(idx));
    idx = logical(idx(:,I));
    K = max(data{4}(idx));
else
    fprintf('retrieveMaxKcat: No matches found!\n')
    K = NaN;
end