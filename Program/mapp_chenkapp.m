function mapp_chenkapp()
%function that exports a xlsx file to Data/Chen/Chen_mapping.xlsx 
%with a fuzzy mapping of reaction names
%from chen  et al to the current yeast model for manual processing. First
%column in the file is the reaction name from chen et al. second column
%gives all reaction names in the current yeast model that have a perfect
%matching substring, third column gives all reaction names that have a substring 
% with minimal levensthein distance to the davide et al reaction name,
% fourht column already contains a final mapping if only 1 reaction with
% perfect matching substring is found
%read apparent catalytic rate data from davidi et al. 
configuration_yeast
kapp=readtable(fullfile(topDir, 'Data', 'Chen', 'chen_kmax.xlsx'));

%read basic model 
model=readGKOmodel(modelFile);
%match reactions using perfect and fuzzy matching of substringst
match=zeros(length(model.rxns), size(kapp, 1));
fzmatch=zeros(length(model.rxns), size(kapp, 1));
for i=1:size(kapp,1)
    query=regexprep(erase(kapp.ReactionIDInYeast8{i}, '_fwd'), '_rvs', '_REV');
    match(:,i)=~cellfun(@isempty, regexp(model.rxns, [query, 'No\d']));
    if sum(match(:,i))==0
        %if no match is found check in reactions without Numbering ( 
        match(:,i)=strcmp(model.rxns, query);
    %if not match with regex is found run fzsearch
        if sum(match(:,i))==0
            for j=1:length(model.rxns)
                min_dist=fzsearch(model.rxns{j},query);
                %min_dist=min_dist{1};
                fzmatch(j,i)=min_dist(1);
            end 
        end
    end
end

%create report cell array first column davidi reaction name, second colunn
%perfect matches of substrings, third column model reaction abbreviations
%with shortes levensthein distance
report=cell(size(kapp,1),4);
report(:,1)=kapp.ReactionIDInYeast8;
for i=1:size(kapp,1)
    %if there is a perfect match assign final match
    if sum(match(:,i))~=0
    report{i,2}=strjoin(model.rxns(logical(match(:,i))),'; ');
    report{i,4}=strjoin(model.rxns(logical(match(:,i))),'; ');
    else
    %report fuzzy matched reactions    
    report{i,3}=strjoin(model.rxns(fzmatch(:,i)==min(fzmatch(:,i))), '; ');
    end
end
if any(~cellfun(@isempty, report(:,3)))
    warning('not all reaction names could be matched automatically manual modifications might be nescessary')
end
report=cell2table(report, 'VariableNames', {'reactionName', 'perfect_match', 'min_levensthein_match', 'final_match'});
writetable(report, fullfile(topDir, 'Data', 'Chen','chen_mapping.xlsx'))

end