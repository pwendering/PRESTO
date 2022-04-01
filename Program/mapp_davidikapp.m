function mapp_davidikapp()
%function that exports a xlsx file to Data/davidi_mapping.xlsx 
%with a fuzzy mapping of reaction names
%from davidi et al to the current e_coli model for manual processing. First
%column in the file is the reaction name from davidi et al. second column
%gives all reaction names in the current e coli model that have a perfect
%matching substring, third column gives all reaction names that have a substring 
% with minimal levensthein distance to the davide et al reaction name,
% fourht column already contains a final mapping if only 1 reaction with
% perfect matching substring is found
%read apparent catalytic rate data from davidi et al. 
configuration_ecoli
dav_kapp=readtable(fullfile(topDir, 'Data', 'Davidi', 'davidi_suppl.xlsx'), 'Sheet', 'kmax 1s');

%read basic model 
model=readGKOmodel(modelFile);
%match reactions using perfect and fuzzy matching of substringst
match=zeros(length(model.rxns), size(dav_kapp, 1));
fzmatch=zeros(length(model.rxns), size(dav_kapp, 1));
for i=1:size(dav_kapp,1)
    match(:,i)=contains(model.rxns, dav_kapp.reaction_modelName_{i}, 'IgnoreCase', true);
    for j=1:length(model.rxns)
        min_dist=fzsearch(model.rxns{j},dav_kapp.reaction_modelName_{i}, 0, 1);
        min_dist=min_dist{1};
        fzmatch(j,i)=min_dist(1);
    end 
end

%create report cell array first column davidi reaction name, second colunn
%perfect matches of substrings, third column model reaction abbreviations
%with shortes levensthein distance
report=cell(size(dav_kapp,1),4);
report(:,1)=dav_kapp.reaction_modelName_;
for i=1:size(dav_kapp,1)
    perfmatch=model.rxns(logical(match(:,i)));
   
    %if there is only one perfect match already assign final match
    if length(perfmatch)==1
        report(i,2)=perfmatch;
        report(i,4)=perfmatch;
    %if there are two perfect matches and they only differ in '_REV' take
    %the forward reaction as final match
    elseif length(perfmatch)==2 & strcmp(regexprep(perfmatch(1), '_REV', ''), ...
            regexprep(perfmatch(2), '_REV', ''))
         report{i,2}=strjoin(model.rxns(logical(match(:,i))), '; ');
         report(i,4)=perfmatch(~contains(perfmatch, '_REV'));
    else
        report{i,2}=strjoin(model.rxns(logical(match(:,i))), '; ');
    end
    report{i,3}=strjoin(model.rxns(fzmatch(:,i)==min(fzmatch(:,i))), '; ');
end
report=cell2table(report, 'VariableNames', {'reactionName', 'perfect_match', 'min_levensthein_match', 'final_match'});
writetable(report, fullfile(topDir, 'Data', 'Davidi','davidi_mapping.xlsx'))
end
