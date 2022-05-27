function[out_tab] = get_kcatori(reportTable, basic_batchmod, filepath)
%Function that fits the origin code from matchkcats.m to the reactions in
%basic_batchmod with the ids given in reportTable.("Reaction ID")
%INPUT
% char filepath:    the path and prefix of the file under which a boxplot
%                   should be saved
if contains(basic_batchmod.name, 'Escherichia coli', 'IgnoreCase', true)
    configuration_ecoli
elseif contains(basic_batchmod.name, 'Yeast', 'IgnoreCase', true)
    configuration_yeast
end
    orig_tab=readtable(kcatoriginFile);
orig_tab.Properties.VariableNames(1)={'rxns'};
origin=nan(size(reportTable, 1), 1);
for i=1:size(reportTable, 1)
    recn=basic_batchmod.rxns(reportTable.("REACTION ID")(i));
    react=extractBefore(recn, 'No');
    isoz_id=cellfun(@str2num, extractAfter(recn, 'No'));
    if sum(ismember(orig_tab.rxns, react)>1)
        error('ambiguos name fitting from ecModel to normal model, rewrite code')
    elseif sum(ismember(orig_tab.rxns, react))==0
        warning(['Reaction ' react{1} ' could not be fitted. setting origin to NA'])
    end
    origin(i)=table2array(orig_tab(ismember(orig_tab.rxns, react), isoz_id+1));
end
out_tab=reportTable;
out_tab.Origin=origin;
if sum(out_tab.Origin==6)>0
    disp([num2str(sum(out_tab.Origin==6)), 'entries belonging to unspecific SA group have not been plotted'])
end
plot_tab=out_tab(out_tab.Origin~=6,:);
boxplot(plot_tab.("FOLD-INCREASE"), plot_tab.Origin)
set(gca, 'YScale', 'log')
tblstat=grpstats(plot_tab(:,{'FOLD-INCREASE', 'DELTA', 'Origin'}), 'Origin');
org_code={'org+subs', 'subs', 'org', 'any', 'SA org'};
xticklabels(strcat(org_code, ' n=', cellstr(string(tblstat.GroupCount))'))
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 4]);
saveas(gcf, [filepath, 'sigmabyorigin.svg'])
end
    