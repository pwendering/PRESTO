function overvTable = comp_kcat2(prst_mod, gko_mod, by, figprefix)
%function to compare presto kcat modifications, with 
%GECKO automatic modifications and GECKO manual modifications and conduct
%pathway enrichment analysis.
%INPUT:
% - struc prst_mod: model structure generatred by PRESTO
% - struc gko_mod: model structure generated by GECKO to compare kcats
% - char gko_bestmod: A string giving the condition name of the best
%                                   performing GECKO model
% - char by: character vector giving the information by which kcats should
%                  be compared. Either 'Protein' or 'both'. both means
%                  protein Id and reaction ID of the catalysed reaction is
%                  taken into account. The kcat pathway enrichment will
%                  always only take into account unique protein terms
%OUTPUT:
% - talbe overvTable: Table with first columns giving protein and 
%                               reaction ID (depending on by) for each changed kcat
%                               remaining columns contain boolean indexes 
%                               indicating wether the given kcat is
%                               correcte by the respective approach

if ~(ismember(by, {'Protein', 'both'}))
    error("Argument by can either be 'Protein' or 'both'")
end

%detect model organism 
if contains(prst_mod.name, 'Escherichia coli', 'IgnoreCase', true)
    configuration_ecoli
    kegg_orgcode='eco';
    %Get enzymes linked to model
    %ATTENTION PRESTO und GECKO mod have differntly ordered enzyme fields!!
    [~, E,~, ~, ~]=readDavidi2016(prst_mod,topDir);
elseif contains(prst_mod.name, 'Yeast', 'IgnoreCase', true)
    configuration_yeast
    kegg_orgcode='sce';
    [~, E,~, ~, ~]=readChenetal(prst_mod, prst_mod, true);
end
%Find enzymes measured in all conditions that where available for
%corrections by PRESTO (for hypergeometric test)
used_E=all(E,2);
%Create directory to save figures if it doesn't allready exist
resdir=fullfile(topDir, 'Results', 'kcats');
if ~isdir(resdir)
    mkdir(resdir)
end


basic_batchmod=readGKOmodel(batchModelFile);
%get GECKO  modifications
 [~, ~ , ~,gkoreportTable] = findKcatCorrections(basic_batchmod, gko_mod, enzRxnPfx);
 %convert to identical format
gkoreportTable.("PROTEIN ID")=string(table2array(gkoreportTable(:, 'PROTEIN ID')));
gkoreportTable.("REACTION ID")=string(table2array(gkoreportTable(:, 'REACTION ID')));
gkoreportTable.Properties.VariableNames(1:2)={'Protein', 'Reaction#'};

 %read manual modifications hard coded in gecko
 %get Imort options
 Ioptions=detectImportOptions(fullfile(topDir, 'Data', [orgBasename '_manmods.txt']));
 %change variable types
 Ioptions=setvartype(Ioptions, {'string', 'string'});
 manmods=readtable(fullfile(topDir, 'Data', [orgBasename '_manmods.txt']), Ioptions);
 manmods.Properties.VariableNames={'Protein', 'Reaction#'};
 manmods=unique(manmods);
 %get PRESTO kcat modifications
 basic_model=readGKOmodel(modelFile);
 [~,~,~,reportTable] = findKcatCorrections(basic_model, prst_mod, enzMetPfx);
% convert reaction and protein ids to strings
reportTable.("PROTEIN ID")=string(table2array(reportTable(:, 'PROTEIN ID')));
reportTable.("REACTION ID")=string(table2array(reportTable(:, 'REACTION ID')));
%adapt identical variable names
reportTable.Properties.VariableNames(1:2)={'Protein', 'Reaction#'};
%Check for paths affected by presto
if  exist('figprefix', 'var')
    KegginfTable=get_KEGGid(prst_mod.enzymes(used_E), kegg_orgcode);
    PWenrich(KegginfTable, unique(reportTable.Protein), fullfile(resdir, figprefix))
end
%report raw numbers of changed kcats
disp(['The total number of manual adjustments in S is ' string(size(manmods, 1))])
disp(['The total number of GECKO adjustments in S is ' string(size(gkoreportTable, 1))])
disp(['The total number of PRESTO adjustments in S is ' string(size(reportTable, 1))])
%add reaction names for visualization
%compile a overview table
%set column ids according to by argument
switch by
    case 'both'
        col_idx=1:2;
    case 'Protein'
        col_idx=1;
end
%only keep unique entries round to avoid doubled entries due to numerical
%precision issues
gkoreportTable=[gkoreportTable(:,1:2), array2table(round(table2array(gkoreportTable(:, 3:end)), 6,'significant'), 'VariableNames', gkoreportTable.Properties.VariableNames(3:end))];
gkoreportTable=unique(gkoreportTable(:, [col_idx 3:end]));
reportTable=[reportTable(:,1:2), array2table(round(table2array(reportTable(:, 3:end)), 6,'significant'), 'VariableNames', reportTable.Properties.VariableNames(3:end))];
reportTable=unique(reportTable(:, [col_idx 3:end]));
manmods=unique(manmods(:,col_idx));
row_names=gkoreportTable(:,col_idx);
row_names=[row_names; manmods(~(ismember(manmods(:,col_idx), row_names)),col_idx)];
row_names=[row_names; reportTable(~(ismember(reportTable(:,col_idx), row_names)), col_idx)];

%compile a table of boolean values indicating if a certain kcat was
%modified in GECKO, PRESTO or manual modification
overlap=[ismember(row_names, manmods(:, col_idx)), ismember(row_names, gkoreportTable(:,col_idx)), ismember(row_names, reportTable(:,col_idx))];
overvTable=[row_names, array2table(overlap, 'VariableNames', {'Manual_Modifications', 'GECKO_Modifications', 'PRESTO_Modifications'})];
[~,Idx]=sort(sum(table2array(overvTable(:, (col_idx(end)+1):end)), 2));
overvTable=overvTable(flip(Idx), :);

%compile a table with updated kcat values from PRESTO and GECKO
kcat_tab=zeros(size(overlap,1),2);
kcat_tab(logical(overlap(:, 2)),1)=gkoreportTable.("KCAT UPDATED [s^-1]")(ismember(gkoreportTable(:,col_idx), row_names(logical(overlap(:,2)),:)));
kcat_tab(logical(overlap(:, 3)),2)=reportTable.("KCAT UPDATED [s^-1]")(ismember(reportTable(:, col_idx), row_names(logical(overlap(:,3)),:)));
kcat_tab=[row_names, array2table(kcat_tab, 'VariableNames', {'GECKO_KCAT[s^-1]', 'PRESTO_KCAT[s^-1]'})];
kcat_tab=kcat_tab(any(table2array(kcat_tab(:, ismember(kcat_tab.Properties.VariableNames, {'GECKO_KCAT[s^-1]', 'PRESTO_KCAT[s^-1]'}))), 2),:);

if exist('figprefix', 'var')
%export results
writetable(overvTable, fullfile(resdir, [figprefix, '_', by, 'kcattab.tsv']), 'FileType', 'text', 'Delimiter', '\t')
writetable(kcat_tab, fullfile(resdir, [figprefix, '_', by,'kcattabvalues.tsv']), 'FileType', 'text', 'Delimiter', '\t')

%% plot log plot
plot_x=reportTable.("KCAT UPDATED [s^-1]")(ismember(reportTable(:,col_idx),...
    overvTable(overvTable.GECKO_Modifications & overvTable.PRESTO_Modifications, col_idx)));
plot_y=gkoreportTable.("KCAT UPDATED [s^-1]")(ismember(gkoreportTable(:,col_idx),...
    overvTable(overvTable.GECKO_Modifications & overvTable.PRESTO_Modifications, col_idx)));
[rho, pval]=corr(log10(plot_x), log10(plot_y), 'Type', 'Spearman');
figure
colo=colormap;
loglog(plot_x, plot_y,'o', 'MarkerFaceColor', colo(1,:))
refline(1,0)
text(min(plot_x)+0.5, max(plot_y)-0.5, {['Spearman \rho: ', num2str(round(rho, 3 , 'significant')), ' '], ...
    ['p-value: ', num2str(round(pval, 3, 'significant'))]})
xlabel('kcat PRESTO [s^{-1}]')
ylabel('kcat GECKO [s^{-1}]')
saveas(gcf, fullfile(resdir, [figprefix, 'kcatcomp.svg']))

% %% Compare with davidi data
% if contains(prst_mod.name, 'Escherichia coli', 'IgnoreCase', true)
%     %read Davidi data 
%     dav_kapp=readtable(fullfile(topDir, 'Data', 'Davidi', 'davidi_suppl.xlsx'), 'Sheet', 'kmax 1s');
%     %get reaction mapping information
%     dav_map=readtable(fullfile(topDir, 'Data', 'Davidi', 'davidi_mapping_manual.xlsx'));
%     %filter unmatched
%     dav_map=dav_map(~cellfun(@isempty, dav_map.final_match),[1,4]);
%     %duplicate isozyme entries. 
%     i=1;
%     while i<=size(dav_map,1)
%         if contains(dav_map.final_match(i), ';')
%             isozms=strsplit(dav_map.final_match{i}, '; ')';
%             tmp_map=dav_map(repelem(i, length(isozms)), :);
%             tmp_map.final_match=isozms;
%             dav_map=[dav_map(1:(i-1), :); tmp_map; dav_map((i+1):end,:)];
%             i=i+length(isozms);
%         else
%             i=i+1;
%         end
%     end
%     %compile table with data 
%     enzRxnIdx = find(contains(prst_mod.rxns,enzMetPfx)&~ismember(prst_mod.rxns, {'prot_pool_exchange'}));
%     enzMetIdx = find(contains(prst_mod.mets, enzMetPfx) &~ismember(prst_mod.mets, {'prot_pool'}));
%     kapp_comp=cell(size(dav_map, 1), 4);
%     kapp_comp(:,1)=dav_map.final_match;
%     prs_kcatsubmat=prst_mod.S(enzMetIdx,:);
%     gko_kcatsubmat=gko_mod.S(enzMetIdx,:);
%     for i=1:size(kapp_comp,1)
%         %retrieve davidi kapp
%         kapp_comp{i,2}=dav_kapp.kmaxPerPolypeptideChain_s_1_(ismember(dav_kapp.reaction_modelName_, dav_map.reactionName(i)));
%         %retrieve max presto kcat 
%         enz_idx=find(prs_kcatsubmat(:,ismember(prst_mod.rxns, dav_map.final_match(i))));
%         if length(enz_idx)~=1
%             if length(unique(prs_kcatsubmat(enz_idx,ismember(prst_mod.rxns, dav_map.final_match(i)))))==1
%                 %complexes should not be found but accept them if kcat is
%                 %identical
%                 enz_idx=enz_idx(1);
%             else
%                 warning(['reaction with missing enzyme or enzyme complex detected comparison impossible. skipping reaction ', ...
%                     prst_mod.rxns{ismember(prst_mod.rxns, dav_map.final_match(i))}])
%                 continue
%             end
%         end
%         kapp_comp{i,3}=-1./prs_kcatsubmat(enz_idx, ismember(prst_mod.rxns, dav_map.final_match(i)))/3600;
%         kapp_comp{i,4}=-1./gko_kcatsubmat(enz_idx, ismember(prst_mod.rxns, dav_map.final_match(i)))/3600;
%     end
%     %remove omitted reactions
%     kapp_comp=kapp_comp(~cellfun(@isempty, kapp_comp(:,3)),:);
%     kapp_comp=cell2table(kapp_comp, 'VariableNames', {'Reaction', 'Davidi_kapp', 'PRESTO_kcat', 'GECKO_kcat'})
%     %plot figure
%     figure
%     tiledlayout (1,2)
%     nexttile
%     scatter(log(kapp_comp.Davidi_kapp), log(kapp_comp.PRESTO_kcat))
%     refline(1,0)
%     title(['log_{10} s⁻1]; r=', num2str(round(corr(log10(kapp_comp.Davidi_kapp), log10(kapp_comp.PRESTO_kcat)),2))])
%     xlabel('Davidi et al.')
%     ylabel('PRESTO')
%     nexttile
%     scatter(log(kapp_comp.Davidi_kapp), log(kapp_comp.GECKO_kcat))
%     refline(1,0)
%     title(['log_{10} s⁻1]; r=', num2str(round(corr(log10(kapp_comp.Davidi_kapp), log10(kapp_comp.GECKO_kcat)),2))])
%     xlabel('Davidi et al.')
%     ylabel('GECKO')
%     saveas(gcf, fullfile(resdir, [figprefix, 'davidicomp.svg']))
% end
end






