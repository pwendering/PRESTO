function plotFIGKcatOvlp(ovlp, condNames, fileprfx, topDir)
%Function to plot the overlap betweeen kcat predictions in the different
%models as heatmap
%INPUT:
% double ovlp:  acutally boolean matrix :D kcat_corrections x conditions as
%               returned by comp_GKOcondmod.m
% cell condNames:   A cell array of character vectors giving the condition
%                   names for each column of ovlp
% char fileprfx:    A character vector giving the prefix for the saved svg
%                   files
% char topDir:  path to the working directory
resdir=fullfile(topDir, 'Results', 'kcats');

if ~isdir(resdir)
    mkdir(resdir)
end

%Order the presence absence matrix rows (specific corrected kcats)
[~, rorder] = sort(sum(ovlp, 2));


%Order the presence absence matrix columns (conditions by sum 
[~, corder] = sort(sum(ovlp,1));

%plot heatmap
h=heatmap(condNames(flip(corder)), 1:size(ovlp,1),ovlp(flip(rorder), flip(corder)), 'CellLabelColor', 'none', 'ColorbarVisible',false);
h.NodeChildren(3).XAxis.TickLabelInterpreter ='none';
set(0, 'DefaultAxesTickLabelInterpreter', 'none')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 10 6]);
saveas(gcf, fullfile(resdir, [fileprfx '_gkokcat_hm.svg']))
end