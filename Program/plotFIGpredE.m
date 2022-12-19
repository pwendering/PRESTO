function plotFIGpredE(prspredE, gkopredE, condNames, fileprfx, topDir)
% Function to plots 
%-  comparison of predicted growth in PRESTO and all GECKO models for each
%   constraint szenario seperately (FIG 2; FIG 4)
% 
%INPUT:
% -cell prspredE:   cell array of length two containing  the
%                   spearman correlation between measured and predicted
%                   enzyme abundance based on the PRESTO model with
%                   no default uptake bounds, fixed
%                   measured growth rate and minimization of the total
%                   enzyme pool
%                       1: only the proteins for which measurements are
%                           available in all conditions
%                       2: all enzymes in the model (unmeasured enzymes in
%                       the condition are assumed to have abundance 0)
% -cell gkorelE:    cell array of length two of spearman correlations for the 
%                   GECKO models, same content
%                   as om prspredE but each cell contains a matrix
%                   1modelx27conditions (the diagonal corresponds to
%                   the input of prspredE)
% -cell condNames: cell of character vectors giving the condition names
% -char fileprfx: character vector giving the prefix for saved svg files
% -char topDir: directory to save results in

resdir=fullfile(topDir, 'Results', 'relE');
if ~isdir(resdir)
    mkdir(resdir)
end
%close any open figures
close all

%% plot a point diagram of all condition specific models of GECKO in all conditions and 1 point for GECKO
%obtain condition names sorting according to PRESTO performance
[~, order]=sort(abs(prspredE{1}));
% x_condNames=categorical(condNames);
% x_condNames=reordercats(x_condNames, condNames(order))
figure

tiledlayout(2,1)
set(0, 'DefaultAxesTickLabelInterpreter', 'none')
for d=1:2
    nexttile
    hold on
    boxplot(abs(gkopredE{d}(:, order)))
    %change line style
    all_lines= findobj(gca, 'Type', 'Line');
    arrayfun( @(x) set(x, 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1), all_lines)
    %delete outliers 
    outliers=findobj(gca, 'Tag', 'Outliers');
    delete(outliers)
    
    %To change box aesthetics
    %myboxes = findobj(ax,'Tag','Box')
    %arrayfun( @(box) patch( box.XData, box.YData, 'm', 'FaceAlpha', 0.5), myboxes(1:5) )
    for j=1:length(condNames)
        %scatter(repelem(j, size(gkorelE{2},1)), gkorelE{2}(:,j), 30, 'Jitter', 'on', 'JitterAmount', 0.1, 'MarkerEdgeColor', 'k')
        pl(j)=scatter(repelem(j, size(gkopredE{d},1)), abs(gkopredE{d}(:,order(j))), 15, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.6, 'Jitter', 'on', 'JitterAmount', 0.2);
    end
%     for i=1:size(gkorelE{d},1)
%         pl(i)=scatter(1:numel(condNames), abs(gkorelE{d}(i,order)), [],'SizeData', 10,...
%             'MarkerEdgeColor', 'none',  'MarkerFaceColor', 'black', 'Marker', 'o',  'MarkerFaceAlpha', 0.35);
%     end
    pl(length(condNames)+1)=scatter(1:numel(condNames),abs(prspredE{d}(order)), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor','none','Marker', 'd', 'SizeData', 30, 'MarkerFaceAlpha', 0.8);
    ylabel('relative error')
    %create padding in y axis
    pad=0.025;
    maxylim=max(abs([gkopredE{d}, prspredE{d}]), [], 'all', 'omitnan');
    maxylim=maxylim*(1+pad);
    minylim=min(abs([gkopredE{d}, prspredE{d}]), [], 'all', 'omitnan')-maxylim*pad;
    ylim([minylim, maxylim]);
    set(gca, 'Box', 'off')
    if (d~=2)
        set(gca,'xticklabel',[])
    else
        xticks(1:numel(condNames))
        xticklabels(condNames(order))
        xtickangle(75)
        l=legend([pl(1) pl(numel(condNames)+1)], 'GECKO', 'PRESTO');
        l.Location='southeast';
    end
    hold off
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 8]);
saveas(gcf, fullfile(resdir, [fileprfx '_condpredE.svg']))
end