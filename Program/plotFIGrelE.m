function plotFIGrelE(prsrelE, gkorelE, condNames, fileprfx, topDir, maxmin_relE)
%Function to plots 
%-  comparison of the predicted growth rate for PRESTO and 
%   for the condition specific GECKo model in the differen constraint
%   scenarios
%-  comparison of predicted growth in PRESTO and all GECKO models for each
%   constraint szenario seperately (FIG 2; FIG4)
%-  comparison of predicted growth in PRESTO and condition specific GECKO
%   models. 
% conditions 
%INPUT:
% -cell prsrelE: cell array of length four containing the perfomance of 
%                       PRESTO models as vector 1modelx1condition the
%                       each cell contains solution for a set of constrains
%                       1: only  pool constrain
%                       2: pool constrain + experimental uptake rates
%                       3: pc + exp. upt. rates + proein abundance
%                       4: only experimental uptake rates (pool=1000)
%                       (non-sparse)
% -cell gkorelE: cell array of of relative errors for GECKO, same content
%                        as om prsrelE but each cell contains a matrix
%                        1modelx27conditions (the diagonal corresponds to
%                        the input of prsrelE)
% -cell condNames: cell of character vectors giving the condition names
% -char fileprfx: character vector giving the prefix for saved svg files
% -char topDir: directory to save results in
% -cell maxmin_relE:    (optional, default {})cell array of the relative errors for a PRESTO
%                       model with additionally down corrected kcats same
%                       format as prsrelE.

if nargin<6
    maxmin_relE={};
end
resdir=fullfile(topDir, 'Results', 'relE');
if ~isdir(resdir)
    mkdir(resdir)
end
%close any open figures
close all
%% Plot a point diagram of relative errors only comparing condition specific
%%modelling results
%For PRESTO
%marker sympols
mk={'.', '+', 's', 'd'};
hold on
for d=1:4
    scatter(categorical(condNames), prsrelE{d}, [], mk{d})
end
hold off
yline(0)
legend({'p', 'p;exp.u', 'p;eu;abun', 'exp.u.', 'perfect fit'})
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 5]);
saveas(gcf, fullfile(resdir, [fileprfx '_constcomp_prst.svg']))
% For GECKO
figure
hold on
for d=1:4
    scatter(categorical(condNames), diag(gkorelE{d}), [], mk{d})
end
hold off
yline(0)
legend({'p', 'p;exp.u', 'p;eu;abun', 'exp.u.', 'perfect fit'})
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 5]);
saveas(gcf, fullfile(resdir, [fileprfx '_constcomp_gecko.svg']))
%% plot a point diagram of all condition specific models of GECKO in all conditions and 1 point for GECKO
%obtain condition names sorting according to PRESTO performance
[~, order]=sort(abs(prsrelE{3}));
% x_condNames=categorical(condNames);
% x_condNames=reordercats(x_condNames, condNames(order))
figure

tiledlayout(3,1)
set(0, 'DefaultAxesTickLabelInterpreter', 'none')
for d=1:3
    nexttile
    hold on
    boxplot(abs(gkorelE{d}(:, order)))
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
        pl(j)=scatter(repelem(j, size(gkorelE{d},1)), abs(gkorelE{d}(:,order(j))), 15, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.6, 'Jitter', 'on', 'JitterAmount', 0.2);
    end
%     for i=1:size(gkorelE{d},1)
%         pl(i)=scatter(1:numel(condNames), abs(gkorelE{d}(i,order)), [],'SizeData', 10,...
%             'MarkerEdgeColor', 'none',  'MarkerFaceColor', 'black', 'Marker', 'o',  'MarkerFaceAlpha', 0.35);
%     end
    pl(length(condNames)+1)=scatter(1:numel(condNames),abs(prsrelE{d}(order)), [], 'MarkerFaceColor', 'red', 'MarkerEdgeColor','none','Marker', 'd', 'SizeData', 30, 'MarkerFaceAlpha', 0.8);
    if ~isempty(maxmin_relE)
        pl(length(condNames)+2)=scatter(1:numel(condNames),abs(maxmin_relE{d}(order)), [], 'MarkerFaceColor', 'cyan', 'MarkerEdgeColor','none','Marker', 'd', 'SizeData', 30, 'MarkerFaceAlpha', 0.8);
        vals=abs([gkorelE{d}, prsrelE{d}, maxmin_relE{d}]);
    else
        vals=abs([gkorelE{d}, prsrelE{d}]);
    end
    ylabel('relative error')
    %create padding in y axis
    pad=0.025;
    maxylim=max(vals, [], 'all', 'omitnan');
    maxylim=maxylim*(1+pad);
    minylim=min(vals, [], 'all', 'omitnan')-maxylim*pad;
    ylim([minylim, maxylim]);
    set(gca, 'Box', 'off')
    if (d~=3)
        set(gca,'xticklabel',[])
    else
        xticks(1:numel(condNames))
        xticklabels(condNames(order))
        xtickangle(75)
        if isempty(maxmin_relE)
            l=legend([pl(1) pl(numel(condNames)+1)], 'GECKO', 'PRESTO');
        else
            l=legend([pl(1) pl(numel(condNames)+1) pl(numel(condNames)+2)], 'GECKO', 'PRESTO','PRESTO + min');
        end
        l.Location='southeast';
    end
    hold off
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 8]);
saveas(gcf, fullfile(resdir, [fileprfx '_condp.svg']))
end