function plotFIGrelE(prsrelE, gkorelE, condNames, fileprfx, topDir)
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
[~, order]=sort(prsrelE{3});
% x_condNames=categorical(condNames);
% x_condNames=reordercats(x_condNames, condNames(order))
figure
tiledlayout(3,1)
set(0, 'DefaultAxesTickLabelInterpreter', 'none')
for d=1:3
    nexttile
    hold on
    for i=1:size(gkorelE{d},1)
        pl(i)=scatter(1:numel(condNames), abs(gkorelE{d}(i,order)), [],'SizeData', 10,...
            'MarkerEdgeColor', 'none',  'MarkerFaceColor', 'black', 'Marker', 'o',  'MarkerFaceAlpha', 0.35);
    end
    pl(numel(condNames)+1)=scatter(1:numel(condNames),abs(prsrelE{d}(order)), [], 'red', '+');
    ylabel('relative error')
    if (d~=3)
        set(gca,'xticklabel',[])
    else
        xticks(1:numel(condNames))
        xticklabels(condNames(order))
        xtickangle(75)
        legend([pl(1) pl(numel(condNames)+1)], 'GECKO', 'PRESTO')
    end
    hold off
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 5 8]);
saveas(gcf, fullfile(resdir, [fileprfx '_condp.svg']))
%% plot scatterplot per condition, condition specific PRESTO agains
%condidtion specific GECKO
plot_cond=cellfun(@(x) strsplit(x, '_'), condNames, 'UniformOutput', false);
plot_cond=cellfun(@(x) x{1}, plot_cond, 'UniformOutput', false);
%[~, plot_cond]=ismember(plot_cond, unique(plot_cond));
un_plotcond=unique(plot_cond);
%titles  - not used in publication
%titles={'pool', 'pool+exp.upt.', 'pool+exp.upt.+prot.abun.'};
figure
tiledlayout (1,3);
%t.Units='inches';
%t.OuterPosition=[0.25, 0.25, 15,5];
if length(un_plotcond)<=7
    colos=lines(length(un_plotcond));
else
    colos=jet(length(un_plotcond));
end
for d=1:3
    gkorelE_plot=abs(diag(gkorelE{d}));
    
    axis_min=Inf;
    axis_max=-Inf;
    nexttile
    hold on
    for i=1:length(un_plotcond)
        scatter(abs(prsrelE{d}(ismember(plot_cond, un_plotcond(i)))), ...
            gkorelE_plot(ismember(plot_cond, un_plotcond(i))),[],colos(i,:), 'fill')
        axis_min=min(axis_min, min([abs(prsrelE{d}(ismember(plot_cond, un_plotcond(i)))); gkorelE_plot(ismember(plot_cond, un_plotcond(i)))]));
        axis_max=max(axis_max, max([abs(prsrelE{d}(ismember(plot_cond, un_plotcond(i)))); gkorelE_plot(ismember(plot_cond, un_plotcond(i)))]));
    end
    hold off
    %set equal axis ranges
    axis([axis_min, axis_max, axis_min, axis_max])
    ticks=0:round(axis_max/10,1):axis_max;
    yticks(ticks);
    xticks(ticks);
    refline(1,0);
    if d==3 %only add to the last plot
        [leg, matlabsucks]=legend(un_plotcond, 'Location', 'southeast');
        set(leg, 'Box', 'off');
    end
    xlabel('relative error PRESTO')
    ylabel('relative error GECKO')
    %title(titles(d))
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 9 3]);
saveas(gcf, fullfile(resdir, [fileprfx '_scat.svg']))
%% Plot with errorbars
figure
tiledlayout (1,3);
%t.Units='inches';
%t.OuterPosition=[0.25, 0.25, 15,5];
%colos={'blue', 'green', 'yellow', 'magenta'};
for d=1:3
    axis_min=Inf;
    axis_max=-Inf;
    nexttile
    hold on
    for i=1:length(un_plotcond)
        eb=errorbar(abs(prsrelE{d}(ismember(plot_cond, un_plotcond(i)))),mean(abs(gkorelE{d}(ismember(plot_cond, un_plotcond(i)),:)), 2, 'omitnan'), ...
            mean(abs(gkorelE{d}(ismember(plot_cond, un_plotcond(i)),:)), 2, 'omitnan')-min(abs(gkorelE{d}(ismember(plot_cond, un_plotcond(i)),:)),[], 2), ...
            max(abs(gkorelE{d}(ismember(plot_cond, un_plotcond(i)),:)), [],2)-mean(abs(gkorelE{d}(ismember(plot_cond, un_plotcond(i)),:)), 2, 'omitnan'),  ...
            'vertical', '.');
        eb.Color=colos(i,:);
        axis_min=min(axis_min, min([abs(prsrelE{d}(ismember(plot_cond, un_plotcond(i)))), abs(gkorelE{d}(ismember(plot_cond, un_plotcond(i)),:))],[], 'all'));
        axis_max=max(axis_max, max([abs(prsrelE{d}(ismember(plot_cond, un_plotcond(i)))), abs(gkorelE{d}(ismember(plot_cond, un_plotcond(i)),:))],[], 'all'));
    end
    hold off
    %set equal axis ranges
    axis([axis_min, axis_max, axis_min, axis_max])
    ticks=0:round(axis_max/10,1):axis_max;
    yticks(ticks);
    xticks(ticks);
    refline(1,0);
    if d==3 %only add to the last plot
        [leg, matlabsucks]=legend(un_plotcond, 'Location', 'southeast');
        set(leg, 'Box', 'off');
    end
    xlabel('relative error PRESTO')
    ylabel('relative error GECKO')
    %title(titles(d))
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 9 3]);
saveas(gcf, fullfile(resdir, [fileprfx '_eb.svg']))
end