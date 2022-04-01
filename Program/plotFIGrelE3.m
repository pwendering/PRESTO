function plotFIGrelE3(prsrelE, max_gkorelE, pfbarelE, filenm, topDir, pfbalab)
%Function to plot a scatter comparison between the predictions of the
%PRESTO, the maximum GECKO and the pFBA supplemented model FIG S8
%INPUT:
% -cell prsrelE: cell array of length four containing the perfomance of 
%                       PRESTO models as vector 1modelx1condition the
%                       each cell contains solution for a set of constrains
%                       1: only  pool constrain
%                       2: pool constrain + experimental uptake rates
%                       3: pc + exp. upt. rates + proein abundance
% -cell max_gkorelE: cell array of of relative errors for max GECKO 
%                       model, same format
%                        as om prsrelE 
% -cell pfbarelE:cell array of of relative errors for pFBA corrected 
%                       model, same format
%                        as om prsrelE 
% -cell condNames: cell of character vectors giving the condition names
% -char filenm: character vector giving the file name under which figure
%               will be saved to
%               topDir/Results/Manuscript/Fig1/<filenm>.svg
% -char topDir: directory to save results in
% -char pfbalab: label for pfba model predictions

resdir=fullfile(topDir, 'Results', 'relE');
if ~isdir(resdir)
    mkdir(resdir)
end
%close any open figure
close all
%% plot scatterplot per model prediction presto vs geckomax and presto vs pFBA 

figure
tiledlayout (1,3);
%t.Units='inches';
%t.OuterPosition=[0.25, 0.25, 15,5];
for d=1:3
    gkorelE_plot=abs(max_gkorelE{d});
    pfbarelE_plot=abs(pfbarelE{d});
    colos=lines(2);
    axis_min=Inf;
    axis_max=-Inf;
    nexttile
    hold on
    scatter(abs(prsrelE{d}), ...
            gkorelE_plot,[],colos(1,:), 'fill')
    scatter(abs(prsrelE{d}), pfbarelE_plot, [], colos(2,:), 'fill')
    axis_min=min(axis_min, min([abs(prsrelE{d}); gkorelE_plot; pfbarelE_plot]));
    axis_max=max(axis_max, max([abs(prsrelE{d}); gkorelE_plot; pfbarelE_plot]));
    hold off
    %set equal axis ranges
    axis([axis_min, axis_max, axis_min, axis_max])
    ticks=0:round(axis_max/10,1):axis_max;
    yticks(ticks);
    xticks(ticks);
    refline(1,0);
    if d==3 %only add to the last plot
        leg=legend({'GECKO', pfbalab} , 'Location', 'southeast');
        set(leg, 'Box', 'off');
    end
    xlabel('relative error PRESTO')
    ylabel('relative error compared model')
    %title(titles(d))
end
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 9 3]);
saveas(gcf, fullfile(resdir, [filenm '.svg']))
