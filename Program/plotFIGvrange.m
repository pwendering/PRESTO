function plotFIGvrange(fluxvar, maxgkofluxvar, batch_model, enzRxnPfx, condNames, fileprfx, topDir)
%Function to plot a boxplot of the distribution of flux ranges for all
%metabolic raections in PRESTO and GECKO models obtained from FBA
    %test if flux ranges are comparable
    resdir=fullfile(topDir, 'Results', 'relE');
    if ~isdir(resdir)
        mkdir(resdir)
    end

    if ~isequal(size(fluxvar), size(maxgkofluxvar))
        error('FVA results from GECKO and PRESTO models have different dimensions. Aborting...')
    end
    %reshape to boxchart input
    plotdat=table(reshape(fluxvar(:,~contains(batch_model.rxns, enzRxnPfx))',[],1),...
        repelem(condNames, sum(~contains(batch_model.rxns, enzRxnPfx)))',...
        repmat({'PRESTO'}, sum(~contains(batch_model.rxns, enzRxnPfx))*size(fluxvar,1),1));
    plotdat=[plotdat;
        table(reshape(maxgkofluxvar(:,~contains(batch_model.rxns, enzRxnPfx))',[],1),...
        repelem(condNames, sum(~contains(batch_model.rxns, enzRxnPfx)))',...
        repmat({'GECKO'}, sum(~contains(batch_model.rxns, enzRxnPfx))*size(fluxvar,1),1))];

    plotdat.Properties.VariableNames={'Range', 'Condition', 'Modeltype'};
    %Transform ranges to log10 scale
    plotdat.Range=log10(plotdat.Range);
    close all
    boxchart(categorical(plotdat.Condition), plotdat.Range, 'GroupByColor', categorical(plotdat.Modeltype))
    legend('Location', 'southeast')
    ylabel('log_{10}(Flux range)')
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 7 5]);
    saveas(gcf, fullfile(resdir, [fileprfx '_fluxranges.svg']))
end