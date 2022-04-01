clear;clc;close all
configuration_yeast
ecYeastGEM = readGKOmodel(modelFile);
% find indices of enzyme usage pseudoreactions
enzRxnIdx = find(contains(ecYeastGEM.rxns,enzMetPfx));
enzMetIdx = find(any(ecYeastGEM.S(:,enzRxnIdx),2));
ecYeastGEM.enzymes = cellfun(@(x)regexp(x,'[A-Z0-9]{6}','match'),ecYeastGEM.metNames(enzMetIdx));
if exist(mwFile,'file')
    ecYeastGEM.protMW = table2array(readtable(mwFile,'ReadVariableNames',false));
else
    ecYeastGEM.protMW = 0;
end

if ~exist(mwFile,'file') || numel(ecYeastGEM.protMW)~=numel(ecYeastGEM.enzymes)
    ecYeastGEM.protMW = getMWfromUniProtID(ecYeastGEM.enzymes);
    writetable(array2table(ecYeastGEM.protMW),mwFile,'WriteVariableNames',false)
end

[condNames, E, expVal, nutrExch, P]=readChenetal(topDir, ecYeastGEM);

idxExclude = contains(condNames,{'Lahtvee2017_Temp33',...
    'Lahtvee2017_Temp36','Lahtvee2017_Temp38'});
expVal = expVal(~idxExclude);
P = P(~idxExclude);
E = E(:,~idxExclude);
nutrExch = nutrExch(:,~idxExclude);
condNames = condNames(~idxExclude);

%% PCA growth rates by experimental conditions
subplot(1,2,1)
[~,score,~,~,explained] = pca(...
    [P' table2array(nutrExch)']...
);

rxnNames = ecYeastGEM.rxnNames(findRxnIDs(ecYeastGEM,nutrExch.Properties.RowNames));

studies = strtok(condNames,'_');
uniqStudies = unique(studies,'stable');
uniqColors = {'#cb6b47','#79a251','#be5ba7','#7b98df'};

plotColors = cellfun(@(x)uniqColors(ismember(uniqStudies,x)),studies);

for i=1:size(score,1)
    scatter(score(i,1),score(i,2),80,'MarkerFaceColor',plotColors{i},...
        'MarkerEdgeColor',plotColors{i})
    hold on
end

ylim([min(score(:,2))-abs(0.2*min(score(:,2))) max(score(:,2))+0.2*abs(max(score(:,2)))]);
xlim([min(score(:,1))-abs(0.2*min(score(:,1))) max(score(:,1))+0.2*abs(max(score(:,1)))]);

xlabel(sprintf('PC1 (%.2f %%)',explained(1)))
ylabel(sprintf('PC2 (%.2f %%)',explained(2)))

text(0.03,0.95,'a','units','normalized','FontSize',16,'FontWeight','bold',...
    'FontName','Arial')

h = zeros(numel(uniqStudies),1);
for i=1:numel(uniqStudies)
    h(i) = plot(NaN,NaN,'.','MarkerSize', 30, 'Color',uniqColors{i});
end
legend(h,uniqStudies,'box','off','Position',[0.16 0.68 0.18 0.22])

box on 
set(gca,'FontSize',14,'LineWidth',1.3)


%% PCA with protein abundances
subplot(1,2,2)

[~,score,~,~,explained] = pca(...
    E'...
);

studies = strtok(condNames,'_');
uniqStudies = unique(studies,'stable');
uniqColors = {'#cb6b47','#79a251','#be5ba7','#7b98df'};

plotColors = cellfun(@(x)uniqColors(ismember(uniqStudies,x)),studies);

for i=1:size(score,1)
    scatter(score(i,1),score(i,2),80,'filled',...
        'MarkerFaceColor',plotColors{i})
    hold on
end
hold off

xlabel(sprintf('PC1 (%.2f %%)',explained(1)))
ylabel(sprintf('PC2 (%.2f %%)',explained(2)))

box on 
set(gca,'FontSize',14,'LineWidth',1.3)
text(0.03,0.95,'b','units','normalized','FontSize',16,'FontWeight','bold',...
    'FontName','Arial')

set(gcf,'OuterPosition',[0.2157    0.2150    1.0573    0.4940]*1000)

exportgraphics(gcf,fullfile('Results/PCA_conditions.tiff'),'Resolution',300)