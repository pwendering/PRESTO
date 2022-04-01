close all
clear;clc

% S. cerevisiae
ecYeastGEM = readGKOmodel(fullfile('Data','rawecYeast.mat'));
enzMetIdx = startsWith(ecYeastGEM.mets,'prot_');
kcat_matrix = -1./ecYeastGEM.S(enzMetIdx,:)./3600;
kcat_matrix(isinf(kcat_matrix)) = NaN;

rxns_per_enzyme_sce = sum(kcat_matrix>0,2);
uniq_k_per_enzyme_sce = arrayfun(@(i)numel(unique(kcat_matrix(i,kcat_matrix(i,:)>0))),...
    1:size(kcat_matrix,1));

% E. coli
ecEcoli = readGKOmodel(fullfile('Data','rawecEcoli.mat'));
enzMetIdx = startsWith(ecEcoli.mets,'prot_');
kcat_matrix = -1./ecEcoli.S(enzMetIdx,:)./3600;
kcat_matrix(isinf(kcat_matrix)) = NaN;
protein_ids = strrep(ecEcoli.mets(enzMetIdx),'prot_','');

rxns_per_enzyme_eco = sum(kcat_matrix>0,2);
uniq_k_per_enzyme_eco = arrayfun(@(i)numel(unique(kcat_matrix(i,kcat_matrix(i,:)>0))),...
    1:size(kcat_matrix,1));

% plot histograms
tl = tiledlayout(2,2,'TileSpacing', 'compact');

nexttile
histogram(rxns_per_enzyme_sce,'BinMethod','integers','NumBins',20,...
    'FaceColor','#6999e6','BinLimits',[1 20])
fprintf('S. cerevisiae: number of reaction per enzyme > 20: %d\n',sum(rxns_per_enzyme_sce>20))
box on
set(gca,'FontSize',10,'FontName','Arial','LineWidth',1.5)
text(.9,.9,'a','FontSize',14,'FontName','Arial','Units','normalized')
% text(.4,.9,sprintf('n(> 20) = %d',sum(rxns_per_enzyme_sce>20)),...
%     'FontSize',10,'FontName','Arial','Units','normalized');

nexttile
histogram(uniq_k_per_enzyme_sce,'BinMethod','integers',...
    'FaceColor','#6999e6')
box on
set(gca,'FontSize',10,'FontName','Arial','LineWidth',1.5)
text(.9,.9,'b','FontSize',14,'FontName','Arial','Units','normalized')

nexttile
histogram(rxns_per_enzyme_eco,'BinMethod','integers','NumBins',20,...
    'FaceColor','#6999e6','BinLimits',[1 20])
fprintf('E. coli number of reaction per enzyme > 20: %d\n',sum(rxns_per_enzyme_eco>20))
xlabel('Number of reactions per protein')
box on
set(gca,'FontSize',10,'FontName','Arial','LineWidth',1.5)
text(.9,.9,'c','FontSize',14,'FontName','Arial','Units','normalized')
% text(.4,.9,sprintf('n(> 20) = %d',sum(rxns_per_enzyme_eco>20)),...
%     'FontSize',10,'FontName','Arial','Units','normalized');

nexttile
histogram(uniq_k_per_enzyme_eco,'BinMethod','integers',...
    'FaceColor','#6999e6')
xlabel('Number of unique k_{cat} values per protein')
box on
set(gca,'FontSize',10,'FontName','Arial','LineWidth',1.5)
text(.9,.9,'d','FontSize',14,'FontName','Arial','Units','normalized')

ylabel(tl,'Count','FontSize',14,'FontName','Arial')

set(gcf,'OuterPosition',[-859.0000  164.3333  806.0000  473.3333])

exportgraphics(gcf,'Results/Fig_S8.tiff','Resolution',300)