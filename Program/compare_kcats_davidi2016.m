%% Compare corrected kcat values with estimated values from Davidi et al. 2016
clear;clc

%% Read data
% read S1 Dataset from Davidi et al. 2016
davidi2016_s1_dataset = readtable(...
    fullfile('Data','Davidi','pnas.1514240113.sd01.xlsx'),...
    'NumHeaderLines',2,'Sheet','kmax 1s');
kmax = davidi2016_s1_dataset.kmaxPerPolypeptideChain_s_1_;
rxns = davidi2016_s1_dataset.reaction_modelName_;

% modify Davidi2016 reaction IDs to enable comparison to GECKO model
rxns = strrep(rxns,'_reverse','_REV');
rxns(ismember(rxns,'ACOAD1f_REV')) = {'ACOAD1fr'};

% read eciML1515 model
eciML1515 = readGKOmodel(fullfile('Data','rawecEcoli.mat'));

% read iJO1366 model
load('Data/Xu/iJO1366.mat');

% load corrected models (PRESTO)
load(fullfile('Results','corrected_models_e_coli_1_e_minus_5'))

enzMetIdx = find(startsWith(eciML1515.mets,'prot_'));
kcat_max_brenda = nan(numel(rxns),1);
kcat_max_presto = nan(numel(rxns),1);

eciML1515.rxns = regexprep(eciML1515.rxns,'No\d+','');
ishomomeric = false(numel(rxns),1);
for i=1:numel(rxns)
    % find reaction IDs in the model that contain the current ID
    matchIdxML1515 = cellfun(@(x)strcmp(x,rxns{i}),eciML1515.rxns);
    matchIdxiJO1366 = cellfun(@(x)strcmp(x,rxns{i}),iJO1366.rxns);
    ishomomeric(i) = all(~cellfun(@isempty,iJO1366.grRules(matchIdxiJO1366)) & ~contains(iJO1366.grRules(matchIdxiJO1366),{'and','or'}));
    if sum(matchIdxML1515) > 0
        % determine max kcat in orgiginal model (from BRENDA)
        kcat_max_brenda(i) = max(max(-1./eciML1515.S(enzMetIdx,matchIdxML1515)/3600));
        % determine max kcat in corrected model (PRESTO)
        kcat_max_presto(i) = max(max(-1./corrModels{1}.S(enzMetIdx,matchIdxML1515)/3600));
    else
        fprintf('Reaction not found: %s\n', rxns{i})
    end
end
kcat_max_brenda(isinf(kcat_max_brenda)) = NaN;
kcat_max_presto(isinf(kcat_max_presto)) = NaN;

% filter for reactions catalyzed by homomeric enzyme
kcat_max_brenda(~ishomomeric)=NaN;
kcat_max_presto(~ishomomeric)=NaN;
kmax(~ishomomeric)=NaN;

%% Plot results
% fig_davidi2016_kcat_comp = figure;
% tiledlayout('flow')

% PRESTO vs Davidi
nexttile
line([-5 7],[-5 7], 'Color', 'k', 'LineWidth', 1)

hold on
scatter(log10(kmax),log10(kcat_max_presto),'filled')
xlabel('log_{10} k_{cat}^{max} (Davidi2016) [s^{-1}]')
ylabel('log_{10} k_{cat}^{max} (PRESTO) [s^{-1}]')
xlim([-5 7]); ylim([-5 7])
box on
[rho,Pval]=corr(kmax,kcat_max_presto,'rows','complete','type','Spearman');
if Pval<0.05
    text(.65,.1,sprintf('\\rho_S^* = %.2f',rho),...
        'Units','normalized', 'HorizontalAlignment', 'left')
else
        text(.65,.1,sprintf('\\rho_S = %.2f',rho),...
        'Units','normalized', 'HorizontalAlignment', 'left')
end
text(.05,.92,'c','FontWeight','bold','Units','normalized','FontSize',12)
set(gca,'FontName','Arial')

% Davidi vs BRENDA
nexttile
line([-5 7],[-5 7], 'Color', 'k', 'LineWidth', 1)
hold on
scatter(log10(kmax),log10(kcat_max_brenda),'filled')
xlabel('log_{10} k_{cat}^{max} (Davidi2016) [s^{-1}]')
ylabel('log_{10} k_{cat}^{max} (BRENDA) [s^{-1}]')
xlim([-5 7]); ylim([-5 7])
box on
[rho,Pval]=corr(kmax,kcat_max_brenda,'rows','complete','type','Spearman');
if Pval<0.05
    text(.65,.1,sprintf('\\rho_S^* = %.2f',rho),...
        'Units','normalized', 'HorizontalAlignment', 'left')
else
        text(.65,.1,sprintf('\\rho_S = %.2f',rho),...
        'Units','normalized', 'HorizontalAlignment', 'left')
end
text(.05,.92,'d','FontWeight','bold','Units','normalized','FontSize',12)
set(gca,'FontName','Arial')

