%% Compare corrected kcat values with estimated values from Chen et al. 2021
clear;clc

%% Read data
% read S4 Dataset from Chen et al. 2021
chen2021_s4_dataset = readtable(...
    fullfile('Data','Chen','pnas.2108391118.sd04.xlsx'),...
    'NumHeaderLines',1);
kmax = chen2021_s4_dataset.kmax__s_;
rxns = chen2021_s4_dataset.ReactionIDInYeast8;

% modify Chen2021 reaction IDs to enable comparison to GECKO model
rxns = erase(rxns,'_fwd');
rxns = strrep(rxns,'_rvs','_REV');

% read yeast 8 model
ecYeastGEM = readGKOmodel(fullfile('Data','rawecYeast.mat'));

% load corrected models (PRESTO)
load(fullfile('Results','corrected_models_s_cerevisiae_1_e_minus_7'))

enzMetIdx = find(startsWith(ecYeastGEM.mets,'prot_'));
kcat_max_brenda = nan(numel(rxns),1);
kcat_max_presto = nan(numel(rxns),1);

ecYeastGEM.rxns = regexprep(ecYeastGEM.rxns,'No\d+','');
ishomomeric = false(numel(rxns),1);
for i=1:numel(rxns)
    % find reaction IDs in the model that contain the current ID
    matchIdx = ~cellfun(@isempty,regexp(ecYeastGEM.rxns,['^' rxns{i} '[^_]*$']));
    ishomomeric(i) = all(~cellfun(@isempty,ecYeastGEM.grRules(matchIdx)) & ~contains(ecYeastGEM.grRules(matchIdx),{'and','or'}));
    if sum(matchIdx) > 0
    % determine max kcat in orgiginal model (from BRENDA)
    kcat_max_brenda(i) = max(max(-1./ecYeastGEM.S(enzMetIdx,matchIdx)/3600));
    % determine max kcat in corrected model (PRESTO)
    kcat_max_presto(i) = max(max(-1./corrModels{1}.S(enzMetIdx,matchIdx)/3600));
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
fig_chen2021_kcat_comp = figure;
tiledlayout('flow')
% Chen vs PRESTO
nexttile
line([-5 7],[-5 7], 'Color', 'k', 'LineWidth', 1)
box on
hold on
scatter(log10(kmax),log10(kcat_max_presto),'filled')
xlabel('log_{10} k_{cat}^{max} (Chen2021) [s^{-1}]')
ylabel('log_{10} k_{cat}^{max} (PRESTO) [s^{-1}]')
xlim([-5 7]); ylim([-5 7])
box on
[rho,Pval] = corr(kmax,kcat_max_presto,'rows','complete','type','Spearman');
if Pval<0.05
    text(.65,.1,sprintf('\\rho_S^* = %.2f',rho),...
        'Units','normalized', 'HorizontalAlignment', 'left')
else
        text(.65,.1,sprintf('\\rho_S = %.2f',rho),...
        'Units','normalized', 'HorizontalAlignment', 'left')
end
text(.05,.92,'a','FontWeight','bold','Units','normalized','FontSize',12)
set(gca,'FontName','Arial')

% Chen vs BRENDA
nexttile
line([-5 7],[-5 7], 'Color', 'k', 'LineWidth', 1)
hold on
scatter(log10(kmax),log10(kcat_max_brenda),'filled')
xlabel('log_{10} k_{cat}^{max} (Chen2021) [s^{-1}]')
ylabel('log_{10} k_{cat}^{max} (BRENDA) [s^{-1}]')
xlim([-5 7]); ylim([-5 7])
box on
[rho,Pval] = corr(kmax,kcat_max_brenda,'rows','complete','type','Spearman');
if Pval<0.05
    text(.65,.1,sprintf('\\rho_S^* = %.2f',rho),...
        'Units','normalized', 'HorizontalAlignment', 'left')
else
        text(.65,.1,sprintf('\\rho_S = %.2f',rho),...
        'Units','normalized', 'HorizontalAlignment', 'left')
end
text(.05,.92,'b','FontWeight','bold','Units','normalized','FontSize',12)
set(gca,'FontName','Arial')
