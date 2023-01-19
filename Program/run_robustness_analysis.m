% Run this script after running the correct_kcats_yeast.m script, so all
% required variables are created:
% models:       condition-specific enzyme-constraint models
% expVal:       experimentally measured specific growth rates
% epsilon:      allowed fold change for kcat corrections
% lambda:       parameter tuning the trade-off between introduced 
%               corrections and relative errors between experimentally
%               measured and predicted growth rates 
% theta:        maximum allowed relative error between measured and
%               predicted growth rates

% run PRESTO to obtain the set of corrections with all conditions used
[~,~,~,changeTab] = PRESTO(models,expVal,E,...
    'epsilon',epsilon,'lambda',lambda,'theta',theta);

[pids_all,uniq_idx] = unique(erase(changeTab.(1),'prot_'));
delta_all = changeTab{uniq_idx,3};

% set the number of iterations and set sizes for robustness analysis
n_repeats = 50;
set_sizes = [3 5 10 15];

JI_mats = cell(numel(set_sizes),1);
corr_P_mats = cell(numel(set_sizes),1);
corr_S_mats = cell(numel(set_sizes),1);

plot_mats = cell(numel(set_sizes),1);

ji_presto_full = zeros(numel(set_sizes),n_repeats);
corr_P_presto_full = zeros(numel(set_sizes),n_repeats);
n_shared_presto_full = zeros(numel(set_sizes),n_repeats);

for s = 1:numel(set_sizes)
    
    % run robustness analysis
    [pIds,deltas] = RA_PRESTO(models,expVal,E,epsilon,lambda,theta, n_repeats, set_sizes(s));
    
    JI = zeros(n_repeats);
    correlation_P = zeros(n_repeats);
    correlation_S = zeros(n_repeats);
    plot_mat = ones(n_repeats);
    
    for i=1:n_repeats-1
        
        for j=i+1:n_repeats
            
            % calculate Jaccard index and correlation
            [prot_intersect,idx_i,idx_j] = intersect(pIds{i},pIds{j});
            JI(i,j) = numel(prot_intersect) / ...
                numel(union(pIds{i},pIds{j}));
            correlation_P(i,j) = corr(deltas{i}(idx_i),deltas{j}(idx_j));
            correlation_S(i,j) = corr(deltas{i}(idx_i),deltas{j}(idx_j),...
                'type','Spearman');
            
            % combine JI and correlation into a composite matrix, which
            % contained JI in the upper right triangle and correlation in
            % the lower left triangle
            plot_mat(i,j) = JI(i,j);
            plot_mat(j,i) = correlation_P(i,j);
            
        end
    end
    
    JI_mats{s} = JI;
    corr_P_mats{s} = correlation_P;
    corr_S_mats{s} = correlation_S;
    plot_mats{s} = plot_mat;
    
    % calculate JI and correlation to corrections using the full set
    for i=1:n_repeats
       [prot_intersect,idx_i,idx_full] = intersect(pIds{i},pids_all);
       ji_presto_full(s,i) = numel(prot_intersect) / ...
           numel(union(pIds{i},pids_all));
       corr_P_presto_full(s,i) = corr(deltas{i}(idx_i),delta_all(idx_full));
       n_shared_presto_full(s,i) = numel(prot_intersect);
    end

end

% create figure that compares the iterations for each set size
figure
tiledlayout('flow')
letters = {'a','b','c','d'};
for i=1:numel(set_sizes)
    nexttile
    colormap('bone')
    imagesc(plot_mats{i})
    text(-0.1,1,letters{i},...
        'units','normalized',...
        'fontsize',14,...
        'fontweight','bold')
    title(['Size = ' num2str(set_sizes(i))])
    xlabel('\rho_P','fontsize',14)
    ylabel('JI','fontsize',14)
end
colorbar
set(gcf,'outerposition',[353.6667  115.6667  722.0000  583.3333])
exportgraphics(gcf,'Results/similarity_random_subsets_presto_yeast_n50.tiff','Resolution',300)
% exportgraphics(gcf,'Results/similarity_random_subsets_presto_ecoli_n50.tiff','Resolution',300)
% create figure that compares the results from each iteration and set size
% with the solution obtained using all experimental conditions
figure
colormap('bone')
tiledlayout('flow')

nexttile
imagesc(ji_presto_full)
yticks(1:numel(set_sizes));
yticklabels(set_sizes)
ylabel('Set size','fontsize',14)
xlabel('Iteration','fontsize',14)
text(-0.1,1,letters{1},...
        'units','normalized',...
        'fontsize',14,...
        'fontweight','bold')
colorbar
hcb = colorbar;
hcb.Title.String = "JI";
nexttile
imagesc(corr_P_presto_full)
yticks(1:numel(set_sizes));
yticklabels(set_sizes)
ylabel('Set size','fontsize',14)
xlabel('Iteration','fontsize',14)
text(-0.1,1,letters{2},...
        'units','normalized',...
        'fontsize',14,...
        'fontweight','bold')
hcb = colorbar;
hcb.Title.String = "\rho_P";

set(gcf,'Outerposition',[353.6667  357.0000  817.3333  342.0000])
exportgraphics(gcf,'Results/robustness_presto_yeast_n50.tiff','Resolution',300)
% exportgraphics(gcf,'Results/robustness_presto_ecoli_n50.tiff','Resolution',300)

%% boxplot
tiledlayout('flow')
nexttile
boxchart(...
    ji_presto_full',...
    'LineWidth',1.3,...
    'JitterOutliers','on'...
);
ylim([0.9*min(min(ji_presto_full)) 1.05*max(max(ji_presto_full))])
xticklabels(set_sizes)
y_ticks = get(gca,'YTick');
set(gca,'YTick',y_ticks(y_ticks<=1))
xlabel('Set size')
ylabel('Jaccard Index')
text(-0.15,1,letters{1},...
        'units','normalized',...
        'fontsize',14,...
        'fontweight','bold')
set(gca,'FontSize',14,'FontName','Arial','Box','on','LineWidth',1.3)

    
nexttile
boxchart(...
    corr_P_presto_full',...
    'LineWidth',1.3,...
    'JitterOutliers','on'...
);
ylim([0.9*min(min(corr_P_presto_full)) 1.05*max(max(corr_P_presto_full))])
xticklabels(set_sizes)
y_ticks = get(gca,'YTick');
set(gca,'YTick',y_ticks(y_ticks<=1))
xlabel('Set size')
ylabel('\rho_P')
text(-0.15,1,letters{2},...
        'units','normalized',...
        'fontsize',14,...
        'fontweight','bold')
    
set(gca,'FontSize',14,'FontName','Arial','Box','on','LineWidth',1.3)
set(gcf,'OuterPosition',1000*[-1.4103    0.1557    0.7320    0.3200])

exportgraphics(gcf,'Results/robustness_presto_yeast_boxplot.tiff','Resolution',300)
% exportgraphics(gcf,'Results/robustness_presto_ecoli_boxplot.tiff','Resolution',300)

