
[~,~,~,changeTab] = PRESTO(models,expVal,E,...
    'epsilon',epsilon,'lambda',lambda,'theta',theta);

[pids_all,uniq_idx] = unique(erase(changeTab.(1),'prot_'));
delta_all = changeTab{uniq_idx,3};

n_repeats = 10;
set_sizes = [3 5 10 15];

JI_mats = cell(numel(set_sizes),1);
corr_P_mats = cell(numel(set_sizes),1);
corr_S_mats = cell(numel(set_sizes),1);
plot_mats = cell(numel(set_sizes),1);
ji_presto_full = zeros(numel(set_sizes),n_repeats);
corr_P_presto_full = zeros(numel(set_sizes),n_repeats);
n_shared_presto_full = zeros(numel(set_sizes),n_repeats);
for s = 1:numel(set_sizes)
    [pIds,deltas] = RA_PRESTO(models,expVal,E,epsilon,lambda,theta, n_repeats, set_sizes(s));
    JI = zeros(n_repeats);
    correlation_P = zeros(n_repeats);
    correlation_S = zeros(n_repeats);
    plot_mat = ones(n_repeats);
    for i=1:n_repeats-1
        for j=i+1:n_repeats
            
            [prot_intersect,idx_i,idx_j] = intersect(pIds{i},pIds{j});
            JI(i,j) = numel(prot_intersect) / ...
                numel(union(pIds{i},pIds{j}));
            correlation_P(i,j) = corr(deltas{i}(idx_i),deltas{j}(idx_j));
            correlation_S(i,j) = corr(deltas{i}(idx_i),deltas{j}(idx_j),...
                'type','Spearman');
            plot_mat(i,j) = JI(i,j);
            plot_mat(j,i) = correlation_P(i,j);
            
        end
    end
    JI_mats{s} = JI;
    corr_P_mats{s} = correlation_P;
    corr_S_mats{s} = correlation_S;
    plot_mats{s} = plot_mat;
  
    for i=1:n_repeats
       [prot_intersect,idx_i,idx_full] = intersect(pIds{i},pids_all);
       ji_presto_full(s,i) = numel(prot_intersect) / ...
           numel(union(pIds{i},pids_all));
       corr_P_presto_full(s,i) = corr(deltas{i}(idx_i),delta_all(idx_full));
       n_shared_presto_full(s,i) = numel(prot_intersect);
    end

end

figure
tiledlayout('flow')
letters = {'a','b','c','d'};
for i=1:numel(set_sizes)
    nexttile
    colormap('summer')
    imagesc(plot_mats{i})
    text(-0.1,1,letters{i},...
        'units','normalized',...
        'fontsize',14,...
        'fontweight','bold')
    title(['Size = ' num2str(set_sizes(i))])
    xlabel('\rho_S','fontsize',14)
    ylabel('JI','fontsize',14)
end
colorbar

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
exportgraphics(gcf,'../Results/robustness_presto.tiff','Resolution',300)

save('robustness_ws')