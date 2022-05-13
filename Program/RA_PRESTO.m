function [pIds,deltas] = RA_PRESTO(models,expVal,E,epsilon,lambda,theta, n_repeats, set_size)
% Robustness analysis for PRESTO

% reset random number generator
rng('default')
pIds = cell(n_repeats,1);
deltas = cell(n_repeats,1);

for i=1:n_repeats
    rand_idx = randsample(numel(models),set_size,false);
    % run PRESTO with the selected random subset
    [~,~,~,changeTab] = PRESTO(...
        models(rand_idx),...
        expVal(rand_idx),...
        E(:,rand_idx),...
        'epsilon',epsilon,'lambda',lambda,'theta',theta);
    
    [pIds{i},uniq_idx] = unique(erase(changeTab.(1),'prot_'));
    deltas{i} = changeTab{uniq_idx,3};
    
end