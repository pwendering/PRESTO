function [pIds,deltas] = RA_PRESTO(models,expVal,E,epsilon,lambda,theta, n_repeats, set_size)
%% [pIds,deltas] = RA_PRESTO(models,expVal,E,epsilon,lambda,theta, n_repeats, set_size)
% Robustness analysis for PRESTO using randomly-selected subsets of
% experimental conditions
% Input:
%   cell models:            enzyme-constraint metabolic model(s) generated
%                           using GECKO and with adjusted biomass
%                           composition dependent on total protein content
%   double expVal:          scalar or vector containing experimentally
%                           measured growth rate(s)
%   double E:               column vector or #proteins x #conditions matrix
%                           that contains absolute proteomics data for the
%                           proteins in the model (NaN or zero if unmeasured)
%   double epsilon:         allowed fold change of a k_cat value
%   double lambda:          weight for the deviation of
%                           predicted from experimental growth
%   double theta:           allowed deviation of the predicted growth rate
%                           from the given experimental value
%   double n_repeats:       number of iterations for testing of random
%                           subsets
%   double set_size:        size of randomly-selected subsets
% 
% Output:
%   cellstr pIds:           unique protein IDs associated to corrected kcat
%                           values
%   double deltas:          kcat corrections assocated to unique protein
%                           IDs

% reset random number generator
rng('default')

% initialize output variables
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