function [score,mu_pred] = scoreKcatCorrection(model,mu,E,enzRxnIdx,signed)
%% score = scoreKcatCorrection(model,mu,E)
% Input:
%   struct model:       GECKO metabolic model (kcats corrected using
%                       correctKcats)
%   double mu:          experimental growth rate [h^-1]
%   double E:           array of enzyme abundances per enzyme in the model
%                       (zero or NaN of unmeasured)
%   logical signed:     (optional) whether to return a signed relative
%                       error; otherwise, the absolute value of the
%                       difference between experimental and predicted value
%                       will be used; default: false
% Output:
%   double scores:      relative errors of predicted growth rate to the
%                       experimental growth rate

if nargin < 5
    signed = false;
end

% set upper bounds of enzyme usage reactions, leave the remaining reactions
% at their default upper bound
model.ub(enzRxnIdx(E>0)) = E(E>0);
% run FBA
try
    solution = optimizeCbModel(model);
    mu_pred = solution.f;
catch ME
    warning(ME.message)
    disp("Setting solution toNaN")
    mu_pred=NaN;
end
% calculate relative error to expected growth rate
score = (mu_pred-mu)/mu;

if ~signed
    score = abs(score);
end

end