function [score,mu_pred, varargout] = scoreKcatCorrection(model,mu,E,enzRxnIdx,signed, optout)
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
%   cell optout:        Cell array of character vector indicating which 
%                       additional solutions should be output can contain 
%                       'predE' or 'FVA'}. FVA - varargout contains
%                       a vector of flux ranges when fixing the growth to
%                       the maximum predicted. If options are given the
%                       order of varargout is {FVA, predE}predE - varargout
%                       contains a vector of predicted enzyme usage coefficients
%                       (according to enzRxnIdx) when growth rate is fixed
%                       to mu and pool is minized. 
%                       Example:
%                       [score, mu_pred, FVA, predE] =
%                       scoreKcatCorrection(model, mu, E, enzRxnIdx, true,
%                       {'FVA', 'predE'}
%                       
% Output:
%   double score:      relative error of predicted growth rate to the
%                       experimental growth rate
%   double mu_pred:     Predicted growth rate
%   double FVA:         (see. INPUT optout) Vector of flux ranges for FVA
%   doulbe E_pred;      (see. INPUT optout) Vector of predicted enzyme 
%                       usage coefficients (according to enzRxnIdx) when 
%                       growth rate is fixed to mu and pool is minized

if nargin < 5
    signed = false;
end
if nargout>2
    if nargin <6
        error("If more than two output arguments are required optout has to be specified")
    elseif length(optout)>2 || ~all(ismember(optout, {'FVA' ,'predE'}))
        error("Argument optout can only contain character vectors 'FVA' or 'predE'")
    end
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
    disp("Modelling error in FBA. Setting solution to NaN")
    mu_pred=NaN;
end
if nargout>2
    varargout={};
    if ismember('FVA', optout)
        if isempty(gcp('nocreate'))
            warning('No parallel pool detected, FVA with GECKO models will likely take considerable runtime (days)')
        end
        %Do flux variability analysis with solutions at least 99% of
        %optimal objective value
        try
            [vmin, vmax]= fluxVariability(model, 97);
            varargout=[varargout, {vmax-vmin}];
        catch ME
            warning(ME.message)
            disp("Modelling error in FVA. Setting solution to NaN")
            varargout=[varargout, {nan(length(model.rxns),1)}];
        end
    end
    if ismember('predE', optout)
        if any(E>0)
            error("Enzyme usage can not be predicted if proteomics data is used as constrain")
        end
        %release enzyme pool contrain
        model.ub(end)=inf;
        %fix growth rate to measured
        model.ub(logical(model.c))=mu;
        model.lb(logical(model.c))=0.99*mu;
        model.c(logical(model.c))=0;
        %change optimization target to minimize enzyme pool
        model.c(end)=-1;
        try
            solution2=optimizeCbModel(model);
            varargout=[varargout, {solution2.v(enzRxnIdx)}];
        catch ME
            warning(ME.message)
            disp("Modelling error in minimize enzyme pool FBA. Setting solution to NaN")
            varargout=[varargout, {nan(length(enzRxnIdx),1)}];
        end
    end
end
% calculate relative error to expected growth rate
score = (mu_pred-mu)/mu;

if ~signed
    score = abs(score);
end

end