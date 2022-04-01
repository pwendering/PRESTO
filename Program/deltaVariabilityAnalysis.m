function [minDelta,maxDelta] = deltaVariabilityAnalysis(varargin)
%% [minDelta,maxDelta] = deltaVariabilityAnalysis(varargin)
% Perform a variability analysis for deltas in kcat correction. The sum of
% deltas and relative errors per condition are fixed with a tolerance of 0.001
% and one percent, respectively. The deltas are then maximized and minimized
% in subsequent correction runs.
%
% Input:
%   struct corrLP:          inear program that has been solved to
%                           correct kcat values (output from correctKcats)
%   struct corrSol:         solution to corrLP (output from correctKcats)
%   char enzMetPfx:         (optional) prefix for enzyme mass balance constraints
%                           (i.e. values in model.metNames); default: 'prot_'
% Output:
%   double minDelta:        minimal values for deltas
%   double maxDelta:        maximal values for deltas
% 
% 22.03.2022 Philipp Wendering, University of Potsdam, philipp.wendering@gmail.com

% parse input arguments
p = parseInput(varargin);

corrLP = p.Results.corrLP;
corrSol = p.Results.corrSol;
enzMetPfx = p.Results.enzMetPfx;

% add sum constraint for deltas
deltaIdx = contains(corrLP.rxns,['delta_' enzMetPfx]);
sumDelta = sum(corrSol.x(deltaIdx));

% get the indices of non-zero deltas at do not allow changes in the remaining deltas
nzDeltaIdx = corrSol.x>0 & deltaIdx;
corrLP.ub(deltaIdx & ~nzDeltaIdx) = 0;
nzDeltaIdx = find(nzDeltaIdx);

deltaSumLHS = zeros(2,size(corrLP.S,2));
deltaSumLHS(1,nzDeltaIdx) = 1;
deltaSumLHS(2,nzDeltaIdx) = -1;
deltaSumRHS = [sumDelta+1e-3; 1e-3-sumDelta];
fprintf('\nSum of deltas: %.2g\n',sumDelta)

deltaSumLHS = deltaSumLHS / 1e6;
deltaSumRHS = deltaSumRHS / 1e6;

corrLP.S = [corrLP.S; deltaSumLHS];
corrLP.b = [corrLP.b; deltaSumRHS];
corrLP.csense = [corrLP.csense; ['L';'L']];
corrLP.mets = [corrLP.mets; {'delta_sum_constraint_max';'delta_sum_constraint_min'}];
corrLP.metNames = [corrLP.metNames; {'delta_sum_constraint_min';'delta_sum_constraint_min'}];

% add constraints on omegas (relative errors per condition)
omegaIdx = contains(corrLP.rxns,'omega_cond');
omegaVal = corrSol.x(omegaIdx);

corrLP.lb(omegaIdx) = omegaVal - 0.01*omegaVal;
corrLP.ub(omegaIdx) = omegaVal + 0.01*omegaVal;

% initialize output arrays
minDelta = nan(numel(nzDeltaIdx),1);
maxDelta = nan(numel(nzDeltaIdx),1);

for i=1:numel(nzDeltaIdx)
    fprintf('Current Enzyme: %s\n',corrLP.rxnNames{nzDeltaIdx(i)})
    % change objective so current delta
    tmpLP = changeObjective(corrLP,corrLP.rxns(nzDeltaIdx(i)));
    % minimization
    tmpLP.osenseStr = 'min';
    try
        solution = optimizeCbModel(tmpLP);
    catch ME
        disp(ME.message)
        solution.stat = 3;
    end
    if solution.stat == 1
        minDelta(i) = solution.x(tmpLP.c==1);
        disp(minDelta(i))
    end
    % maximization
    tmpLP.osenseStr = 'max';
    try
        solution = optimizeCbModel(tmpLP);
    catch ME
        disp(ME.message)
        solution.stat = 3;
    end
    if solution.stat == 1
        maxDelta(i) = solution.x(tmpLP.c==1);
        disp(maxDelta(i))
    end
end

    function p = parseInput(arguments)
        % set default values
        ENZ_MET_PFX_DEFAULT = 'prot_';
        
        validateLP = @(M) all(isfield(M,{'S','b','c','lb','ub','rxns',...
            'rxnNames','mets','metNames'}));
        validateSol= @(S) all(isfield(S,{'x','f','stat'}));
        
        
        p = inputParser;
        p.FunctionName = 'deltaVariabilityAnalysis';
        
        addRequired(p,'corrLP',validateLP)
        addRequired(p,'corrSol',validateSol)
        addParameter(p,'enzMetPfx',ENZ_MET_PFX_DEFAULT,@ischar)
        
        parse(p,arguments{:})
    end
end
