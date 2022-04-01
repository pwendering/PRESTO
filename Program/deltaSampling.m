function deltaSamplingMat = deltaSampling(varargin)
%% deltaSamplingMat = deltaSampling(varargin)
% Performs random sampling of delta values that satisfy the constraints
% defined in kcat correction and the relative error per condition. Sampling
% is done by minimizing the absolute distance of delta values to a uniformly 
% sampled random vector of deltas within the bound obtained from delta
% variability analysis.
% 
% Input:
%   struct corrLP:                  linear program that has been solved to
%                                   correct kcat values (output from
%                                   correctKcats)
%   struct corrSol:                 solution to corrLP (output from
%                                   correctKcats)
%   double minDelta:                array of minimum delta values
%   double maxDelta:                array of maximum delta values
%   scalar int nSamples:            (optional) number of samples that should be
%                                   returned; default: 100
%   char enzMetPfx:                 (optional) prefix for enzyme mass balance
%                                   constraints (i.e. values in model.metNames);
%                                   default: 'prot_'
% 
% Output:
%   double deltaSamplingMat:        matrix that contains the sampled delta
%                                   values for those deltas that have been
%                                   non-zero in corrSol across nSamples
%                                   sampling iterations
% 22.03.2022 Philipp Wendering, University of Potsdam, philipp.wendering@gmail.com

% parse input arguments
p = parseInput(varargin);
corrLP = p.Results.corrLP;
corrSol = p.Results.corrSol;
minDelta = p.Results.minDelta;
maxDelta = p.Results.maxDelta;
nSamples = p.Results.nSamples;
enzMetPfx = p.Results.enzMetPfx;

% get indices of enzyme mass-balance constraints
enzMetIdx = find(contains(corrLP.mets,enzMetPfx));

% set seed for random numbers
rng('default');

% get indices and numbers of all deltas and those that were non-zero in corrSol
deltaIdx = contains(corrLP.rxns,['delta_' enzMetPfx]);
numAllDelta = sum(deltaIdx);
nzDeltaIdx = corrSol.x > 0 & deltaIdx;
numNzDelta = sum(nzDeltaIdx);

% check dimension of delta ranges
if numel(minDelta) ~= sum(numNzDelta)
    error('Column dimension of minDelta does not match number of non-zero deltas in solution')
elseif numel(maxDelta) ~= sum(numNzDelta)
    error('Column dimension of maxDelta does not match number of non-zero deltas in solution')
end

% add constraints on the sum of deltas and on omegas
corrLP.ub(deltaIdx & ~nzDeltaIdx) = 0;

% sum of deltas from correction
sumDelta = sum(corrSol.x(deltaIdx));

% delta sum constraint matrix
deltaSumLHS = zeros(2,size(corrLP.S,2));
deltaSumLHS(1,nzDeltaIdx) = 1;
deltaSumLHS(2,nzDeltaIdx) = -1;
% delta sum right-hand side
deltaSumRHS = [sumDelta+1e-3; 1e-3-sumDelta];
deltaSumLHS = deltaSumLHS / 1e6;
deltaSumRHS = deltaSumRHS / 1e6;
% update LP fields
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

% add additional variables tao- and tao+ for minimizing the absolute distance to
% a random vector
taoPlusMat = speye(numAllDelta);
taoPlusMat(~nzDeltaIdx(deltaIdx),:) = [];
taoPlusMat(:,~nzDeltaIdx(deltaIdx)) = [];
deltaTaoMat = -speye(numAllDelta);
deltaTaoMat(~nzDeltaIdx(deltaIdx),:) = [];
deltaTaoMat(:,~nzDeltaIdx(deltaIdx)) = [];
taoMatrix = [sparse(numNzDelta, size(corrLP.S,2)) taoPlusMat -taoPlusMat];
taoMatrix(:,nzDeltaIdx) = deltaTaoMat;
clear deltaTaoMat taoPlusMat

% update LP fields
corrLP.c = [sparse(size(corrLP.S,2),1); ones(2*numNzDelta,1)];
corrLP.S = [corrLP.S sparse(size(corrLP.S,1),2*numNzDelta); taoMatrix]; clear taoMatrix
corrLP.b = [corrLP.b; zeros(numNzDelta,1)];
corrLP.lb = [corrLP.lb; zeros(2*numNzDelta,1)];
corrLP.ub = [corrLP.ub; inf(2*numNzDelta,1)];
protIds = strrep(corrLP.mets(enzMetIdx(nzDeltaIdx(deltaIdx))),'_cond_1','');
corrLP.rxns = [corrLP.rxns;...
    cellfun(@(x)sprintf('tao_plus_%s',x),protIds,'un',0);...
    cellfun(@(x)sprintf('tao_minus_%s',x),protIds,'un',0)];
corrLP.rxnNames = [corrLP.rxnNames;...
    cellfun(@(x)sprintf('tao_plus_%s',x),protIds,'un',0);...
    cellfun(@(x)sprintf('tao_minus_%s',x),protIds,'un',0)];
corrLP.mets = [corrLP.mets; cellfun(@(x)sprintf('tao_const_%s',x),protIds,'un',0)];
corrLP.metNames = [corrLP.metNames; cellfun(@(x)sprintf('tao_const_%s',x),protIds,'un',0)];
corrLP.osenseStr = 'min';

% toa constraint indices
taoRHSIdx = contains(corrLP.mets,'tao_const_');

% initialize matrix for delta sampling results
deltaSamplingMat = nan(numNzDelta,nSamples);

for i=1:nSamples
    tmpLP = corrLP;
    
    % random vector of deltas, uniformly sampled between minima and maxima
    randDeltaVector = (maxDelta-minDelta) .* rand(numNzDelta,1) + minDelta;
    tmpLP.b(taoRHSIdx) = -randDeltaVector;
    
    % run sampling by minimization
    sol.stat = 0;
    try
        sol = optimizeCbModel(tmpLP);
    catch ME
        disp(ME.message)
    end
    
    if sol.stat == 1
        deltaSamplingMat(:,i) = sol.x(nzDeltaIdx);
    end
    
end

    function p = parseInput(arguments)
        % set default values
        ENZ_MET_PFX_DEFAULT = 'prot_';
        SAMPLES_DEFAULT = 100;
        
        % validation functions
        validateLP = @(M) all(isfield(M,{'S','b','c','lb','ub','rxns',...
            'rxnNames','mets','metNames'}));
        validateSol= @(S) all(isfield(S,{'x','f','stat'}));
        validScalarInt = @(x)~ischar(x)&isscalar(x)&floor(x)==x;

        p = inputParser;
        p.FunctionName = 'deltaSampling';
        
        addRequired(p,'corrLP',validateLP)
        addRequired(p,'corrSol',validateSol)
        addRequired(p,'minDelta',@isnumeric)
        addRequired(p,'maxDelta',@isnumeric)
        addParameter(p,'nSamples',SAMPLES_DEFAULT,validScalarInt)
        addParameter(p,'enzMetPfx',ENZ_MET_PFX_DEFAULT,@ischar)
        
        parse(p,arguments{:})
    end


end
