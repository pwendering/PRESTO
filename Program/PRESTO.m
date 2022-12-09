function [solution,models,relError,changeTab,LP] = PRESTO(varargin)
%% [solution,models,relError,changeTab,LP] = PRESTO(varargin)
% Function to correct turnover values in an enzyme-constraint metabolic
% model generated using the GECKO (v2.0) pipeline (DOI:10.15252/msb.20167411).
% The minimum kcat of every measured enzyme across all reactions is used
% to impose updated upper bounds on the reactions, introducing a variable
% delta, which is represents the increase in the minimum kcat that is required
% to achieve the given experimental growth rate:
%                       v <= (k_cat + delta) * E .
% Input:
%   cell models:            enzyme-constraint metabolic model(s) generated
%                           using GECKO and with adjusted biomass
%                           composition dependent on total protein content
%   double expVal:          scalar or vector containing experimentally
%                           measured growth rate(s)
%   double E:               column vector or #proteins x #conditions matrix
%                           that contains absolute proteomics data for the
%                           proteins in the model (NaN or zero if unmeasured)
%   char enzMetPfx:         (optional) prefix for enzyme mass balance constraints
%                           (i.e. values in model.metNames); default: 'prot_'
%   char enzRxnPfx:         (optional) prefix for enzyme usage pseudoreactions
%                           (i.e. values in model.rxns); default: 'prot_'
%   double epsilon:         (optional) allowed fold change of a k_cat value
%                           default: 1e6
%   double theta:           (optional) allowed deviation of the predicted growth rate
%                           from the given experimental value; default: 0.5
%   double lambda:          (optional) weight for the deviation of
%                           predicted from experimental growth (default:
%                           1e-7)
%   cellstr enzBlackList:   (optional) names of enzyme metabolites that
%                           should be excluded from the correction
%                           default: {''}
%   double K:               (optional) maximum kcat that is allowed to be
%                           reached using delta values [s^-1];
%                           default: 57500000 s^-1
%   logical includeUM:      (optional) whether to include unmeasured
%                           enzymes in the correction procedure; if
%                           specified, the model must contain a field
%                           'protMW'; default: false
%   double f_n:             (optional) mass fraction of all unmeasured proteins
%                           included in the model (see GECKO documentation);
%                           This is only relevant when includeUM=true;
%                           default: 0.5 for all conditions
%   double sigma:           (optional) average in vivo enzyme saturation
%                           (can be fitted using GECKO sigmaFitter);
%                           default: 0.5
%   double pCorrFactor:     (optional) correction factor that has been used
%                           to increase the enzyme abundances (to avoid
%                           very low f factors when enzyme pool is used)
%                           default: all equal to one
%   logical negCorrFlag:    (optional) add second step that finds negative
%                           corrections for kcat values;
%                           default: false
% Output:
%   struct solution:        solution to the linear program as returned by
%                           solveCobraLP
%   struct model:           metabolic models with updated kcat values
%   double relError:        relative error of predicted growth rate with
%                           respect to the experimental growth rate
%   table changeTab:        table containing information on introduced
%                           changes (enzyme name, original kcat, updated
%                           kcat, and fold change)
%   struct LP:              linear optimization program that was used to
%                           obtain corrected kcat values
% 
% 22.03.2022 Philipp Wendering, University of Potsdam, philipp.wendering@gmail.com
% 09.12.2022 Philipp Wendering
% * added option for negative adjustments (i.e. decrease of kcat values) in
%   a second optimization step with fixed positive corrections and a lower
%   bound on the sum of relative errors; the minimum possible kcat after
%   the second correction step is the minimum kcat seen in the model before
%   applying any correction
% * fixed error with dealing with unmeasured proteins by changing into
%   organism-specific GECKO directory before calling 'sumProtein'
% * changed FOLD_INCREASE to FOLD_CHANGE in reporting table and added
%   DELTA_LB columns for delta lower bounds
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parse input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
p = parseInput(varargin);

models = p.Results.models;
expVal = p.Results.expVal;
E = p.Results.E;
enzMetPfx = p.Results.enzMetPfx;
enzRxnPfx = p.Results.enzRxnPfx;
epsilon = p.Results.epsilon;
theta = p.Results.theta;
lambda = p.Results.lambda;
enzBlackList = p.Results.enzBlackList;
K = 3600*p.Results.K;
includeUM = p.Results.includeUM;
f_n = ones(size(E,2),1).*p.Results.f_n;
sigma = p.Results.sigma;
pCorrFactor = p.Results.pCorrFactor;
negCorrFlag = p.Results.negCorrFlag;

if includeUM && ~isfield(models{1},'protMW')
    error('The provided models lack a field ''protMW'' for the inclusion of unmeasured proteins')
elseif includeUM
    MW = models{1}.protMW/1000; % [g/mmol]
end

% number of conditions
NCOND = size(E,2);
if numel(expVal) ~= NCOND
    error('Number of provided growth rates does not match the number of conditions')
elseif numel(models) ~= NCOND
    error('Number of provided models does not match the number of conditions')
end

if isstruct(models)
    models = {models};
end

% if pCorrFactor is a scalar, transform to vector with as many entries as
% models
if isscalar(pCorrFactor)
    pCorrFactor = pCorrFactor .* ones(numel(models),1);
end

% find commmon indices using first model as representative
model = models{1};

% find enzyme indices in model metabolites
[enzMetIdx,enzRxnIdx] = findEnzIdx(model,enzMetPfx,enzRxnPfx);
% get row and column dimensions of stoichiometric matrix
[RDIM,CDIM] = size(model.S);
% determine the number of deltas from the number of enzyme usage reactions
NDELTA = numel(enzRxnIdx);
nBiochemRxns = CDIM-NDELTA;
% find biomass reaction indices in models
bioIdx = findBioIdx(models,nBiochemRxns);
% get blacklist indices (enzymes, which should be excluded from correction)
blackListIdx = ismember(model.mets(enzMetIdx),enzBlackList);

if includeUM
    idxExclude = false(size(E,1),1);
    idxUM = ~any(E>0,2);
    
    % divide E by protein abundance correction factor
    E = (E'./pCorrFactor)';
    
    % set initial abundances for each enzyme by uniformly distributing the
    % remaining protein mass
    currDir = pwd;
    geckoPath = regexp(which('sumProtein'), '.*GECKO_[^\\]+', 'match');
    cd(fullfile(char(geckoPath), 'geckomat', 'limit_proteins'))
    pTot = cellfun(@(M)sumProtein(M),models)';
    pRemain = pTot-sum(bsxfun(@times,E,MW))';
    f = zeros(size(E,2),1);
    for i=1:size(E,2)
        tmpIdx = ~E(:,i)>0;
        E(tmpIdx,i) = min(E(E(:,i)>0,i));
        pRemain(i) = pRemain(i)-E(tmpIdx,i)'*MW(tmpIdx);
        % consider f and sigma
        f_m = (pTot(i)-pRemain(i)) / pTot(i);
        f(i) = f_n(i)./(1-f_m);
        pRemain(i) = sigma*f(i)*pRemain(i);
    end
    
    % again multiply E by protein abundance correction factor
    E = (pCorrFactor.*E')';
    cd(currDir)
else
    % find enzymes for which there are no abundances available in at least one
    % condition
    idxExclude = ~all(E>0,2);
end

% store original and minimum kcat values per enzyme (identical across all models)
[~,minKcat] = getKcatFromEcModel(model, enzMetPfx, enzRxnPfx);

% initialize matrix and rhs for the final linear optimization program
nz = NCOND * sum(sum(model.S~=0)) + sum(any(E,2)) - sum(blackListIdx);
lpMatrix = spalloc(NCOND*RDIM, NCOND*nBiochemRxns+NDELTA, nz);
% right hand side
lpRHS = zeros(size(lpMatrix,1),1);
% constraint sense (COBRA style)
tmpSense = model.csense;
tmpSense(enzMetIdx) = 'L';
lpSense = repmat(tmpSense,NCOND,1);
clear model tmpSense
% LP objective (minimize deltas)
deltaColStartIdx = NCOND*nBiochemRxns+1;
lpObj = zeros(size(lpMatrix,2),1);

% weigh deltas by lambda, divided by the number of kcats that can be corrected
lpObj(deltaColStartIdx:end) =  lambda / (sum(any(E,2)) - sum(blackListIdx));

for i=1:numel(models)
    % find appropriate row and column ranges for current stoichiometric
    % matrix
    tmpModel = models{i};
    
    % set entries of enzyme usage reactions to zero
    tmpModel.S(enzMetIdx,enzRxnIdx) = 0;
    
    % loop over enzyme metabolites and change entries for kinetic constants
    % and enzyme pseudoreactions
    for j=1:numel(enzMetIdx)
        % find 1/kcat entries
        enzKcatIdx = tmpModel.S(enzMetIdx(j),:)<0;
        if sum(enzKcatIdx) > 0
            if ismember(enzMetIdx(j),enzMetIdx(idxExclude))
                tmpModel.S(enzMetIdx(j),:) = 0;
            else
                % if entries were found, and a measurement for the enzyme
                % is available, change the respective kcat entry to one
                tmpModel.S(enzMetIdx(j),enzKcatIdx) = 1;
                % add the delta index
                tmpModel.S(enzMetIdx(j),nBiochemRxns+j) = -E(j,i);
            end
        end
    end
    clear enzKcatIdx
    
    % separate stoichiometrix matrix into a part of biochemical reactions
    % and deltas
    tmpSMat = tmpModel.S(:,1:nBiochemRxns);
    tmpDeltaMat = tmpModel.S(:,enzRxnIdx);
    
    % row and column indices for the current block
    rowIdx = (i-1)*RDIM+1:i*RDIM;
    colIdx = (i-1)*nBiochemRxns+1:i*nBiochemRxns;
    
    % update LP components
    lpMatrix(rowIdx,colIdx) = tmpSMat;
    lpMatrix(rowIdx,deltaColStartIdx:end) = tmpDeltaMat;
    tmpRHS = tmpModel.b;
    tmpRHS(enzMetIdx) = minKcat.*E(:,i);
    lpRHS(rowIdx) = tmpRHS;
end
lpMatrix = sparse(lpMatrix);
clear rowIdx colIdx enzKcatIdx tmpSMat tmpDeltaMat tmpModel tmpRHS

% lower bounds
lpLB = cellfun(@(M)M.lb(1:nBiochemRxns)',models,'UniformOutput',false);
deltaLB = zeros(NDELTA,1);
lpLB = [[lpLB{:}]'; deltaLB];

% upper bounds
lpUB = cellfun(@(M)M.ub(1:nBiochemRxns)',models,'UniformOutput',false);
% minimum of allowed fold-change and cap value
deltaUB = min((epsilon-1)*minKcat,K-minKcat);
% set deltas for blacklist enzymes and unmeasured enzymes to zero
deltaUB(blackListIdx) = 0;
deltaUB(idxExclude) = 0;

lpUB = [[lpUB{:}]'; deltaUB];

% construct omega constraints
omega = sparse(2*NCOND,size(lpMatrix,2)+NCOND);
for i=1:2:2*NCOND
    condIdx = (i+1)/2;
    omega(i,[bioIdx(condIdx) size(lpMatrix,2)+condIdx]) = -[1 expVal(condIdx)];
    omega(i+1,[bioIdx(condIdx) size(lpMatrix,2)+condIdx]) = [1 -expVal(condIdx)];
end
clear condIdx

if includeUM
    % create matrix for uncertainties for unmeasured proteins (i.e., similar to
    % the 'protein Pool' in the GECKO model
    
    % create block for rhos
    rhoMatSingle = sparse(RDIM,NDELTA);
    rhoMatSingle(enzMetIdx,:) = -diag(idxUM.*minKcat);
    % initialize rho matrix and additional summation constraint to respect
    % the available protein content
    rhoMatFull = sparse(RDIM*NCOND,NDELTA*NCOND+NCOND);
    rhoSumConst = sparse(NCOND,NCOND*NDELTA);
    for i=1:NCOND
        rowIdx = (i-1)*RDIM+1:i*RDIM;
        colIdx = (i-1)*NDELTA+1:i*NDELTA;
        rhoMatFull(rowIdx,colIdx) = rhoMatSingle;
        rhoSumConst(i,colIdx) = MW.*idxUM;
    end
    rhoSumConst = [sparse(NCOND,NCOND*nBiochemRxns+NDELTA+NCOND) rhoSumConst -speye(NCOND)];
    rhoLB = zeros(NCOND*NDELTA+NCOND,1);
    rhoUB = [ones(NCOND*NDELTA,1); 0.1*pRemain];
    rhoUB(repmat(~idxUM,1,NCOND)) = 0;
    rhoObj = [zeros(NCOND*NDELTA,1); ones(NCOND,1)];
    rhoSumRHS = pRemain;
    rhoSumSense = repmat('L',NCOND,1);
else
    rhoMatFull = sparse(1,0);
    rhoLB = zeros(0,1);
    rhoUB = zeros(0,1);
    rhoSumConst = sparse(0,1);
    rhoSumRHS = zeros(0,1);
    rhoSumSense = char;
    rhoObj = zeros(0,1);
end

% create linear optimization problem
LP.S = [
    [lpMatrix sparse(size(lpMatrix,1),NCOND) rhoMatFull];
    [omega sparse(2*NCOND, size(rhoMatFull,2))];
    rhoSumConst];
LP.b = [
    lpRHS;
    repelem(columnVector(expVal),2,1).*repmat([-1;1],NCOND,1);
    rhoSumRHS];
LP.lb = [lpLB; zeros(NCOND,1); rhoLB];
LP.ub = [lpUB; theta.*ones(NCOND,1); rhoUB];
LP.csense = [lpSense; repmat('L',2*NCOND,1); rhoSumSense];
LP.c = [lpObj; ones(NCOND,1) / numel(models); rhoObj];

LP.osenseStr = 'min';
clear lpMatrix lpRHS lpLB lpUB lpObj omega rhoMatFull rhoLB rhoUB rhoObj ...
    rhoSumConst rhoSumRHS rhoSumSense

% update reaction and metabolite names and IDs
rxnCondNrStr = strtrim(cellstr(num2str(repelem((1:NCOND)',nBiochemRxns,1))));
LP.rxns = [
    strcat(strcat(repmat(models{1}.rxns(1:nBiochemRxns),NCOND,1),'_cond_'), rxnCondNrStr);...
    strcat('delta_', strtok(models{1}.metNames(enzMetIdx),'['));...
    strcat({'omega_cond_'}, num2str((1:NCOND)'))
    ];

LP.rxnNames = [
    strcat(strcat(repmat(models{1}.rxnNames(1:nBiochemRxns),NCOND,1),{' - condition '}), rxnCondNrStr);...
    strcat('delta_', strtok(models{1}.metNames(enzMetIdx),'['));...
    strcat({'omega_cond_'}, num2str((1:NCOND)'))
    ];
if includeUM
    rxnCondNrStr = strtrim(cellstr(num2str(repelem((1:NCOND)',NDELTA,1))));
    LP.rxns = [LP.rxns;...
        strcat(repmat(strcat('rho_', strtok(models{1}.metNames(enzMetIdx),'[')),NCOND,1),...
        strcat('_cond_', rxnCondNrStr));...
        arrayfun(@(i)sprintf('delta_ptot_%d',i),1:NCOND,'un',0)'];
    LP.rxnNames = [LP.rxnNames; strcat(repmat(strcat('rho_',...
        strtok(models{1}.metNames(enzMetIdx),'[')),NCOND,1),strcat('_cond_', rxnCondNrStr));...
        arrayfun(@(i)sprintf('delta_ptot_%d',i),1:NCOND,'un',0)'];
end

metCondNrStr = strtrim(cellstr(num2str(repelem((1:NCOND)',numel(models{1}.mets),1))));
LP.mets = [strcat(strcat(repmat(models{1}.mets,NCOND,1),'_cond_'), metCondNrStr);...
    strcat(repmat({'omega_1_cond_'; 'omega_2_cond_'},NCOND,1),num2str(repelem((1:NCOND)',2,1)))] ;
LP.metNames = [strcat(strcat(repmat(models{1}.metNames,NCOND,1),{' - condition '}), metCondNrStr);...
    strcat(repmat({'omega_1_cond_'; 'omega_2_cond_'},NCOND,1),num2str(repelem((1:NCOND)',2,1)))];

if includeUM
    LP.mets = [LP.mets; strcat({'rho_sum_const_'}, num2str((1:NCOND)'))];
    LP.metNames = [LP.metNames; strcat({'rho_sum_const_'}, num2str((1:NCOND)'))];
end
clear rxnCondNrStr metCondNrStr

% solve the LP
solution = optimizeCbModel(LP);

if negCorrFlag
    
    LP_min = LP;
    
    % fix non-zero deltas and test feasibility
    deltaIdx = deltaColStartIdx:deltaColStartIdx+NDELTA-1;
    deltaVal = solution.x(deltaIdx);
    nzDeltaIdx = deltaVal~=0;
    feasTol = getCobraSolverParams('LP', 'feasTol');
    LP_min.lb(deltaIdx(nzDeltaIdx)) = deltaVal(nzDeltaIdx) - feasTol;
    LP_min.ub(deltaIdx(nzDeltaIdx)) = deltaVal(nzDeltaIdx) + feasTol;
    
    % set upper bound for sum of relative errors
    omegaIdx = startsWith(LP_min.rxns, 'omega_cond_');
    omegaSumRow = sparse(1, size(LP_min,2));
    omegaSumRow(omegaIdx) = 1;
    omegaSumRHS = sum(solution.x(omegaIdx)) + feasTol;
    omegaSumSense = 'L';
    LP_min.S = [LP_min.S; omegaSumRow];
    LP_min.b = [LP_min.b; omegaSumRHS];
    LP_min.csense = [LP_min.csense; omegaSumSense];
    LP_min.mets = [LP_min.mets; 'sum_omega_constraint'];
    LP_min.metNames = [LP_min.metNames; 'sum_omega_constraint'];
    
    % allow only for negative deltas if not changed in first step
    deltaLB = max((1/epsilon-1)*minKcat,min(minKcat)-minKcat);
    deltaLB(blackListIdx) = 0;
    deltaLB(idxExclude) = 0;
    
    % also exclude deltas for reactions that are inactive across all
    % conditions
    idxAllZeroFlux = false(NDELTA, 1);
    for i = 1:NDELTA
        % find reactions associated to current delta
        protID = erase(LP.rxns(deltaIdx(i)), 'delta_');
        rxnIDs = models{1}.rxns(models{1}.S(findMetIDs(models{1},protID),:)<0);
        rxnIdx = arrayfun(@(i)i+nBiochemRxns*(0:NCOND-1),...
            findRxnIDs(models{1}, rxnIDs), 'UniformOutput', false);
        rxnIdx = [rxnIdx{:}];
        idxAllZeroFlux(i) = ~any(solution.x(rxnIdx));
    end
    deltaLB(idxAllZeroFlux) = 0;
    
    LP_min.lb(deltaIdx(~nzDeltaIdx)) = deltaLB(~nzDeltaIdx);
    LP_min.ub(deltaIdx(~nzDeltaIdx)) = zeros(sum(~nzDeltaIdx),1);
    
    % minimize the sum of allowed deltas
    LP_min.c(deltaIdx) = 0;
    LP_min.c(deltaIdx(~nzDeltaIdx)) = lambda / (sum(any(E,2)) - sum(blackListIdx) - sum(nzDeltaIdx));
    
    solution_min = optimizeCbModel(LP_min);
    
    solution = solution_min;
end

if solution.stat == 1
    % get delta values
    deltaVal = solution.x(deltaColStartIdx:deltaColStartIdx+NDELTA-1);
    % find indices of enzymes with corrected kcats
    enzMetChangedIdx = find(deltaVal~=0);
    % reduce delta values to non-zero values
    deltaVal(deltaVal==0) = [];
    enzMetNames = regexprep(LP.metNames(enzMetIdx(enzMetChangedIdx)),' - condition \d+','');
    changeTab = cell2table(repmat({'',0,0,0,0,0,0,true},numel(enzMetChangedIdx),1),'VariableNames',...
        {'ENZYME_MET_NAME','KCAT_ORIG [s^-1]','KCAT_UPDATED [s^-1]',...
        'FOLD_CHANGE','DELTA','DELTA_LB','DELTA_UB', 'MEASURED'});
    for i=1:numel(enzMetChangedIdx)
        % calculate updated kcat value
        updtKcat = minKcat(enzMetChangedIdx(i))+deltaVal(i);
        % update result table
        changeTab.(1)(i) = enzMetNames(i);
        changeTab.(2)(i) = minKcat(enzMetChangedIdx(i))/3600;
        changeTab.(3)(i) = updtKcat/3600;
        changeTab.(4)(i) = updtKcat/minKcat(enzMetChangedIdx(i));
        changeTab.(5)(i) = deltaVal(i)/3600;
        changeTab.(6)(i) = deltaLB(enzMetChangedIdx(i))/3600;
        changeTab.(7)(i) = deltaUB(enzMetChangedIdx(i))/3600;
        changeTab.(8)(i) = true;
    end
    
    if includeUM
        isMeasured = ~any(rhoMatSingle);
        changeTab.(7) = full(isMeasured(enzMetChangedIdx))';
    end
    % update kcats in models
    for i=1:numel(models)
        tmpModel = models{i};
        for j=1:numel(enzMetChangedIdx)
            rowIdx = enzMetIdx(enzMetChangedIdx(j));
            % find kcat entries
            enzKcatIdx = find(tmpModel.S(rowIdx,:)<0);
            origKcats = -1./tmpModel.S(rowIdx,enzKcatIdx);
            % find indices if minimum
            idxMin = origKcats==minKcat(enzMetChangedIdx(j));
            updtKcat = minKcat(enzMetChangedIdx(j)) + deltaVal(j);
            tmpModel.S(rowIdx,enzKcatIdx(idxMin)) = -1./updtKcat;
        end
        models{i} = tmpModel;
    end
    
    % get relative error
    relError = solution.x(contains(LP.rxns,'omega_cond_'));
    
else
    changeTab = table;
    relError = zeros(sum(contains(LP.rxns,'omega_cond_')),1);
end

    function [origKcat,minKcat,nUniqPerEnzyme] = getKcatFromEcModel(model,enzMetPfx,enzRxnPfx)
        % retrieve kcat values from GECKO stoichiometric matrix
        [enzMetIdx,enzRxnIdx] = findEnzIdx(model,enzMetPfx,enzRxnPfx);
        % find kcats
        minKcat = zeros(numel(enzMetIdx),1);
        nUniqPerEnzyme = zeros(numel(enzMetIdx),1);
        nTotKcat = sum(sum(model.S(enzMetIdx,~enzRxnIdx)<0));
        origKcat = zeros(nTotKcat,1);
        c = 0;
        for e=1:numel(enzMetIdx)
            % find reactions associated with the current enzyme
            enzKcatIdx = find(model.S(enzMetIdx(e),:)<0);
            tmpCoeff = zeros(numel(enzKcatIdx),1);
            
            for k=1:numel(enzKcatIdx)
                c = c + 1;
                % find unique kcats per enzyme
                tmpCoeff(k) = model.S(enzMetIdx(e),enzKcatIdx(k));
                origKcat(c) = -1/tmpCoeff(k);
            end
            nUniqPerEnzyme(e) = numel(unique(tmpCoeff));
            % find maximum kcat for this enzyme
            minKcat(e) = min(-1./tmpCoeff);
            
        end
    end
    function bioIdx = findBioIdx(models,nBiochemRxns)
        % find biomass reaction indices from model objectives and return a
        % double array that contains the positions within the complete linear
        % optimization problem
        bioIdx = nan(numel(models),1);
        for m=1:numel(models)
            bioIdx(m) = (m-1)*nBiochemRxns+find(models{m}.c);
        end
        
    end
    function [enzMetIdx,enzRxnIdx] = findEnzIdx(model,enzMetPfx,enzRxnPfx)
        % Find indices of enzyme metabolites and associated enzyme usage
        % pseudoreaction in the GECKO stoichiometric matrix
        
        enzMetIdx = find(contains(model.metNames,enzMetPfx));
        if isempty(enzMetIdx)
            error('No enzyme metabolites found in model with given prefix')
        end
        
        enzRxnIdx = find(contains(model.rxns,enzRxnPfx));
        if sum(enzRxnIdx) == 0
            error('No enzyme usage pseudoreactions found in model with given prefix')
        end
    end
    function p = parseInput(arguments)
        % set default values
        ENZ_MET_PFX_DEFAULT = 'prot_';
        ENZ_RXN_pfx_DEFAULT = 'prot_';
        EPSILON_DEFAULT = 1e6;
        THETA_DEFAULT = 0.5;
        LAMBDA_DEFAULT = 1e-7;
        K_DEFAULT = 57500000; % Pyrococcus furiosus; 5.3.1.1; D-glyceraldehyde 3-phosphate
        INCL_UM_DEFAULT = false;
        F_N_DEFAULT = 0.5;
        SIGMA_DEFAULT = 0.5;
        NEG_CORR_F_DEFAULT = false;
        
        % validation functions
        validateModel = @(M) iscell(M)||isstruct(M);
        validScalarDouble = @(v)~ischar(v)&isscalar(v);
        
        p = inputParser;
        p.FunctionName = 'PRESTO';
        
        addRequired(p,'models',validateModel)
        addRequired(p,'expVal',@isnumeric)
        addRequired(p,'E',@isnumeric)
        addParameter(p,'enzMetPfx',ENZ_MET_PFX_DEFAULT,@ischar)
        addParameter(p,'enzRxnPfx',ENZ_RXN_pfx_DEFAULT,@ischar)
        addParameter(p,'epsilon',EPSILON_DEFAULT,validScalarDouble)
        addParameter(p,'theta',THETA_DEFAULT,@iscolumn)
        addParameter(p,'lambda',LAMBDA_DEFAULT,validScalarDouble)
        addParameter(p,'enzBlackList',{''},@iscellstr)
        addParameter(p,'K',K_DEFAULT,validScalarDouble)
        addParameter(p,'includeUM',INCL_UM_DEFAULT,@islogical)
        addParameter(p,'f_n',F_N_DEFAULT,@isnumeric)
        addParameter(p,'sigma',SIGMA_DEFAULT,validScalarDouble)
        addParameter(p,'pCorrFactor',1,@isnumeric)
        addParameter(p,'negCorrFlag',NEG_CORR_F_DEFAULT,@islogical)
        
        parse(p,arguments{:})
    end
end