function [relErr,errVar,sumDelta,objVal,avJD,corrKcatProts] = cvLambdaFitting(varargin)
%% [relErr,errVar,sumDelta,objVal,avJD,corrKcatProts] = cvLambdaFitting(varargin)
% Input:
%   struct model:       base metabolic model (GECKO2.0 format)
%   double expGrowth:   experimental growth rates for all conditions
%   double pTot:        total protein contents for all conditions
%   double E:           enzyme abundance matrix (#model proteins x #conditions)
%   double lambdaParams:array of weighting parameters between minimization 
%                       between relative error to experimental growth rates
%                       and the sum of corrections
%   table nutrExch:     nutrient uptake information (first
%                       column contains reaction IDs, the subsequent
%                       columns contain values for every condition)
%   double kfold:       (optional) number for k-fold cross validation;
%                       default: 3
%   double nIter:       (optional) number of iterations for cross-valiation;
%                       default: 10
%   logical randomize:  (optional) whether sets should be chosen at random
%                       or consecutively into groups of size setSize
%   double epsilon:     (optional) maximum fold change of kcats
%                       (default: 1e6)
%   double theta:       (optional) allowed deviation of the predicted growth rate
%                       from the given experimental value; default: 0.5
%   char enzMetPfx:     (optional) prefix for enzyme mass balance constraints
%                       (i.e. values in model.metNames);default: 'prot_'
%   char enzRxnPfx:     (optional) prefix for enzyme usage pseudoreactions
%                       (i.e. values in model.rxns); default: 'prot_'
%   cellstr enzBlackList:(optional) names of enzyme metabolites that
%                       should be excluded from the correction
%                       default: {''}
%   logical runParallel:(optional) whether or not the function should be
%                       run on multiple threads (default: false)
%   double K:           (optional) maximum kcat that is allowed to be
%                       reached using delta values [s^-1]; default: 57500000 s^-1
%   double GAM:         (optional) growth-associated maintenance:
%                       either as scalar (applied to all conditions) or
%                       as an array (for each condition separately);
%                       default: 50
%   double f:           (optional) protein mass fraction accounted for by
%                       proteins in the model (see GECKO documentation);
%                       default: 0.5
%   double sigma:       (optional) average in vitro enzyme saturation
%                       (can be fitted using GECKO sigmaFitter);
%                       default: 0.5
%   logical negCorrFlag:(optional) add second step that finds negative
%                       corrections for kcat values;
%                       default: false
% Output:
%   double relErr:      average relative error for each lambda in
%                       lambdaParams across all cross-validations in each
%                       iteration (dimension: |lambdaParams| x nIter)
%   double errVar:      average variance of the relative errors across
%                       the cross-validation folds within each iteration
%                       (dimension: |lambdaParams| x nIter)
%   double sumDelta:    average sum of deltas for each lambda in
%                       lambdaParams across all cross-validations in each
%                       iteration (dimension: |lambdaParams| x nIter)
%   double objVal:      average objective values each lambda in
%                       lambdaParams across all cross-validations in each
%                       iteration (dimension: |lambdaParams| x nIter)
%   double avJD:        average Jaccard distance between corrected enzyme
%                       sets between CV folds (dimension: |lambdaParams| x nIter)
%   cell corrKcatProts: unique set of corrected kcats over the CV folds per
%                       iteration (dimension: |lambdaParams| x nIter)
% 22.03.2022 Philipp Wendering, University of Potsdam, philipp.wendering@gmail.com
% 10.12.2022 Philipp Wendering
% * pass optional argument that allows for introdition of negative deltas
%   to the PRESTO function

% parse and check input arguments and set required variables
p = parseInput(varargin);

model = p.Results.model;
expGrowth = p.Results.expGrowth;
pTot = p.Results.pTot;
E = p.Results.E;
nutrExch = p.Results.nutrExch;
lambdaParams = p.Results.lambdaParams;
kfold = p.Results.kfold;
nIter = p.Results.nIter;
epsilon = p.Results.epsilon;
theta = p.Results.theta;
enzMetPfx = p.Results.enzMetPfx;
enzRxnPfx = p.Results.enzRxnPfx;
enzBlackList = p.Results.enzBlackList;
runParallel = p.Results.runParallel;
K = p.Results.K;
GAM = ones(1,size(E,2)).*p.Results.GAM;
f = p.Results.f;
sigma = p.Results.sigma;
negCorrFlag = p.Results.negCorrFlag;

% check dimensions of condition-specific variables
if all([numel(expGrowth) numel(pTot)]==size(E,2))
    nCond = size(E,2);
else
    error('Dimensions of growth rates, protein contents and enzyme abundances do not match')
end

if ~isempty(nutrExch) && size(nutrExch,2)~=size(E,2)
    error('Dimensions nutrient exchanges and enzyme abundances do not match')
end

% find indices of enzyme metabolites that correspond to each of the enzyme
% usage reactions
enzMetIdx = find(contains(model.mets,enzMetPfx));

cmdsz = matlab.desktop.commandwindow.size;
fprintf([repmat('-',1,cmdsz(1)) '\n'])
fprintf('\tRunning %ix %i-fold cross validation scheme\n',nIter,kfold)
fprintf([repmat('-',1,cmdsz(1)) '\n'])
clear cmdsz

relErr   = nan(numel(lambdaParams),nIter);
errVar   = nan(numel(lambdaParams),nIter);
sumDelta = zeros(numel(lambdaParams), nIter);
avJD = zeros(numel(lambdaParams), nIter);
corrKcatProts = cell(numel(lambdaParams), nIter);

objVal   = nan(numel(lambdaParams),nIter);
for l=1:numel(lambdaParams)
    fprintf('EPSILON\t\t%.2g\nTHETA\t\t%.2g\nLAMBDA\t\t%.2g\n\n',epsilon,theta,lambdaParams(l))
    for i=1:nIter
        % create cross-validation partition object
        c = cvpartition(nCond,'KFold',kfold);
        
        fprintf('--- Iteration #%i ---\n',i)
        
        % initialize error array for current iteration
        iterErr      = nan(kfold,1);
        iterSumDelta = zeros(kfold,1);
        iterObj      = nan(kfold,1);
        iterVar      = nan(kfold,1);
        iterUniqCorrKcatProts = cell(kfold,1);
        
        if runParallel
            environment = getEnvironment;
            solver = getCobraSolver('LP');
            parfor s=1:kfold
                log = [];
                % so solver options are still set
                restoreEnvironment(environment);
                changeCobraSolver(solver, 'LP', 0, -1);
                log = [log sprintf('Set #%i...',s)];
                
                % get training and set indices
                trainIdx = training(c,s);
                testIdx = find(test(c,s));
                
                % adjust base model to current conditions
                models = adjBaseModel(model,pTot(trainIdx),nutrExch(:,trainIdx),GAM(trainIdx));
                
                % run kcat correction
                solution.stat = 0;
                tmpModels = models;
                try
                    [solution,tmpModels,~,~,LP] = PRESTO(models,expGrowth(trainIdx),E(:,trainIdx),...
                        'epsilon',epsilon,'lambda',lambdaParams(l),'theta',theta,...
                        'K',K,'enzBlackList',enzBlackList,'enzMetPfx',enzMetPfx,...
                        'enzRxnPfx',enzRxnPfx,'negCorrFlag',negCorrFlag);
                catch ME
                    log = [log ME.stack(1).name ME.message];
                end
                
                if solution.stat == 1
                    iterObj(s) = solution.f;
                    
                    log = [log 'validating...'];
                    
                    % score correction using the remaining conditions
                    tmpErr = nan(numel(testIdx),1);
                    
                    % add pool draw reactions to corrected model
                    corrBatchModel = tmpModels{1};
                    
                    % this is similar to the content of the constrainPool
                    % function in GECKO
                    for j=1:numel(enzMetIdx)
                        corrBatchModel = addReaction(corrBatchModel,['draw_' corrBatchModel.mets{enzMetIdx(j)}],...
                            'metaboliteList', {'prot_pool' corrBatchModel.mets{enzMetIdx(j)}},...
                            'stoichCoeffList', [-corrBatchModel.protMW(j)/1000 1],...
                            'reversible', false);
                        corrBatchModel = removeRxns(corrBatchModel,{[corrBatchModel.mets{enzMetIdx(j)} '_exchange']});
                    end
                    
                    models = adjBaseModel(corrBatchModel,pTot(testIdx),...
                        nutrExch(:,testIdx),GAM(testIdx));
                    
                    % add condition-specific pool reaction
                    for j=1:numel(models)
                        models{j} = addReaction(models{j},'prot_pool_exchange',...
                            'metaboliteList', {'prot_pool'},...
                            'stoichCoeffList', 1,...
                            'upperBound',pTot(testIdx(j))*f*sigma,...
                            'reversible', false);
                        sol = optimizeCbModel(models{j});
                        tmpErr(j) = abs(sol.f-expGrowth(testIdx(j)))/expGrowth(testIdx(j));
                    end
                    
                    iterVar(s) = var(tmpErr,'omitnan');
                    iterErr(s) = mean(tmpErr,'omitnan');
                    
                    log = [log 'DONE'];
                    
                    % get the sum of seltas
                    deltaIdx = contains(LP.rxns,['delta_' enzMetPfx]);
                    iterSumDelta(s) = sum(solution.x(deltaIdx));
                    
                    % store corrected kcats for the CV folds
                    [~,~,pids] = findKcatCorrections(model,tmpModels{1},enzMetPfx);
                    iterUniqCorrKcatProts{s} = unique(pids);
                else
                    log = [log 'INFEASIBLE']
                end
                disp(log)
            end
        else
            for s=1:kfold
                fprintf('Set #%i...',s)
                
                % get training and set indices
                trainIdx = training(c,s);
                testIdx = find(test(c,s));
                
                % adjust base model to current conditions
                models = adjBaseModel(model,pTot(trainIdx),nutrExch(:,trainIdx),GAM(trainIdx));
                
                % run kcat correction
                solution.stat = 0;
                tmpModels = models;
                try
                    [solution,tmpModels,~,~,LP] = PRESTO(models,expGrowth(trainIdx),E(:,trainIdx),...
                        'epsilon',epsilon,'lambda',lambdaParams(l),'theta',theta,...
                        'K',K,'enzBlackList',enzBlackList,'enzMetPfx',enzMetPfx,...
                        'enzRxnPfx',enzRxnPfx,'negCorrFlag',negCorrFlag);

                catch ME
                    if contains(ME.message,'Dual optimality')
                        fprintf([repmat('\b',1,length('OPTIMAL')+1) 'INFEASIBLE: Error: '...
                            ME.stack(1).name ': ' ME.message])
                    else
                        disp(ME)
                    end
                end
                
                if solution.stat == 1
                    iterObj(s) = solution.f;
                    
                    fprintf('validating...')
                    
                    % score correction using the remaining conditions
                    tmpErr = nan(numel(testIdx),1);
                    
                    % add pool draw reactions to corrected model
                    corrBatchModel = tmpModels{1};
                    
                    % this is similar to the content of the constrainPool
                    % function in GECKO
                    for j=1:numel(enzMetIdx)
                        corrBatchModel = addReaction(corrBatchModel,['draw_' corrBatchModel.mets{enzMetIdx(j)}],...
                            'metaboliteList', {'prot_pool' corrBatchModel.mets{enzMetIdx(j)}},...
                            'stoichCoeffList', [-corrBatchModel.protMW(j)/1000 1],...
                            'reversible', false);
                        corrBatchModel = removeRxns(corrBatchModel,{[corrBatchModel.mets{enzMetIdx(j)} '_exchange']});
                    end
                    
                    models = adjBaseModel(corrBatchModel,pTot(testIdx),...
                        nutrExch(:,testIdx),GAM(testIdx));
                    
                    % add condition-specific pool reaction
                    for j=1:numel(models)
                        models{j} = addReaction(models{j},'prot_pool_exchange',...
                            'metaboliteList', {'prot_pool'},...
                            'stoichCoeffList', 1,...
                            'upperBound',pTot(testIdx(j))*f*sigma,...
                            'reversible', false);
                        sol = optimizeCbModel(models{j});
                        tmpErr(j) = abs(sol.f-expGrowth(testIdx(j)))/expGrowth(testIdx(j));
                    end
                    clear corrBatchModel sol
                    
                    iterVar(s) = var(tmpErr,'omitnan');
                    iterErr(s) = mean(tmpErr,'omitnan');
                    
                    fprintf('DONE\n')
                    
                    % get the sum of seltas
                    deltaIdx = contains(LP.rxns,['delta_' enzMetPfx]);
                    iterSumDelta(s) = sum(solution.x(deltaIdx));
                    
                    % store corrected kcats for the CV folds
                    [~,~,pids] = findKcatCorrections(model,tmpModels{1},enzMetPfx);
                    iterUniqCorrKcatProts{s} = unique(pids);
                    
                else
                    fprintf('INFEASIBLE\n')
                end
            end
        end
        errVar(l,i)   = mean(iterVar,'omitnan');
        relErr(l,i)   = mean(iterErr,'omitnan');
        sumDelta(l,i) = mean(iterSumDelta,'omitnan');
        objVal(l,i)   = mean(iterObj,'omitnan');
        corrKcatProts{l,i} = unique(vertcat(iterUniqCorrKcatProts{:}));
        
        jd_mat = nan(kfold);
        for idx1=1:kfold-1
            for idx2=idx1+1:kfold
                jd_mat(idx1,idx2) = 1 - numel(intersect(iterUniqCorrKcatProts{idx1},...
                    iterUniqCorrKcatProts{idx2})) / numel(union(iterUniqCorrKcatProts{idx1},...
                    iterUniqCorrKcatProts{idx2}));
            end
        end
        avJD(l,i) = mean(reshape(jd_mat,1,numel(jd_mat)),'omitnan');
    end
end

    function p = parseInput(arguments)
        % define default variables
        FOLD_DEFAULT = 3;
        EPSILON_DEFAULT = 1e6;
        THETA_DEFAULT = 0.5;
        ENZ_MET_PFX_DEFAULT = 'prot_';
        ENZ_RXN_PFX_DEFAULT = 'prot_';
        NUTR_UPT_DEFAULT = table({''},'RowNames',{'empty'});
        PAR_DEFAULT = false;
        K_DEFAULT = 57500000; % Pyrococcus furiosus; 5.3.1.1; D-glyceraldehyde 3-phosphate
        GAM_DEFAULT = 50;
        ITER_DEFAULT = 10;
        F_DEFAULT = 0.5;
        SIGMA_DEFAULT = 0.5;
        NEG_CORR_F_DEFAULT = false;
        
        % validation for scalar integer inputs
        validScalarDouble = @(v)~ischar(v)&isscalar(v);
        validScalarInt = @(v)validScalarDouble(v)&(floor(v)==ceil(v));
        validateModel = @(m)isstruct(m)&isfield(m,'protMW');
        
        % parse input argument array
        p = inputParser;
        p.FunctionName = 'findBestPerformingSetsByCV';
        
        addRequired(p,'model',validateModel)
        addRequired(p,'expGrowth',@isnumeric)
        addRequired(p,'pTot',@isnumeric)
        addRequired(p,'E',@isnumeric)
        addRequired(p,'lambdaParams',@isnumeric)
        addOptional(p,'nutrExch',NUTR_UPT_DEFAULT,@istable)
        addOptional(p,'kfold',FOLD_DEFAULT,validScalarInt)
        addOptional(p,'nIter',ITER_DEFAULT,validScalarInt)
        addOptional(p,'epsilon',EPSILON_DEFAULT,validScalarDouble)
        addOptional(p,'theta',THETA_DEFAULT,validScalarDouble)
        addOptional(p,'enzMetPfx',ENZ_MET_PFX_DEFAULT,@ischar)
        addOptional(p,'enzRxnPfx',ENZ_RXN_PFX_DEFAULT,@ischar)
        addOptional(p,'enzBlackList',{''},@iscellstr)
        addOptional(p,'runParallel',PAR_DEFAULT,@islogical)
        addOptional(p,'K',K_DEFAULT,validScalarDouble)
        addOptional(p,'GAM',GAM_DEFAULT,@isnumeric)
        addOptional(p,'f',F_DEFAULT,@isnumeric)
        addOptional(p,'sigma',SIGMA_DEFAULT,@isnumeric)
        addOptional(p,'negCorrFlag',NEG_CORR_F_DEFAULT,@islogical)
        
        parse(p,arguments{:})
    end

end