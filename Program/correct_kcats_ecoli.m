%% Configuration / Input
clear;clc
configuration_ecoli

% start diary
% diary(logFileName)

% set up parallel pool if cross-validation should be run parallelized
if runParallel
    runParallel = setup_parallel(ncpu);
end

%% read data from file and prepare for kcat correction
% read model from file
fprintf('Reading model files...\n')

% load model
model = readGKOmodel(fullfile(modelFile));
model.ub(isinf(model.ub)&model.ub>0) = 1000;
model.lb(isinf(model.lb)&model.lb<0) = -1000;

% block all uptake reactions
model.ub(contains(model.rxns,'_REV')&model.ub>0&sum(model.S~=0)'==1) = 0;

% find indices of enzyme usage pseudoreactions
enzRxnIdx = find(contains(model.rxns,enzRxnPfx));
enzMetIdx = find(any(model.S(:,enzRxnIdx),2));
model.enzymes = cellfun(@(x)regexp(x,'[A-Z0-9]{6}','match'),model.metNames(enzMetIdx));
if exist(mwFile,'file')
    model.protMW = table2array(readtable(mwFile,'ReadVariableNames',false));
else
    model.protMW = 0;
end

if ~exist(mwFile,'file') || numel(model.protMW)~=numel(model.enzymes)
    model.protMW = getMWfromUniProtID(model.enzymes)';
    writetable(array2table(model.protMW),mwFile,'WriteVariableNames',false)
end

% read experimental rata from file
[condNames, E, expVal, nutrExch, P] = readDavidi2016(model);

% set allowed nutrient uptakes to default flux bound
defaultNutrExch = nutrExch;
defaultNutrExch{:,:} = cell2mat(arrayfun(@(i)1000*double(sign(nutrExch.(i))),1:size(nutrExch,2),'un',0));

% use a single GAM value
if any(isnan(GAM), 'all')
cd(fullfile(geckoDir, 'geckomat', 'limit_proteins'))
[model,GAM] = scaleBioMass(model,P(1),[]);
GAM = repelem(GAM,numel(condNames),1);
cd(topDir)
end

% use only the subset of abundances for those proteins that have been
% measured across all conditions
E(~all(E>0,2),:) = 0;

% add COBRA constaint sense field
model.csense = repelem('E',size(model.S,1),1);

%% Find the best lambda
lambdaParams = logspace(-14,-1,14);
% run cross-validation
[relErr,errVar,sumsDelta,objVal,avJD,corrKcatProts] = ...
    cvLambdaFitting(model,expVal,P,E,lambdaParams,...
    defaultNutrExch,'epsilon',epsilon,'theta',theta,'runParallel',runParallel,...
    'GAM',GAM, 'f', f, 'sigma', sigma);
save(fullfile('Results','lambda_cv','lambda_fitting_e_coli'),'relErr','errVar','sumsDelta','objVal',...
    'lambdaParams','corrKcatProts','avJD')

%% Run kcat correction with all conditions
models = adjBaseModel(model,P,defaultNutrExch,GAM);
[solution,corrModels,relErrCorr,changeTab,LP] = PRESTO(models,expVal,E,...
    'epsilon',epsilon,'lambda',lambda,'theta',theta);

[o,c,proteinIds] = findKcatCorrections(model,corrModels{1},enzMetPfx);
[proteinIds,ia] = unique(proteinIds,'stable');
deltaVal = c(ia)-o(ia);
corrKcats = c(ia);
clear c o

%% Run variability analysis and sampling for deltas
[minDelta,maxDelta] = deltaVariabilityAnalysis(LP,solution);

deltaSamplingMat = deltaSampling(LP,solution,minDelta,maxDelta,'nSamples',1e4);

save(fullfile('Results','delta_va',['deltaVariability','_e_coli','_lambda_',num2str(lambda),'.mat']),...
    'minDelta','maxDelta','deltaVal','corrKcats','proteinIds','deltaSamplingMat')

delete(gcp('nocreate'));
diary off
%% Run validation comparison with pFBA and GECKO models; Run corrected kcat analysis

[relE, Mu]=comp_PRTcondmod(corrModels, false, 'escherichia coli', 'lambda1e5');

function runParallel = setup_parallel(ncpu)
runParallel = true;
p = gcp('nocreate');
if isempty(p)
    try
        parpool(ncpu);
    catch ME
        disp([ME.identifier ': ' ME.message]); clear ME
        fprintf('Starting analysis without parallelization.\n\n')
        runParallel = false;
    end
elseif p.NumWorkers ~= ncpu
    fprintf('Parallel pool with %d workers exist',p.NumWorkers)
    try
        delete(p);
        parpool(ncpu);
    catch ME
        disp([ME.identifier ': ' ME.message]); clear ME
        fprintf('Starting analysis without parallelization.\n\n')
        runParallel = false;
    end
end
end
