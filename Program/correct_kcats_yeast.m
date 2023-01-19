%% Configuration / Input
clear;clc
run(fullfile('Program', 'configuration_yeast.m'))

% start diary
diary(logFileName)

% set up parallel pool if cross-validation should be run parallelized
if runParallel
    runParallel = setup_parallel(ncpu);
end

%% read data from file and prepare for kcat correction
% read model from file
fprintf('Reading model files...\n')
ecYeastGEM = readGKOmodel(modelFile);

% find indices of enzyme usage pseudoreactions
enzRxnIdx = find(contains(ecYeastGEM.rxns,enzRxnPfx));
enzMetIdx = find(any(ecYeastGEM.S(:,enzRxnIdx),2));
enzMetRxnMatch = arrayfun(@(i)find(any(ecYeastGEM.S(enzMetIdx,i)==1,2)),...
    enzRxnIdx);
ecYeastGEM.enzymes = cellfun(@(x)regexp(x,'[A-Z0-9]{6}','match'),ecYeastGEM.metNames(enzMetIdx));
if exist(mwFile,'file')
    ecYeastGEM.protMW = table2array(readtable(mwFile,'ReadVariableNames',false));
else
    ecYeastGEM.protMW = 0;
end

    if ~exist(mwFile,'file') || numel(ecYeastGEM.protMW)~=numel(ecYeastGEM.enzymes)
        ecYeastGEM.protMW = getMWfromUniProtID(ecYeastGEM.enzymes);
        writetable(array2table(ecYeastGEM.protMW),mwFile,'WriteVariableNames',false)
    end

% read experimental data
[condNames, E, expVal, nutrExch, P] = readChenetal(topDir, ecYeastGEM);

%% Exclude inappropriate datasets
% exclude Lathvee2017_Temp33,Lathvee2017_Temp36,Lathvee2017_Temp38
idxExclude = contains(condNames,{'Lahtvee2017_Temp33',...
    'Lahtvee2017_Temp36','Lahtvee2017_Temp38'});
expVal = expVal(~idxExclude);
P = P(~idxExclude);
E = E(:,~idxExclude);
nutrExch = nutrExch(:,~idxExclude);
condNames = condNames(~idxExclude);

%% fit growth-associated maintenance
% adjust GAM
if isnan(GAM)
    % read batch model
    batchModel = readGKOmodel(batchModelFile);
    % parse growth data into chemostatData.tsv (GECKO toolbox)
    header = {'Drate' 'GlucoseUptake' 'O2uptake' 'CO2production'};
    searchStrings = {'r_1714','r_1992','r_1672'};
    rowIdx = cellfun(@(x)find(ismember(nutrExch.Row,x)),searchStrings);
    chemostatMatrix = abs(table2array(nutrExch(rowIdx,:))');
    chemostatMatrix = [expVal' chemostatMatrix];
    
    GAM = nan(1,size(chemostatMatrix,1));
    for i=1:size(chemostatMatrix,1)
        writetable(array2table(chemostatMatrix(i,:),'VariableNames',header),...
            fullfile(geckoDir, 'Databases', 'chemostatData.tsv'),'FileType','text',...
            'Delimiter','\t')
        fprintf('fitting GAM for condition #%i\n',i)
        % change into GECKO directory to be able to call fitGAM function
        cd(fullfile(geckoDir, 'geckomat', 'limit_proteins'))
        
        try
            [~,GAM(i)] = scaleBioMass(batchModel,P(i));
        catch ME
            disp(ME.message)
        end
        cd(topDir)
    end
    
    clear header searchStrings rowIdx chemostatMatrix batchModel
end
GAM(isnan(GAM)) = mean(GAM,'omitnan');

% use only the subset of abundances for those proteins that have been
% measured across all conditions
E(~all(E>0,2),:) = 0;

%% Find the best lambda
lambdaParams = logspace(-14,-1,14);
% run cross-validation
[relErr,errVar,sumsDelta,objVal,avJD,corrKcatProts] = ...
    cvLambdaFitting(ecYeastGEM,expVal,P,E,lambdaParams,...
    nutrExch,'epsilon',epsilon,'theta',theta,'runParallel',runParallel,...
    'GAM',GAM, 'f', f, 'sigma', sigma);
save(fullfile('Results','lambda_cv','lambda_fitting_s_cerevisiae'),'relErr','errVar','sumsDelta',...
    'objVal','lambdaParams','corrKcatProts','avJD')
% find best lambda value
lambda = selectBestLambda(lambdaParams, relErr, sumsDelta);

%% Run kcat correction with all conditions
models = adjBaseModel(ecYeastGEM,P,nutrExch,GAM);
[solution,corrModels,relError,~,LP] = PRESTO(models,expVal,E,...
    'epsilon',epsilon,'lambda',lambda,'theta',theta);

% find kcat corrections
[o,c,proteinIds] = findKcatCorrections(ecYeastGEM,corrModels{1},enzMetPfx);
[proteinIds,ia] = unique(proteinIds,'stable');
deltaVal = c(ia)-o(ia);
corrKcats = c(ia);
clear c o

%% Run variability analysis and sampling for deltas
[minDelta,maxDelta] = deltaVariabilityAnalysis(LP,solution);

save(fullfile('Results','delta_va',['deltaVariability' '_s_cerevisiae' '_lambda_' num2str(lambda)]),...
    'minDelta','maxDelta','deltaVal','corrKcats','proteinIds','deltaSamplingMat')

delete(gcp('nocreate'));
%% Run validation comparison with pFBA and GECKO models; Run corrected kcat analysis
[relE, Mu]=comp_PRTcondmod(corrModels, false, 'saccharomyces cerevisiae', 'realprot_lambda1e_7');
diary off

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
