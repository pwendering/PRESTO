%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%  
%%% Configuration file for PRESTO                                       %%%
%%%  -execute this script in the top-level directory of the approach    %%%
%%%                                                                     %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify organism name
orgName = 'Saccharomyces cerevisiae';
%Model basename used
orgBasename = 'ecYeast';
% Growth-associated maintenance (GAM); set to NaN if it should be fitted
GAM = [23.4 9 9 67.5 24.8 21.5 73.3 47.6 93.8 95.2 53.8 161 9 87.2 13.2 63.4 ...
    83.3 66.5 69.5 98.1 12.3 9 70.4 124.4 136.7 105.9 91.2];
% GAM = NaN;
% mass fraction of all proteins included in the model (see GECKO documentation)
f = 0.4423;
% mass fraction of unmeasured proteins according to PAX DB (also discounting
% proteins that do not have a measured abundance across all conditions)
f_n = 0.1124;
% average in vitro enzyme saturation (fitted in GECKO)
sigma = 0.5;
% correction factor for protein abudances
protCorrFact = [1.8898,2.0621,1.6025,1.7599,1.3023,1.4009, 1.382,1.9242,...
    2.5238,3.5113,3.4852, 1.385,1.3633,2.2578,2.4316,3.2938,2.8668, 2.909,...
    2.191, 1.877,1.5624,2.2166,2.4049,1.9211,1.9156,2.4875,1.7164]';
% define top-level directory
topDir = char(regexp(pwd, '^.*PRESTO', 'match'));
% solver for linear optimization
cobraSolver = 'gurobi';
% define log file name (defaults to current date and time)
logFileName = fullfile('Logs', ['presto_' regexprep(datestr(datetime),'[:-\s]','_') '.log']);
% prefixes for enzyme metabolites and enzyme usage reaction in the GECKO
% model
enzMetPfx = 'prot_';
enzRxnPfx = 'prot_';
% Specify whether the approach should be run parallelized
runParallel = true;
ncpu = 20;
% set the number of iterations of k-fold cross-validation
nIter = 50;

% model and data files
mwFile = fullfile(topDir, 'Data','MW_yeast.txt');
protAbcFile = fullfile(topDir, 'Data','Chen','abs_proteomics_yeast.tsv');
growthDataFile = fullfile(topDir, 'Data','Chen','growth_rates_yeast.tsv');
ptotFile = fullfile(topDir, 'Data', 'Chen', 'total_protein_yeast.tsv');
gkologFile=fullfile(topDir, 'Logs', 'yeast_getcondec.log');

% correction parameters
epsilon = 1e5;
lambda = 1e-7;
theta = 0.6; % relative error
% add functions to MATLAB path
addpath(genpath(fullfile(topDir, 'Program')))
% check if COBRA toolbox is installed
try
    changeCobraSolver(cobraSolver,'LP',0);
    changeCobraSolverParams('LP','feasTol',1e-9);
catch
    error('COBRA toolbox is not installed')
end

% get maximum kcat value
geckoDir = fullfile(topDir, 'GECKO_S_cerevisiae');
maxKcatFile = fullfile(geckoDir, 'databases', 'max_KCAT.txt');
K = retrieveMaxKcat(maxKcatFile,orgName);
clear maxKcatFile
modelFile = fullfile(geckoDir, 'models', 'ecYeast', 'rawecYeast.mat');
batchModelFile = fullfile(geckoDir,  'models', 'ecYeast', 'rawecYeast_batch.mat');
kcatoriginFile=fullfile(geckoDir, 'models', 'ecYeast', 'ecYeast_kcatOrigins.txt');
