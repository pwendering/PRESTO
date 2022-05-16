%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%  
%%% Configuration file for PRESTO                                       %%%
%%%  -execute this script in the top-level directory of the approach    %%%
%%%                                                                     %%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify organism name
orgName = 'Escherichia coli';
%Model basename used
orgBasename = 'ecEcoli';
% Growth-associated maintenance (GAM); set to NaN if it should be fitted
GAM = repelem(75.5522,31);
% mass fraction of all proteins included in the model (see GECKO documentation)
f = 0.3827;
% mass fraction of unmeasured proteins according to PAX DB (also discounting
% proteins that do not have a measured abundance across all conditions)
f_n = 0.0633;
% average in vitro enzyme saturation (fitted in GECKO)
sigma = 0.5;
% correction factor for protein abudances
protCorrFact = NaN;
% define top-level directory
topDir = pwd;
% solver for linear optimization
cobraSolver = 'gurobi';
% define log file name )
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

mwFile = fullfile(topDir, 'Data', 'MW_ecoli.txt');
%log from GECKO toolbox run
gkologFile=fullfile(topDir, 'Logs', 'ecoli_getcondec.log');
% correction parameters
epsilon = 1e5;
lambda = 1e-5;
theta = 0.6; % relative error
% add functions to MATLAB path
addpath(genpath(fullfile(topDir, 'Program')))
% check if COBRA toolbox is installed
try
    changeCobraSolver(cobraSolver,'LP',0);
    changeCobraSolverParams('LP', 'feasTol', 1e-9);
catch
    error('COBRA toolbox is not installed')
end

% get maximum kcat value
geckoDir = fullfile(topDir, 'GECKO_E_coli');
maxKcatFile = fullfile(geckoDir, 'databases', 'max_KCAT.txt');
% model and data files
modelFile = fullfile(geckoDir, 'models', 'ecEcoli', 'rawecEcoli.mat');
batchModelFile = fullfile(geckoDir, 'models', 'ecEcoli', 'rawecEcoli_batch.mat');
kcatoriginFile=fullfile(geckoDir, 'models', 'ecEcoli', 'ecEcoli_kcatOrigins.txt');
K = retrieveMaxKcat(maxKcatFile,orgName);
clear maxKcatFile
