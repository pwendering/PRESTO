function get_GKOmod_ecoli(unmod, start_i)
%Function to obtain EC models from GECKO branch
%INPUT:
%  - log unmod: (optional, defaul: true) logcial indicating if models should
%                   be build excluding GECKO manual modifications
%  - int start_i: (optional, default:1)giving the index of model condition to
%            start with, to avoid looping over all models after error etc
%  
%Import environment variables
disp('ATTENTION: GECKO is restricted to create identical files for unmod and manual modified models.')
disp('SAVE THE OUTPUT IN GECKO/models/ BEFORE CALLING THIS FUNCTION WITH A DIFFERENT unmod argument')
if nargin<1
    unmod=true;
end
if nargin<2
    start_i=1;
end
configuration_ecoli

if ~isfolder('Logs')
    mkdir Logs
end
diary('Logs/ecoli_getcondec.log')
%Load GEM
load('Data/ecoli-GEM.mat')
%read in chen data
[condNames, E, expVal, nutrExch, P] = readDavidi2016([], topDir);
cd([geckoDir '/geckomat'])
load('../databases/parameters.mat')
 %check for correct organism 
 if ~strcmp(parameters.org_name, 'escherichia coli')
     error('GECKO parameters do not belong to Ecoli model, check if scripts have been copied from ecModels folder')
 end
  %save initial parameters
 sav_Ptot=parameters.Ptot;
 sav_gR_exp=parameters.gR_exp;
 clear parameters
%select conditions to build model from 
rng(2021)
cond_idx=randsample(1:length(condNames), length(condNames));
%For each selected condition
for i=start_i:length(cond_idx)
    tmp_idx=cond_idx(i)
    %Update the parameters object
    load('../databases/parameters.mat')
    parameters.Ptot=P(tmp_idx);
    parameters.gR_exp=expVal(tmp_idx);
    save('../databases/parameters.mat', 'parameters')
    clear parameters
    
    enhanceGEM(model, 'COBRA', unmod, ['ecEcoli_', condNames{tmp_idx}], '1')
end

 load('../databases/parameters.mat')
 parameters.Ptot=sav_Ptot;
 parameters.gR_exp=sav_gR_exp;
 clear parameters
cd(topDir)
diary off
end



