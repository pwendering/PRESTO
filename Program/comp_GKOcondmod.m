function [relError, Mu, predE, ovlp, maxrelError, maxMu, maxpredE, maxfluxvar, max_gkomod] = comp_GKOcondmod (moddir, logfile, org_name, prot_cor)
% Function to compare produced GECKO models for Yeast model and 
%Chen et al data
%INPUT: 
% - char moddir: absolute path to directory containing one directory per 
%                         produce gecko model named ecYeast_<cond> 
%                         eg. ecYeast_Lahtvee2017_Osmo06
% - char logfile:   absolute path to the logfile produce by get get_GKOmod.m
% - char org_name:  The model organism
% - logic prot_cor: Passed to comp_condmod indicates 
%                   if protein abundancies should be relaxed using the correction factor
%OUTPUT:
% - cell relError:  relative errors from comp_condmod for all models in all
%                   conditions
% - cell Mu:    predicted growth rates from comp_condmod for all models in 
%               all conditions
% - cell predE: spearman correlation betwen predicted and measured protein
%               abundances from comp_condmod fro all models in all
%               conditions
% - logic ovlp:     A matrix of logical values rows indicated kcat
%                   corrections columns give the models ordered according
%                   to condNames.                
% - cell maxrelError:   rel errors from comp_condmod for model containing the
%                       maximumg combination of kcat adjustments combining
%                       all condition specific adjustments
% - cell maxMu: predicted growht rates from maximum model in all conditiosn
% - double maxpredE:    spearman correlation betwen predicted and measured protein
%                       abundances from comp_condmodfor model containing the
%                       maximumg combination of kcat adjustments combining
%                       all condition specific adjustments
% - struct max_gkomod:  GECKO model with pool constrains with maximum of
%                       all kcat adaptions (combining all condition
%                       specific adaptions)


switch org_name
    case 'saccharomyces cerevisiae'
        configuration_yeast
        condNames=readChenetal(topDir, [], true);
    case 'escherichia coli'
        configuration_ecoli
        condNames=readDavidi2016([], topDir);
end

loginf= getloginf(logfile, org_name);
basic_batchmod=readGKOmodel(batchModelFile);


%% read in models
models=cell(size(loginf, 1),1);


%save gecko cat changes per model
kcat_tabs=cell(size(loginf,1),1);
for i=1:size(loginf, 1)
    files=dir(fullfile(moddir, [orgBasename, '_', loginf.condNames{i}]));
    files=files(~[files.isdir]);
    files={files.name};
    match=~cellfun(@isempty, regexp(files, '^adp.*_batch.mat'));
    if sum(match, 'all')~=1
    error('Multiple modelfiles for load detected adapt regex to find model file or check moddir argument')
    end
    models{i}=readGKOmodel(fullfile(moddir, [orgBasename, '_', loginf.condNames{i}], files{match}));
    %get GECKO  modifications
    %% create cumulative table of kcat corrections

    [~, ~ , ~,gkoreportTable] = findKcatCorrections(basic_batchmod, models{i}, enzRxnPfx);
    kcat_tabs{i}=gkoreportTable;
    if i==1
        kcat_tab=gkoreportTable;
    else
        kcat_tab=[kcat_tab; gkoreportTable(~(ismember(gkoreportTable, kcat_tab)),:)];
    end
end

%generate a overlap matrix
ovlp=nan(size(kcat_tab,1),length(kcat_tabs));
for i=1:length(kcat_tabs)
    temptab=kcat_tabs{i};
    ovlp(:,i)=ismember(kcat_tab(:,1:2), temptab(:,1:2));
end

%generate the union of 
clear gkoreportTable

%create a maximum correction model
max_gkomod=models{1};
[~,~, ~, gkoreportTable]= findKcatCorrections(basic_batchmod, models{1}, enzRxnPfx);
%add corrections not yet present in the model 1

 for i=find(~ismember(kcat_tab, gkoreportTable))'
     %if the reported intial kcat is consistent
     if -1/max_gkomod.S(contains(max_gkomod.mets, kcat_tab.("PROTEIN ID")(i)), kcat_tab.("REACTION ID")(i))/3600 - ...
             kcat_tab.("KCAT ORIG [s^-1]")(i)<1e-6
         %update kcat
         max_gkomod.S(contains(max_gkomod.mets, kcat_tab.("PROTEIN ID")(i)), kcat_tab.("REACTION ID")(i))=-1/(kcat_tab.("KCAT UPDATED [s^-1]")(i)*3600);
     else
         error('Stochiometric matrix entries are not consistent, check raw batch model and code...')
     end
 end

%if there is a model for each modelling condition reorder them  and the overlap table previous to
%optimization 
if length(loginf.condNames)==length(condNames)
    [match, idx]=ismember(condNames, loginf.condNames);
    if any(~match, 'all')
        error("Duplicated modelling conditions, check log file")
    end
    models=models(idx);
    ovlp=ovlp(:,idx);
else
    warning(["less models then experimental conditions, The returned matrix has rows", ...
        "reordered so that the upper part matches the conditions in which models are available"])
end
    
%report performance for all gecko models
[relError, Mu, predE]=comp_condmod(models,org_name, true, prot_cor);
%report performance for maximum correction model 
[maxrelError, maxMu, maxpredE, maxfluxvar]=comp_condmod(repelem({max_gkomod}, length(condNames)), org_name, false, prot_cor);
%if there are less models then conditions reorder them previous now
if length(loginf.condNames)~=length(condNames)%order results
     [match, idx]=ismember( loginf.condNames, condNames);
     newidx=[idx(match)', setdiff(1:length(condNames), idx(match))];
    for i=1:length(relError)
        relError{i}=relError{i}(newidx, :);
        Mu{i}=Mu{i}(newidx, :);
        maxrelError{i}(newidx,:);
        maxMu{i}(newidx,:);
    end
end

%old example code for heat map creation
% heatmap(loginf.condNames, condNames(newidx), relError{4}(newidx,:), ...
%      'Title', 'rel. error predicted growth unbound + corrected proteomics nopool')
%  saveas(gcf, 'Results/gko_unbcorprotbatch.svg')
end

