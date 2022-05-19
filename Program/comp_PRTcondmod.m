function [relE, Mu] = comp_PRTcondmod(models,  prot_cor, org_name, figprefix)
%Function to compare the final PRESTO models
%with GECKO and pFBA models
%INPUT:
% - cell models: A cell array containing models in GECKO format as produced
%                by the kcat correction code
% - logic prot_cor: Passed to comp_condmod indicates if protein abundancies
%                   should be relaxed using the correction factor
% - char org_name:  The model organism
% - char figprefix: (optional, default: []) A character vector given the
%                     file prefix with which the plots are saved to
%                     "Results/Manuscript/Fig1/figprefix..." and
%                     "Results/Manuscript/Fig2/figprefix..."

%OUTPUT:
% - unbError: matrix giving for each model the rel error in every
%                   experimental condition when no nutrient uptake rates
%                   are used as additional constrains
% - bError: matrix of rel. error with uptake rates as additional model 
%               constrains
if nargin<4
    figprefix=[];
end
%set the colormap
switch org_name
    case 'saccharomyces cerevisiae'
        configuration_yeast
        condNames=readChenetal(topDir, [], true);
        %set label for pfba dataset used for comparison 
        pfbalab='Chen';
        %Create alternative cond_name prefixes according to carbon source
        altcondNames=repelem({'GLUC_'}, length(condNames));
        altcondNames(contains(condNames,{'EtOH', 'Etoh'}))={'ETOH_'};
        altcondNames(contains(condNames,{'Gln'}))={'GLUC+GLN_'};
        altcondNames(contains(condNames,{'Phe'}))={'GLUC+PHE_'};
        altcondNames(contains(condNames,{'Ile'}))={'GLUC+ILE_'};
        altcondNames=strcat(altcondNames', condNames');
    case 'escherichia coli'
        configuration_ecoli
        condNames=readDavidi2016([],topDir);
        %set lebel for pfba dataset used for comparison 
        pfbalab='Davidi';
        %Create alternative cond Names based on publication
        altcondNames=cell(1, length(condNames));
        for i=1:length(condNames)
            tmpSplit=strsplit(condNames{i}, '_');
            switch(tmpSplit{end})
                case 'V'
                    altcondNames{i}='Valgepea2013_';
                case 'S'
                    altcondNames{i}='Schmidt2015_';
                case 'P'
                    altcondNames{i}='Peebo2015_';
            end
        end
        altcondNames=strcat(altcondNames', extractBefore(condNames, '_mu')');
end

%import the sigma values used in GECKO models to set equal constrains
[loginf, ~]= getloginf(gkologFile, org_name);
if length(loginf.condNames)~=length(condNames)
    error("GECKO model creation log does not contain model for each condition. Matchin of sigma values impossible")
elseif length(models)~=length(condNames)
    error("Number of models does not match number of conditions. Matching of Sigma ismpossible")
end
%order GECKO model creation logfile according to conditions
[match, idx]=ismember(condNames, loginf.condNames);
if any(~match, 'all')
    error("Conditions in GECKO model creation log does not match experimental conditions")
end
loginf=loginf(idx,:);


%convert modeld to batch models
cd([geckoDir, '/geckomat/limit_proteins'])

batch_models=cell(1,length(models));
%get the amount of protein mass accounted for 
%by enzymes in the model
[f,~] = measureAbundance(models{1}.enzymes);
disp(['Fraction of model enzymes in total protein content (f) = ', string(f)])
%build initial model
load('../../databases/parameters.mat')
%check if COBRA toolbox is setup correctly
if ~strcmp(parameters.org_name, org_name)
    error(['matlab toolbox is not configured for organism ' org_name])
end
parameters.sigma=loginf.Sigma(1);
models{1}.MWs=models{1}.protMW/1000; %GECKO uses kg per mol 
batch_models{1}= constrainEnzymes(models{1}, f, true, [], sumProtein(models{1}));
batch_models{1}.csense=repmat('E', length(batch_models{1}.mets),1);

% restore default sigma
load('../../databases/parameters.mat')
parameters.sigma=0.5;
save('../../databases/parameters.mat', 'parameters')


%Adapt Biomass and protein pool constrain for all other conditions
for i=2:length(models)
    batch_models{i}=scaleBioMass(batch_models{1}, sumProtein(models{i}), GAM(i));
    batch_models{i}.ub(end)=f*loginf.Sigma(i)*sumProtein(models{i});
end
cd(topDir)
%since PRESTO models contain only one set of kcat adaptions only test
%differenc conditions
[relE, Mu]=comp_condmod(batch_models, org_name, false, prot_cor);

if ~isempty(figprefix)
    %add figure prefix to file 
    switch org_name
    case 'saccharomyces cerevisiae'
        figprefix=['yeast' figprefix];
    case 'escherichia coli'
        figprefix=['ecoli' figprefix];
    end
    
    %get gekco model perfomance for comparison
[gkorelE, ~, maxgkorelE, ~, max_gkomod]=comp_GKOcondmod(fullfile(geckoDir, 'models'), gkologFile, org_name, prot_cor);

%generate result path
if ~isdir(fullfile(topDir, 'Results', 'relE'))
    mkdir(fullfile(topDir, 'Results', 'relE'))
end
%Generate a report table of the errors witout protein pool constrain
%(Support Table S1)
gko_mean=mean(gkorelE{4}, 2);
gko_sd=std(gkorelE{4},1, 2);
nopool_relEtab=table(condNames', relE{4}, gko_mean, gko_sd, 'VariableNames', {'Condition Chen', 'PRESTO relative Error', 'mean GECKO relative Error', 'sd GECKO relative Error'});
writetable(nopool_relEtab, fullfile('Results', 'relE', [figprefix, 'nopool_relE.tsv']), 'FileType', 'text')

%Get model performance for pFBA fitted kcats
[pFBArelE, pFBAmu]=comp_condmod(repelem({create_pFBAmod(readGKOmodel(batchModelFile), enzMetPfx)}, length(condNames)), org_name, false, prot_cor);

%compare kcats and plot enrichment analysis and scatter of GECKO PRESTO
%(Fig3, Fig 12)
%overlap, save a boolean table indicating whhich kcats are changed in the
%manual modifications, GECKO and PRESTO and save a table with the actual
%corrected values of PRESTO and GECKO (Supplementary Table S2,3)
kcat_comptab=comp_kcat2(models{1}, max_gkomod, 'Protein',figprefix);
%plot for no uptake and proteomics constrains for bestmod
%for all

%% rel Error all GECKO and PRESTO (Figure 2 among others)
plotFIGrelE(relE([1 2 3 5]), gkorelE([1 2 3 5]), condNames, figprefix, topDir)
%same plots with alternative condNames (Figure 4 among others
plotFIGrelE(relE([1 2 3 5]), gkorelE([1 2 3 5]), altcondNames, [figprefix '_alt'], topDir)
%plot relative Error of 
plotFIGrelE3(relE(1:3), maxgkorelE(1:3), pFBArelE(1:3), [figprefix '_pFBAscat'], topDir, pfbalab);

end
end
