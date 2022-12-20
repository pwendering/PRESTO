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
        altcondNames=strrep(condNames, 'Yu2021', 'Y21');
        altcondNames=strrep(altcondNames, 'Yu2020', 'Y20');
        altcondNames=strrep(altcondNames, 'Lahtvee2017', 'L');
        altcondNames=strrep(altcondNames, 'DiBartolomeo2020', 'D');
        maxmin_relE={};
        maxmin_predE={};
    case 'escherichia coli'
        configuration_ecoli
        condNames=readDavidi2016([],topDir);
        %Put the study abbrevation at the beginning of the strind
        altcondNames=cell(1,length(condNames));

        for i=1:length(condNames)
            tmpSplit=strsplit(condNames{i}, '_');
            altcondNames{i}=strjoin([tmpSplit(end), tmpSplit(1:(end-1))], '_');
        end
        %set lebel for pfba dataset used for comparison
        pfbalab='Davidi';
        %For now just load model with minimzed kcats in second step
        maxmin_mod=readGKOmodel(fullfile('Data', 'model_max_min_deltas.mat'));
        maxmin_batch_models=addGKOconstrains(maxmin_mod, org_name);
        %calculate performacnce for these models
        [maxmin_relE, ~, maxmin_predE] = comp_condmod(maxmin_batch_models, org_name, false, prot_cor);
end

batch_models=addGKOconstrains(models, org_name);
%since PRESTO models contain only one set of kcat adaptions only test
%differenc conditions
[relE, Mu, predE, fluxvar]=comp_condmod(batch_models, org_name, false, prot_cor);


if ~isempty(figprefix)
    %add figure prefix to file 
    switch org_name
    case 'saccharomyces cerevisiae'
        figprefix=['yeast' figprefix];
    case 'escherichia coli'
        figprefix=['ecoli' figprefix];
    end
    
    %get gekco model perfomance for comparison
[gkorelE, ~, gkopredE, maxgkorelE, ~, maxgkopredE, maxgkofluxvar, max_gkomod]=comp_GKOcondmod(fullfile(geckoDir, 'models'), gkologFile, org_name, prot_cor);

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
%Use shortened condition names for plot labels
plotFIGrelE(relE([1 2 3 5]), gkorelE([1 2 3 5]), altcondNames, figprefix, topDir)
%plot figure of spearman correlation between measured and corrected 
plotFIGpredE(predE, gkopredE, altcondNames, figprefix, topDir)
if ~isempty(maxmin_relE)
    %plot a supplementary plot with the min model performance
    plotFIGrelE(relE([1 2 3 5]), gkorelE([1 2 3 5]), altcondNames, [figprefix '_maxmin'], topDir, maxmin_relE([1 2 3 5]))
    %plot figure of spearman correlation between measured and corrected 
    plotFIGpredE(predE, gkopredE, altcondNames, figprefix, topDir, maxmin_predE)
end
%plot relative Error of 
plotFIGrelE3(relE(1:3), maxgkorelE(1:3), pFBArelE(1:3), [figprefix '_pFBAscat'], topDir, pfbalab);

plotFIGvrange(fluxvar, maxgkofluxvar, batch_models{1}, enzRxnPfx, altcondNames, figprefix, topDir)
end
end
