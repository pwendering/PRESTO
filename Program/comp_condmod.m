function[relError, Mu, pred_E, fluxVar] = comp_condmod(models, org_name, test_unseen, prot_cor)
%Function to run the benchmarking using the respective Ptot and
%biomass function for each condition.
%INPUT:
% - struc models;: a cell array of model structures.
% - char org_name:  The model organism
% - logical test_unseen: optional logical if ture models will be
%                        validated in all conditions, if false only in 
%                        the one they are build for. If true relError and Mu 
%                        cell arrays contain a matrix (conditions x models) else
%                        they only contain a vector where each model is validated in 
%                        respective condition. If true pred_E and fluxVar
%                        return a matrix (condition x models x values) 
%                        else a matrix 2D matrix where each row corresponds to one
%                        condition specific model prediction/variance 
%                        defaul: true
% - logical prot_cor: optional logical if true total abundancy values are
%                       corrected by the correction factor. default: false
%OUTPUT:
% - cell relError:  1x5 cell giving the relative error for models
%                   constrained by: 1 - total protein pool constrain,
%                   2-pool + experimental uptake rates (if available else 1=2), 3
%                   pool, uptake + absolute protein abundancies, 4- only
%                   absolute protein abundancies 5- only experimental
%                   uptake rates (if available, else large overprediction)
% - cell Mu:    1x5 cell giving the respective predicted growth rates (s.
%               relError)
% - double pred_E:     matrix (conditions x models x enzRxns giving the predicted enyzme usage 
%                      coefficients with  no constrains on uptake 
% - double fluxVar: matrix (conditions x models x rxns) giving the flux
%                   ranges calculated using FVA in the szenario constrained
%                   with pool, uptake and absolute protein abundance
if nargin<3
    test_unseen=true;
elseif ~test_unseen
    disp("test_unssen is false, IT IS ASSUMED MODELS ARE ORDERED ACCORDING TO CondNames")
end
if nargin<4
    prot_cor=false;
    disp('prot_cor not specified, protein abundances are not corrected')
end

switch org_name
    case 'saccharomyces cerevisiae'
        configuration_yeast
    case 'escherichia coli'
        configuration_ecoli
end
%Check if all models have the same reaction ordering
%Here the pool reaction can be included
enz_mappable=true;
for i=2:length(models)
    if ~isequal(models{1}.rxns, models{i}.rxns)
        if length(models{1}.rxns) ~=length(models{i}.rxns)
            error("Different number of model reactions detected. Aborting")
        else
            warning("models have different ordering of proteins. Enzyme abundance constrains not supported")
            enz_mappable=false;
            switch org_name
                case 'saccharomyces cerevisiae'
                    [condNames, ~, expVal, nutrExch, P]=readChenetal(topDir, [], true);
                case 'escherichia coli'
                    [condNames, ~, expVal, nutrExch, P]=readDavidi2016([], topDir);
            end
            break
        end
    end
end
enzRxnIdx = find(contains(models{1}.rxns,enzMetPfx)&~ismember(models{1}.rxns, {'prot_pool_exchange'}));

if enz_mappable
    %Create a mapping for proteomics to enzyme constrains
    enzMetIdx = find(any(models{1}.S(:,enzRxnIdx),2)&~ismember(models{1}.mets, {'prot_pool'}));
    enzMetRxnMatch = arrayfun(@(i)find(any(models{1}.S(enzMetIdx,i)==1,2)),...
        enzRxnIdx);
                switch org_name
                case 'saccharomyces cerevisiae'
                    [condNames, E, expVal, nutrExch, P]=readChenetal(topDir, models{1}, true);
                case 'escherichia coli'
                    [condNames, E, expVal, nutrExch, P]=readDavidi2016(models{1},topDir);
            end

    %create relaxed enzyme abundancy constrains using the correction
    %factor
    if prot_cor
        if isnan(protCorrFact)
            error('protCorrFact is not declared, corrected protein abundancies not available')
        end
        nsparseE=(protCorrFact.*E')';
    else
        nsparseE=E;
    end
    %create enzyme constrains only for proteins measured in all conditions
     nsparseE(~all(nsparseE,2),:) = 0;
end

%create default maximum condition specific uptake rates
unboundNutrExch = nutrExch;
unboundNutrExch{:,:} = cell2mat(arrayfun(@(i)1000*double(sign(nutrExch.(i))),1:size(nutrExch,2),'un',0));

%in case of ecoli only default uptake rates are used

if strcmp(org_name, 'escherichia coli')
    nutrExch=unboundNutrExch;
    warning('For Ecoli no experimental uptade rates are used')
end

%create dummy enzyme concentration matrix for ScorKcatCorrection
E_0=zeros(length(models{1}.enzymes), length(condNames));


%Intialize array of unbound model results
unbError=nan(length(condNames), length(models));
unbMu=nan(length(condNames), length(models));
%Initialize a 3D arrazy to store enzyme usage predictions
unbpredE=nan(length(condNames), length(models), length(enzRxnIdx));
%Initialize array of experimental uptake rate bounds results
bError=nan(length(condNames), length(models));
bMu=nan(length(condNames), length(models));
if enz_mappable
    %intialize array for for corrected proteomic constrains
    pbcorError=nan(length(condNames), length(models));
    pbcorMu=nan(length(condNames), length(models));
    if nargout>3
        fluxVar=nan(length(condNames), length(models), length(models{1}.rxns));
    end
    %intitalize array with experimental + proteomics constrains
    %but without the enzyme pool constrain
    pbnopoolError=nan(length(condNames), length(models));
    pbnopoolMu=nan(length(condNames), length(models));
    %initialize arra for model only constrained by experimental uptake
    %rates
    bnopoolError=nan(length(condNames), length(models));
    bnopoolMu=nan(length(condNames), length(models));
end
for i=1:length(models)
    
    %Simulate for all conditions
    for j=1:length(condNames)
        %if test_unseen is false only test the condition for which the
        %model was build
        if ~test_unseen
            j=i;
        end
        tmpmod=models{i};
        %adjust protein batch constrain
        cd(fullfile(geckoDir, 'geckomat', 'limit_proteins'))
        tmpmod.ub(end)=tmpmod.ub(end)*(P(j)/sumProtein(tmpmod));
        cd(topDir)
        %adjust model for unbounded growth in specific culture conditions;get gam from
        %configuration script
        unbmod = adjBaseModel(tmpmod,P(j),unboundNutrExch(:,j), GAM(j));
        unbmod=unbmod{:};
        %adjust model with experimentally measured uptake rates;
        bmod= adjBaseModel(tmpmod,P(j),nutrExch(:,j), GAM(j));
        bmod=bmod{:};
        %create a model with experimental measured uptakte rates and
        %liftet pool constrain
        bnopoolmod=changeRxnBounds(bmod, bmod.rxns(end), inf, 'u');
        %switch constrain on enzyme usage rates back to 1000
        bnopoolmod.ub(enzRxnIdx)=1000;
        %Calculate prediction error of grwoth rate for all conditions
        %without uptake fluxes
        [unbError(j,i),unbMu(j,i), unbpredE(j,i,:)] = scoreKcatCorrection(unbmod,expVal(j),...
            E_0(:,j),enzRxnIdx, true, {'predE'});
        %Calculate prediction error of grwoth rate for all conditions
        [bError(j,i),bMu(j,i)] = scoreKcatCorrection(bmod,expVal(j),...
            E_0(:,j),enzRxnIdx, true);
        if enz_mappable
            %solve model with relaxed proteomics constrains
            if nargout>3
            [pbcorError(j,i),pbcorMu(j,i), fluxVar(j,i,:)] = scoreKcatCorrection(bmod,expVal(j),...
                nsparseE(enzMetRxnMatch,j),enzRxnIdx, true, {'FVA'});
            else
                [pbcorError(j,i),pbcorMu(j,i)] = scoreKcatCorrection(bmod,expVal(j),...
                nsparseE(enzMetRxnMatch,j),enzRxnIdx, true);
            end
            %solve model with relaxed proteomics constrains only for enzymes 
            %quantified in all conditions and experimental constrains
            %but without protein pool constrain to compare with presto
            %scores
             [pbnopoolError(j,i),pbnopoolMu(j,i)] = scoreKcatCorrection(bnopoolmod,expVal(j),...
                nsparseE(enzMetRxnMatch,j),enzRxnIdx, true);
            %Solve model only constrained by uptake
             [bnopoolError(j,i),bnopoolMu(j,i)] = scoreKcatCorrection(bnopoolmod,expVal(j),...
                E_0(enzMetRxnMatch,j),enzRxnIdx, true);
            
        end
        %if test_unseen is false only test the condition in which model was
        %build
        if ~test_unseen
            break
        end
    end
end
relError={unbError, bError};
Mu={unbMu, bMu};
pred_E=unbpredE;
if enz_mappable
    relError=[relError, {pbcorError, pbnopoolError, bnopoolError}];
    Mu=[Mu, {pbcorMu, pbnopoolMu, bnopoolMu}];
end
if ~test_unseen
    %if test_unseen is false only the diagonal has values
    for i=1:length(relError)
        relError{i}=diag(relError{i});
        Mu{i}=diag(Mu{i});
    end
    for i=1:size(pred_E,3)
        pred_E(:,1,i)=diag(pred_E(:,:,i));
    end
    pred_E=squeeze(pred_E(:,1,:));
    if nargout>3
        for i=1:size(fluxVar,3)
            fluxVar(:,1,i)=diag(fluxVar(:,:,i));
        end
        fluxVar=squeeze(fluxVar(:,1,:));
    end
end
