# Correction of turnover numbers in enzyme-constraint metabolic models

## Publication

## OS
* Code was tested on Ubuntu 20.04.3 LTS and Windows 10

## Dependencies
* Matlab (tested with 2020a/b)
* [COBRA toolbox 3.0](https://github.com/opencobra/cobratoolbox)
* [RAVEN toolbox 2.0](https://github.com/SysBioChalmers/RAVEN)
* [Gurobi solver](https://www.gurobi.com/) (tested with version 9.1.1)
* 
## Setup
* Set up the COBRA toolbox following the [installation instructions](https://opencobra.github.io/cobratoolbox/stable/installation.html)
* Set up the RAVEN toolbox following the [installation instructions](https://github.com/SysBioChalmers/RAVEN/wiki/Installation#installation)
* Set up the Gurobi solver and connect it with the COBRA toolbox using install instructions that can be found [here](https://opencobra.github.io/cobratoolbox/stable/installation.html#gurobi)

## Run k<sub>cat</sub> correction with PRESTO

### To reproduce the results presented in the paper run the following scripts
* _S. cerevisiae_: `correct_kcats_yeast`
* _E. coli_: `correct_kcats_ecoli`

### Generate raw and adapted ecModels
1. `get_rawecmod()` will generate among others raw enzyme constrained models (no kcat adaption, no manual modifications, no protein pool constraint) in the model folder of the respective GECKO path e.g`GECKO_S_cerevisiae/model/ecYeast/rawYeast.mat`
2. call `get_GKOmod()` or `get_GKOmod_ecoli()` to generate all condition specific adapted GECKO models using the experimentally determinded protein content, growth rate and uptake rates where available and also save them to the `model/` folder .

### To apply PRESTO to another model
1. create ecModel using the [GECKO toolbox](https://github.com/SysBioChalmers/GECKO)
2. make adjustments to parameters and input file names in the configuration script

| Parameter | Explanation |
| :---:         | --- |
| _orgName_     | name of the organism |
| _modelFile_   | file of the ecModel as Matlab workspace (.mat) |
| _cobraSolver_ | preferred solver for linear optimization problems (default: gurobi) |
| _runParallel_ | whether the correction should be run on multiple threads (default: true) |
| _ncpu_        | number of threads (default: 20) |
| _epsilon_     | upper limit for fold-change of k<sub>cat</sub> values |
| _lambda_      | weight for the minimization of absolute difference between measured and predicted growth rate(s) |
| _theta_       | upper limit for the difference between measured and predicted growth rate(s) |
| _GAM_         | value of growth-associated maintencance (put NaN if unknown, will be fitted using all provided conditions) |
| _f_           | mass fraction of all proteins accounted for by the model (see [GECKO publication](https://doi.org/10.15252/msb.20167411)) |
| _f\_n_        | mass fraction of unmeasured proteins in the ecModel (for inclusion of unmeasured proteins; put NaN if unknown, will be fitted for each condition separately) (see [GECKO publication](https://doi.org/10.15252/msb.20167411)) |
| _sigma_       | average saturation of enzymes in the model (see [GECKO publication](https://doi.org/10.15252/msb.20167411)); can be fitted using GECKO _sigmaFitter_) |

3. run `cvLambdaFitting` to estimate the optimal weighting parameter $\lambda$
```
[relErr,errVar,sumsDelta,objVal,avJD,corrKcatProts] = cvLambdaFitting(...
    model,...                       % GECKO ecModel
    expGrowth,...                   % experimental growth rates for all conditions
    P,...                           % total protein contents for all conditions
    E,...                           % enzyme abundance matrix (#model proteins x #conditions)
    lambdaParams,...                % array of lambda parameters to be explored
    nutrExch,...                    % nutrient exchange rates
    'epsilon', epsilon,...          % (optional) maximum allowed fold change of kcats
    'theta', theta,...              % (optional) maximum allowed relative error
    'runParallel', runParallel,...  % (optional) whether to run the cross-validation on multiple workers
    'GAM', GAM,...                  % (optional) growth associated maintenance
    'f', f,...                      % f factor for protein pool (see GECKO paper or description above)
    'sigma', sigma...               % sigma factor for protein pool (see GECKO paper or description above)
    );
```

4. adjust ecModel to experimental conditions using `adjBaseModel`
```
adj_models = adjBaseModel(...
    model,...       % GECKO ecModel
    P,...           % total protein contents for all conditions
    nutrExch,...    % nutrient exchange rates
    GAM...          % growth associated maintenance
    );
```

5. run `PRESTO` to obtain k<sub>cat</sub> corrections
```
[solution,corr_models,relError,changeTab,LP] = PRESTO(...
    adj_models,...     % enzyme-constraint metabolic model(s)
    expGrowth,...      % experimental growth rates for all conditions
    E...               % enzyme abundance matrix (#model proteins x #conditions)
    );
```
