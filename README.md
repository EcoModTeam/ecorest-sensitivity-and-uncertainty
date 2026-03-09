# ecorest-sensitivity-and-uncertainty
Sensitivity and uncertainty analysis for habitat suitability index (HSI) models in the ecorest R package

## Purpose of this repository
The code in this repository allows users to conduct a batch sensitivity and uncertainty analysis of the habitat suitability index models (HSI) available in the 'ecorest' package in R (https://CRAN.R-project.org/package=ecorest), both using their original HSI equation structures and changing the structure to an arithmetic mean, geometric mean, limiting factor, or multiplicative structure.

## Files in this repository

### ecorest_SA_UA.R:
This file can be used to conduct a variance-based sensitivity and uncertainty analysis in R using the 'sensobol' package (Puy et al., 2022). The script iteratively runs through all of the models in ecorest and performs a sensitivity and uncertainty analysis for every model with more than one input parameter. The analysis adheres to best practices defined in Puy et al. (2022) and Saltelli et al. (2010). The script runs the analyses using the model equation stored in HSImetadata in ecorest, but also tests each model using an arithmetic mean, geometric mean, limiting factor, and multiplicative structure.

The outputs of this script include RDS files tracking model progress ('progress_') and containing the results of the sensitivity and uncertainty analysis at the model and the parameter level for each model structure ('model_').

Each output file with a name beginning with 'model' contains the following information for every model in ecorest:
- name: (character) the name of the model, as defined in the ecorest 'HSImetadata' file
- runs: (integer) the number of model evaluations conducted during the sensitivity/uncertainty analysis
- inputs: (integer) the number of input parameters in the model (SIVs in ecorest)
- components: (integer) the number of input components in the model (see HSImetadata in ecorest)
- cat_inputs: (integer) the total number of categorical inputs in the model
- n: (double) the base sample size required for model convergence (default is 60,000)
- ind: (list) the results of the sensitivity analysis, including estimates, bias, standard error, and bootstrapped confidence intervals for the first (Si) and total (Ti) order indices of each input parameter for all model formats (original, arithmetic mean, geometric mean, limiting factor, multiplicative) 
- HSI: (double) a vector containing all of the outputs (HSI scores) from every model evaluation during the sensitivity and uncertainty analysis for every model format (original, arithmetic mean, geometric mean, limiting factor, multiplicative) 
- quantiles: (double) a vector containing the 1, 2.5, 5, 25, 50, 75, 95, 97.5, 99, and 100% quantiles obtained from the output HSI score distribution
- dummy: (list) estimates of numerical approximation error, in the form of dummy parameter estimates, following the same format as results for every model format (original, arithmetic mean, geometric mean, limiting factor, multiplicative)
- converged: (logical) a vector that reports whether the model's sensitivity indices converged for every model format (original, arithmetic mean, geometric mean, limiting factor, multiplicative). Convergence was determined based on whether the range of the 95% confidence interval for each sensitivity index was less than or equal to 0.05 (Sarrazin et al., 2016)

### SA_UA_stats.R

The purpose of this script is to assess patterns of sensitivity and uncertainty in ecorest models and compare how these patterns change when the model structure is altered. 

The script requires the following inputs:
- the RDS files with names beginning with 'model' generated using ecorest_SA_UA.R
- HSImetadata.RData (or HSImetadata in the ecorest package)
- HSImodels.RData (or HSImodels in the ecorest package)
- equation_types.csv

### HSImetadata.RData
This is equivalent to HSImetadata from ecorest v2.0.1.

### HSImodels.RData
This is equivalent to HSImodels from ecorest v2.0.1.

### equation_types.csv
This is a csv file containing a list of model structures for the models in ecorest. Models could be categorized as: single value, limiting factor, multiplicative, arithmetic mean, geometric mean, weighted arithmetic mean, weighted geometric mean, or author-specified.

## References
Puy, A., Piano, S. L., Saltelli, A., & Levin, S. A. (2022). sensobol: An R package to compute variance-based sensitivity indices. Journal of Statistical Software, 102(5). https://doi.org/10.18637/jss.v102.i05

Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M., & Tarantola, S. (2010). Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index. Computer Physics Communications, 181(2), 259–270. 

Sarrazin, F., Pianosi, F., & Wagener, T. (2016). Global sensitivity analysis of environmental models: Convergence and validation. Environmental Modelling & Software, 79, 135–152. https://doi.org/10.1016/j.envsoft.2016.02.005
