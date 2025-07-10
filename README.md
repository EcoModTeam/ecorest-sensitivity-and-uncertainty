# ecorest-sensitivity-and-uncertainty
Sensitivity and uncertainty analysis for habitat suitability index (HSI) models in the ecorest R package

## Purpose of this repository
The code in this repository allows users to conduct a batch sensitivity and uncertainty analysis of the habitat suitability index models (HSI) available in the 'ecorest' package in R (https://CRAN.R-project.org/package=ecorest), both using their original HSI equation structures and changing the structure to an arithmetic mean, geometric mean, limiting factor, or multiplicative structure.

## Files in this repository

### ecorest_sensitivity_uncertainty_analysis.R:
This file can be used to conduct a variance-based sensitivity and uncertainty analysis in R using the 'sensobol' package (Puy et al., 2022). The script iteratively runs through all of the models in ecorest and performs a sensitivity and uncertainty analysis for every model with more than one input parameter. The analysis adheres to best practices defined in Puy et al. (2022) and Saltelli et al. (2010). The script runs the analyses using the model equation stored in HSImetadata in ecorest, but also tests each model using an arithmetic mean, geometric mean, limiting factor, and multiplicative structure.

The output of this script are ten RDS files containing the results of the sensitivity and uncertainty analysis at the model and the parameter level for each model structure, stored in an S4 class.

Each output file with a name beginning with 'model' contains the following information for every model in ecorest:
- name: (character) the name of the model, as defined in the ecorest 'HSImetadata' file
- runs: (integer) the number of model evaluations conducted during the sensitivity/uncertainty analysis
- first: (character) the estimator used to calculate the first  order Sobol' index in sensobol. Default is 'saltelli'
- total: (character) the estimator used to calculate the total order Sobol' index in sensobol. Default is 'jansen'
- inputs: (integer) the number of input parameters in the model (SIVs in ecorest)
- components: (integer) the number of input components in the model (see HSImetadata in ecorest)
- cat_inputs: (double) the total number of categorical inputs in the model
- breakpoints: (list) a data frame containing the data from 'HSImodels' in ecorest, used to convert input parameters to SIV values
- HSI: (double) a vector containing all of the outputs (HSI scores) from every model evaluation during the sensitivity and uncertainty analysis
- quantiles: (double) a vector containing the 1, 2.5, 5, 25, 50, 75, 95, 97.5, 99, and 100% quantiles obtained from the output HSI score distribution
- sum: (double) the total amount of uncertainty in HSI scores that could be accounted for by all first order Sobol' indices
- results: (list) results of the sensitivity analysis, including estimates, bias, standard error, and bootstrapped confidence intervals for the first (Si) and total (Ti) order indices of each input parameter
- dummy: (list) estimates of numerical approximation error, in the form of dummy parameter estimates, following the same format as results

Each output file with a name beginning with 'params' contains the following information for every parameter in each ecorest model:
  - model: (character) the name of the model in ecorest that the parameter belongs to
  - name: (character) the name of the parameter
  - type: (character) the type of parameter (numeric or categorical)
  - num_breaks: (integer) the total number of breakpoints in the parameter's SIV graph
  - Si: (double) the first order Sobol' index estimate for the parameter
  - Si_std_error: (double) the standard error of the first order Sobol' index estimate for the parameter
  - Si_low_ci: (double) the lower bound of the 95% confidence interval for the first order Sobol' index estimate, obtained by bootstrapping 10^3 times
  - Si_high_ci: (double) the upper bound of the 95% confidence interval for the first order Sobol' index estimate, obtained by bootstrapping 10^3 times
  - Ti: (double) the total order Sobol' index estimate for the parameter
  - Ti_std_error: (double) the standard error of the total order Sobol' index estimate for the parameter
  - Ti_low_ci: (double) the lower bound of the 95% confidence interval for the total order Sobol' index estimate, obtained by bootstrapping 10^3 times
  - Ti_high_ci: (double) the upper bound of the 95% confidence interval for the total order Sobol' index estimate, obtained by bootstrapping 10^3 times
  - influential: (double) value indicating whether the parameter is influential (1) or not (0) based on whether the confidence intervals of both the Si and Ti estimates overlap their corresponding dummy parameter confidence intervals, indicating that the parameter's influence on uncertainty in HSI scores is indistinguishable from numerical approximation error

### ecorest_statistics_sensitivity_uncertainty_analysis.R

The purpose of this script is to assess patterns of sensitivity and uncertainty in ecorest models and compare how these patterns change when the model structure is altered. 

The script requires the following inputs:
- the 10 RDS files generated using ecorest_sensitivity_uncertainty_analysis.R
- HSImetadata.RData (or HSImetadata in the ecorest package)
- HSImodels.RData (or HSImodels in the ecorest package)
- equation_types.csv

### HSImetadata.RData
This is an updated version of the HSImetadata available in ecorest v2.0.0.

### HSImodels.RData
This is an updated version of HSImodels available in ecorest v2.0.0.

### equation_types.csv
This is a csv file containing a list of model structures for the models in ecorest. Models could be categorized as: single value, limiting factor, multiplicative, arithmetic mean, geometric mean, weighted arithmetic mean, weighted geometric mean, or author-specified.

## References
Puy, A., Piano, S. L., Saltelli, A., & Levin, S. A. (2022). sensobol: An R package to compute variance-based sensitivity indices. Journal of Statistical Software, 102(5). https://doi.org/10.18637/jss.v102.i05

Saltelli, A., Annoni, P., Azzini, I., Campolongo, F., Ratto, M., & Tarantola, S. (2010). Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index. Computer Physics Communications, 181(2), 259–270. 
