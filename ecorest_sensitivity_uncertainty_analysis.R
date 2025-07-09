# Batch code for running sensitivity analysis using ecorest models
# Authors: Kiara Cushway, COlton Shaw, Kyle McKay, Todd Swannack

## Clear workspace
rm(list = ls()) ## remove stored files and objects
gc(T) ## garbage collection
graphics.off() ## Turn off graphics

## Install packages if needed
### install.packages("ecorest")
### install.packages("data.table")
### install.packages("ggplot2")
### install.packages("sensobol")
### install.packages("tidyr")

## Load libraries
library(here) ## used for file paths
library(ecorest) ## contains HSI models and function to calculate HSI
library(sensobol) ## used to run sensitivity analysis and plot results
library(ggplot2) ## used for plotting
library(stringr)## used to manipulate strings
library(extraDistr) ## used to generate discrete uniform distribution
library(tidyr) ## used for data manipulation
library(dplyr) ## used to manage model inclusion
library(patchwork) ## used to combine plots

load("HSImetadata.RData") ## load updated HSImetadata
load("HSImodels.RData") ## load updated HSImodels

################################################################################

## Define function to calculate multiplicative HSI
HSImult <- function(x){
  HSI <- min(1, prod(x, na.rm=TRUE))
  
  if(HSI< 0 | HSI > 1){
    HSIout <- "Habitat suitability index not within 0 to 1 range."
  } else {
    HSIout <- HSI
  }
  
  # Return HSI outcome
  return(HSIout)
}

## Define function to calculate limiting factor HSI
HSIlimit <- function(x){
  HSI <- min(x, na.rm=TRUE)
  
  if(HSI< 0 | HSI > 1){
    HSIout <- "Habitat suitability index not within 0 to 1 range."
  } else {
    HSIout <- HSI
  }
  
  # Return HSI outcome
  return(HSIout)
}

################################################################################

## Parameters needed for sensitivity analysis (using best practices from literature)
M = NULL ## number of input variables, will vary with model used
n = 10000 ## base sample size; n should be 500 or higher (Saltelli et al. 2010)
N = n*(M+2) ## appropriate number of model evaluations for a given model
params = NULL ## parameters in equation (inputs of HSIeqtn())
matrices = c("A", "B", "AB") ## defines matrix to outline sampling space (columns= model inputs; rows=sampling points)
first = "saltelli" ## best practice first order estimator (Puy et al. 2022b)
total = "jansen" ## best practice total order estimator (Puy et al. 2022b)
order = "first" ## type of estimator you are interested in
R = 10^3 ## number of bootstrap replicas
type = "QRN" ## type of sampling scheme; Quasi-Random (QRN) is best when functions have unknown behavior and tends to have smaller unexplored volumes (Puy et al. 2022b)
conf = 0.95 ## size of confidence interval

cat_vars = list() ## list to store categorical variables in
not_na = function(x) sum(!is.na(x)) ## function to determine whether a value is NA
num_cat = 0 ## initial number of categorical variables
cat_names = c() ## initial list of categorical names

################################################################################

## Set seed for reproducible results
set.seed(349)

## Load metadata into data frame, excluding the two models with only one variable
HSImetadata = HSImetadata %>% filter(!model %in% c("redwingedblackbirdB", "woodduckWinter"))

## Remove the two models with one variable from HSImodels as well
HSImodels$redwingedblackbirdB = NULL
HSImodels$woodduckWinter = NULL

## Create a class "Model_Info" where we will store the model results from the SA
setClass("Model_Info", slots = list(name = "character",  ## name of model
                                  runs = "numeric",    ## number of runs
                                  first = "character", ## first order estimator
                                  total = "character", ## total order estimator
                                  inputs = "numeric",  ## number of inputs
                                  components = "numeric", ## number of components
                                  cat_inputs = "numeric", ## number of categorical inputs
                                  breakpoints = "data.frame", ## SIV data frame
                                  HSI = "numeric", ## HSI output
                                  quantiles = "numeric", ## uncertainty quantiles
                                  sum = "numeric", ## sum of first order indices; when 1, model is additive with no interactions
                                  results = "list", ## results from sensitivity analysis
                                  dummy = "list" ## dummy parameter to estimate approximation error
                                  ))

## Create a class "Model_Info_Arithmetic" where we will store the model results from the SA using arithmetic means
setClass("Model_Info_Arithmetic", slots = list(name = "character",  ## name of model
                                  runs = "numeric",    ## number of runs
                                  first = "character", ## first order estimator
                                  total = "character", ## total order estimator
                                  inputs = "numeric",  ## number of inputs
                                  components = "numeric", ## number of components
                                  cat_inputs = "numeric", ## number of categorical inputs
                                  breakpoints = "data.frame", ## SIV data frame
                                  HSI = "numeric", ## HSI output
                                  quantiles = "numeric", ## uncertainty quantiles
                                  sum = "numeric", ## sum of first order indices; when 1, model is additive with no interactions
                                  results = "list", ## results from sensitivity analysis
                                  dummy = "list" ## dummy parameter to estimate approximation error
))

## Create a class "Model_Info_Geometric" where we will store the model results from the SA using geometric means
setClass("Model_Info_Geometric", slots = list(name = "character",  ## name of model
                                  runs = "numeric",    ## number of runs
                                  first = "character", ## first order estimator
                                  total = "character", ## total order estimator
                                  inputs = "numeric",  ## number of inputs
                                  components = "numeric", ## number of components
                                  cat_inputs = "numeric", ## number of categorical inputs
                                  breakpoints = "data.frame", ## SIV data frame
                                  HSI = "numeric", ## average HSI output
                                  quantiles = "numeric", ## uncertainty quantiles
                                  sum = "numeric", ## sum of first order indices; when 1, model is additive with no interactions
                                  results = "list", ## results from sensitivity analysis
                                  dummy = "list" ## dummy parameter to estimate approximation error
))

## Create a class "Model_Info_Multiplicative" where we will store the model results from the SA using geometric means
setClass("Model_Info_Multiplicative", slots = list(name = "character",  ## name of model
                                              runs = "numeric",    ## number of runs
                                              first = "character", ## first order estimator
                                              total = "character", ## total order estimator
                                              inputs = "numeric",  ## number of inputs
                                              components = "numeric", ## number of components
                                              cat_inputs = "numeric", ## number of categorical inputs
                                              breakpoints = "data.frame", ## SIV data frame
                                              HSI = "numeric", ## average HSI output
                                              quantiles = "numeric", ## uncertainty quantiles
                                              sum = "numeric", ## sum of first order indices; when 1, model is additive with no interactions
                                              results = "list", ## results from sensitivity analysis
                                              dummy = "list" ## dummy parameter to estimate approximation error
))

## Create a class "Model_Info_Limiting" where we will store the model results from the SA using geometric means
setClass("Model_Info_Limiting", slots = list(name = "character",  ## name of model
                                              runs = "numeric",    ## number of runs
                                              first = "character", ## first order estimator
                                              total = "character", ## total order estimator
                                              inputs = "numeric",  ## number of inputs
                                              components = "numeric", ## number of components
                                              cat_inputs = "numeric", ## number of categorical inputs
                                              breakpoints = "data.frame", ## SIV data frame
                                              HSI = "numeric", ## average HSI output
                                              quantiles = "numeric", ## uncertainty quantiles
                                              sum = "numeric", ## sum of first order indices; when 1, model is additive with no interactions
                                              results = "list", ## results from sensitivity analysis
                                              dummy = "list" ## dummy parameter to estimate approximation error
))

## Create a class "Parameter_Info" where we will store the parameter results from the SA
setClass("Parameter_Info", slots=list(model = "character", ## name of model
                                  name = "character",  ## name of parameter
                                  type = "character",  ## categorical or numeric
                                  num_breaks = "numeric", ## number of breakpoints
                                  Si = "numeric",  ## first order index
                                  Si_std_error = "numeric", ## first order standard error
                                  Si_low_ci = "numeric", ## first order low confidence interval
                                  Si_high_ci = "numeric", ## first order high confidence interval
                                  Ti = "numeric", ## total order index
                                  Ti_std_error = "numeric", ## total order standard error
                                  Ti_low_ci = "numeric", ## total order low confidence interval
                                  Ti_high_ci = "numeric", ## total order high confidence interval
                                  influential = "numeric" ## whether both Si and Ti low CI fall below dummy parameter high CI (i.e., not influential; 1= influential, 0= uninfluential)
                                  ))

## Create a class "Parameter_Info_Arithmetic" where we will store the parameter results from the SA using arithmetic means
setClass("Parameter_Info_Arithmetic", slots=list(model = "character", ## name of model
                                  name = "character",  ## name of parameter
                                  type = "character",  ## categorical or numeric
                                  num_breaks = "numeric", ## number of breakpoints
                                  Si = "numeric",  ## first order index
                                  Si_std_error = "numeric", ## first order standard error
                                  Si_low_ci = "numeric", ## first order low confidence interval
                                  Si_high_ci = "numeric", ## first order high confidence interval
                                  Ti = "numeric", ## total order index
                                  Ti_std_error = "numeric", ## total order standard error
                                  Ti_low_ci = "numeric", ## total order low confidence interval
                                  Ti_high_ci = "numeric", ## total order high confidence interval
                                  influential = "numeric" ## whether both Si and Ti low CI fall below dummy parameter high CI (i.e., not influential; 1= influential, 0= uninfluential)
))

## Create a class "Parameter_Info_Geometric" where we will store the parameter results from the SA using geometric means
setClass("Parameter_Info_Geometric", slots=list(model = "character", ## name of model
                                  name = "character",  ## name of parameter
                                  type = "character",  ## categorical or numeric
                                  num_breaks = "numeric", ## number of breakpoints
                                  Si = "numeric",  ## first order index
                                  Si_std_error = "numeric", ## first order standard error
                                  Si_low_ci = "numeric", ## first order low confidence interval
                                  Si_high_ci = "numeric", ## first order high confidence interval
                                  Ti = "numeric", ## total order index
                                  Ti_std_error = "numeric", ## total order standard error
                                  Ti_low_ci = "numeric", ## total order low confidence interval
                                  Ti_high_ci = "numeric", ## total order high confidence interval
                                  influential = "numeric" ## whether both Si and Ti low CI fall below dummy parameter high CI (i.e., not influential; 1= influential, 0= uninfluential)
))

## Create a class "Parameter_Info_Multiplicative" where we will store the parameter results from the SA using geometric means
setClass("Parameter_Info_Multiplicative", slots=list(model = "character", ## name of model
                                                name = "character",  ## name of parameter
                                                type = "character",  ## categorical or numeric
                                                num_breaks = "numeric", ## number of breakpoints
                                                Si = "numeric",  ## first order index
                                                Si_std_error = "numeric", ## first order standard error
                                                Si_low_ci = "numeric", ## first order low confidence interval
                                                Si_high_ci = "numeric", ## first order high confidence interval
                                                Ti = "numeric", ## total order index
                                                Ti_std_error = "numeric", ## total order standard error
                                                Ti_low_ci = "numeric", ## total order low confidence interval
                                                Ti_high_ci = "numeric", ## total order high confidence interval
                                                influential = "numeric" ## whether both Si and Ti low CI fall below dummy parameter high CI (i.e., not influential; 1= influential, 0= uninfluential)
))

## Create a class "Parameter_Info_Limiting" where we will store the parameter results from the SA using geometric means
setClass("Parameter_Info_Limiting", slots=list(model = "character", ## name of model
                                                name = "character",  ## name of parameter
                                                type = "character",  ## categorical or numeric
                                                num_breaks = "numeric", ## number of breakpoints
                                                Si = "numeric",  ## first order index
                                                Si_std_error = "numeric", ## first order standard error
                                                Si_low_ci = "numeric", ## first order low confidence interval
                                                Si_high_ci = "numeric", ## first order high confidence interval
                                                Ti = "numeric", ## total order index
                                                Ti_std_error = "numeric", ## total order standard error
                                                Ti_low_ci = "numeric", ## total order low confidence interval
                                                Ti_high_ci = "numeric", ## total order high confidence interval
                                                influential = "numeric" ## whether both Si and Ti low CI fall below dummy parameter high CI (i.e., not influential; 1= influential, 0= uninfluential)
))

## Create an empty object for each class of models
model_object = c()
model_object_arithmetic = c()
model_object_geometric = c()
model_object_multiplicative = c()
model_object_limiting = c()

## Create an empty object for parameters
parameter_object = c()
parameter_object_arithmetic = c()
parameter_object_geometric = c()
parameter_object_multiplicative = c()
parameter_object_limiting = c()

################################################################################

## Run sensitivity analysis and populate model class with information from models

## Pull pertinent data from HSImetadata and HSImodels
for(i in 1:nrow(HSImetadata)){
  breakpoints = HSImodels[[paste(HSImetadata[i,1])]] ## grab breakpoints of model
  Comp_exist = HSImetadata[i,41:(length(HSImetadata)-1)] ## store components from component columns in HSImetadata
  comp = Comp_exist[,colSums(is.na(Comp_exist)) == 0] ## subset components that aren't NA
  SIV_exist = HSImetadata[i,9:40] ## store SIV data from all SIV cols HSImetadata
  siv = SIV_exist[,colSums(is.na(SIV_exist)) == 0] ## subset SIVs that aren't NA
  params = unlist(siv,use.names = FALSE) ## Set parameters for sensitivity analysis equal to SIV names
  M = length(params) ## set number of input variables equal to number of SIVs in model
  N = n*(M+2) ## set number of model evaluations for model based on M and n
  breakpoints = breakpoints[,colSums(is.na(breakpoints)) < nrow(breakpoints)] ## select only columns with breakpoint data
  mat = sobol_matrices(N = n,params = params,order = order,type = type) ## create matrix of probabilities to compute Sobol' first and total order indices
  for (col in seq(from = 1,to = length(breakpoints),by = 2)) {
    if (is.factor(breakpoints[1,col]) == TRUE) {
      has_cat = TRUE ## set to true if parameter is categorical
      num_cat = num_cat + 1 ## add 1 to count of categorical variables in model
      cat_names = append(cat_names,colnames(breakpoints[col + 1])) ## create list of categorical variables
      }
    }
  
  ## Identify categorical versus continuous variables and store sample data in mat
  for(s in 1:length(params)){
    ## if parameter is categorical, sample from discrete uniform distribution
    if (!is.numeric(breakpoints[,s*2-1])) {
      name = colnames(breakpoints[2*s-1]) ## store column name
      cat_vars[name] = not_na(breakpoints[2*s-1]) ## count number of non-NA breakpoints
      length = 1:get(name,cat_vars) ## create list of numbers equal to number of categories
      ## convert mat probabilities into sample values using discrete uniform distribution
      mat[,s] = qdunif(mat[,s],min=min(length),max=max(length))
    } else {
        if (breakpoints[1,s*2] == 0 && breakpoints[2,s*2] == 0) {
          ## Store the minimum and maximum values within the non-zero range of the SIV
          is_zero = na.omit(breakpoints[,s*2]) == 0 # test whether each value is zero
          keep_min <- is_zero & (c(FALSE, is_zero[-length(is_zero)])) # for min: keep second zero of adjacent pairs
          # Find minimum value
          min = min(breakpoints[,s*2-1][keep_min], na.rm = T)
        } else {
          min = min(breakpoints[,s*2-1],na.rm = TRUE) ## store min value for each SIV
        }
        if (breakpoints[length(na.omit(breakpoints[,s*2])),s*2] == 0 && breakpoints[(length(na.omit(breakpoints[,s*2])) - 1),s*2] == 0) {
          ## Store the minimum and maximum values within the non-zero range of the SIV
          is_zero = na.omit(breakpoints[,s*2]) == 0 # test whether each value is zero
          keep_max <- is_zero & (c(is_zero[-length(is_zero)], FALSE)) # for max: keep first zero of adjacent pairs
          # Find  maximum value
          max = max(breakpoints[,s*2-1][keep_max], na.rm = T)
        } else {
          max = max(breakpoints[,s*2-1],na.rm = TRUE) ## store max value for each SIV
        }
        ## If a portion of the data is zero, select only the non-zero range
        mat[, s] = qunif(mat[, s], min = min, max = max) ## convert probabilities in mat to SIV values between min and max using uniform distribution
      }
  }
  
  ## Create a matrix to store SIV calculations for each variable in the model
  SIV_scores = matrix(rep(NA,length(mat[,1])),nrow = length(mat[,1]),ncol = ncol(breakpoints)/2)
  colnames(SIV_scores) = params
  
  ## Create a matrix to store final HSI scores
  HSI = matrix(rep(NA,length(mat[,1])),nrow = length(mat[,1]),ncol = 5)
  colnames(HSI) = c("HSI","HSI_Arithmetic","HSI_Geometric", "HSI_Multiplicative", "HSI_Limiting")
  
  ## Replace mat contents with letters for categorical variables
  for (c in 1:ncol(mat)) {
    if (colnames(mat)[c] %in% cat_names) {
      mat[,c] = letters[as.numeric(mat[,c])]
    }
  }
  
  ## Calculate SIV and HSI scores for each quasi-random combination of values from mat
  for(ii in 1:length(mat[,1])){
    SIV_scores[ii,] = SIcalc(breakpoints, cbind(mat[ii,1:ncol(mat)])) ## calculate SIV score
    stopifnot((SIV_scores[ii,] <= 1)) ## throw error if SIV scores are not 1 or below
    HSI[ii,1] = HSIeqtn(paste(HSImetadata[i,1]),cbind(SIV_scores[ii,1:ncol(SIV_scores)]), HSImetadata, exclude = NULL) ## calculate final HSI
    HSI[ii,2] = HSIarimean(cbind(SIV_scores[ii,1:ncol(SIV_scores)])) ## calculate final arithmetic HSI
    HSI[ii,3] = HSIgeomean(cbind(SIV_scores[ii,1:ncol(SIV_scores)])) ## calculate final geometric HSI
    HSI[ii,4] = HSImult(cbind(SIV_scores[ii,1:ncol(SIV_scores)])) ## calculate final multiplicative HSI
    HSI[ii,5] = HSIlimit(cbind(SIV_scores[ii,1:ncol(SIV_scores)])) ## calculate final limiting factor HSI
    stopifnot((HSI[ii,] <= 1)) ## throw error if HSI values are not 1 or below
  }
  
  ## Plot uncertainty in HSI values
  uncertainty_graphs = plot_uncertainty(Y = as.vector(HSI[,1]), N = N) + 
    labs(y = "Frequency", x = "HSI score", title = "Original") + plot_uncertainty(Y = as.vector(HSI[,2]), N = N) + 
    labs(y = "Frequency", x = "HSI score", title = "Arithmetic mean") + plot_uncertainty(Y = as.vector(HSI[,3]), N = N) + 
    labs(y = "Frequency", x = "HSI score", title = "Geometric mean")  + plot_uncertainty(Y = as.vector(HSI[,4]), N = N) + 
    labs(y = "Frequency", x = "HSI score", title = "Multiplicative")  + plot_uncertainty(Y = as.vector(HSI[,5]), N = N) + 
    labs(y = "Frequency", x = "HSI score", title = "Limiting factor") + plot_layout(ncol = 1)
  
  print(uncertainty_graphs)
  
  ## Get quantiles from uncertainty plot
  quantiles = quantile(HSI[,1], probs = c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99, 1))
  quantiles_ari = quantile(HSI[,2], probs = c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99, 1))
  quantiles_geo = quantile(HSI[,3], probs = c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99, 1))
  quantiles_mult = quantile(HSI[,4], probs = c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99, 1))
  quantiles_limit = quantile(HSI[,5], probs = c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99, 1))
  
  ## Store parameter names
  params = lapply(params,as.character) ## make SIV names characters
  call_params = do.call(rbind.data.frame, params) ## paste SIV names in dataframe
  
  ## Compute Sobol' indices for original, arithemtic, and geometric models
  ## Original results
  ind = sobol_indices(Y = as.vector(HSI[,1]), N = n, params = call_params[,1], 
                      first = first, total = total, boot = TRUE, 
                      R = R,type = "percent",conf = conf)
  cols = colnames(ind$results)[1:5] ## column names
  ind$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
  ind
  
  ## Arithmetic results
  ind_ari = sobol_indices(Y = as.vector(HSI[,2]), N = n, params = call_params[,1], 
                          first = first, total = total, boot = TRUE, 
                      R = R,type = "percent",conf = conf)
  cols = colnames(ind_ari$results)[1:5] ## column names
  ind_ari$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
  ind_ari
  
  ## Geometric results
  ind_geo = sobol_indices(Y = as.vector(HSI[,3]), N = n, params = call_params[,1], 
                          first = first, total = total, boot = TRUE, 
                      R = R,type = "percent",conf = conf)
  cols = colnames(ind_geo$results)[1:5] ## column names
  ind_geo$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
  ind_geo
  
  ## Multiplicative results
  ind_mult = sobol_indices(Y = as.vector(HSI[,4]), N = n, params = call_params[,1],
                           first = first, total = total, boot = TRUE, 
                          R = R,type = "percent",conf = conf)
  cols = colnames(ind_mult$results)[1:5] ## column names
  ind_mult$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
  ind_mult
  
  ## Limiting factor results
  ind_limit = sobol_indices(Y = as.vector(HSI[,5]), N = n, params = call_params[,1],
                            first = first, total = total, boot = TRUE, 
                          R = R,type = "percent",conf = conf)
  cols = colnames(ind_limit$results)[1:5] ## column names
  ind_limit$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
  ind_limit
  
  ## Calculate dummy parameter for model (estimates approximation error)
  ## Only parameters whose lower confidence intervals are not below the Si and Ti indices of the dummy parameter can be considered truly influential
  ind.dummy = sobol_dummy(Y = as.vector(HSI[,1]),N = n,params = params,boot = TRUE,R = R)
  ind.dummy.ari = sobol_dummy(Y = as.vector(HSI[,2]),N = n,params = params,boot = TRUE,R = R)
  ind.dummy.geo = sobol_dummy(Y = as.vector(HSI[,3]),N = n,params = params,boot = TRUE,R = R)
  ind.dummy.mult = sobol_dummy(Y = as.vector(HSI[,4]),N = n,params = params,boot = TRUE,R = R)
  ind.dummy.limit = sobol_dummy(Y = as.vector(HSI[,5]),N = n,params = params,boot = TRUE,R = R)
  
  ## Graph Sobol' indices and dummy parameter
  sobol_graphs = plot(ind, dummy = ind.dummy) + labs(title = "Original") + 
    plot(ind_ari,dummy = ind.dummy.ari) + labs(title = "Arithmetic mean") + 
    plot(ind_geo,dummy = ind.dummy.geo) + labs(title = "Geometric mean")  + 
    plot(ind_mult, dummy = ind.dummy.mult) + labs(title = "Multiplicative")  + 
    plot(ind_limit,dummy = ind.dummy.limit) + labs(title = "Limiting factor") + plot_layout(ncol = 1)
  
  print(sobol_graphs)
  
  ## Store original model info
  current_model_object = new("Model_Info", 
                              name = paste(HSImetadata[i,1]),
                              runs = ind$C,
                              first = ind$first,
                              total = ind$total,
                              inputs = M,
                              components = length(comp),
                              cat_inputs = num_cat,
                              breakpoints = HSImodels[[paste(HSImetadata[i,1])]],
                              HSI = HSI[,1],
                              quantiles = quantiles,
                              sum = ind$si.sum,
                              results = ind$results,
                              dummy = ind.dummy)
  
  ## Store arithmetic model info
  current_model_object_ari = new("Model_Info_Arithmetic", 
                             name = paste(HSImetadata[i,1]),
                             runs = ind_ari$C,
                             first = ind_ari$first,
                             total = ind_ari$total,
                             inputs = M,
                             components = length(comp),
                             cat_inputs = num_cat,
                             breakpoints = HSImodels[[paste(HSImetadata[i,1])]],
                             HSI = HSI[,2],
                             quantiles = quantiles_ari,
                             sum = ind_ari$si.sum,
                             results = ind_ari$results,
                             dummy = ind.dummy.ari)
  
  ## Store geometric model info
  current_model_object_geo = new("Model_Info_Geometric", 
                             name = paste(HSImetadata[i,1]),
                             runs = ind_geo$C,
                             first = ind_geo$first,
                             total = ind_geo$total,
                             inputs = M,
                             components = length(comp),
                             cat_inputs = num_cat,
                             breakpoints = HSImodels[[paste(HSImetadata[i,1])]],
                             HSI = HSI[,3],
                             quantiles = quantiles_geo,
                             sum = ind_geo$si.sum,
                             results = ind_geo$results,
                             dummy = ind.dummy.geo)
  
  ## Store multiplicative model info
  current_model_object_mult = new("Model_Info_Multiplicative", 
                                 name = paste(HSImetadata[i,1]),
                                 runs = ind_mult$C,
                                 first = ind_mult$first,
                                 total = ind_mult$total,
                                 inputs = M,
                                 components = length(comp),
                                 cat_inputs = num_cat,
                                 breakpoints = HSImodels[[paste(HSImetadata[i,1])]],
                                 HSI = HSI[,4],
                                 quantiles = quantiles_mult,
                                 sum = ind_mult$si.sum,
                                 results = ind_mult$results,
                                 dummy = ind.dummy.mult)
  
  ## Store geometric model info
  current_model_object_limit = new("Model_Info_Limiting", 
                                 name = paste(HSImetadata[i,1]),
                                 runs = ind_limit$C,
                                 first = ind_limit$first,
                                 total = ind_limit$total,
                                 inputs = M,
                                 components = length(comp),
                                 cat_inputs = num_cat,
                                 breakpoints = HSImodels[[paste(HSImetadata[i,1])]],
                                 HSI = HSI[,5],
                                 quantiles = quantiles_limit,
                                 sum = ind_limit$si.sum,
                                 results = ind_limit$results,
                                 dummy = ind.dummy.limit)
  

  ## Append current model info to model objects
  model_object = append(model_object, current_model_object)
  model_object_arithmetic = append(model_object_arithmetic, current_model_object_ari)
  model_object_geometric = append(model_object_geometric, current_model_object_geo)
  model_object_multiplicative = append(model_object_multiplicative, current_model_object_mult)
  model_object_limiting = append(model_object_limiting, current_model_object_limit)
  
  ## Reset loop variables
  cat_names = c() ## reset cat_names to be empty
  num_cat = 0 ## reset num_cat to be 0
  print(i) ## track loop progress
  }


## Name components of list based on model name
names(model_object)=HSImetadata[,1]
names(model_object_arithmetic)=HSImetadata[,1]
names(model_object_geometric)=HSImetadata[,1]
names(model_object_multiplicative)=HSImetadata[,1]
names(model_object_limiting)=HSImetadata[,1]

################################################################################

## Populate parameter class with information about parameters
for (model in 1:length(model_object)) {
  mod = model_object[model] ## grab model from list
  mod_ari = model_object_arithmetic[model] ## grab arithmetic version
  mod_geo = model_object_geometric[model] ## grab geometric version
  mod_mult = model_object_multiplicative[model] ## grab multiplicative version
  mod_limit = model_object_limiting[model] ## grab limiting factor version
  model_name = names(mod[1]) ## store model name
  breakpoints = data.frame(get(model_name,mod)@breakpoints) ## store breakpoints
  breakpoints = breakpoints[,colSums(is.na(breakpoints)) < nrow(breakpoints)] ## select only columns with breakpoint data
  SIV_exist = HSImetadata[model,9:40] ## store SIV data from all SIV cols HSImetadata
  siv = SIV_exist[,colSums(is.na(SIV_exist)) == 0] ## store SIVs that aren't NA
  params = unlist(siv,use.names = FALSE) ## Store parameter names
  results = data.frame(get(model_name,mod)@results) ## Store results of sensitivity analysis
  dummy = data.frame(get(model_name,mod)@dummy) ## Store dummy variables for model Si and Ti
  results_ari = data.frame(get(model_name,mod_ari)@results) ## Store results of sensitivity analysis
  dummy_ari = data.frame(get(model_name,mod_ari)@dummy) ## Store dummy variables for model Si and Ti
  results_geo = data.frame(get(model_name,mod_geo)@results) ## Store results of sensitivity analysis
  dummy_geo = data.frame(get(model_name,mod_geo)@dummy) ## Store dummy variables for model Si and Ti
  results_mult = data.frame(get(model_name,mod_mult)@results) ## Store results of sensitivity analysis
  dummy_mult = data.frame(get(model_name,mod_mult)@dummy) ## Store dummy variables for model Si and Ti
  results_limit = data.frame(get(model_name,mod_limit)@results) ## Store results of sensitivity analysis
  dummy_limit = data.frame(get(model_name,mod_limit)@dummy) ## Store dummy variables for model Si and Ti
  
  
  
  for (param in 1:length(params)) {
    name = params[param]
    ## If parameter is categorical, set type to "categorical"; else set to numeric
    if (is.factor(breakpoints[1,which(colnames(breakpoints) == params[param])-1])) {
      type = "categorical" ## set type as categorical
      num_breaks = not_na(breakpoints[,which(colnames(breakpoints) == params[param])]) ## count non-NA breakpoints
    } else {
      type = "numeric" ## set type as numeric
      num_breaks = not_na(breakpoints[,which(colnames(breakpoints) == params[param])]) ## count non-NA breakpoints
    }
    ## Test whether BOTH Si and Ti low CIs are indistinguishable from approximation error (i.e., dummy parameter) in original model
    if ((results$low.ci[which(results$parameters == params[param] & results$sensitivity == "Si")] <= dummy$high.ci[1]) & 
        (results$low.ci[which(results$parameters == params[param] & results$sensitivity == "Ti")]<= dummy$high.ci[2])) {
      influential = 0
    } else {
      influential = 1
    }
    ## Test whether BOTH Si and Ti low CIs are indistinguishable from approximation error (i.e., dummy parameter) in arithmetic model
    if ((results_ari$low.ci[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Si")] <= dummy_ari$high.ci[1]) & 
        (results_ari$low.ci[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Ti")]<= dummy_ari$high.ci[2])) {
      influential_ari = 0
    } else {
      influential_ari = 1
    }
    ## Test whether BOTH Si and Ti low CIs are indistinguishable from approximation error (i.e., dummy parameter) in geometric model
    if ((results_geo$low.ci[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Si")] <= dummy_geo$high.ci[1]) & 
        (results_geo$low.ci[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Ti")]<= dummy_geo$high.ci[2])) {
      influential_geo = 0
    } else {
      influential_geo = 1
    }
    ## Test whether BOTH Si and Ti low CIs are indistinguishable from approximation error (i.e., dummy parameter) in multiplicative model
    if ((results_mult$low.ci[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Si")] <= dummy_mult$high.ci[1]) & 
        (results_mult$low.ci[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Ti")]<= dummy_mult$high.ci[2])) {
      influential_mult = 0
    } else {
      influential_mult = 1
    }
    ## Test whether BOTH Si and Ti low CIs are indistinguishable from approximation error (i.e., dummy parameter) in limiting factor model
    if ((results_limit$low.ci[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Si")] <= dummy_limit$high.ci[1]) & 
        (results_limit$low.ci[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Ti")]<= dummy_limit$high.ci[2])) {
      influential_limit = 0
    } else {
      influential_limit = 1
    }
    
    ## Store parameter info for original models
    current_parameter_object = new("Parameter_Info", 
                                   model = model_name,
                                   name = name,
                                   type = type,  
                                   num_breaks = num_breaks, 
                                   Si = results$original[which(results$parameters == params[param] & results$sensitivity == "Si")],  
                                   Si_std_error = results$std.error[which(results$parameters == params[param] & results$sensitivity == "Si")], 
                                   Si_low_ci = results$low.ci[which(results$parameters == params[param] & results$sensitivity == "Si")], 
                                   Si_high_ci = results$high.ci[which(results$parameters == params[param] & results$sensitivity == "Si")], 
                                   Ti = results$original[which(results$parameters == params[param] & results$sensitivity == "Ti")], 
                                   Ti_std_error = results$std.error[which(results$parameters == params[param] & results$sensitivity == "Ti")], 
                                   Ti_low_ci = results$low.ci[which(results$parameters == params[param] & results$sensitivity == "Ti")], 
                                   Ti_high_ci = results$high.ci[which(results$parameters == params[param] & results$sensitivity == "Ti")],
                                   influential = influential)
    
    ## Store parameter info for arithmetic models
    current_parameter_object_ari = new("Parameter_Info_Arithmetic", 
                                   model = model_name,
                                   name = name,
                                   type = type,  
                                   num_breaks = num_breaks, 
                                   Si = results_ari$original[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Si")],  
                                   Si_std_error = results_ari$std.error[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Si")], 
                                   Si_low_ci = results_ari$low.ci[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Si")], 
                                   Si_high_ci = results_ari$high.ci[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Si")], 
                                   Ti = results_ari$original[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Ti")], 
                                   Ti_std_error = results_ari$std.error[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Ti")], 
                                   Ti_low_ci = results_ari$low.ci[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Ti")], 
                                   Ti_high_ci = results_ari$high.ci[which(results_ari$parameters == params[param] & results_ari$sensitivity == "Ti")],
                                   influential = influential_ari)
    
    ## Store parameter info for geometric models
    current_parameter_object_geo = new("Parameter_Info_Geometric", 
                                   model = model_name,
                                   name = name,
                                   type = type,  
                                   num_breaks = num_breaks, 
                                   Si = results_geo$original[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Si")],  
                                   Si_std_error = results_geo$std.error[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Si")], 
                                   Si_low_ci = results_geo$low.ci[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Si")], 
                                   Si_high_ci = results_geo$high.ci[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Si")], 
                                   Ti = results_geo$original[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Ti")], 
                                   Ti_std_error = results_geo$std.error[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Ti")], 
                                   Ti_low_ci = results_geo$low.ci[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Ti")], 
                                   Ti_high_ci = results_geo$high.ci[which(results_geo$parameters == params[param] & results_geo$sensitivity == "Ti")],
                                   influential = influential_geo)
    
    ## Store parameter info for multiplicative models
    current_parameter_object_mult = new("Parameter_Info_Multiplicative", 
                                       model = model_name,
                                       name = name,
                                       type = type,  
                                       num_breaks = num_breaks, 
                                       Si = results_mult$original[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Si")],  
                                       Si_std_error = results_mult$std.error[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Si")], 
                                       Si_low_ci = results_mult$low.ci[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Si")], 
                                       Si_high_ci = results_mult$high.ci[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Si")], 
                                       Ti = results_mult$original[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Ti")], 
                                       Ti_std_error = results_mult$std.error[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Ti")], 
                                       Ti_low_ci = results_mult$low.ci[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Ti")], 
                                       Ti_high_ci = results_mult$high.ci[which(results_mult$parameters == params[param] & results_mult$sensitivity == "Ti")],
                                       influential = influential_mult)
    
    ## Store parameter info for geometric model
    current_parameter_object_limit = new("Parameter_Info_Limiting", 
                                       model = model_name,
                                       name = name,
                                       type = type,  
                                       num_breaks = num_breaks, 
                                       Si = results_limit$original[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Si")],  
                                       Si_std_error = results_limit$std.error[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Si")], 
                                       Si_low_ci = results_limit$low.ci[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Si")], 
                                       Si_high_ci = results_limit$high.ci[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Si")], 
                                       Ti = results_limit$original[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Ti")], 
                                       Ti_std_error = results_limit$std.error[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Ti")], 
                                       Ti_low_ci = results_limit$low.ci[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Ti")], 
                                       Ti_high_ci = results_limit$high.ci[which(results_limit$parameters == params[param] & results_limit$sensitivity == "Ti")],
                                       influential = influential_limit)
    
    
    ## Append model info to model objects
    parameter_object = append(parameter_object,current_parameter_object)
    parameter_object_arithmetic = append(parameter_object_arithmetic,current_parameter_object_ari)
    parameter_object_geometric = append(parameter_object_geometric,current_parameter_object_geo)
    parameter_object_multiplicative = append(parameter_object_multiplicative,current_parameter_object_mult)
    parameter_object_limiting = append(parameter_object_limiting,current_parameter_object_limit)
  }
  
}



################################################################################

## Assign names to each class object in the list

par_names=c() ## Create an empty vector to store parameter names

## Store names for each parameter in the list
for (param in 1:length(parameter_object)) {
  par_names[param] = paste0(parameter_object[[param]]@model,".",parameter_object[[param]]@name,sep="")
}

## Set names of parameter_object equal to par_names
names(parameter_object) = par_names
names(parameter_object_arithmetic) = par_names
names(parameter_object_geometric) = par_names
names(parameter_object_multiplicative) = par_names
names(parameter_object_limiting) = par_names

## Save data
saveRDS(model_object, file = "model_object.rds")
saveRDS(model_object_arithmetic, file = "model_object_arithmetic.rds")
saveRDS(model_object_geometric, file = "model_object_geometric.rds")
saveRDS(model_object_multiplicative, file = "model_object_multiplicative.rds")
saveRDS(model_object_limiting, file = "model_object_limiting.rds")
saveRDS(parameter_object, file = "parameter_object.rds")
saveRDS(parameter_object_arithmetic, file = "parameter_object_arithmetic.rds")
saveRDS(parameter_object_geometric, file = "parameter_object_geometric.rds")
saveRDS(parameter_object_multiplicative, file = "parameter_object_multiplicative.rds")
saveRDS(parameter_object_limiting, file = "parameter_object_limiting.rds")


