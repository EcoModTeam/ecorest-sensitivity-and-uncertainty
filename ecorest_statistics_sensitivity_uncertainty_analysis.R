# The purpose of this document is to analyze the results of our ecorest sensitivity and uncertainty analyses
# We will examine patterns of uncertainty and sensitivity across ecorest models
# Date: July 9th, 2025
# Author: Kiara Cushway


################################################################################
# File set-up

# Clear workspace
rm(list = ls()) ## remove stored files and objects
gc(T) ## garbage collection
graphics.off() ## turn off graphics

# Load libraries
library(ecorest) ## contains HSI models and function to calculate HSI
library(ggplot2) ## used for plotting
library(stringr)## used to manipulate strings
library(tidyr) ## used for data manipulation
library(dplyr) ## used for data manipulation
library(ggpubr) ## used for graphing
library(here) ## used for file paths
library(ggrain) ## used for raincloud plots
library(scales)## for graph text
library(moments) ## for assessing skewness

load("HSImetadata.RData") ## load updated HSImetadata
load("HSImodels.RData") ## load updated HSImodels

## Load HSImetadata into data frame, excluding the two models with only one variable
HSImetadata = HSImetadata %>% filter(!model %in% c("redwingedblackbirdB", "woodduckWinter"))

## Remove the two models with one variable from HSImodels as well
HSImodels$redwingedblackbirdB = NULL
HSImodels$woodduckWinter = NULL

# Colorblind friendly color palette for graphs
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Bring in data for each model structure from sensitivity/uncertainty analysis
models_orig = readRDS("model_object.rds") ## original model
models_ari = readRDS("model_object_arithmetic.rds") ## arithmetic mean
models_geo = readRDS("model_object_geometric.rds") ## geometric mean
models_mult = readRDS("model_object_multiplicative.rds") ## multiplicative
models_limit = readRDS("model_object_limiting.rds") ## limiting factor

# Bring in parameter data for each model structure
params_orig = readRDS("parameter_object.rds") ## original model parameters
params_ari = readRDS("parameter_object_arithmetic.rds") ## arithmetic mean model parameters
params_geo = readRDS("parameter_object_geometric.rds") ## geometric mean model parameters
params_mult = readRDS("parameter_object_multiplicative.rds") ## multiplicative model parameters
params_limit = readRDS("parameter_object_limiting.rds") ## limiting factor model parameters

################################################################################

# Pull data from lists and place in data frames for ease of use

# Original model data
uncertainty_mat = matrix(NA,nrow = 347,ncol = 29) ## create empty matrix to store data
uncertainty_df = data.frame(uncertainty_mat) ## convert matrix to data frame to store multiple data types
colnames(uncertainty_df) = c("model", "perc_1", "perc_2.5", "perc_5", "perc_25", 
                             "perc_50", "perc_75", "perc_95","perc_97.5", "perc_99", 
                             "min", "max", "mean", "variance", "std", 
                             "cv", "skewness", "kurtosis", "num_params", "num_components", 
                             "num_categorical", "sum", "dummy_low_Si", "dummy_Si", 
                             "dummy_high_Si", "dummy_low_Ti", "dummy_Ti", "dummy_high_Ti", 
                             "structure")


# Use for loop to fill data frame with data from class
for (model in 1:length(models_orig)) {
  uncertainty_df$model[model] = models_orig[[model]]@name ## store model name
  uncertainty_df$perc_1[model] = quantile(models_orig[[model]]@HSI, 0.01)  ## store 0.01 quantile for HSI values
  uncertainty_df$perc_2.5[model] = quantile(models_orig[[model]]@HSI, 0.025) ## store 0.025 quantile for HSI values
  uncertainty_df$perc_5[model] = quantile(models_orig[[model]]@HSI, 0.05) ## store 0.05 quantile for HSI values
  uncertainty_df$perc_25[model] = quantile(models_orig[[model]]@HSI, 0.25) ## store 0.25 quantile for HSI values
  uncertainty_df$perc_50[model] = quantile(models_orig[[model]]@HSI, 0.5) ## store 0.5 quantile for HSI values
  uncertainty_df$perc_75[model] = quantile(models_orig[[model]]@HSI, 0.75) ## store 0.75 quantile for HSI values
  uncertainty_df$perc_95[model] = quantile(models_orig[[model]]@HSI, 0.95) ## store 0.95 quantile for HSI values
  uncertainty_df$perc_97.5[model] = quantile(models_orig[[model]]@HSI, 0.975) ## store 0.975 quantile for HSI values
  uncertainty_df$perc_99[model] = quantile(models_orig[[model]]@HSI, 0.99) ## store 0.99 quantile for HSI values
  uncertainty_df$mean[model] = mean(models_orig[[model]]@HSI) ## store average HSI value
  uncertainty_df$variance[model] = var(models_orig[[model]]@HSI) ## store variance in HSI values
  uncertainty_df$std[model] = sd(models_orig[[model]]@HSI) ## store standard deviation in HSI values
  uncertainty_df$cv[model] = sd(models_orig[[model]]@HSI)/mean(models_orig[[model]]@HSI) ## store coefficient of variation in HSI values
  uncertainty_df$skewness[model] = skewness(models_orig[[model]]@HSI) # store output HSI distribution skewness
  uncertainty_df$kurtosis[model] = kurtosis(models_orig[[model]]@HSI) ## store output HSI distribution kurtosis
  uncertainty_df$min[model] = min(models_orig[[model]]@HSI) ## store minimum HSI value
  uncertainty_df$max[model] = max(models_orig[[model]]@HSI) ## store max HSI value
  uncertainty_df$num_params[model] = models_orig[[model]]@inputs ## store number of parameters in model
  uncertainty_df$num_components[model] = models_orig[[model]]@components ## store number of components in model
  uncertainty_df$num_categorical[model] = models_orig[[model]]@cat_inputs ## store number of categorical inputs in model
  uncertainty_df$sum[model] = models_orig[[model]]@sum ## store total percent of variance explained by first order indices
  uncertainty_df$dummy_low_Si[model] = models_orig[[model]]@dummy$low.ci[1] ## store lower bound of 95% CI for Si dummy parameter
  uncertainty_df$dummy_Si[model] = models_orig[[model]]@dummy$original[1] ## store dummy parameter estimate for Si
  uncertainty_df$dummy_high_Si[model] = models_orig[[model]]@dummy$high.ci[1] ## store upper bound of 95% CI for Si dummy parameter
  uncertainty_df$dummy_low_Ti[model] = models_orig[[model]]@dummy$low.ci[2] ## store lower bound of 95% CI for Ti dummy parameter 
  uncertainty_df$dummy_Ti[model] = models_orig[[model]]@dummy$original[2] ## store dummy parameter estimate for Ti
  uncertainty_df$dummy_high_Ti[model] = models_orig[[model]]@dummy$high.ci[2] ## store upper bound of 95% CI for Ti dummy parameter 
}

# Bring in data with information about model structures
eqt_types = read.csv("equation_types.csv")

# Assign structures to models
uncertainty_df$structure = eqt_types$eqtn_type_final

# Store parameter data from original models in data frame for ease of use
par_mat = matrix(NA, nrow = 2687, ncol = 13) ## create empty matrix to store parameter data
par_df = data.frame(par_mat) ## convert matrix to data frame to store multiple data types
colnames(par_df) = c("model","parameter","type","num_breakpoints","Si","Si_std_error",
                     "Si_low_ci","Si_high_ci","Ti","Ti_std_error","Ti_low_ci",
                     "Ti_high_ci","influential")

# Use for loop to fill data frame with class data
for (par in 1:length(params_orig)) {
  par_df$model[par] = params_orig[[par]]@model ## store model name
  par_df$parameter[par] = params_orig[[par]]@name ## store parameter name
  par_df$type[par] = params_orig[[par]]@type ## store whether parameter is numeric or categorical
  par_df$num_breakpoints[par] = params_orig[[par]]@num_breaks ## store number of breakpoints in SIV graph
  par_df$Si[par] = params_orig[[par]]@Si ## store Si estimate for parameter
  par_df$Si_std_error[par] = params_orig[[par]]@Si_std_error ## store standard error of Si estimate
  par_df$Si_low_ci[par] = params_orig[[par]]@Si_low_ci ## store lower bound of 95% CI for Si estimate
  par_df$Si_high_ci[par] = params_orig[[par]]@Si_high_ci ## store upper bound of 95% CI for Si estimate
  par_df$Ti[par] = params_orig[[par]]@Ti ## store Ti estimate for parameter
  par_df$Ti_std_error[par] = params_orig[[par]]@Ti_std_error ## store standard error of Ti estimate
  par_df$Ti_low_ci[par] = params_orig[[par]]@Ti_low_ci ## store lower bound for 95% CI for Ti estimate
  par_df$Ti_high_ci[par] = params_orig[[par]]@Ti_high_ci ## store upper bound for 95% CI for Ti estimate
  par_df$influential[par] = params_orig[[par]]@influential ## store whether parameter is influential (1) or not (0)
}


# Find number of influential parameters in each model and add to model summary
num_influential = par_df %>% group_by(model) %>% summarise(important = sum(influential > 0))
num_influential = num_influential[1:347,]
uncertainty_df$influential_params = num_influential$important

# Find proportion of influential parameters
uncertainty_df$prop_influential = uncertainty_df$influential_params / uncertainty_df$num_params

# Find the proportion of input parameters that are categorical
uncertainty_df = uncertainty_df %>% mutate(prop_categorical = num_categorical / num_params)

# Find the percent of input parameters that are influential
uncertainty_df = uncertainty_df %>% mutate(percent_influential = prop_influential * 100)


# Arithmetic mean model data
uncertainty_df_ari = data.frame(uncertainty_mat) ## convert matrix to data frame to store multiple data types
colnames(uncertainty_df_ari) = c("model", "perc_1", "perc_2.5", "perc_5", "perc_25", 
                             "perc_50", "perc_75", "perc_95","perc_97.5", "perc_99", 
                             "min", "max", "mean", "variance", "std", 
                             "cv", "skewness", "kurtosis", "num_params", "num_components", 
                             "num_categorical", "sum", "dummy_low_Si", "dummy_Si", 
                             "dummy_high_Si", "dummy_low_Ti", "dummy_Ti", "dummy_high_Ti", 
                             "structure")


# Use for loop to fill data frame with data from class
for (model in 1:length(models_ari)) {
  uncertainty_df_ari$model[model] = models_ari[[model]]@name ## store model name
  uncertainty_df_ari$perc_1[model] = quantile(models_ari[[model]]@HSI, 0.01)  ## store 0.01 quantile for HSI values
  uncertainty_df_ari$perc_2.5[model] = quantile(models_ari[[model]]@HSI, 0.025) ## store 0.025 quantile for HSI values
  uncertainty_df_ari$perc_5[model] = quantile(models_ari[[model]]@HSI, 0.05) ## store 0.05 quantile for HSI values
  uncertainty_df_ari$perc_25[model] = quantile(models_ari[[model]]@HSI, 0.25) ## store 0.25 quantile for HSI values
  uncertainty_df_ari$perc_50[model] = quantile(models_ari[[model]]@HSI, 0.5) ## store 0.5 quantile for HSI values
  uncertainty_df_ari$perc_75[model] = quantile(models_ari[[model]]@HSI, 0.75) ## store 0.75 quantile for HSI values
  uncertainty_df_ari$perc_95[model] = quantile(models_ari[[model]]@HSI, 0.95) ## store 0.95 quantile for HSI values
  uncertainty_df_ari$perc_97.5[model] = quantile(models_ari[[model]]@HSI, 0.975) ## store 0.975 quantile for HSI values
  uncertainty_df_ari$perc_99[model] = quantile(models_ari[[model]]@HSI, 0.99) ## store 0.99 quantile for HSI values
  uncertainty_df_ari$mean[model] = mean(models_ari[[model]]@HSI) ## store average HSI value
  uncertainty_df_ari$variance[model] = var(models_ari[[model]]@HSI) ## store variance in HSI values
  uncertainty_df_ari$std[model] = sd(models_ari[[model]]@HSI) ## store standard deviation in HSI values
  uncertainty_df_ari$cv[model] = sd(models_ari[[model]]@HSI)/mean(models_ari[[model]]@HSI) ## store coefficient of variation in HSI values
  uncertainty_df_ari$skewness[model] = skewness(models_ari[[model]]@HSI) # store output HSI distribution skewness
  uncertainty_df_ari$kurtosis[model] = kurtosis(models_ari[[model]]@HSI) ## store output HSI distribution kurtosis
  uncertainty_df_ari$min[model] = min(models_ari[[model]]@HSI) ## store minimum HSI value
  uncertainty_df_ari$max[model] = max(models_ari[[model]]@HSI) ## store max HSI value
  uncertainty_df_ari$num_params[model] = models_ari[[model]]@inputs ## store number of parameters in model
  uncertainty_df_ari$num_components[model] = models_ari[[model]]@components ## store number of components in model
  uncertainty_df_ari$num_categorical[model] = models_ari[[model]]@cat_inputs ## store number of categorical inputs in model
  uncertainty_df_ari$sum[model] = models_ari[[model]]@sum ## store total percent of variance explained by first order indices
  uncertainty_df_ari$dummy_low_Si[model] = models_ari[[model]]@dummy$low.ci[1] ## store lower bound of 95% CI for Si dummy parameter
  uncertainty_df_ari$dummy_Si[model] = models_ari[[model]]@dummy$original[1] ## store dummy parameter estimate for Si
  uncertainty_df_ari$dummy_high_Si[model] = models_ari[[model]]@dummy$high.ci[1] ## store upper bound of 95% CI for Si dummy parameter
  uncertainty_df_ari$dummy_low_Ti[model] = models_ari[[model]]@dummy$low.ci[2] ## store lower bound of 95% CI for Ti dummy parameter 
  uncertainty_df_ari$dummy_Ti[model] = models_ari[[model]]@dummy$original[2] ## store dummy parameter estimate for Ti
  uncertainty_df_ari$dummy_high_Ti[model] = models_ari[[model]]@dummy$high.ci[2] ## store upper bound of 95% CI for Ti dummy parameter 
}

# Assign structure to models
uncertainty_df_ari$structure = "Arithmetic mean"

# Store parameter data from arithmetic mean models in data frame for ease of use
par_df_ari = data.frame(par_mat) ## convert matrix to data frame to store multiple data types
colnames(par_df_ari) = c("model","parameter","type","num_breakpoints","Si","Si_std_error",
                     "Si_low_ci","Si_high_ci","Ti","Ti_std_error","Ti_low_ci",
                     "Ti_high_ci","influential")

# Use for loop to fill data frame with class data
for (par in 1:length(params_ari)) {
  par_df_ari$model[par] = params_ari[[par]]@model ## store model name
  par_df_ari$parameter[par] = params_ari[[par]]@name ## store parameter name
  par_df_ari$type[par] = params_ari[[par]]@type ## store whether parameter is numeric or categorical
  par_df_ari$num_breakpoints[par] = params_ari[[par]]@num_breaks ## store number of breakpoints in SIV graph
  par_df_ari$Si[par] = params_ari[[par]]@Si ## store Si estimate for parameter
  par_df_ari$Si_std_error[par] = params_ari[[par]]@Si_std_error ## store standard error of Si estimate
  par_df_ari$Si_low_ci[par] = params_ari[[par]]@Si_low_ci ## store lower bound of 95% CI for Si estimate
  par_df_ari$Si_high_ci[par] = params_ari[[par]]@Si_high_ci ## store upper bound of 95% CI for Si estimate
  par_df_ari$Ti[par] = params_ari[[par]]@Ti ## store Ti estimate for parameter
  par_df_ari$Ti_std_error[par] = params_ari[[par]]@Ti_std_error ## store standard error of Ti estimate
  par_df_ari$Ti_low_ci[par] = params_ari[[par]]@Ti_low_ci ## store lower bound for 95% CI for Ti estimate
  par_df_ari$Ti_high_ci[par] = params_ari[[par]]@Ti_high_ci ## store upper bound for 95% CI for Ti estimate
  par_df_ari$influential[par] = params_ari[[par]]@influential ## store whether parameter is influential (1) or not (0)
}


# Find number of influential parameters in each model and add to model summary
num_influential_ari = par_df_ari %>% group_by(model) %>% summarise(important = sum(influential > 0))
num_influential_ari = num_influential_ari[1:347,]
uncertainty_df_ari$influential_params = num_influential_ari$important

# Find proportion of influential parameters
uncertainty_df_ari$prop_influential = uncertainty_df_ari$influential_params / uncertainty_df_ari$num_params

# Find the proportion of input parameters that are categorical
uncertainty_df_ari = uncertainty_df_ari %>% mutate(prop_categorical = num_categorical / num_params)

# Find the percent of input parameters that are influential
uncertainty_df_ari = uncertainty_df_ari %>% mutate(percent_influential = prop_influential * 100)



# Geometric mean model data
uncertainty_df_geo = data.frame(uncertainty_mat) ## convert matrix to data frame to store multiple data types
colnames(uncertainty_df_geo) = c("model", "perc_1", "perc_2.5", "perc_5", "perc_25", 
                             "perc_50", "perc_75", "perc_95","perc_97.5", "perc_99", 
                             "min", "max", "mean", "variance", "std", 
                             "cv", "skewness", "kurtosis", "num_params", "num_components", 
                             "num_categorical", "sum", "dummy_low_Si", "dummy_Si", 
                             "dummy_high_Si", "dummy_low_Ti", "dummy_Ti", "dummy_high_Ti", 
                             "structure")


# Use for loop to fill data frame with data from class
for (model in 1:length(models_geo)) {
  uncertainty_df_geo$model[model] = models_geo[[model]]@name ## store model name
  uncertainty_df_geo$perc_1[model] = quantile(models_geo[[model]]@HSI, 0.01)  ## store 0.01 quantile for HSI values
  uncertainty_df_geo$perc_2.5[model] = quantile(models_geo[[model]]@HSI, 0.025) ## store 0.025 quantile for HSI values
  uncertainty_df_geo$perc_5[model] = quantile(models_geo[[model]]@HSI, 0.05) ## store 0.05 quantile for HSI values
  uncertainty_df_geo$perc_25[model] = quantile(models_geo[[model]]@HSI, 0.25) ## store 0.25 quantile for HSI values
  uncertainty_df_geo$perc_50[model] = quantile(models_geo[[model]]@HSI, 0.5) ## store 0.5 quantile for HSI values
  uncertainty_df_geo$perc_75[model] = quantile(models_geo[[model]]@HSI, 0.75) ## store 0.75 quantile for HSI values
  uncertainty_df_geo$perc_95[model] = quantile(models_geo[[model]]@HSI, 0.95) ## store 0.95 quantile for HSI values
  uncertainty_df_geo$perc_97.5[model] = quantile(models_geo[[model]]@HSI, 0.975) ## store 0.975 quantile for HSI values
  uncertainty_df_geo$perc_99[model] = quantile(models_geo[[model]]@HSI, 0.99) ## store 0.99 quantile for HSI values
  uncertainty_df_geo$mean[model] = mean(models_geo[[model]]@HSI) ## store average HSI value
  uncertainty_df_geo$variance[model] = var(models_geo[[model]]@HSI) ## store variance in HSI values
  uncertainty_df_geo$std[model] = sd(models_geo[[model]]@HSI) ## store standard deviation in HSI values
  uncertainty_df_geo$cv[model] = sd(models_geo[[model]]@HSI)/mean(models_geo[[model]]@HSI) ## store coefficient of variation in HSI values
  uncertainty_df_geo$skewness[model] = skewness(models_geo[[model]]@HSI) # store output HSI distribution skewness
  uncertainty_df_geo$kurtosis[model] = kurtosis(models_geo[[model]]@HSI) ## store output HSI distribution kurtosis
  uncertainty_df_geo$min[model] = min(models_geo[[model]]@HSI) ## store minimum HSI value
  uncertainty_df_geo$max[model] = max(models_geo[[model]]@HSI) ## store max HSI value
  uncertainty_df_geo$num_params[model] = models_geo[[model]]@inputs ## store number of parameters in model
  uncertainty_df_geo$num_components[model] = models_geo[[model]]@components ## store number of components in model
  uncertainty_df_geo$num_categorical[model] = models_geo[[model]]@cat_inputs ## store number of categorical inputs in model
  uncertainty_df_geo$sum[model] = models_geo[[model]]@sum ## store total percent of variance explained by first order indices
  uncertainty_df_geo$dummy_low_Si[model] = models_geo[[model]]@dummy$low.ci[1] ## store lower bound of 95% CI for Si dummy parameter
  uncertainty_df_geo$dummy_Si[model] = models_geo[[model]]@dummy$original[1] ## store dummy parameter estimate for Si
  uncertainty_df_geo$dummy_high_Si[model] = models_geo[[model]]@dummy$high.ci[1] ## store upper bound of 95% CI for Si dummy parameter
  uncertainty_df_geo$dummy_low_Ti[model] = models_geo[[model]]@dummy$low.ci[2] ## store lower bound of 95% CI for Ti dummy parameter 
  uncertainty_df_geo$dummy_Ti[model] = models_geo[[model]]@dummy$original[2] ## store dummy parameter estimate for Ti
  uncertainty_df_geo$dummy_high_Ti[model] = models_geo[[model]]@dummy$high.ci[2] ## store upper bound of 95% CI for Ti dummy parameter 
}

# Assign structures to models
uncertainty_df_geo$structure = "Geometric mean"

# Store parameter data from geometric mean models in data frame for ease of use
par_df_geo = data.frame(par_mat) ## convert matrix to data frame to store multiple data types
colnames(par_df_geo) = c("model","parameter","type","num_breakpoints","Si","Si_std_error",
                     "Si_low_ci","Si_high_ci","Ti","Ti_std_error","Ti_low_ci",
                     "Ti_high_ci","influential")

# Use for loop to fill data frame with class data
for (par in 1:length(params_geo)) {
  par_df_geo$model[par] = params_geo[[par]]@model ## store model name
  par_df_geo$parameter[par] = params_geo[[par]]@name ## store parameter name
  par_df_geo$type[par] = params_geo[[par]]@type ## store whether parameter is numeric or categorical
  par_df_geo$num_breakpoints[par] = params_geo[[par]]@num_breaks ## store number of breakpoints in SIV graph
  par_df_geo$Si[par] = params_geo[[par]]@Si ## store Si estimate for parameter
  par_df_geo$Si_std_error[par] = params_geo[[par]]@Si_std_error ## store standard error of Si estimate
  par_df_geo$Si_low_ci[par] = params_geo[[par]]@Si_low_ci ## store lower bound of 95% CI for Si estimate
  par_df_geo$Si_high_ci[par] = params_geo[[par]]@Si_high_ci ## store upper bound of 95% CI for Si estimate
  par_df_geo$Ti[par] = params_geo[[par]]@Ti ## store Ti estimate for parameter
  par_df_geo$Ti_std_error[par] = params_geo[[par]]@Ti_std_error ## store standard error of Ti estimate
  par_df_geo$Ti_low_ci[par] = params_geo[[par]]@Ti_low_ci ## store lower bound for 95% CI for Ti estimate
  par_df_geo$Ti_high_ci[par] = params_geo[[par]]@Ti_high_ci ## store upper bound for 95% CI for Ti estimate
  par_df_geo$influential[par] = params_geo[[par]]@influential ## store whether parameter is influential (1) or not (0)
}


# Find number of influential parameters in each model and add to model summary
num_influential_geo = par_df_geo %>% group_by(model) %>% summarise(important = sum(influential > 0))
num_influential_geo = num_influential_geo[1:347,]
uncertainty_df_geo$influential_params = num_influential_geo$important

# Find proportion of influential parameters
uncertainty_df_geo$prop_influential = uncertainty_df_geo$influential_params / uncertainty_df_geo$num_params

# Find the proportion of input parameters that are categorical
uncertainty_df_geo = uncertainty_df_geo %>% mutate(prop_categorical = num_categorical / num_params)

# Find the percent of input parameters that are influential
uncertainty_df_geo = uncertainty_df_geo %>% mutate(percent_influential = prop_influential * 100)



# Limiting factor model data
uncertainty_df_limit = data.frame(uncertainty_mat) ## convert matrix to data frame to store multiple data types
colnames(uncertainty_df_limit) = c("model", "perc_1", "perc_2.5", "perc_5", "perc_25", 
                             "perc_50", "perc_75", "perc_95","perc_97.5", "perc_99", 
                             "min", "max", "mean", "variance", "std", 
                             "cv", "skewness", "kurtosis", "num_params", "num_components", 
                             "num_categorical", "sum", "dummy_low_Si", "dummy_Si", 
                             "dummy_high_Si", "dummy_low_Ti", "dummy_Ti", "dummy_high_Ti", 
                             "structure")


# Use for loop to fill data frame with data from class
for (model in 1:length(models_limit)) {
  uncertainty_df_limit$model[model] = models_limit[[model]]@name ## store model name
  uncertainty_df_limit$perc_1[model] = quantile(models_limit[[model]]@HSI, 0.01)  ## store 0.01 quantile for HSI values
  uncertainty_df_limit$perc_2.5[model] = quantile(models_limit[[model]]@HSI, 0.025) ## store 0.025 quantile for HSI values
  uncertainty_df_limit$perc_5[model] = quantile(models_limit[[model]]@HSI, 0.05) ## store 0.05 quantile for HSI values
  uncertainty_df_limit$perc_25[model] = quantile(models_limit[[model]]@HSI, 0.25) ## store 0.25 quantile for HSI values
  uncertainty_df_limit$perc_50[model] = quantile(models_limit[[model]]@HSI, 0.5) ## store 0.5 quantile for HSI values
  uncertainty_df_limit$perc_75[model] = quantile(models_limit[[model]]@HSI, 0.75) ## store 0.75 quantile for HSI values
  uncertainty_df_limit$perc_95[model] = quantile(models_limit[[model]]@HSI, 0.95) ## store 0.95 quantile for HSI values
  uncertainty_df_limit$perc_97.5[model] = quantile(models_limit[[model]]@HSI, 0.975) ## store 0.975 quantile for HSI values
  uncertainty_df_limit$perc_99[model] = quantile(models_limit[[model]]@HSI, 0.99) ## store 0.99 quantile for HSI values
  uncertainty_df_limit$mean[model] = mean(models_limit[[model]]@HSI) ## store average HSI value
  uncertainty_df_limit$variance[model] = var(models_limit[[model]]@HSI) ## store variance in HSI values
  uncertainty_df_limit$std[model] = sd(models_limit[[model]]@HSI) ## store standard deviation in HSI values
  uncertainty_df_limit$cv[model] = sd(models_limit[[model]]@HSI)/mean(models_limit[[model]]@HSI) ## store coefficient of variation in HSI values
  uncertainty_df_limit$skewness[model] = skewness(models_limit[[model]]@HSI) # store output HSI distribution skewness
  uncertainty_df_limit$kurtosis[model] = kurtosis(models_limit[[model]]@HSI) ## store output HSI distribution kurtosis
  uncertainty_df_limit$min[model] = min(models_limit[[model]]@HSI) ## store minimum HSI value
  uncertainty_df_limit$max[model] = max(models_limit[[model]]@HSI) ## store max HSI value
  uncertainty_df_limit$num_params[model] = models_limit[[model]]@inputs ## store number of parameters in model
  uncertainty_df_limit$num_components[model] = models_limit[[model]]@components ## store number of components in model
  uncertainty_df_limit$num_categorical[model] = models_limit[[model]]@cat_inputs ## store number of categorical inputs in model
  uncertainty_df$sum[model] = models_limit[[model]]@sum ## store total percent of variance explained by first order indices
  uncertainty_df_limit$dummy_low_Si[model] = models_limit[[model]]@dummy$low.ci[1] ## store lower bound of 95% CI for Si dummy parameter
  uncertainty_df_limit$dummy_Si[model] = models_limit[[model]]@dummy$original[1] ## store dummy parameter estimate for Si
  uncertainty_df_limit$dummy_high_Si[model] = models_limit[[model]]@dummy$high.ci[1] ## store upper bound of 95% CI for Si dummy parameter
  uncertainty_df_limit$dummy_low_Ti[model] = models_limit[[model]]@dummy$low.ci[2] ## store lower bound of 95% CI for Ti dummy parameter 
  uncertainty_df_limit$dummy_Ti[model] = models_limit[[model]]@dummy$original[2] ## store dummy parameter estimate for Ti
  uncertainty_df_limit$dummy_high_Ti[model] = models_limit[[model]]@dummy$high.ci[2] ## store upper bound of 95% CI for Ti dummy parameter 
}


# Assign structures to models
uncertainty_df_limit$structure = "Limiting factor"

# Store parameter data from limiting factors models in data frame for ease of use
par_df_limit = data.frame(par_mat) ## convert matrix to data frame to store multiple data types
colnames(par_df_limit) = c("model","parameter","type","num_breakpoints","Si","Si_std_error",
                     "Si_low_ci","Si_high_ci","Ti","Ti_std_error","Ti_low_ci",
                     "Ti_high_ci","influential")

# Use for loop to fill data frame with class data
for (par in 1:length(params_limit)) {
  par_df_limit$model[par] = params_limit[[par]]@model ## store model name
  par_df_limit$parameter[par] = params_limit[[par]]@name ## store parameter name
  par_df_limit$type[par] = params_limit[[par]]@type ## store whether parameter is numeric or categorical
  par_df_limit$num_breakpoints[par] = params_limit[[par]]@num_breaks ## store number of breakpoints in SIV graph
  par_df_limit$Si[par] = params_limit[[par]]@Si ## store Si estimate for parameter
  par_df_limit$Si_std_error[par] = params_limit[[par]]@Si_std_error ## store standard error of Si estimate
  par_df_limit$Si_low_ci[par] = params_limit[[par]]@Si_low_ci ## store lower bound of 95% CI for Si estimate
  par_df_limit$Si_high_ci[par] = params_limit[[par]]@Si_high_ci ## store upper bound of 95% CI for Si estimate
  par_df_limit$Ti[par] = params_limit[[par]]@Ti ## store Ti estimate for parameter
  par_df_limit$Ti_std_error[par] = params_limit[[par]]@Ti_std_error ## store standard error of Ti estimate
  par_df_limit$Ti_low_ci[par] = params_limit[[par]]@Ti_low_ci ## store lower bound for 95% CI for Ti estimate
  par_df_limit$Ti_high_ci[par] = params_limit[[par]]@Ti_high_ci ## store upper bound for 95% CI for Ti estimate
  par_df_limit$influential[par] = params_limit[[par]]@influential ## store whether parameter is influential (1) or not (0)
}


# Find number of influential parameters in each model and add to model summary
num_influential_limit = par_df_limit %>% group_by(model) %>% summarise(important = sum(influential > 0))
num_influential_limit = num_influential_limit[1:347,]
uncertainty_df_limit$influential_params = num_influential_limit$important

# Find proportion of influential parameters
uncertainty_df_limit$prop_influential = uncertainty_df_limit$influential_params / uncertainty_df_limit$num_params

# Find the proportion of input parameters that are categorical
uncertainty_df_limit = uncertainty_df_limit %>% mutate(prop_categorical = num_categorical / num_params)

# Find the percent of input parameters that are influential
uncertainty_df_limit = uncertainty_df_limit %>% mutate(percent_influential = prop_influential * 100)



# Multiplicative model data
uncertainty_df_mult = data.frame(uncertainty_mat) ## convert matrix to data frame to store multiple data types
colnames(uncertainty_df_mult) = c("model", "perc_1", "perc_2.5", "perc_5", "perc_25", 
                             "perc_50", "perc_75", "perc_95","perc_97.5", "perc_99", 
                             "min", "max", "mean", "variance", "std", 
                             "cv", "skewness", "kurtosis", "num_params", "num_components", 
                             "num_categorical", "sum", "dummy_low_Si", "dummy_Si", 
                             "dummy_high_Si", "dummy_low_Ti", "dummy_Ti", "dummy_high_Ti", 
                             "structure")


# Use for loop to fill data frame with data from class
for (model in 1:length(models_mult)) {
  uncertainty_df_mult$model[model] = models_mult[[model]]@name ## store model name
  uncertainty_df_mult$perc_1[model] = quantile(models_mult[[model]]@HSI, 0.01)  ## store 0.01 quantile for HSI values
  uncertainty_df_mult$perc_2.5[model] = quantile(models_mult[[model]]@HSI, 0.025) ## store 0.025 quantile for HSI values
  uncertainty_df_mult$perc_5[model] = quantile(models_mult[[model]]@HSI, 0.05) ## store 0.05 quantile for HSI values
  uncertainty_df_mult$perc_25[model] = quantile(models_mult[[model]]@HSI, 0.25) ## store 0.25 quantile for HSI values
  uncertainty_df_mult$perc_50[model] = quantile(models_mult[[model]]@HSI, 0.5) ## store 0.5 quantile for HSI values
  uncertainty_df_mult$perc_75[model] = quantile(models_mult[[model]]@HSI, 0.75) ## store 0.75 quantile for HSI values
  uncertainty_df_mult$perc_95[model] = quantile(models_mult[[model]]@HSI, 0.95) ## store 0.95 quantile for HSI values
  uncertainty_df_mult$perc_97.5[model] = quantile(models_mult[[model]]@HSI, 0.975) ## store 0.975 quantile for HSI values
  uncertainty_df_mult$perc_99[model] = quantile(models_mult[[model]]@HSI, 0.99) ## store 0.99 quantile for HSI values
  uncertainty_df_mult$mean[model] = mean(models_mult[[model]]@HSI) ## store average HSI value
  uncertainty_df_mult$variance[model] = var(models_mult[[model]]@HSI) ## store variance in HSI values
  uncertainty_df_mult$std[model] = sd(models_mult[[model]]@HSI) ## store standard deviation in HSI values
  uncertainty_df_mult$cv[model] = sd(models_mult[[model]]@HSI)/mean(models_mult[[model]]@HSI) ## store coefficient of variation in HSI values
  uncertainty_df_mult$skewness[model] = skewness(models_mult[[model]]@HSI) # store output HSI distribution skewness
  uncertainty_df_mult$kurtosis[model] = kurtosis(models_mult[[model]]@HSI) ## store output HSI distribution kurtosis
  uncertainty_df_mult$min[model] = min(models_mult[[model]]@HSI) ## store minimum HSI value
  uncertainty_df_mult$max[model] = max(models_mult[[model]]@HSI) ## store max HSI value
  uncertainty_df_mult$num_params[model] = models_mult[[model]]@inputs ## store number of parameters in model
  uncertainty_df_mult$num_components[model] = models_mult[[model]]@components ## store number of components in model
  uncertainty_df_mult$num_categorical[model] = models_mult[[model]]@cat_inputs ## store number of categorical inputs in model
  uncertainty_df_mult$sum[model] = models_mult[[model]]@sum ## store total percent of variance explained by first order indices
  uncertainty_df_mult$dummy_low_Si[model] = models_mult[[model]]@dummy$low.ci[1] ## store lower bound of 95% CI for Si dummy parameter
  uncertainty_df_mult$dummy_Si[model] = models_mult[[model]]@dummy$original[1] ## store dummy parameter estimate for Si
  uncertainty_df_mult$dummy_high_Si[model] = models_mult[[model]]@dummy$high.ci[1] ## store upper bound of 95% CI for Si dummy parameter
  uncertainty_df_mult$dummy_low_Ti[model] = models_mult[[model]]@dummy$low.ci[2] ## store lower bound of 95% CI for Ti dummy parameter 
  uncertainty_df_mult$dummy_Ti[model] = models_mult[[model]]@dummy$original[2] ## store dummy parameter estimate for Ti
  uncertainty_df_mult$dummy_high_Ti[model] = models_mult[[model]]@dummy$high.ci[2] ## store upper bound of 95% CI for Ti dummy parameter 
}

# Assign structures to models
uncertainty_df_mult$structure = "Multiplicative"

# Store parameter data from multiplicative models in data frame for ease of use
par_df_mult = data.frame(par_mat) ## convert matrix to data frame to store multiple data types
colnames(par_df_mult) = c("model","parameter","type","num_breakpoints","Si","Si_std_error",
                     "Si_low_ci","Si_high_ci","Ti","Ti_std_error","Ti_low_ci",
                     "Ti_high_ci","influential")

# Use for loop to fill data frame with class data
for (par in 1:length(params_mult)) {
  par_df_mult$model[par] = params_mult[[par]]@model ## store model name
  par_df_mult$parameter[par] = params_mult[[par]]@name ## store parameter name
  par_df_mult$type[par] = params_mult[[par]]@type ## store whether parameter is numeric or categorical
  par_df_mult$num_breakpoints[par] = params_mult[[par]]@num_breaks ## store number of breakpoints in SIV graph
  par_df_mult$Si[par] = params_mult[[par]]@Si ## store Si estimate for parameter
  par_df_mult$Si_std_error[par] = params_mult[[par]]@Si_std_error ## store standard error of Si estimate
  par_df_mult$Si_low_ci[par] = params_mult[[par]]@Si_low_ci ## store lower bound of 95% CI for Si estimate
  par_df_mult$Si_high_ci[par] = params_mult[[par]]@Si_high_ci ## store upper bound of 95% CI for Si estimate
  par_df_mult$Ti[par] = params_mult[[par]]@Ti ## store Ti estimate for parameter
  par_df_mult$Ti_std_error[par] = params_mult[[par]]@Ti_std_error ## store standard error of Ti estimate
  par_df_mult$Ti_low_ci[par] = params_mult[[par]]@Ti_low_ci ## store lower bound for 95% CI for Ti estimate
  par_df_mult$Ti_high_ci[par] = params_mult[[par]]@Ti_high_ci ## store upper bound for 95% CI for Ti estimate
  par_df_mult$influential[par] = params_mult[[par]]@influential ## store whether parameter is influential (1) or not (0)
}


# Find number of influential parameters in each model and add to model summary
num_influential_mult = par_df_mult %>% group_by(model) %>% summarise(important = sum(influential > 0))
num_influential_mult = num_influential_mult[1:347,]
uncertainty_df_mult$influential_params = num_influential_mult$important

# Find proportion of influential parameters
uncertainty_df_mult$prop_influential = uncertainty_df_mult$influential_params / uncertainty_df_mult$num_params

# Find the proportion of input parameters that are categorical
uncertainty_df_mult = uncertainty_df_mult %>% mutate(prop_categorical = num_categorical / num_params)

# Find the percent of input parameters that are influential
uncertainty_df_mult = uncertainty_df_mult %>% mutate(percent_influential = prop_influential * 100)

################################################################################

# Create data frame for evaluating model structures that contains all model versions

# Remove original models that were geometric mean, arithmetic mean, limiting factor, or multiplicative formats
uncertainty_df_orig = uncertainty_df %>% filter(structure != "geometric mean" & 
                                                  structure != "arithmetic mean" & 
                                                  structure != "multiplicative" & 
                                                  structure != "limiting factor")
uncertainty_df_orig$structure = "Original"

# Filter out same equations from original model that were already geometric or arithmetic means
uncertainty_df_ari = uncertainty_df_ari %>% filter(model %in% uncertainty_df_orig$model)
uncertainty_df_geo = uncertainty_df_geo %>% filter(model %in% uncertainty_df_orig$model)
uncertainty_df_mult = uncertainty_df_mult %>% filter(model %in% uncertainty_df_orig$model)
uncertainty_df_limit = uncertainty_df_limit %>% filter(model %in% uncertainty_df_orig$model)

# Combine all versions into a single data frame
all_versions = rbind(uncertainty_df_orig, uncertainty_df_ari, uncertainty_df_geo, 
                     uncertainty_df_mult, uncertainty_df_limit)



################################################################################
# Perform some basic quality control measures

# Check that any negative Si values have confidence intervals that overlap zero
for (param in 1:nrow(par_df)) {
  if (par_df$Si[param] < 0) {
    if(par_df$Si_high_ci[param] < 0) {
      print(par_df$model[param])
      print(par_df$parameter[param])
      print("Check this one!")
    }
  }
}

# Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Si_summary = par_df %>% group_by(model) %>% summarise(sum_low_ci = sum(Si_low_ci), 
                                                      sum_estimate = sum(Si), 
                                                      sum_high_ci = sum(Si_high_ci))

# Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Ti_summary = par_df %>% group_by(model) %>% summarise(sum_low_ci = sum(Ti_low_ci), 
                                                      sum_estimate = sum(Ti), 
                                                      sum_high_ci = sum(Ti_high_ci))

## Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Si_summary_ari = par_df_ari %>% group_by(model) %>% summarise(sum_low_ci = sum(Si_low_ci), 
                                                              sum_estimate = sum(Si), 
                                                              sum_high_ci = sum(Si_high_ci))

# Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Ti_summary_ari = par_df_ari %>% group_by(model) %>% summarise(sum_low_ci = sum(Ti_low_ci), 
                                                              sum_estimate = sum(Ti), 
                                                              sum_high_ci = sum(Ti_high_ci))

# Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Si_summary_geo = par_df_geo %>% group_by(model) %>% summarise(sum_low_ci = sum(Si_low_ci), 
                                                              sum_estimate = sum(Si), 
                                                              sum_high_ci = sum(Si_high_ci))

# Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Ti_summary_geo = par_df_geo %>% group_by(model) %>% summarise(sum_low_ci = sum(Ti_low_ci), 
                                                              sum_estimate = sum(Ti), 
                                                              sum_high_ci = sum(Ti_high_ci))

# Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Si_summary_limit = par_df_limit %>% group_by(model) %>% summarise(sum_low_ci = sum(Si_low_ci), 
                                                                  sum_estimate = sum(Si), 
                                                                  sum_high_ci = sum(Si_high_ci))

# Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Ti_summary_limit = par_df_limit %>% group_by(model) %>% summarise(sum_low_ci = sum(Ti_low_ci), 
                                                                  sum_estimate = sum(Ti), 
                                                                  sum_high_ci = sum(Ti_high_ci))

# Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Si_summary_mult = par_df_mult %>% group_by(model) %>% summarise(sum_low_ci = sum(Si_low_ci), 
                                                                sum_estimate = sum(Si), 
                                                                sum_high_ci = sum(Si_high_ci))

# Check whether the sum of Si values isn't greater than 1 (or much greater at least)
Ti_summary_mult = par_df_mult %>% group_by(model) %>% summarise(sum_low_ci = sum(Ti_low_ci), 
                                                                sum_estimate = sum(Ti), 
                                                                sum_high_ci = sum(Ti_high_ci))

################################################################################

# Examining model structure in original USFWS models

# Summarize the number and percent of models exhibiting each structure
uncertainty_df %>% group_by(structure) %>% summarise(count = n(), percent = n()/347*100)

# Determine range of num_params
range(uncertainty_df$num_params)

# Determine range of num_components
range(uncertainty_df$num_components)

# Determine the number of models that have at least one categorical parameter
uncertainty_df %>% filter(num_categorical > 0) %>% summarise(count = n(), percent = (n() / 347) *100)

# Determine the number of models that have at least one numeric parameter
uncertainty_df %>% filter((num_params - num_categorical) > 0) %>% summarise(count = n(), percent = (n() / 347) *100)

################################################################################

# Patterns of HSI model sensitivity in USFWS ecorest models

# Count the number of models that have at least one uninfluential parameter
uncertainty_df %>% filter(percent_influential < 100) %>% summarise(count = n(), percent = (n() / 347) * 100)

# Find average, sd, and median number of influential parameters for ecorest models
uncertainty_df %>% summarise(mean = mean(percent_influential), 
                             sd = sd(percent_influential), 
                             median = median(percent_influential))

# Find number of models with more than half of all parameters that are uninfluential
uncertainty_df %>% filter(percent_influential < 50) %>% summarise(count = n(), percent = (n() / 347) * 100)

# Find the maximum number of input parameters in a model that are all influential
uncertainty_df %>% filter(percent_influential == 100) %>% summarise(max = max(num_params))

# Find the correlation between the percent of influential parameters and the total number of parameters in a model
## Plot against each other
plot(x = uncertainty_df$num_params, y = uncertainty_df$percent_influential, 
     xlab = "Number of input parameters", ylab = "Percent of total parameters that are influential")

## Run Spearman rank correlation
cor_test_num_params = cor.test(uncertainty_df$percent_influential, uncertainty_df$num_params, method = "spearman")

## Display results of correlation test
cor_test_num_params


# Find the correlation between the percent of influential parameters and the total number of parameters in a model
## Plot against each other
plot(x = uncertainty_df$num_components, y = uncertainty_df$percent_influential, 
     xlab = "Number of input components", ylab = "Percent of total parameters that are influential")

## Run Spearman rank correlation
cor_test_num_comps = cor.test(uncertainty_df$percent_influential, uncertainty_df$num_components, method = "spearman")

## Display results of correlation test
cor_test_num_comps


# Find the correlation between the percent of influential parameters and the proportion of categorical parameters
## Plot against each other
plot(x = uncertainty_df$prop_categorical, y = uncertainty_df$percent_influential, 
     xlab = "Proportion of categorical input parameters", ylab = "Percent of total parameters that are influential")

## Run Spearman rank correlation
cor_test_prop_cat = cor.test(uncertainty_df$percent_influential, uncertainty_df$prop_categorical, method = "spearman")

## Display results of correlation test
cor_test_prop_cat

# FDR adjustment of p values for correlation tests using percent_influential
p_vals_influential = c(cor_test_num_params$p.value, cor_test_num_comps$p.value, cor_test_prop_cat$p.value)
adjusted_p_vals_influential = p.adjust(p_vals_influential, method = "fdr")
adjusted_p_vals_influential

# Create figure 4: scatter plot of percent influential versus num_params, num_components, and prop_categorical
# Num_params graph
fig_4a = uncertainty_df %>% ggplot(aes(x = num_params, y = percent_influential)
                  ) + geom_point(size = 3, color = "#E69F00") + labs (y = "Percent of influential input parameters", 
  x = "Total number of input parameters", title = "A") + theme(text = element_text(color = "black", 
  size = 18)) + theme(axis.text = element_text(color = "black")
  ) + theme(axis.title.y = element_text(hjust=0.75)) + labs(color =
     "Structure") + theme(legend.position = "none")  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# Num_components graph
fig_4b = uncertainty_df %>% ggplot(aes(x = num_components, y = percent_influential)
                  ) + geom_point(size = 3, color = "#56B4E9") + labs (y = "Percent of influential input parameters", 
  x = "Total number of input components", title = "B") + theme(text = element_text(color = "black", 
  size = 18)) + theme(axis.text = element_text(color = "black")
  ) + theme(axis.title.y = element_text(hjust=0.75)) + labs(color =
     "Structure") + theme(legend.position = "none")  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# Prop_categorical graph
fig_4c = uncertainty_df %>% ggplot(aes(x = prop_categorical, y = percent_influential)
                  ) + geom_point(size = 3, color = "#009E73") + labs (y = "Percent of influential input parameters", 
  x = "Proportion of categorical inputs", title = "C") + theme(text = element_text(color = "black", 
  size = 18)) + theme(axis.text = element_text(color = "black")
  ) + theme(axis.title.y = element_text(hjust=0.75)) + labs(color =
     "Structure") + theme(legend.position = "none")  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# Combine graphs
fig_4 = ggarrange(fig_4a, fig_4b, fig_4c, ncol = 3)

# Display graph
fig_4

# Save figure as a png file
### ggsave("Figure 4.png", dpi = 600, height = 5, width = 15)


################################################################################

# Influence of parameter characteristics on sensitivity

# Run Spearman correlations to test whether Si is correlated with parameter characteristics

# Plot relationship between breakpoints and Si
plot(x = par_df$num_breakpoints, y = par_df$Si, 
     xlab = "Number of parameter breakpoints", 
     ylab = "First order Sobol' index estimate (Si)")

# Run correlation test
cor_test_breakpoints_Si = cor.test(par_df$Si, par_df$num_breakpoints, method = "spearman")

# Display results
cor_test_breakpoints_Si

# Plot Si for categorical versus numeric parameters
plot(x = as.factor(par_df$type), y = par_df$Si, xlab = "Parameter type", 
     ylab = "First order Sobol' index estimate (Si)")

# Isolate numeric and categorical parameter Sis
numeric = par_df %>% filter(type == "numeric")
categorical = par_df %>% filter(type == "categorical")

# Find average and standard deviation of Si for numeric and categorical parameters
par_df %>% group_by(type) %>% summarise(mean_Si = mean(Si), 
                                        sd_Si = sd(Si), 
                                        median_Si = median(Si))

# Run Wilcoxon test to compare Si across parameter types
wilcox_type_Si = wilcox.test(x = numeric$Si, y = categorical$Si, alternative = "two.sided")

# Display results
wilcox_type_Si


# Test whether the shape of SIV curves influences Si

## Determine whether parameters are monotonically increasing, decreasing, or neither

### Correct spelling errors of a few parameter names
for (par in 1:nrow(par_df)) {
  if (par_df$parameter[par] == "avg.pct.veg.strmbank.allochtonous.SIV") {
    par_df$parameter[par] = "avg.pct.veg.strmbank.allochthonous.SIV"
  }
}

## Create column to store information about curve
par_df$curve_shape = NA

## Assign curve type to model using for loop
for (model in 1:length(HSImodels)) {
  breakpoints = HSImodels[[paste(HSImetadata[model,1])]] ## isolate breakpoint data for model
  breakpoints = breakpoints[,colSums(is.na(breakpoints)) < nrow(breakpoints)] ## select only columns with breakpoint data
  
  ### Loop through columns to determine the direction of curve
  for (col in seq(from = 2,to = length(breakpoints),by = 2)) {
    if (is.numeric(breakpoints[,col])) {
      if (all(diff(breakpoints[,col]) >= 0, na.rm = T)) {
        par_df[which(par_df$parameter == colnames(breakpoints)[col] & par_df$model 
                     == HSImetadata$model[model]), "curve_shape"] = "monotonically increasing"
      } else if (all(diff(breakpoints[,col]) <= 0, na.rm = T)) {
        par_df[which(par_df$parameter == colnames(breakpoints)[col] & par_df$model 
                     == HSImetadata$model[model]), "curve_shape"] = "monotonically decreasing"
      } else {
        par_df[which(par_df$parameter == colnames(breakpoints)[col] & par_df$model 
                     == HSImetadata$model[model]), "curve_shape"] = "variable"
      }
    } else {
      par_df[which(par_df$parameter == colnames(breakpoints)[col] & par_df$model 
                   == HSImetadata$model[model]), "curve_shape"] = "categorical"
    }
  }
}

par_df$curve_shape = as.factor(par_df$curve_shape) ## set as factor


## Determine whether curve shape of parameters influences sensitivity

## Constrain analysis to numeric variables only
numeric_par_df = par_df %>% filter(type == "numeric")

## Randomized anova for first order sensitivity (Si)
set.seed(3) ## Set seed for repeatable results

## Step 1: get and store observed F value
obs.aov = aov(Si ~ curve_shape, data = numeric_par_df)
summary(obs.aov)

obs.curve.F = summary(obs.aov)[[1]]$F[1]

## Step 2: set up variables for permutation loop
perms = 10000 ### number of permutations
perm.curve.F = rep(NA,perms) ### create empty file to place permuted F ratios in

## Step 3: run permutation
for (i in 1:perms) {
  Randomized.Resp = sample(numeric_par_df$Si,length(numeric_par_df$Si),replace = FALSE) ### randomly reorder data
  perm.aov = aov(Randomized.Resp ~ numeric_par_df$curve_shape) ### run anova on randomly ordered data
  perm.curve.F[i] = summary(perm.aov)[[1]]$F[1] ### store f value in empty file
}

## Step 4: compute p-values
p.curve.Si = mean(perm.curve.F > obs.curve.F) ### prob of obtaining F value as large or larger than observed F with randomly ordered data

## Step 5: print p-values
print(p.curve.Si)

## Post-hoc test (Table S1)
pairwise.wilcox.test(numeric_par_df$Si, numeric_par_df$curve_shape, p.adjust.method = "fdr")

# FDR adjustment of p values for statistical tests using Si
p_vals_Si = c(cor_test_breakpoints_Si$p.value, wilcox_type_Si$p.value, p.curve.Si)
adjusted_p_vals_Si = p.adjust(p_vals_Si, method = "fdr")
adjusted_p_vals_Si




# Run Spearman correlations to test whether Ti is correlated with parameter characteristics

# Plot relationship between breakpoints and Ti
plot(x = par_df$num_breakpoints, y = par_df$Ti, 
     xlab = "Number of parameter breakpoints", 
     ylab = "Total order Sobol' index estimate (Ti)")

# Run correlation test
cor_test_breakpoints_Ti = cor.test(par_df$Ti, par_df$num_breakpoints, method = "spearman")

# Display results
cor_test_breakpoints_Ti

# Plot Si for categorical versus numeric parameters
plot(x = as.factor(par_df$type), y = par_df$Ti, xlab = "Parameter type", 
     ylab = "Total order Sobol' index estimate (Ti)")

# Find average and standard deviation of Ti for numeric and categorical parameters
par_df %>% group_by(type) %>% summarise(mean_Ti = mean(Ti), 
                                        sd_Ti = sd(Ti), 
                                        median_Ti = median(Ti))

# Run Wilcoxon test to compare Si across parameter types
wilcox_type_Ti = wilcox.test(x = numeric$Ti, y = categorical$Ti, alternative = "two.sided")

# Display results
wilcox_type_Ti


# Test if curve shape matters for total order sensitivity (Ti)
## Randomized anova for total order sensitivity (Ti)
set.seed(3) ## Set seed for repeatable results

## Step 1: get and store observed F value
obs.aov = aov(Ti ~ curve_shape, data = numeric_par_df)
summary(obs.aov)

obs.curve.F = summary(obs.aov)[[1]]$F[1]

## Step 2: set up variables for permutation loop
perms = 10000 ### number of permutations
perm.curve.F = rep(NA,perms) ### create empty file to place permuted F ratios in

## Step 3: run permutation
for (i in 1:perms) {
  Randomized.Resp = sample(numeric_par_df$Ti,length(numeric_par_df$Ti),replace = FALSE) ### randomly reorder data
  perm.aov = aov(Randomized.Resp ~ numeric_par_df$curve_shape) ### run anova on randomly ordered data
  perm.curve.F[i] = summary(perm.aov)[[1]]$F[1] ### store f value in empty file
}

## Step 4: compute p-values
p.curve.Ti = mean(perm.curve.F > obs.curve.F) ### prob of obtaining F value as large or larger than observed F with randomly ordered data

## Step 5: print p-values
print(p.curve.Ti)

# Post-hoc test (Table S1)
pairwise.wilcox.test(numeric_par_df$Ti, numeric_par_df$curve_shape, p.adjust.method = "fdr")

# FDR adjustment of p values for statistical tests using Si
p_vals_Ti = c(cor_test_breakpoints_Ti$p.value, wilcox_type_Ti$p.value, p.curve.Ti)
adjusted_p_vals_Ti = p.adjust(p_vals_Ti, method = "fdr")
adjusted_p_vals_Ti



################################################################################

# Patterns of HSI score distribution and uncertainty in USFWS ecorest models

# Determine range, mean, and sd of median HSI scores for ecorest models
range(uncertainty_df$perc_50)
mean(uncertainty_df$perc_50)
sd(uncertainty_df$perc_50)

# Determine how many models have median scores falling below each USFWS rating
uncertainty_df %>% filter(perc_50 < 0.25) %>% summarise(count = n(), 
                                                        percent = (n() / 347) * 100)

uncertainty_df %>% filter(perc_50 < 0.5) %>% summarise(count = n(), 
                                                        percent = (n() / 347) * 100)

uncertainty_df %>% filter(perc_50 < 0.75) %>% summarise(count = n(), 
                                                        percent = (n() / 347) * 100)

uncertainty_df %>% filter(perc_50 < 1) %>% summarise(count = n(), 
                                                        percent = (n() / 347) * 100)


# Determine how many models have maximum scores falling below each USFWS rating
uncertainty_df %>% filter(max < 0.25) %>% summarise(count = n(), 
                                                        percent = (n() / 347) * 100)

uncertainty_df %>% filter(max < 0.5) %>% summarise(count = n(), 
                                                       percent = (n() / 347) * 100)

uncertainty_df %>% filter(max < 0.75) %>% summarise(count = n(), 
                                                        percent = (n() / 347) * 100)

uncertainty_df %>% filter(max < 1) %>% summarise(count = n(), 
                                                     percent = (n() / 347) * 100)


# Explore skewness of models

## Find range, mean, sd, and median of skewness values
range(uncertainty_df$skewness)
mean(uncertainty_df$skewness)
sd(uncertainty_df$skewness)
median(uncertainty_df$skewness)

## Find how many models were positively versus negatively skewed
uncertainty_df %>% filter(skewness > 0) %>% summarise(count = n(), percent = (n() / 347) * 100)
uncertainty_df %>% filter(skewness < 0) %>% summarise(count = n(), percent = (n() / 347) * 100)

# Explore kurtosis of models

## Find range, mean, sd, and median of kurtosis values
range(uncertainty_df$kurtosis)
mean(uncertainty_df$kurtosis)
sd(uncertainty_df$kurtosis)
median(uncertainty_df$kurtosis)

## Find how many models were platykurtic versus leptokurtic
uncertainty_df %>% filter(kurtosis < 3) %>% summarise(count = n(), percent = (n() / 347) * 100) ## platykurtic
uncertainty_df %>% filter(kurtosis > 3) %>% summarise(count = n(), percent = (n() / 347) * 100) ## leptokurtic



# Run Spearman correlation tests to investigate median HSI scores

## Plot relationship between number of parameters and median HSI scores
plot(x = uncertainty_df$num_params, y = uncertainty_df$perc_50, 
     xlab = "Number of input parameters", ylab = "Median HSI score")

## Run Spearman correlation
corr_test_num_params_HSI = cor.test(uncertainty_df$perc_50, uncertainty_df$num_params, method = "spearman")

## Display results
corr_test_num_params_HSI


## Plot relationship between number of components and median HSI scores
plot(x = uncertainty_df$num_components, y = uncertainty_df$perc_50, 
     xlab = "Number of input components", ylab = "Median HSI score")

## Run Spearman correlation
corr_test_num_comps_HSI = cor.test(uncertainty_df$perc_50, uncertainty_df$num_components, method = "spearman")

## Display results
corr_test_num_comps_HSI


## Plot relationship between proportion of categorical parameters and median HSI scores
plot(x = uncertainty_df$prop_categorical, y = uncertainty_df$perc_50, 
     xlab = "Proportion of input parameters that are categorical", ylab = "Median HSI score")

## Run Spearman correlation
corr_test_prop_cat_HSI = cor.test(uncertainty_df$perc_50, uncertainty_df$prop_categorical, method = "spearman")

## Display results
corr_test_prop_cat_HSI

# FDR adjustment of p values for statistical tests using median HSI
p_vals_median = c(corr_test_num_params_HSI$p.value, corr_test_num_comps_HSI$p.value, corr_test_prop_cat_HSI$p.value)
adjusted_p_vals_median = p.adjust(p_vals_median, method = "fdr")
adjusted_p_vals_median




# Run Spearman correlation tests to investigate skewness of HSI scores

## Plot relationship between number of parameters and skewness
plot(x = uncertainty_df$num_params, y = uncertainty_df$skewness, 
     xlab = "Number of input parameters", ylab = "Skewness of HSI score distributions")

## Run Spearman correlation
corr_test_num_params_skewness = cor.test(uncertainty_df$skewness, uncertainty_df$num_params, method = "spearman")

## Display results
corr_test_num_params_skewness


## Plot relationship between number of components and skewness of HSI scores
plot(x = uncertainty_df$num_components, y = uncertainty_df$skewness, 
     xlab = "Number of input components", ylab = "Skewness of HSI score distributions")

## Run Spearman correlation
corr_test_num_comps_skewness = cor.test(uncertainty_df$skewness, uncertainty_df$num_components, method = "spearman")

## Display results
corr_test_num_comps_skewness


## Plot relationship between proportion of categorical parameters and skewness of HSI scores
plot(x = uncertainty_df$prop_categorical, y = uncertainty_df$skewness, 
     xlab = "Proportion of input parameters that are categorical", ylab = "Skewness of HSI score distributions")

## Run Spearman correlation
corr_test_prop_cat_skewness = cor.test(uncertainty_df$skewness, uncertainty_df$prop_categorical, method = "spearman")

## Display results
corr_test_prop_cat_skewness

# FDR adjustment of p values for statistical tests using skewness
p_vals_skewness = c(corr_test_num_params_skewness$p.value, corr_test_num_comps_skewness$p.value, corr_test_prop_cat_skewness$p.value)
adjusted_p_vals_skewness = p.adjust(p_vals_skewness, method = "fdr")
adjusted_p_vals_skewness



# Run Spearman correlation tests to investigate kurtosis of HSI scores

## Plot relationship between number of parameters and kurtosis
plot(x = uncertainty_df$num_params, y = uncertainty_df$kurtosis, 
     xlab = "Number of input parameters", ylab = "Kurtosis of HSI score distributions")

## Run Spearman correlation
corr_test_num_params_kurtosis = cor.test(uncertainty_df$kurtosis, uncertainty_df$num_params, method = "spearman")

## Display results
corr_test_num_params_kurtosis


## Plot relationship between number of components and kurtosis of HSI scores
plot(x = uncertainty_df$num_components, y = uncertainty_df$kurtosis, 
     xlab = "Number of input components", ylab = "Kurtosis of HSI score distributions")

## Run Spearman correlation
corr_test_num_comps_kurtosis = cor.test(uncertainty_df$kurtosis, uncertainty_df$num_components, method = "spearman")

## Display results
corr_test_num_comps_kurtosis


## Plot relationship between proportion of categorical parameters and kurtosis of HSI scores
plot(x = uncertainty_df$prop_categorical, y = uncertainty_df$kurtosis, 
     xlab = "Proportion of input parameters that are categorical", ylab = "Kurtosis of HSI score distributions")

## Run Spearman correlation
corr_test_prop_cat_kurtosis = cor.test(uncertainty_df$kurtosis, uncertainty_df$prop_categorical, method = "spearman")

## Display results
corr_test_prop_cat_kurtosis

# FDR adjustment of p values for statistical tests using kurtosis
p_vals_kurtosis = c(corr_test_num_params_kurtosis$p.value, corr_test_num_comps_kurtosis$p.value, corr_test_prop_cat_kurtosis$p.value)
adjusted_p_vals_kurtosis = p.adjust(p_vals_kurtosis, method = "fdr")
adjusted_p_vals_kurtosis



# Run Spearman correlation tests to investigate coefficient of variation (cv) of HSI scores

## Plot relationship between number of parameters and cv
plot(x = uncertainty_df$num_params, y = uncertainty_df$cv, 
     xlab = "Number of input parameters", ylab = "Coefficient of variation of HSI scores")

## Run Spearman correlation
corr_test_num_params_cv = cor.test(uncertainty_df$cv, uncertainty_df$num_params, method = "spearman")

## Display results
corr_test_num_params_cv


## Plot relationship between number of components and cv of HSI scores
plot(x = uncertainty_df$num_components, y = uncertainty_df$cv, 
     xlab = "Number of input components", ylab = "Coefficient of variation of HSI score distributions")

## Run Spearman correlation
corr_test_num_comps_cv = cor.test(uncertainty_df$cv, uncertainty_df$num_components, method = "spearman")

## Display results
corr_test_num_comps_cv


## Plot relationship between proportion of categorical parameters and cv of HSI scores
plot(x = uncertainty_df$prop_categorical, y = uncertainty_df$cv, 
     xlab = "Proportion of input parameters that are categorical", ylab = "Coefficient of variation of HSI score distributions")

## Run Spearman correlation
corr_test_prop_cat_cv = cor.test(uncertainty_df$cv, uncertainty_df$prop_categorical, method = "spearman")

## Display results
corr_test_prop_cat_cv

# FDR adjustment of p values for statistical tests using cv
p_vals_cv = c(corr_test_num_params_cv$p.value, corr_test_num_comps_cv$p.value, corr_test_prop_cat_cv$p.value)
adjusted_p_vals_cv = p.adjust(p_vals_cv, method = "fdr")
adjusted_p_vals_cv

################################################################################

# Patterns of numerical approximation error in sensitivity analyses for ecorest models

# Find range, mean, sd, and median of upper bounds of 95% CIs for Si 
range(uncertainty_df$dummy_high_Si * 100)
mean(uncertainty_df$dummy_high_Si * 100)
sd(uncertainty_df$dummy_high_Si * 100)
median(uncertainty_df$dummy_high_Si * 100)

# Find range, mean, sd, and median of upper bounds of 95% CIs for Ti 
range(uncertainty_df$dummy_high_Ti * 100)
mean(uncertainty_df$dummy_high_Ti * 100)
sd(uncertainty_df$dummy_high_Ti * 100)
median(uncertainty_df$dummy_high_Ti * 100)

# Run Spearman correlations to examine relationship between UB of 95% CI for Si

# Plot relationship between Si UB and num_params
plot(x = uncertainty_df$num_params, y = uncertainty_df$dummy_high_Si,
     xlab = "Number of input parameters", ylab = "Upper bounds of 95% CI for Si estimates")

# Run Spearman correlation
cor_test_num_params_Sidummy = cor.test(uncertainty_df$num_params, uncertainty_df$dummy_high_Si, method = "spearman")

# Display results
cor_test_num_params_Sidummy


# Plot relationship between Si UB and num_components
plot(x = uncertainty_df$num_components, y = uncertainty_df$dummy_high_Si,
     xlab = "Number of input components", ylab = "Upper bounds of 95% CI for Si estimates")

# Run Spearman correlation
cor_test_num_comps_Sidummy = cor.test(uncertainty_df$num_components, uncertainty_df$dummy_high_Si, method = "spearman")

# Display results
cor_test_num_comps_Sidummy


# Plot relationship between Si UB and proportion of categorical parameters
plot(x = uncertainty_df$prop_categorical, y = uncertainty_df$dummy_high_Si,
     xlab = "Proportion of categorical parameters", ylab = "Upper bounds of 95% CI for Si estimates")

# Run Spearman correlation
cor_test_prop_cat_Sidummy = cor.test(uncertainty_df$prop_categorical, uncertainty_df$dummy_high_Si, method = "spearman")

# Display results
cor_test_prop_cat_Sidummy


# FDR adjustment of p values for statistical tests using Si UB
p_vals_Sidummy = c(cor_test_num_params_Sidummy$p.value, cor_test_num_comps_Sidummy$p.value, cor_test_prop_cat_Sidummy$p.value)
adjusted_p_vals_Sidummy = p.adjust(p_vals_Sidummy, method = "fdr")
adjusted_p_vals_Sidummy



# Run Spearman correlations to examine relationship between UB of 95% CI for Ti

# Plot relationship between Ti UB and num_params
plot(x = uncertainty_df$num_params, y = uncertainty_df$dummy_high_Ti,
     xlab = "Number of input parameters", ylab = "Upper bounds of 95% CI for Ti estimates")

# Run Spearman correlation
cor_test_num_params_Tidummy = cor.test(uncertainty_df$num_params, uncertainty_df$dummy_high_Ti, method = "spearman")

# Display results
cor_test_num_params_Tidummy


# Plot relationship between Ti UB and num_components
plot(x = uncertainty_df$num_components, y = uncertainty_df$dummy_high_Ti,
     xlab = "Number of input components", ylab = "Upper bounds of 95% CI for Ti estimates")

# Run Spearman correlation
cor_test_num_comps_Tidummy = cor.test(uncertainty_df$num_components, uncertainty_df$dummy_high_Ti, method = "spearman")

# Display results
cor_test_num_comps_Tidummy


# Plot relationship between Ti UB and proportion of categorical parameters
plot(x = uncertainty_df$prop_categorical, y = uncertainty_df$dummy_high_Ti,
     xlab = "Proportion of categorical parameters", ylab = "Upper bounds of 95% CI for Ti estimates")

# Run Spearman correlation
cor_test_prop_cat_Tidummy = cor.test(uncertainty_df$prop_categorical, uncertainty_df$dummy_high_Ti, method = "spearman")

# Display results
cor_test_prop_cat_Tidummy


# FDR adjustment of p values for statistical tests using Si UB
p_vals_Tidummy = c(cor_test_num_params_Tidummy$p.value, cor_test_num_comps_Tidummy$p.value, cor_test_prop_cat_Tidummy$p.value)
adjusted_p_vals_Tidummy = p.adjust(p_vals_Tidummy, method = "fdr")
adjusted_p_vals_Tidummy

################################################################################

# Influence of HSI model structure on model sensitivity

# Make structure an ordered factor
all_versions$structure = factor(all_versions$structure, 
                                levels = c("Original", "Arithmetic mean", 
                                           "Geometric mean", "Limiting factor", 
                                           "Multiplicative"))

# Summary statistics for percent of influential parameters across model structures (Table 4)
all_versions %>% group_by(structure) %>% summarise(min = min(percent_influential),
                                                   mean = mean(percent_influential),
                                                   median = median(percent_influential),
                                                   max = max(percent_influential),
                                                   sd = sd(percent_influential))


# Find number of models of each structure with zero influential parameters
all_versions %>% group_by(structure) %>% filter(percent_influential == 0) %>% summarise(count = n(), percent = (n() / 347) * 100)

# Friedman test to test whether structures differ in average percent of influential variables
friedman.test(y = all_versions$percent_influential, groups = all_versions$structure, blocks = all_versions$model)

# Post-hoc pairwise wilcoxon test with fdr correction
pairwise.wilcox.test(x = all_versions$percent_influential, all_versions$structure, p.adj = "fdr")

################################################################################

# Influence of HSI model structure on HSI score distribution and uncertainty

# Summary statistics for median HSI scores, skewness, and kurtosis of different model structures (Table 5)
## Median
all_versions %>% group_by(structure) %>% summarise(average = mean(perc_50),
                                                   sd = sd(perc_50),
                                                   median = median(perc_50),
                                                   min = min(perc_50),
                                                   max = max(perc_50))

## Skewness
all_versions %>% group_by(structure) %>% summarise(average = mean(skewness),
                                                   sd = sd(skewness),
                                                   median = median(skewness),
                                                   min = min(skewness),
                                                   max = max(skewness))

# Kurtosis
all_versions %>% group_by(structure) %>% summarise(average = mean(kurtosis),
                                                   sd = sd(kurtosis),
                                                   median = median(kurtosis),
                                                   min = min(kurtosis),
                                                   max = max(kurtosis))

# Friedman test to test whether structures differ in median HSI score
friedman.test(y = all_versions$perc_50, groups = all_versions$structure, blocks = all_versions$model)

# Post-hoc pairwise wilcoxon test with fdr correction (Table S3)
pairwise.wilcox.test(x = all_versions$perc_50, all_versions$structure, p.adj = "fdr")

# Graph median HSI scores across structures (Figure 5)
all_versions %>% group_by(structure) %>% ggplot(aes(x = structure, 
     y = perc_50, fill = structure, color = structure)) + geom_rain(alpha = .6, 
      boxplot.args = list(color = "black", outlier.shape = NA),
      boxplot.args.pos = list(
        width = 0.2, position = ggpp::position_dodgenudge(x = 0.2, width = 0.25)),
      violin.args.pos = list(
        side = "r",
        width = 1, position = position_nudge(x = 0.35))) + theme_classic() + scale_fill_manual(values = cbPalette) + scale_color_manual(values = cbPalette) + theme(
      legend.position="bottom") + labs(fill="") + guides(
        color = 'none') + scale_x_discrete(labels = label_wrap(10))+ ylab("Median HSI score") + xlab(
          "Model structure") + theme(text = element_text(color = "black", 
        size = 20)) + theme(axis.text = element_text(color = "black"))

# ggsave("Figure 5.png", dpi = 600, height = 5, width = 10)


# Friedman test to test whether structures differ in skewness
friedman.test(y = all_versions$skewness, groups = all_versions$structure, blocks = all_versions$model)

# Post-hoc pairwise wilcoxon test with fdr correction (Table S4)
pairwise.wilcox.test(x = all_versions$skewness, all_versions$structure, p.adj = "fdr")


# Friedman test to test whether structures differ in kurtosis
friedman.test(y = all_versions$kurtosis, groups = all_versions$structure, blocks = all_versions$model)

# Post-hoc pairwise wilcoxon test with fdr correction (Table S5)
pairwise.wilcox.test(x = all_versions$kurtosis, all_versions$structure, p.adj = "fdr")


# Examine patterns in  model uncertainty across structures using cv

# Summary statistics for cv (Table 6)
all_versions %>% group_by(structure) %>% summarise(min = min(cv),
                                                   mean = mean(cv),
                                                   max = max(cv),
                                                   sd = sd(cv))

# Friedman test to test whether structures differ in kurtosis
friedman.test(y = all_versions$cv, groups = all_versions$structure, blocks = all_versions$model)

# Post-hoc pairwise wilcoxon test with fdr correction (Table S6)
pairwise.wilcox.test(x = all_versions$cv, all_versions$structure, p.adj = "fdr")

################################################################################

# Patterns of numerical approximation error across model structures

# Summary statistics for upper bounds of 95% CI across model structures for Si
all_versions %>% group_by(structure) %>% summarise(min = min(dummy_high_Si * 100),
                                                   mean = mean(dummy_high_Si * 100),
                                                   sd = sd(dummy_high_Si * 100),
                                                   median = median(dummy_high_Si * 100),
                                                   max = max(dummy_high_Si * 100))


# Summary statistics for upper bounds of 95% CI across model structures for Ti
all_versions %>% group_by(structure) %>% summarise(min = min(dummy_high_Ti * 100),
                                                   mean = mean(dummy_high_Ti * 100),
                                                   sd = sd(dummy_high_Ti * 100),
                                                   median = median(dummy_high_Ti * 100),
                                                   max = max(dummy_high_Ti * 100))


# Friedman test to test whether structures differ in Si error
friedman.test(y = all_versions$dummy_high_Si, groups = all_versions$structure, blocks = all_versions$model)

# Post-hoc pairwise wilcoxon test with fdr correction (Table S7)
pairwise.wilcox.test(x = all_versions$dummy_high_Si, all_versions$structure, p.adj = "fdr")



# Friedman test to test whether structures differ in Ti error
friedman.test(y = all_versions$dummy_high_Ti, groups = all_versions$structure, blocks = all_versions$model)

# Post-hoc pairwise wilcoxon test with fdr correction (Table S8)
pairwise.wilcox.test(x = all_versions$dummy_high_Ti, all_versions$structure, p.adj = "fdr")



# Table S9 correlations between Si and different variables across structures

# Isolate each structure
original = all_versions %>% filter(structure == "Original")
arithmetic = all_versions %>% filter(structure == "Arithmetic mean")
geometric = all_versions %>% filter(structure == "Geometric mean")
limiting = all_versions %>% filter(structure == "Limiting factor")
multiplicative = all_versions %>% filter(structure == "Multiplicative")

## Original structure
params_confidence_Si_orig = cor.test(x = original$num_params, y = original$dummy_high_Si, method = "spearman")

comps_confidence_Si_orig = cor.test(x = original$num_components, y = original$dummy_high_Si, method = "spearman")

cat_confidence_Si_orig = cor.test(x = original$prop_categorical, y = original$dummy_high_Si, method = "spearman")

### FDR adjustment of p values for Ti tests
p_vals_Si_orig = c(params_confidence_Si_orig$p.value, comps_confidence_Si_orig$p.value, cat_confidence_Si_orig$p.value)
adjusted_p_vals_Si_orig = p.adjust(p_vals_Si_orig, method = "fdr")

### Summary of results
params_confidence_Si_orig
comps_confidence_Si_orig
cat_confidence_Si_orig

adjusted_p_vals_Si_orig

## Arithmetic mean structure
params_confidence_Si_ari = cor.test(x = arithmetic$num_params, y = arithmetic$dummy_high_Si, method = "spearman")

comps_confidence_Si_ari = cor.test(x = arithmetic$num_components, y = arithmetic$dummy_high_Si, method = "spearman")

cat_confidence_Si_ari = cor.test(x = arithmetic$prop_categorical, y = arithmetic$dummy_high_Si, method = "spearman")

### FDR adjustment of p values for Ti tests
p_vals_Si_ari = c(params_confidence_Si_ari$p.value, comps_confidence_Si_ari$p.value, cat_confidence_Si_ari$p.value)
adjusted_p_vals_Si_ari = p.adjust(p_vals_Si_ari, method = "fdr")

### Summary of results
params_confidence_Si_ari
comps_confidence_Si_ari
cat_confidence_Si_ari

adjusted_p_vals_Si_ari



## Geometric mean structure
params_confidence_Si_geo = cor.test(x = geometric$num_params, y = geometric$dummy_high_Si, method = "spearman")

comps_confidence_Si_geo = cor.test(x = geometric$num_components, y = geometric$dummy_high_Si, method = "spearman")

cat_confidence_Si_geo = cor.test(x = geometric$prop_categorical, y = geometric$dummy_high_Si, method = "spearman")

### FDR adjustment of p values for Ti tests
p_vals_Si_geo = c(params_confidence_Si_geo$p.value, comps_confidence_Si_geo$p.value, cat_confidence_Si_geo$p.value)
adjusted_p_vals_Si_geo = p.adjust(p_vals_Si_geo, method = "fdr")

### Summary of results
params_confidence_Si_geo
comps_confidence_Si_geo
cat_confidence_Si_geo

adjusted_p_vals_Si_geo



## Limiting factor structure
params_confidence_Si_limiting = cor.test(x = limiting$num_params, y = limiting$dummy_high_Si, method = "spearman")

comps_confidence_Si_limiting = cor.test(x = limiting$num_components, y = limiting$dummy_high_Si, method = "spearman")

cat_confidence_Si_limiting = cor.test(x = limiting$prop_categorical, y = limiting$dummy_high_Si, method = "spearman")

# FDR adjustment of p values for Ti tests
p_vals_Si_limiting = c(params_confidence_Si_limiting$p.value, comps_confidence_Si_limiting$p.value, cat_confidence_Si_limiting$p.value)
adjusted_p_vals_Si_limiting = p.adjust(p_vals_Si_limiting, method = "fdr")

# Summary of results
params_confidence_Si_limiting
comps_confidence_Si_limiting
cat_confidence_Si_limiting

adjusted_p_vals_Si_limiting



## Multiplicative model structure
params_confidence_Si_multiplicative = cor.test(x = multiplicative$num_params, y = multiplicative$dummy_high_Si, method = "spearman")

comps_confidence_Si_multiplicative = cor.test(x = multiplicative$num_components, y = multiplicative$dummy_high_Si, method = "spearman")

cat_confidence_Si_multiplicative = cor.test(x = multiplicative$prop_categorical, y = multiplicative$dummy_high_Si, method = "spearman")

# FDR adjustment of p values for Ti tests
p_vals_Si_multiplicative = c(params_confidence_Si_multiplicative$p.value, comps_confidence_Si_multiplicative$p.value, cat_confidence_Si_multiplicative$p.value)
adjusted_p_vals_Si_multiplicative = p.adjust(p_vals_Si_multiplicative, method = "fdr")

# Summary of results
params_confidence_Si_multiplicative
comps_confidence_Si_multiplicative
cat_confidence_Si_multiplicative

adjusted_p_vals_Si_multiplicative





# Spearman correlations for Table S9 for Ti

## Original structure
params_confidence_Ti_orig = cor.test(x = original$num_params, y = original$dummy_high_Ti, method = "spearman")

comps_confidence_Ti_orig = cor.test(x = original$num_components, y = original$dummy_high_Ti, method = "spearman")

cat_confidence_Ti_orig = cor.test(x = original$prop_categorical, y = original$dummy_high_Ti, method = "spearman")

### FDR adjustment of p values for Ti tests
p_vals_Ti_orig = c(params_confidence_Ti_orig$p.value, comps_confidence_Ti_orig$p.value, cat_confidence_Ti_orig$p.value)
adjusted_p_vals_Ti_orig = p.adjust(p_vals_Ti_orig, method = "fdr")

### Summary of results
params_confidence_Ti_orig
comps_confidence_Ti_orig
cat_confidence_Ti_orig

adjusted_p_vals_Ti_orig


## Arithmetic mean structure
params_confidence_Ti_ari = cor.test(x = arithmetic$num_params, y = arithmetic$dummy_high_Ti, method = "spearman")

comps_confidence_Ti_ari = cor.test(x = arithmetic$num_components, y = arithmetic$dummy_high_Ti, method = "spearman")

cat_confidence_Ti_ari = cor.test(x = arithmetic$prop_categorical, y = arithmetic$dummy_high_Ti, method = "spearman")

### FDR adjustment of p values for Ti tests
p_vals_Ti_ari = c(params_confidence_Ti_ari$p.value, comps_confidence_Ti_ari$p.value, cat_confidence_Ti_ari$p.value)
adjusted_p_vals_Ti_ari = p.adjust(p_vals_Ti_ari, method = "fdr")

### Summary of results
params_confidence_Ti_ari
comps_confidence_Ti_ari
cat_confidence_Ti_ari

adjusted_p_vals_Ti_ari



## Geometric mean structures
params_confidence_Ti_geo = cor.test(x = geometric$num_params, y = geometric$dummy_high_Ti, method = "spearman")

comps_confidence_Ti_geo = cor.test(x = geometric$num_components, y = geometric$dummy_high_Ti, method = "spearman")

cat_confidence_Ti_geo = cor.test(x = geometric$prop_categorical, y = geometric$dummy_high_Ti, method = "spearman")

### FDR adjustment of p values for Ti tests
p_vals_Ti_geo = c(params_confidence_Ti_geo$p.value, comps_confidence_Ti_geo$p.value, cat_confidence_Ti_geo$p.value)
adjusted_p_vals_Ti_geo = p.adjust(p_vals_Ti_geo, method = "fdr")

### Summary of results
params_confidence_Ti_geo
comps_confidence_Ti_geo
cat_confidence_Ti_geo

adjusted_p_vals_Ti_geo



## Limiting factor structure
params_confidence_Ti_limiting = cor.test(x = limiting$num_params, y = limiting$dummy_high_Ti, method = "spearman")

comps_confidence_Ti_limiting = cor.test(x = limiting$num_components, y = limiting$dummy_high_Ti, method = "spearman")

cat_confidence_Ti_limiting = cor.test(x = limiting$prop_categorical, y = limiting$dummy_high_Ti, method = "spearman")

### FDR adjustment of p values for Ti tests
p_vals_Ti_limiting = c(params_confidence_Ti_limiting$p.value, comps_confidence_Ti_limiting$p.value, cat_confidence_Ti_limiting$p.value)
adjusted_p_vals_Ti_limiting = p.adjust(p_vals_Ti_limiting, method = "fdr")

### Summary of results
params_confidence_Ti_limiting
comps_confidence_Ti_limiting
cat_confidence_Ti_limiting

adjusted_p_vals_Si_limiting



## Multiplicative model structure
params_confidence_Ti_multiplicative = cor.test(x = multiplicative$num_params, y = multiplicative$dummy_high_Ti, method = "spearman")

comps_confidence_Ti_multiplicative = cor.test(x = multiplicative$num_components, y = multiplicative$dummy_high_Ti, method = "spearman")

cat_confidence_Ti_multiplicative = cor.test(x = multiplicative$prop_categorical, y = multiplicative$dummy_high_Ti, method = "spearman")

### FDR adjustment of p values for Ti tests
p_vals_Ti_multiplicative = c(params_confidence_Ti_multiplicative$p.value, comps_confidence_Ti_multiplicative$p.value, cat_confidence_Ti_multiplicative$p.value)
adjusted_p_vals_Ti_multiplicative = p.adjust(p_vals_Ti_multiplicative, method = "fdr")

### Summary of results
params_confidence_Ti_multiplicative
comps_confidence_Ti_multiplicative
cat_confidence_Ti_multiplicative

adjusted_p_vals_Ti_multiplicative

# Create Figure 6 (multipanel scatterplot of relationship between Si/Ti and num_params)
fig_6a = all_versions %>% group_by(structure) %>% ggplot(aes(x = num_params, y = dummy_high_Si, 
 color = structure)) + geom_point() + facet_wrap(~structure, 
  ncol = 5, scale = "free_y") + labs( x = "Number of input parameters", y = "Upper bounds of 95% 
  CI for dummy parameter estimates", title = "A") + theme(text = element_text(color = "black", 
  size = 18, face = "bold")) + theme(axis.text = element_text(color = "black")
  ) + theme(axis.title.y = element_text(hjust=0.75)) + scale_fill_manual(
    values = cbPalette) + scale_color_manual(values = cbPalette)  + labs(color =
     "Structure") + theme(legend.position = "none") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + theme(strip.background =
        element_rect(colour="black", fill = "white"))  + theme(strip.text.x = element_text(size = 18, face = "bold"))

fig_6b = all_versions %>% group_by(structure) %>% ggplot(aes(x = num_params, y = dummy_high_Ti, 
 color = structure)) + geom_point() + facet_wrap(~structure, 
  ncol = 5, scale = "free_y") + labs( x = "Number of input parameters", y = "Upper bounds of 95% 
  CI for dummy parameter estimates", title = "B") + theme(text = element_text(color = "black", 
  size = 18, face = "bold")) + theme(axis.text = element_text(color = "black")
  ) + theme(axis.title.y = element_text(hjust=0.75)) + scale_fill_manual(
    values = cbPalette) + scale_color_manual(values = cbPalette)  + labs(color =
     "Structure") + theme(legend.position = "none") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ theme(strip.background =
        element_rect(colour="black", fill = "white")) + theme(strip.text.x = element_text(size = 18, face = "bold"))

Fig_6 = ggarrange(fig_6a, fig_6b, nrow = 2, ncol = 1)

# Display figure
Fig_6

# ggsave("Figure 6.png", Fig_6, dpi = 600, width=15.5, height=7, limitsize = F)


# Create Figure S1 (CV versus structure)
Fig_S1 = all_versions %>% group_by(structure) %>% ggplot(aes(x = structure, 
     y = cv, fill = structure, color = structure)) + geom_rain(alpha = .6, 
      boxplot.args = list(color = "black", outlier.shape = NA),
      boxplot.args.pos = list(
        width = 0.2, position = ggpp::position_dodgenudge(x = 0.2, width = 0.25)),
      violin.args.pos = list(
        side = "r",
        width = 1, position = position_nudge(x = 0.35))) + theme_classic() + scale_fill_manual(values = cbPalette) + scale_color_manual(values = cbPalette) + theme(
      legend.position="bottom") + labs(fill="") + guides(
        color = 'none') + scale_x_discrete(labels = label_wrap(10))+ ylab("Coefficient of variation (CV)") + xlab(
          "Model structure") + theme(text = element_text(color = "black", 
        size = 20)) + theme(axis.text = element_text(color = "black"))


# ggsave("Figure S1.png", Fig_S1, dpi = 600, height = 5, width = 10)



# Create Fig S2 (Si/Ti dummy parameters versus num_components)
fig_s2a = all_versions %>% group_by(structure) %>% ggplot(aes(x = num_components, y = dummy_high_Si, 
 color = structure)) + geom_point() + facet_wrap(~structure, 
  ncol = 5, scale = "free_y") + labs( x = "Number of input components", y = "Upper bounds of 95% 
  CI for dummy parameter estimates", title = "A") + theme(text = element_text(color = "black", 
  size = 18)) + theme(axis.text = element_text(color = "black")
  ) + theme(axis.title.y = element_text(hjust=0.75)) + scale_fill_manual(
    values = cbPalette) + scale_color_manual(values = cbPalette)  + labs(color =
     "Structure") + theme(legend.position = "none") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + theme(strip.background =element_rect(colour="black", fill = "white"))

fig_s2b = all_versions %>% group_by(structure) %>% ggplot(aes(x = num_components, y = dummy_high_Ti, 
 color = structure)) + geom_point() + facet_wrap(~structure, 
  ncol = 5, scale = "free_y") + labs( x = "Number of input components", y = "Upper bounds of 95% 
  CI for dummy parameter estimates", title = "B") + theme(text = element_text(color = "black", 
  size = 18)) + theme(axis.text = element_text(color = "black")
  ) + theme(axis.title.y = element_text(hjust=0.75)) + scale_fill_manual(
    values = cbPalette) + scale_color_manual(values = cbPalette)  + labs(color =
     "Structure") + theme(legend.position = "none") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ theme(strip.background =element_rect(colour="black", fill = "white"))

Fig_S2 = ggarrange(fig_s2a, fig_s2b, nrow = 2, ncol = 1)

# Display graph
Fig_S2

###ggsave("Figure S2.png", Fig_s2, dpi = 600, width=15.5, height=7, limitsize = F)

