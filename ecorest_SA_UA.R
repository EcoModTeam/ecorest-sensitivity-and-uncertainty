# Batch code for running sensitivity analysis using ecorest models
# Authors: Kiara Cushway, Colton Shaw, Kyle McKay, Todd Swannack

## Clear workspace
rm(list = ls()) ## remove stored files and objects
gc(T) ## garbage collection
graphics.off() ## Turn off graphics

## Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Install packages if needed
###install.packages("here")
###install.packages("ecorest")
###install.packages("sensobol")
###install.packages("ggplot2")
###install.packages("stringr")
###install.packages("extraDistr")
###install.packages("tidyr")
###install.packages("dplyr")
###install.packages("patchwork")
###install.packages("future.apply")
###install.packages("progressr")


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
library(future.apply) ## used for parallel computation
library(progressr) ## used to track progress

# future.scheduling = structure(TRUE, ordering = "random")

plan(multisession, workers = parallel::detectCores() - 1)

HSImetadata = HSImetadata
HSImodels = HSImodels


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
n = 60000 ## initial base sample size
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

################################################################################

# Create adaptive function to determine base sample size
adaptive_n_step <- function(n,
                            ci_widths,
                            ci_tol,
                            min_step = 1000,
                            max_step = 25000,
                            max_factor = 2.5) {
  
  ratio <- max(ci_widths / ci_tol, na.rm = TRUE)
  
  step <- if (ratio <= 1.5) {
    min_step
  } else if (ratio <= 3) {
    ceiling(n * 0.25)
  } else {
    ceiling(n * 0.5)
  }
  
  step <- min(max(step, min_step), max_step)
  step <- min(step, ceiling(n * max_factor))
  
  return(step)
}

################################################################################

run_model_SA <- function(i,
                         HSImetadata,
                         HSImodels,
                         n_start = 60000,
                         ci_tol = 0.05,
                         R = 1000,
                         checkpoint_dir = "checkpoints",
                         do_plots = FALSE) {
  
  model_name <- HSImetadata[i, 1]
  ckpt_file  <- file.path(checkpoint_dir, paste0("model_", model_name, ".rds"))
  
  ## ---- Skip if already done ----
  if (file.exists(ckpt_file)) {
    message("Skipping completed model: ", model_name)
    return(readRDS(ckpt_file))
  }
  
  dir.create(checkpoint_dir, showWarnings = FALSE)
  
  ## ---- Model-specific setup (done ONCE) ----
  breakpoints <- HSImodels[[model_name]]
  Comp_exist  <- HSImetadata[i,41:(ncol(HSImetadata)-1)]
  comp        <- Comp_exist[, colSums(is.na(Comp_exist)) == 0]
  
  SIV_exist <- HSImetadata[i,9:40]
  siv       <- SIV_exist[, colSums(is.na(SIV_exist)) == 0]
  params    <- unlist(siv, use.names = FALSE)
  M         <- length(params)
  
  breakpoints <- breakpoints[, colSums(is.na(breakpoints)) < nrow(breakpoints)]
  
  ## Detect categorical variables ONCE
  cat_info <- lapply(seq_len(M), function(s) {
    is_cat <- !is.numeric(breakpoints[, s*2 - 1])
    list(
      is_cat = is_cat,
      n_cat  = if (is_cat) sum(!is.na(breakpoints[, s*2 - 1])) else NA
    )
  })
  cat_names <- params[vapply(cat_info, `[[`, logical(1), "is_cat")]
  num_cat   <- length(cat_names)
  
  ## ---- Convergence loop ----
  n <- n_start
  run <- 1
  converged <- rep(FALSE, 5)
  names(converged) <- c("orig", "ari", "geo", "mult", "limit")
  
  while (!all(converged)) {
    N <- n * (M + 2)
    mat <- sobol_matrices(N = n, params = params, order = "first", type = "QRN")
    
    ## Fill matrix
    for (s in seq_len(M)) {
      if (cat_info[[s]]$is_cat) {
        mat[, s] <- qdunif(mat[, s], 1, cat_info[[s]]$n_cat)
      } else {
        if (breakpoints[1,s*2] == 0 && breakpoints[2,s*2] == 0) {
          ## Store the minimum and maximum values within the non-zero range of the SIV
          is_zero = na.omit(breakpoints[,s*2]) == 0 # test whether each value is zero
          keep_min <- is_zero & (c(FALSE, is_zero[-length(is_zero)])) # for min: keep second zero of adjacent pairs
          # Find minimum value
          mins = min(breakpoints[,s*2-1][keep_min], na.rm = T)
        } else {
          mins = min(breakpoints[,s*2-1],na.rm = TRUE) ## store min value for each SIV
        }
        if (breakpoints[length(na.omit(breakpoints[,s*2])),s*2] == 0 && breakpoints[(length(na.omit(breakpoints[,s*2])) - 1),s*2] == 0) {
          ## Store the minimum and maximum values within the non-zero range of the SIV
          is_zero = na.omit(breakpoints[,s*2]) == 0 # test whether each value is zero
          keep_max <- is_zero & (c(is_zero[-length(is_zero)], FALSE)) # for max: keep first zero of adjacent pairs
          # Find  maximum value
          maxs = max(breakpoints[,s*2-1][keep_max], na.rm = T)
        } else {
          maxs = max(breakpoints[,s*2-1],na.rm = TRUE) ## store max value for each SIV
        }
        mat[, s] <- qunif(mat[, s], min = mins, max = maxs)
      }
    }
    
    ## Replace categorical codes with letters
    for (c in seq_along(params)) {
      if (params[c] %in% cat_names) {
        mat[, c] <- letters[as.numeric(mat[, c])]
      }
    }
    
    ## ---- Vectorized SIV + HSI evaluation ----
    SIV_scores <- t(apply(mat, 1, function(x) {
      SIcalc(breakpoints, cbind(x))
    }))
    
    # Rowwise HSI calculation
    HSI <- t(apply(SIV_scores, 1, function(x) {
      c(
        HSIeqtn(model_name, x, HSImetadata),   # main HSI
        HSIarimean(x),                         # arithmetic mean
        HSIgeomean(x),                         # geometric mean
        HSImult(x),                             # multiplicative
        HSIlimit(x)                             # limiting factor
      )
    }))
    
    colnames(HSI) = c("HSI","HSI_Arithmetic","HSI_Geometric", "HSI_Multiplicative", "HSI_Limiting")
    
    # Get quantiles
    quantiles <- t(apply(HSI, 2, function(x) {
      quantile(x, probs = c(0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99, 1))
    }))
    row.names(quantiles) = c("orig", "ari", "geo", "mult", "limit")
    
    ## ---- Sobol indices ----
    params = lapply(params,as.character) ## make SIV names characters
    call_params = do.call(rbind.data.frame, params) ## paste SIV names in dataframe
    
    ind_list <- list(
      orig  = sobol_indices(Y = as.vector(HSI[,1]), N = n, params = call_params[,1], first = first, total = total, boot = TRUE, R = R, type = "percent", conf = conf),
      ari   = sobol_indices(Y = as.vector(HSI[,2]), N = n, params = call_params[,1], first = first, total = total, boot = TRUE, R = R, type = "percent", conf = conf),
      geo   = sobol_indices(Y = as.vector(HSI[,3]), N = n, params = call_params[,1], first = first, total = total, boot = TRUE, R = R, type = "percent", conf = conf),
      mult  = sobol_indices(Y = as.vector(HSI[,4]), N = n, params = call_params[,1], first = first, total = total, boot = TRUE, R = R, type = "percent", conf = conf),
      limit = sobol_indices(Y = as.vector(HSI[,5]), N = n, params = call_params[,1], first = first, total = total, boot = TRUE, R = R, type = "percent", conf = conf)
    )
    
    for (k in names(ind_list)) {
      ci_width <- with(ind_list[[k]]$results, high.ci - low.ci)
      converged[k] <- all(ci_width <= ci_tol)
    }
    
    ind_dummy_list <- list(
      orig  = sobol_dummy(Y = as.vector(HSI[,1]), N = n, params = params, boot = TRUE,R = R),
      ari   = sobol_dummy(Y = as.vector(HSI[,2]), N = n, params = params, boot = TRUE,R = R),
      geo   = sobol_dummy(Y = as.vector(HSI[,3]), N = n, params = params, boot = TRUE,R = R),
      mult  = sobol_dummy(Y = as.vector(HSI[,4]), N = n, params = params, boot = TRUE,R = R),
      limit = sobol_dummy(Y = as.vector(HSI[,5]), N = n, params = params, boot = TRUE,R = R)
    )
    
    ## ---- Save partial progress ----
    saveRDS(
      list(
        model = model_name,
        run = run,
        n = n,
        converged = converged
      ),
      file = file.path(checkpoint_dir, paste0("progress_", model_name, ".rds"))
    )
    
    if (!all(converged)) {
      
      ## Collect all CI widths
      all_ci <- unlist(lapply(ind_list, function(x) {
        with(x$results, high.ci - low.ci)
      }))
      
      step <- adaptive_n_step(
        n = n,
        ci_widths = all_ci,
        ci_tol = ci_tol
      )
      
      n <- n + step
      run <- run + 1
    }
    
    if (n > 300000) {
      warning("Model ", model_name, " failed to converge by n = ", n)
      break
    }
    
  }
  
  ## Store data for model in RDS file
  out <- list(
    name = model_name,
    runs = ind_list$orig$C,
    inputs = M,
    components = length(comp),
    cat_inputs = num_cat,
    n = n,
    ind = ind_list,
    HSI = HSI,
    quantiles = quantiles,
    dummy = ind_dummy_list,
    converged = converged
  )
  
  saveRDS(out, ckpt_file)
  return(out)
}

# Apply function
handlers(global = TRUE)

with_progress({
  p <- progressor(along = seq(nrow(HSImetadata)))
  
  results <- future_lapply(seq(nrow(HSImetadata)), function(i) {
    res <- run_model_SA(i, HSImetadata, HSImodels)
    p(sprintf("Model %s complete", HSImetadata[i,1]))
    res
  }, future.seed = 349, future.scheduling = structure(TRUE, ordering = "random"))
})

