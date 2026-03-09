# Statistical analysis for ecorest SA/UA
# Authors: Kiara Cushway, Colton Shaw, Kyle McKay, Todd Swannack

## Clear workspace
rm(list = ls()) ## remove stored files and objects
gc(T) ## garbage collection
graphics.off() ## Turn off graphics

setwd("E:/ecorest/R Code/Final code")

# Load packages
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
library(purrr) ## for combining RDS files
library(paletteer) ## for heatmap gradient

load("HSImetadata.RData") ## load updated HSImetadata
load("HSImodels.RData") ## load updated HSImodels

# Bring in data with information about model structures
eqt_types = read.csv("equation_types.csv")

## Load HSImetadata into data frame, excluding the two models with only one variable
HSImetadata = HSImetadata %>% filter(!model %in% c("redwingedblackbirdB", "woodduckWinter"))

## Remove the two models with one variable from HSImodels as well
HSImodels$redwingedblackbirdB = NULL
HSImodels$woodduckWinter = NULL

# Colorblind friendly color palette for graphs
cbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Path to folder containing RDS files
folder_path = "E:/ecorest/R Code/Final code/checkpoints"

# Get file paths
file_paths = list.files(path = folder_path,
                        pattern = "^model_.*\\.rds$", # Filter for files ending with .rds
                        full.names = TRUE)    # Get the full path


# Store all RDS files in a list
models = lapply(file_paths, readRDS)

# Name objects in the list based on model name
names(models) <- vapply(
  models,
  function(x) x$name,
  character(1)
)

# Check which models converged
failed_runs <- lapply(models, function(m) {
  which(!m$converged)
})

summary_df <- data.frame(
  model = names(models),
  failed_variations = sapply(failed_runs, function(x)
    if(length(x) == 0) NA else paste(names(x), collapse = ", ")
  )
)

# Find models that have original forms that failed to converge
summary_df %>% filter(str_detect(failed_variations, "orig"))

# Find the range in total evaluations and base sample size needed to converge
n <- vapply(models, function(x) x$n, numeric(1))

min(n) ## minimum number of runs
max(n) ## maximum number of runs

# Find the range in total evaluations and base sample size needed to converge
runs <- vapply(models, function(x) x$runs, numeric(1))

min(runs) ## minimum number of runs
max(runs) ## maximum number of runs

# Exclude models that did not converge
models$woodduckYear = NULL

# Find the range in number of parameters across models
params_count <- vapply(models, function(x) x$inputs, numeric(1))

min(params_count) ## minimum number of parameters
max(params_count) ## maximum number of parameters

# Find the range in number of components across models
comps_count <- vapply(models, function(x) x$components, numeric(1))

min(comps_count) ## minimum number of parameters
max(comps_count) ## maximum number of parameters

# Find the number of models with at least one categorical parameter
cat_count <- vapply(models, function(x) x$cat_inputs, numeric(1))
sum(cat_count > 0) 
sum(cat_count > 0) / length(models) * 100

# Find the number of models with no numeric parameters
inputs = data.frame(params_count, cat_count)
length(models) - sum(inputs$params_count == inputs$cat_count)
(length(models) - sum(inputs$params_count == inputs$cat_count)) / length(models) * 100


# Store parameter information from original models in a data frame for analysis
results_df <- map_dfr(
  models,
  ~ .x$ind$orig$results,
  .id = "model"
)

# Store information about the dummy variables
dummy_df <- map_dfr(
  models,
  ~ .x$dummy$orig,
  .id = "model"
)

# Add 'dummy' to column names
dummy_df <- dummy_df %>%
  rename_with(
    ~ paste0("dummy_", .x),
    .cols = c(original:high.ci)
  )

# Combine results and dummy data frames
params_df <- results_df %>%
  left_join(
    dummy_df[, 1:7],
    by = c("model", "sensitivity")
  )

# Make it so each parameter has a single row
params_wide <- params_df %>%
  pivot_wider(
    id_cols = c(model, parameters),
    names_from = sensitivity,
    values_from = c(
      original, bias, std.error, low.ci, high.ci,
      dummy_original, dummy_bias, dummy_std.error,
      dummy_low.ci, dummy_high.ci
    ),
    names_glue = "{sensitivity}_{.value}"
  )

# Determine whether each parameter is influential or not
params_wide <- params_wide %>%
  rowwise() %>%
  mutate(
    influential = if_else(
      (Si_low.ci <= Si_dummy_high.ci) & (Ti_low.ci <= Ti_dummy_high.ci),
      0L,
      1L
    )
  ) %>%
  ungroup()

# Count number of models with non-influential parameters
length(unique(params_wide$model[params_wide$influential == 0]))

# Count percent of models with non-influential parameters
length(unique(params_wide$model[params_wide$influential == 0])) / length(models)

# Find the average number of non-influential parameters per model
model_params = params_wide %>% 
  group_by(model) %>% 
  summarise(num_inf = sum(influential == 1), total = n()) %>% 
  mutate(perc_noninf = 100 - (num_inf / total * 100)) 

mean(model_params$perc_noninf)
sd(model_params$perc_noninf)

# Count number of models with less than half of parameters influential
model_params %>% summarise(half_noninf = sum(perc_noninf >= 50) / length(models) * 100)

# Determine which model structure has the highest proportion of models with non-influential parameters
model_params = model_params %>% left_join(eqt_types, by = 'model')

model_params$all_influential = NULL ## initialize new column

# Assign yes to models where all parameters are influential
for (row in 1:nrow(model_params)) {
  if (model_params$perc_noninf[row] > 0) {
    model_params$all_influential[row] = "No"
  } else {
    model_params$all_influential[row] = "Yes"
  }
}

# Graph using pie chart
plot_models <- model_params %>%
  count(eqtn_type_final, all_influential) %>%
  group_by(eqtn_type_final) %>%
  mutate(prop = n / sum(n))

# Add up total number of models
totals <- model_params %>%
  count(eqtn_type_final) %>%
  mutate(label = paste0("n = ", n))

labeller_eqtn_type = labeller(
  eqtn_type_final = c(
    "arithmetic mean" = paste0("Arithmetic mean ", "\n", totals[1, "label"]),
    "author-specified" = paste0("Author-specified ", "\n", totals[2, "label"]),
    "geometric mean" = paste0("Geometric mean ", "\n", totals[3, "label"]),
    "limiting factor" = paste0("Limiting factor ", "\n", totals[4, "label"]),
    "multiplicative" = paste0("Multiplicative ", "\n", totals[5, "label"]),
    "weighted arithmetic mean" = paste0("Weighted arithmetic", "\n", "mean ", "\n", totals[6, "label"]),
    "weighted geometric mean" = paste0("Weighted geometric", "\n", "mean ", "\n", totals[7, "label"])
  )
)

plot_models$eqtn_type_final <- factor(plot_models$eqtn_type_final, 
                                      levels = c("author-specified", 
                                                 "limiting factor", 
                                                 "multiplicative",
                                                 "geometric mean",
                                                 "weighted geometric mean",
                                                 "arithmetic mean",
                                                 "weighted arithmetic mean"))

ggplot(plot_models, aes(x = "", y = prop, fill = all_influential)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~eqtn_type_final, nrow = 1, labeller = labeller_eqtn_type) +
  scale_fill_manual(values = cbPalette) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_stack(vjust = 0.5),
            size = 4) +
  theme_void() + theme(text = element_text(size = 16, face = "bold"))  +
  theme(legend.position = "bottom") + labs(fill = "") + scale_fill_discrete(
    palette = cbPalette, labels = c("Models with non-influential parameters", 
                                    "Models with no non-influential parameters"))

ggsave("Figure S1.png", dpi = 600)

# Test whether the percent of non-influential parameters increased as total parameters increased
## Plot against each other
plot(x = model_params$total, y = model_params$perc_noninf, 
     xlab = "Number of input parameters", ylab = "Percent of total parameters that are non-influential")

## Run Spearman rank correlation
cor_test_num_params = cor.test(model_params$perc_noninf, model_params$total, method = "spearman")

## Display results of correlation test
cor_test_num_params

# Test whether the percent of non-influential parameters increased as total components increased
comps_count = data.frame(comps_count)
comps_count$model = row.names(comps_count)

model_params <- model_params %>%
  left_join(
    comps_count,
    by = c("model")
  )

## Plot against each other
plot(x = model_params$comps_count, y = model_params$perc_noninf, 
     xlab = "Number of input components", ylab = "Percent of total parameters that are non-influential")

## Run Spearman rank correlation
cor_test_num_comps = cor.test(model_params$perc_noninf, model_params$comps_count, method = "spearman")

## Display results of correlation test
cor_test_num_comps

# FDR adjustment of p values for correlation tests using percent_influential
p_vals_influential = c(cor_test_num_params$p.value, cor_test_num_comps$p.value)
adjusted_p_vals_influential = p.adjust(p_vals_influential, method = "fdr")
adjusted_p_vals_influential

# Create figure 4: scatter plot of percent influential versus num_params, num_components, and prop_categorical
# Num_params graph
fig_4a = model_params %>% ggplot(aes(x = total, y = perc_noninf)
                  ) + geom_point(size = 3, color = "#E69F00", position = "jitter") + labs (y = "Percent of non-influential input parameters", 
  x = "Total number of input parameters", title = "A") + theme(text = element_text(color = "black", 
  size = 18)) + theme(axis.text = element_text(color = "black")
  ) + theme(axis.title.y = element_text(hjust=0.75)) + theme(legend.position = "none")  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_cartesian(ylim = c(0, 100))

# Num_components graph
fig_4b = model_params %>% ggplot(aes(x = comps_count, y = perc_noninf)
                  ) + geom_point(size = 3, color = "#56B4E9", position = "jitter") + labs (y = "", 
  x = "Total number of input components", title = "B") + theme(text = element_text(color = "black", 
  size = 18)) + theme(axis.text = element_text(color = "black")
  ) + theme(axis.title.y = element_text(hjust=0.75)) + theme(legend.position = "none")  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_cartesian(ylim = c(0, 100))


# Combine graphs
fig_4 = ggarrange(fig_4a, fig_4b, ncol = 2)

# Display graph
fig_4

# Save figure as a png file
### ggsave("Figure_5.png", dpi = 600)

# Find the maximum number of input parameters with no non-influential parameters
model_params %>% group_by(all_influential) %>% summarise(max_inputs = max(total))


# Explore how input parameter complexity influences sensitivity

for (row in 1:nrow(params_wide)) {
  curves = HSImodels[[params_wide$model[row]]]
  param = curves[, match(params_wide$parameters[row], names(curves))]
  if (is.numeric(curves[, match(params_wide$parameters[row], names(curves)) - 1])) {
    params_wide$type[row] = "continuous"
    if (all(diff(param) >= 0, na.rm = T)) {
      params_wide$curve_shape[row] = "monotonic"
    } else if (all(diff(param) <= 0, na.rm = T)) {
      params_wide$curve_shape[row] = "monotonic"
    } else {
      params_wide$curve_shape[row] = "variable"
    }
  } else {
    params_wide$curve_shape[row] = NA
    params_wide$type[row] = "categorical"
  }
  params_wide$num_breakpoints[row] = sum(complete.cases(param))
}

# Test whether Si is influenced by number of breakpoints
## Plot against each other
plot(x = params_wide$num_breakpoints, y = params_wide$Si_original, 
     xlab = "Number of breakpoints", ylab = "First order sensitivity")

## Run Spearman rank correlation
cor_test_num_breaks = cor.test(params_wide$Si_original, params_wide$num_breakpoints, method = "spearman")

## Display results of correlation test
cor_test_num_breaks

# Test whether Si is influenced by curve shape
## Plot against each other
plot(x = as.factor(params_wide$curve_shape), y = params_wide$Si_original, 
     xlab = "Curve shape", ylab = "First order sensitivity")

# Run Wilcoxon test
t.test.curve = t.test(params_wide$Si_original ~ as.factor(params_wide$curve_shape), na.rm = T)

# Display results
t.test.curve

# Test whether categorical parameters are more influential than continuous parameters
## Plot against each other
plot(x = as.factor(params_wide$type), y = params_wide$Si_original, 
     xlab = "Parameter type", ylab = "First order sensitivity")

# Run Wilcoxon test
t.test.type = t.test(params_wide$Si_original ~ as.factor(params_wide$type), na.rm = T)

# Display results
t.test.type

# FDR adjustment of p values for correlation tests using percent_influential
p_vals_influential = c(cor_test_num_breaks$p.value, t.test.curve$p.value, t.test.type$p.value)
adjusted_p_vals_influential = p.adjust(p_vals_influential, method = "fdr")
adjusted_p_vals_influential

################################################################################

# Uncertainty in original models
for (model in names(models)) {
  model_params[match(model, model_params$model), "minimum"] = min(models[[model]]$HSI[, 1])
  model_params[match(model, model_params$model), "maximum"] = max(models[[model]]$HSI[, 1])
  model_params[match(model, model_params$model), "median"] = median(models[[model]]$HSI[, 1])
  model_params[match(model, model_params$model), "box1"] = mean(models[[model]]$HSI[, 1] <= 0.1)
  model_params[match(model, model_params$model), "box2"] = mean(models[[model]]$HSI[, 1] > 0.1 & models[[model]]$HSI[, 1] <= 0.2)
  model_params[match(model, model_params$model), "box3"] = mean(models[[model]]$HSI[, 1] > 0.2 & models[[model]]$HSI[, 1] <= 0.3)
  model_params[match(model, model_params$model), "box4"] = mean(models[[model]]$HSI[, 1] > 0.3 & models[[model]]$HSI[, 1] <= 0.4)
  model_params[match(model, model_params$model), "box5"] = mean(models[[model]]$HSI[, 1] > 0.4 & models[[model]]$HSI[, 1] <= 0.5)
  model_params[match(model, model_params$model), "box6"] = mean(models[[model]]$HSI[, 1] > 0.5 & models[[model]]$HSI[, 1] <= 0.6)
  model_params[match(model, model_params$model), "box7"] = mean(models[[model]]$HSI[, 1] > 0.6 & models[[model]]$HSI[, 1] <= 0.7)
  model_params[match(model, model_params$model), "box8"] = mean(models[[model]]$HSI[, 1] > 0.7 & models[[model]]$HSI[, 1] <= 0.8)
  model_params[match(model, model_params$model), "box9"] = mean(models[[model]]$HSI[, 1] > 0.8 & models[[model]]$HSI[, 1] <= 0.9)
  model_params[match(model, model_params$model), "box10"] = mean(models[[model]]$HSI[, 1] > 0.9 & models[[model]]$HSI[, 1] <= 1)
}

# Test how many models have ranges that span from 0 to 1
for (row in 1:nrow(model_params)) {
  if (isTRUE(all.equal(c(0, 1), unlist(unname(c(model_params[row, "minimum"], model_params[row, "maximum"]))), tolerance = 0.05))) {
    model_params$entire_range[row] = "Yes"
  } else {
    model_params$entire_range[row] = "No"
  }
}
# Test how many models score 1
for (row in 1:nrow(model_params)) {
  if (isTRUE(all.equal(1, unlist(unname(model_params[row, "maximum"])), tolerance = 0.05))) {
    model_params$one[row] = "Yes"
  } else {
    model_params$one[row] = "No"
  }
}

model_params %>% group_by(entire_range) %>% summarise(n = n()) %>% mutate(perc = n / nrow(model_params) * 100)

# Create heatmap of model scores
heatmap_data = model_params %>% select(c(model, box1:box10)) # subset data

box_labels <- paste(seq(0, 0.9, by = 0.1), seq(0.1, 1, by = 0.1), sep = " – ")

colnames(heatmap_data)[-1] = box_labels

heatmap_data_long <- heatmap_data %>% pivot_longer(cols = -model, names_to = "box", values_to = "proportion")

ggplot(heatmap_data_long, aes(x = model, y = box, fill = proportion)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = rev(paletteer_c("grDevices::YlOrRd", 10)),
    limits = c(0, 1)
  ) +
  scale_y_discrete(limits = rev) +
  theme_minimal(base_size = 25) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  labs(
    y = "HSI scores",
    fill = "Proportion of scores",
    x = "Model"
  ) + theme(legend.position = "bottom") + theme(
    legend.key.width = unit(3.5, "cm")
  ) + theme(
    legend.title = element_text(
      size = 25,
      margin = margin(b = 30)  # 5 points below title
    )
  )

ggsave("Fig_6.png", dpi = 600)

# Median scores
range(model_params$median)
mean(model_params$median)
sd(model_params$median)

################################################################################

# Compare models when structured differently

# Remove all models that are already limiting factor, multiplicative, arithmetic, or geometric
include = model_params %>% 
subset(!eqtn_type_final %in% c('limiting factor', 'multiplicative', 
                               'arithmetic mean', 'geometric mean')) %>% select(model)

# Remove models that did not converge
not_converged = summary_df %>% filter(!is.na(failed_variations))

include = include %>% filter(!model %in% not_converged$model)

# Define the versions you want to process
versions <- c("orig", "ari", "geo", "mult", "limit")

# Loop over versions and combine results
all_versions <- map_dfr(versions, function(version) {
  
  # Filter models to only those in 'include'
  models_sub <- models[names(models) %in% include$model]
  
  # 1. Extract parameter results for this version
  results_df <- map_dfr(
    models_sub,
    ~ .x$ind[[version]]$results,
    .id = "model"
  )
  
  # 2. Extract dummy variable info for this version
  dummy_df <- map_dfr(
    models_sub,
    ~ .x$dummy[[version]],
    .id = "model"
  )
  
  # 3. Rename dummy columns
  dummy_df <- dummy_df %>%
    rename_with(
      ~ paste0("dummy_", .x),
      .cols = c(original:high.ci)
    )
  
  # 4. Combine results and dummy data
  params_df <- results_df %>%
    left_join(
      dummy_df[, 1:7],
      by = c("model", "sensitivity")
    )
  
  # 5. Pivot wider to get one row per parameter
  params_wide <- params_df %>%
    pivot_wider(
      id_cols = c(model, parameters),
      names_from = sensitivity,
      values_from = c(
        original, bias, std.error, low.ci, high.ci,
        dummy_original, dummy_bias, dummy_std.error,
        dummy_low.ci, dummy_high.ci
      ),
      names_glue = "{sensitivity}_{.value}"
    )
  
  # 6. Determine influential parameters
  params_wide <- params_wide %>%
    rowwise() %>%
    mutate(
      influential = if_else(
        (Si_low.ci <= Si_dummy_high.ci) & (Ti_low.ci <= Ti_dummy_high.ci),
        0L,
        1L
      )
    ) %>%
    ungroup()
  
  # 7. Compute model-level summary
  model_params <- params_wide %>%
    group_by(model) %>%
    summarise(
      num_inf = sum(influential == 1),
      total = n(),
      perc_noninf = 100 - (num_inf / total * 100),
      .groups = "drop"
    ) %>%
    left_join(eqt_types, by = "model") %>%
    mutate(
      all_influential = if_else(perc_noninf > 0, "No", "Yes"),
      version = version,
      # Replace eqtn_type_final with version name for all included models (optional)
      eqtn_type_final = version
    )
  
  return(model_params)
})

# Summary statistics
all_versions %>% group_by(eqtn_type_final) %>% summarise(min = min(perc_noninf),
                                                   mean = mean(perc_noninf),
                                                   median = median(perc_noninf),
                                                   max = max(perc_noninf),
                                                   sd = sd(perc_noninf))

all_versions %>% group_by(eqtn_type_final, all_influential) %>% summarise(n = n())

# Friedman test to test whether structures differ in average percent of non-influential variables
friedman.test(y = all_versions$perc_noninf, groups = all_versions$eqtn_type_final, blocks = all_versions$model)

# Post-hoc pairwise wilcoxon test with fdr correction
pairwise.wilcox.test(x = all_versions$perc_noninf, all_versions$eqtn_type_final, p.adj = "fdr")

# Graph using pie chart
plot_models_all <- all_versions %>%
  count(eqtn_type_final, all_influential) %>%
  group_by(eqtn_type_final) %>%
  mutate(prop = n / sum(n))

# Add up total number of models
totals_all <- all_versions %>%
  count(eqtn_type_final) %>%
  mutate(label = paste0("n = ", n))

labeller_eqtn_type_all = labeller(
  eqtn_type_final = c(
    "ari" = paste0("Arithmetic mean ", "\n", totals_all[1, "label"]),
    "orig" = paste0("Original ", "\n", totals_all[5, "label"]),
    "geo" = paste0("Geometric mean ", "\n", totals_all[2, "label"]),
    "limit" = paste0("Limiting factor ", "\n", totals_all[3, "label"]),
    "mult" = paste0("Multiplicative ", "\n", totals_all[4, "label"])
  )
)

plot_models_all$eqtn_type_final <- factor(plot_models_all$eqtn_type_final, 
                                      levels = c("orig", 
                                                 "limit", 
                                                 "mult",
                                                 "geo",
                                                 "ari"))

ggplot(plot_models_all, aes(x = "", y = prop, fill = all_influential)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~eqtn_type_final, nrow = 1, labeller = labeller_eqtn_type_all) +
  scale_fill_manual(values = cbPalette) +
  #geom_text(aes(label = scales::percent(prop, accuracy = 1)),
           # position = position_stack(vjust = 0.5),
           # size = 6) +
  theme_void() + theme(text = element_text(size = 25, face = "bold"))  +
  theme(legend.position = "bottom") + labs(fill = "") + scale_fill_discrete(
    palette = cbPalette, labels = c("Models with non-influential parameters", 
                                    "Models with no non-influential parameters"))

ggsave("Figure_7.png", dpi = 600)

# Examine model uncertainty across structures
all_versions$minimum = NA
all_versions$median = NA
all_versions$maximum = NA

version_num = c(1:5)
versions = data.frame(versions, version_num)

for (model in all_versions$model) {
  for (version in versions$versions) {
    row_ind <- which(
      all_versions$model == model &
        all_versions$eqtn_type_final == version
    )
    
    all_versions[row_ind, "minimum"] <- min(models[[model]]$HSI[, versions[match(version, versions$versions), "version_num"]])
    all_versions[row_ind, "maximum"] <- max(models[[model]]$HSI[,  versions[match(version, versions$versions), "version_num"]])
    all_versions[row_ind, "median"]  <- median(models[[model]]$HSI[,  versions[match(version, versions$versions), "version_num"]])
  }
}

# Test how many models have ranges that span from 0 to 1
for (row in 1:nrow(all_versions)) {
  if (isTRUE(all.equal(c(0, 1), unlist(unname(c(all_versions[row, "minimum"], all_versions[row, "maximum"]))), tolerance = 0.05))) {
    all_versions$entire_range[row] = "Yes"
  } else {
    all_versions$entire_range[row] = "No"
  }
}
# Test how many models score 1
for (row in 1:nrow(all_versions)) {
  if (isTRUE(all.equal(1, unlist(unname(all_versions[row, "maximum"])), tolerance = 0.05))) {
    all_versions$one[row] = "Yes"
  } else {
    all_versions$one[row] = "No"
  }
}

# Find the percent of each model type that did not span the entire range of HSI scores
all_versions %>% group_by(eqtn_type_final, entire_range) %>% summarise(perc = n() / (nrow(all_versions) / 5) * 100)

# Friedman test to test whether structures differ in median HSI score
friedman.test(y = all_versions$median, groups = all_versions$eqtn_type_final, blocks = all_versions$model)

# Post-hoc pairwise wilcoxon test with fdr correction (Table S3)
pairwise.wilcox.test(x = all_versions$median, all_versions$eqtn_type_final, p.adj = "fdr")

all_versions %>% group_by(eqtn_type_final) %>% summarise(mean_median = mean(median))

# Make structure an ordered factor
all_versions$eqtn_type_final = factor(all_versions$eqtn_type_final, 
                                levels = c("orig", "ari", 
                                           "geo", "limit", 
                                           "mult"))

letters_df <- all_versions %>%
  group_by(eqtn_type_final) %>%
  summarise(
    y_pos = 1.05  # adjust spacing if needed
  ) %>%
  arrange(eqtn_type_final) %>%
  mutate(sig_letter = c("a", "b", "c", "d", "e"))  # replace with your actual letters

# Graph median HSI scores across structures (Figure 5)
all_versions %>% group_by(eqtn_type_final) %>% ggplot(aes(x = eqtn_type_final, 
     y = median, fill = eqtn_type_final, color = eqtn_type_final)) + geom_rain(alpha = .6, 
      boxplot.args = list(color = "black", outlier.shape = NA),
      boxplot.args.pos = list(
        width = 0.2, position = ggpp::position_dodgenudge(x = 0.2, width = 0.25)),
      violin.args.pos = list(
        side = "r",
        width = 1, position = position_nudge(x = 0.35))) + theme_classic() + scale_fill_manual(values = cbPalette) + scale_color_manual(values = cbPalette) + theme(
      legend.position="none") + labs(fill="") + guides(
        color = 'none') + scale_x_discrete(labels = label_wrap(10))+ ylab("Median HSI score") + xlab(
          "Model structure") + theme(text = element_text(color = "black", 
        size = 20)) + theme(axis.text = element_text(color = "black"))  +
  scale_x_discrete(labels=c( "Original", "Arithmetic \n mean", 
                             "Geometric \n mean", 
                             "Limiting \n factor", 
                             "Multiplicative")) + geom_text(data = letters_df,
                                                            aes(x = eqtn_type_final, y = y_pos, label = sig_letter),
                                                            inherit.aes = FALSE,
                                                            size = 7,
                                                            nudge_x = -0.05)

ggsave("Figure_8.png", dpi = 600, height = 5, width = 10)
