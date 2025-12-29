# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare abundance detected)
# Data:    Modelled and extracted data.
# Task:    Extract the mean and SE of 1000 sampling designs
# Author:  Lise Fournier-Carnoy
# Date:    September 2024

# -----------------------------------------------------------------------------

# Clear memory
rm(list=ls())

# Load libraries
library(raster)         # Dealing with rasters
library(sp)             # Dealing with shapefiles
library(sf)             # Dealing with shapefiles
library(terra)          # Dealing with shapefiles
library(tidyverse)      # Data manipulation
library(ggstatsplot)    # For plotting with stats integrated
library(ggtext)         # for linebreaks in plot text
library(patchwork)      # For plot arranging
library(terra)          # for manipulating rasters to calculate realised abundance
library(sf)             # for manipulating shapefiles to calculate realised abundance
library(geepack)        # for significance testing clustered designs
library(glmmTMB)        # for modelling
library(scales)         # for gradient_n_pal, in plots

set.seed(123)

## Custom plotting parameters -------------------------------------------------

colour_palette <- eval(parse(text = readLines("chapter_colours.txt")))
source("custom_theme.R")

crs <- 7850
crs_rast <- "epsg:7850"


## Load data ------------------------------------------------------------------

study_site <- "waatern" # waatern or waatu

# Define multiple species
spp_list <- c("Sparidae_Chrysophrys_auratus", 
              "Glaucosomatidae_Glaucosoma_hebraicum", 
              "Labridae_Ophthalmolepis_lineolatus")

all_data <- list()
for (spp in spp_list) {
  
  # file paths for extracted values
  file_simp_spabal          <- paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_simple_spatially_balanced_", spp, "_extracted_values_edge_effect.rds")
  
  file_str_spabal_2525_val  <- paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_stratified_spatially_balanced_2525_",  spp, "_extracted_values_edge_effect.rds")
  file_str_spabal_2030_val  <- paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_stratified_spatially_balanced_2030_",  spp, "_extracted_values_edge_effect.rds")
  
  file_pref_2525_val        <- paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_preferential_2525_", spp, "_extracted_values_edge_effect.rds")
  file_pref_2030_val        <- paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_preferential_2030_", spp, "_extracted_values_edge_effect.rds")
  
  file_clump_val            <- paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_clustered_", spp, "_extracted_values_edge_effect.rds")
  
  # file paths for predicted abundance raster
  file_raster               <- paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_", spp, ".rds")
  
  # read the extracted value files
  simp_spabal               <- readRDS(file_simp_spabal) %>% mutate(species = spp)
  
  str_spabal_2525           <- readRDS(file_str_spabal_2525_val) %>% mutate(species = spp)
  str_spabal_2030           <- readRDS(file_str_spabal_2030_val) %>% mutate(species = spp)
  
  pref_2525                 <- readRDS(file_pref_2525_val) %>% mutate(species = spp)
  pref_2030                 <- readRDS(file_pref_2030_val) %>% mutate(species = spp) 
  
  clump                     <- readRDS(file_clump_val) %>% mutate(species = spp)
  
  # transform crs
  target_crs                <- st_crs(str_spabal_2525)
  
  simp_spabal               <- st_transform(simp_spabal, target_crs)
  str_spabal_2525           <- st_transform(str_spabal_2525, target_crs)
  str_spabal_2030           <- st_transform(str_spabal_2030, target_crs)
  pref_2525                 <- st_transform(pref_2525, target_crs)
  pref_2030                 <- st_transform(pref_2030, target_crs)
  clump                     <- st_transform(clump, target_crs)
  
  # clump has an extra column which we want to keep, make NAs into the other dataframes
  str_spabal_2525$cluster_id <- NA
  str_spabal_2030$cluster_id <- NA
  pref_2525$cluster_id       <- NA
  pref_2030$cluster_id       <- NA
  simp_spabal$cluster_id     <- NA
  
  # combine extracted values for the current species
  species_data <- rbind(simp_spabal, str_spabal_2030, str_spabal_2525, pref_2030, pref_2525, clump)
  
  # append only extracted values to the list
  all_data[[spp]] <- species_data
  
  # Load the predicted abundance raster (but don't store it in the list)
  if (file.exists(file_raster)) {  
    assign(paste0("raster_", spp), rast(file_raster), envir = .GlobalEnv)  # Assign raster to global environment
  } else {
    warning(paste("Raster file not found:", file_raster))
  }
}

# Combine all species data into a single dataframe
dat <- do.call(rbind, all_data) %>%
  glimpse()
summary(dat$fit_values)
## Select an abundance value within the normal distribution -------------------

# For each point sampled from the fake abundance map, I want to select a random abundance within the normal distribution of mean +/- SE
# For this I've made a loop that selects a point within the random distribution (mean +/- SE) of each sampling point

# Create a new column to store the randomly selected values of abundance
dat$random_abundance <- NA
set.seed(12345)
# commented out as it takes a bit of time and has already been run. re-run if any models change after 16/09/2025
for(i in 1:nrow(dat)) {
  print(paste(i, "out of", nrow(dat)))
  
  # get the mean and standard deviation for the current row
  mean_value <- dat$fit_value[i]
  sd_value <- dat$se_value[i]

  random_abundance <- rnorm(1, mean = mean_value, sd = sd_value) # generate a random value from the normal distribution
  random_abundance <- pmax(random_abundance, 0) # replace negative values with 0
  dat$random_abundance[i] <- random_abundance # assign the random value to the new column
} # loop to select a random abundance observed for each sampling point.
glimpse(dat)
summary(dat$random_abundance)

saveRDS(dat, paste0('outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_', study_site, '_sampling_design_points_all_spp_edge_effect.rds'))
dat <- readRDS(paste0('outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_', study_site, '_sampling_design_points_all_spp_edge_effect.rds'))


## Calculate 'true realised abundance' ----------------------------------------

# To compare the observed abundance, we need to know what the actual difference
# inside and outside the NTZ is. Even though we uniformly increased the abundance
# By 80%, it might not be exactly 1.8x because of spatial heterogeneity

sim_sz                <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ.shp")) %>% st_transform(crs)
sampling_area         <- readRDS(paste0("QGIS layers/produced_from_code/D04_", study_site, "_sampling_area_deeper_than_7m.rds")) %>% st_transform(crs)

species_lookup <- data.frame( # lookup table so that different dataframes recognise the species names
  species_formatted = c(
    "Yijarup/Pink Snapper",
    "Djubitj/West Australian Dhufish",
    "Southern Maori Wrasse"
  ),
  species_raw = c(
    "Sparidae_Chrysophrys_auratus",
    "Glaucosomatidae_Glaucosoma_hebraicum",
    "Labridae_Ophthalmolepis_lineolatus"
  )
)


# Initialize an empty list to store realised abundance ratios - commented out as it's been run and saved.
realised_ab_ratio_list <- list()
mean_ab_inside_list <- list()
mean_ab_outside_list <- list()

# Loop over each species in the species list
for (spp in spp_list) {

  # Construct the file path for the raster of the current species
  file_raster <- paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_smooth_edge_effect_", spp, ".rds")

  # Check if the raster file exists and load it
  if (file.exists(file_raster)) {
    # Load the raster for the species
    spp_raster <- rast(file_raster) %>% project(crs_rast)

    # Crop the raster to the smapling area and to the SZ
    masked_ras <- mask(spp_raster, sampling_area)
    masked_ras_sz <- mask(spp_raster, sim_sz)
    outside_area <- st_difference(sampling_area, sim_sz)
    masked_ras_out_sz <- mask(spp_raster, outside_area)

    # Calculate the mean realised abundance inside and outside the NTZ
    mean_abundance_inside <- global(masked_ras_sz, fun = "mean", na.rm = TRUE)
    mean_abundance_outside <- global(masked_ras_out_sz, fun = "mean", na.rm = TRUE)

    # Calculate the realised abundance ratio
    realised_abundance_ratio <- mean_abundance_inside[1,] / mean_abundance_outside[1,]
    mean_ab_inside_list[[spp]] <- mean_abundance_inside[1,]
    mean_ab_outside_list[[spp]] <- mean_abundance_outside[1,]

    # Store the result in the list
    realised_ab_ratio_list[[spp]] <- realised_abundance_ratio
    print(paste("Realised abundance ratio for", spp, ":", realised_abundance_ratio))

  } else {
    warning(paste("Raster file not found for", spp, ":", file_raster))
  }
}

realised_abundance_df <- data.frame(
  species_raw = spp_list,
  realised_abundance_ratio = unlist(realised_ab_ratio_list),
  mean_abundance_inside = unlist(mean_ab_inside_list),
  mean_abundance_outside = unlist(mean_ab_outside_list)
)

print(realised_abundance_df)
saveRDS(realised_abundance_df, paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_realised_abundances_smooth_edge_effect.rds"))
realised_abundance_df <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_realised_abundances_smooth_edge_effect.rds"))

## Calculate ratios -----------------------------------------------------------

glimpse(dat)

ratio <- dat %>%
  group_by(species, SD, design_id, in_SZ) %>%
  summarise(mean_abundance = mean(random_abundance, na.rm = TRUE))

ratio2 <- st_drop_geometry(ratio) %>%
  pivot_wider(names_from = in_SZ, values_from = mean_abundance) %>%
  mutate(observed_abundance_ratio = `TRUE` / `FALSE`)
unique(ratio2$SD)
ratio3 <- ratio2 %>% # this takes forever because there are almost a million rows to recode.
  mutate(SD = recode(SD,
                     "simple_spabal" = "Si_SB",
                     "stratified_spabal_2030" = "St_SB_2030",
                     "stratified_spabal_2525" = "St_SB_2525",
                     "pref_2030" = "P_2030",
                     "pref_2525" = "P_2525",
                     "clustered" = "C"),
         SD = factor(SD, levels = c("Si_SB", "St_SB_2030", "St_SB_2525", "P_2030", "P_2525", "C")),  # Reorder here
         ) %>%
  rename(
    mean_abundance_outside_SZ = `FALSE`,
    mean_abundance_inside_SZ = `TRUE`
  ) %>%
  glimpse()

ratio <- ratio3 %>%
  left_join(realised_abundance_df, by = c("species" = "species_raw")) %>%
  glimpse()

names(ratio) <- c("species", "SD", "design_id",
                  "mean_obs_ab_outside", "mean_obs_ab_inside", "obs_ab_ratio",
                  "real_ab_ratio", "real_ab_inside", "real_ab_outside")

saveRDS(ratio, paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_ratios_edge_effect.rds"))
ratio <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_ratios_edge_effect.rds"))


## Calculate performance metrics ----------------------------------------------

mean_bias <- ratio %>%
  group_by(design_id, SD, species) %>%
  summarise(
    mean_bias = mean(obs_ab_ratio - real_ab_ratio, na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  glimpse()

rmse <- ratio %>%
  group_by(design_id, SD, species) %>%
  summarise(
    rmse = sqrt(mean((obs_ab_ratio - real_ab_ratio)^2, na.rm = TRUE))
  ) %>% 
  ungroup() %>% 
  glimpse()

interval_coverage <- ratio %>%
  filter(is.finite(obs_ab_ratio)) %>% 
  group_by(SD, species) %>%
  summarise(
    sd_ratio = sd(obs_ab_ratio, na.rm = TRUE),
    n = sum(!is.na(obs_ab_ratio)),
    int_coverage = mean(abs(obs_ab_ratio - real_ab_ratio) <= (1.96 * sd_ratio), na.rm = TRUE),.groups = "drop") %>% 
  glimpse()

saveRDS(mean_bias, paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_mean_bias_edge_effect.rds"))
saveRDS(rmse, paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_RMSE_edge_effect.rds"))
saveRDS(interval_coverage, paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_interval_coverage_edge_effect.rds"))


## Plot raw data --------------------------------------------------------------

realised_long <- realised_abundance_df %>%
  dplyr::select("species" = "species_raw", #species_formatted = species_formatted.x,
                SZ = mean_abundance_inside, fished = mean_abundance_outside) %>%
  pivot_longer(cols = c(SZ, fished),
               names_to = "status",
               values_to = "real_ab_mean") %>% 
  glimpse()

ratio_long <- ratio %>%
  pivot_longer(cols = c(mean_obs_ab_inside, mean_obs_ab_outside),
               names_to = "status", 
               values_to = "obs_ab_mean") %>% 
  mutate(status = ifelse(status == "mean_obs_ab_inside", "SZ", "fished")) %>% 
  dplyr::select(-c(obs_ab_ratio, real_ab_ratio)) %>% 
  glimpse()

saveRDS(ratio_long, paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_raw_data_edge_effect.rds"))

## Plot raw data --------------------------------------------------------------

ggplot(ratio_long, aes(x = SD, y = obs_ab_mean, fill = status)) +
  geom_boxplot(width = 1, color = "black", alpha = 0.5, position = position_dodge(0.9), outlier.size = 0.5) +
  
  # Add horizontal dashed lines for realised abundance
  geom_hline(data = realised_long,
             aes(yintercept = real_ab_mean,
                 color = status,
                 group = interaction(species, status)),
             linetype = "dashed",
             size = 0.8) +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    fill = "Status",
    x = "Sampling Design",
    y = "Observed Abundance"
  ) +
  scale_fill_manual(values = colour_palette[1:2], labels = c("fished", "NTZ")) +
  scale_color_manual(values = colour_palette[1:2], guide = "none") +  # hide redundant legend
  #custom_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, hjust = 1))


## Plot performance -----------------------------------------------------------

gradient_fun <- scales::gradient_n_pal(colour_palette[4:6])
SD_levels <- levels(mean_bias$SD)
n_levels <- length(SD_levels)
gradient_colors <- gradient_fun(seq(0, 1, length.out = n_levels))
names(gradient_colors) <- SD_levels  # Name the colors by SD level


p1 <- ggplot(mean_bias, aes(x = SD, y = mean_bias, fill = SD)) +
  geom_boxplot(position = position_dodge(0.9), outlier.size = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = colour_palette[1], size = 1) +  
  facet_grid(cols = vars(species), scales = "free_y", switch = "y") +
  labs(y = "Mean Bias") +
  custom_theme +
  theme(
    axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank()  # Removes x-axis title
  ) +
  scale_fill_manual(values = gradient_colors) +
  scale_color_manual(values = c("white" = colour_palette[4], "black" = "black")) +
  coord_cartesian(ylim = c(NA, 2))

p2 <- ggplot(rmse, aes(x = SD, y = rmse, fill = SD)) +
  geom_boxplot(position = position_dodge(0.9), outlier.size = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = colour_palette[1], size = 1) +  
  facet_grid(cols = vars(species)) +
  theme_minimal() +
  labs(y = "RMSE") +
  custom_theme +
  theme(
    axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank()  # Removes x-axis title
  ) +
  scale_fill_manual(values = gradient_colors) +
  scale_color_manual(values = c("white" = colour_palette[4], "black" = "black")) +
  scale_y_continuous(limits = c(0, NA)) +  # Cut off values below 0 - rmse cannot be negative, but the violin plot makes it look like there are neg values even though there aren't
  coord_cartesian(ylim = c(0, 2))

p3 <- ggplot(interval_coverage, aes(x = SD, y = int_coverage, fill = SD)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~species) +
  labs(y = "Coverage") +
  custom_theme +
  theme(
    strip.text.x = element_blank(),  # Removes column facet titles (species names)
    axis.title.x = element_blank()  # Removes x-axis title
  ) +
  geom_hline(data = interval_coverage, 
             aes(yintercept = 0.95), linetype = "dashed", color = colour_palette[1], size = 1) +  
  scale_fill_manual(values = gradient_colors) +
  coord_cartesian(ylim = c(0.5, 1))


# Combine Plots
p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1))  

## 'p-value significance' rate workflow ---------------------------------------

glimpse(dat)

make_formula <- function(SD, include_cluster, include_strata) {

  fixed <- "in_SZ"
  random_str <- if (include_strata) "+ (1|strata)" else ""
  random_clus <- if (include_cluster) "+ (1|cluster_id)" else ""

  as.formula(paste("random_abundance ~", fixed, random_str, random_clus))
} # this function creates a formula depending on the sampling design

fit_glmmtmb_model <- function(data, formula, family) { # glmmTMB
  tryCatch({
    glmmTMB(formula,
            data = data,
            family = family)
  }, error = function(e) NULL)
} # this function fits a glmmTMB model to the data, using the formula specified with the previous function

test_dat <- st_drop_geometry(dat) %>%
  group_by(species, SD, design_id) %>%
  group_split()

models <- list()
results <- list()

for (sim in test_dat) {
  print(paste0(unique(sim$design_id), "__", unique(sim$species)))
  i <- length(results) + 1
  
  # 1. build the formula
  SD <- unique(sim$SD)
  include_cluster <- grepl("clustered", SD)
  include_strata <- grepl("stratified", SD)
  formula <- make_formula(SD, include_cluster, include_strata)

  # 2. fit the model
  model_tmb <- fit_glmmtmb_model(data = sim, formula = formula, family = tweedie(link = "log"))

  # 3. extract information 
  species <- unique(sim$species)
  design_id <- unique(sim$design_id)
  estimate_tmb <- NA
  p_value_tmb <- NA
  message_tmb <- NA

  if (!is.null(model_tmb)) {
    coef_tmb <- summary(model_tmb)$coefficients$cond
    if ("in_SZTRUE" %in% rownames(coef_tmb)) {
      estimate_tmb <- coef_tmb["in_SZTRUE", "Estimate"]
      p_value_tmb <- coef_tmb["in_SZTRUE", "Pr(>|z|)"]
    }
    optinfo <- model_tmb$fit$message
    if (!is.null(optinfo)) {
      message_tmb <- paste("Convergence code:", optinfo)
    }
  } # this extracts estimates and p-values if the model has converged

  models[[i]] <- if (!is.null(model_tmb)) summary(model_tmb) else NULL
  
  results[[i]] <- data.frame(
    species = species,
    SD = SD,
    design_id = design_id,
    model = "glmmTMB",
    estimate = estimate_tmb,
    p_value = p_value_tmb,
    message = message_tmb,
    converged = model_tmb$sdr$pdHess,
    stringsAsFactors = FALSE
  )
}
glmmTMB_results <- do.call(rbind, results)

results_clean <- unique(glmmTMB_results) %>% filter(converged == TRUE) # only keep models that converge well

results_clean$SD <- factor(results_clean$SD, levels = c("simple_spabal",
                                                      "stratified_spabal_2030",
                                                      "stratified_spabal_2525",
                                                      "pref_2030",
                                                      "pref_2525",
                                                      "clustered"))

saveRDS(results_clean, paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_significance_testing_model_results_edge_effect.rds"))
sig <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_significance_testing_model_results_edge_effect.rds"))

# Plotting

summary_stats <- sig %>%
  group_by(species, SD) %>%
  summarise(
    total = n(),
    under_0.05 = sum(p_value < 0.05, na.rm = TRUE),
    perc_under_0.05 = round(100 * under_0.05 / total, 2)
  ) %>%
  arrange(desc(perc_under_0.05))

# p value
ggplot(sig, aes(x = SD, y = p_value)) +
  geom_boxplot(outlier.size = 0.5, varwidth = T) +
  facet_wrap(~ species) +
  geom_text(
    data = summary_stats,
    aes(x = SD, y = 1.05,
        label = paste0("n=", total, "\nsig.=", round(perc_under_0.05, 1), "%")),
    size = 3,
    angle = 0,
    hjust = 1,
    inherit.aes = FALSE
  ) +
  labs(
    title = "",
    x = "",
    y = "P-value for status"
  ) +
  theme_minimal() +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.01, ymax = 0.05), fill = NA, color = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
  coord_flip()


glimpse(realised_abundance_df)

test <- realised_abundance_df %>%
  mutate(ra = mean_abundance_outside/mean_abundance_inside ) %>% 
  rename(species = "species_raw")

test <- sig %>%
  left_join(test, by = c("species"))

# Plot
ggplot(test, aes(x = SD, y = estimate)) +
  geom_boxplot(outlier.shape = NA, varwidth = TRUE) +
  facet_wrap(~ species, scales="free") +
  geom_hline(aes(yintercept = ra), color = "red", linetype = "dashed") +
  labs(
    x = "",
    y = "Estimate for SZ effect"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_cartesian(ylim = c(0.3, 1.25))

### END ###