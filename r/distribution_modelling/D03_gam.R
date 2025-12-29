# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare Lengths detected)
# Data:    2024 MEGlab BRUV habitat data.
# Task:    Fit GAM to the abundance of large fish.
# Author:  Lise Fournier-Carnoy
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  Need to add the final data and select the right models. Reef_lidar was removed because prediction is terrible

# -----------------------------------------------------------------------------

library(tidyverse) # for general data manipulation
library(mgcv) # for GAMs
library(MuMIn) # for GAMs
library(FSSgam) # for full-subset gams
library(CheckEM) # for EventMeasure data
library(patchwork) # for plot layouts
library(terra) # for spatial object manipulation 
library(raster) # for raster manipulation
library(broom) # for model outputs in Rmarkdown
library(stargazer) # for model outputs as table

# Clear memory
rm(list=ls())

study_site <- "waatu" # waatern or waatu

ext_plot <- if(study_site == "waatern") {
  ext <- c(115.32, 115.5426, -33.60944, -33.40019) 
} else { 
  ext <- c(114.92, 115.035, -34.18, -34.02)
} # Extent to crop layers to - this is a crop to the simulated SZ area


## Files used in this script --------------------------------------------------

file_mature_pres      <- paste0("data/tidy/D01_", study_site, "_presence_latlong.rds")
file_all_predictors   <- paste0("data/tidy/D02_", study_site, "_all_predictors.rds") # must be a terra SpatRaster

spp <- c("Glaucosomatidae Glaucosoma hebraicum", "Labridae Ophthalmolepis lineolatus", "Sparidae Chrysophrys auratus" )
model_spp <- spp[2]


## Load data ------------------------------------------------------------------

## Load abundance data #-------------------------------------------------------

# Extracting predictors values to presence points
pres <- readRDS(file_mature_pres) %>% 
  mutate(longitude = as.numeric(longitude),
         latitude = as.numeric(latitude)) %>%
  glimpse()

predictors <- readRDS(file_all_predictors); predictors <- predictors[[names(predictors) != "aspect_250m"]]
names(predictors) <- c("bathy", #"aspect",
                       "roughness", "detrended",
                       "reef", "sand", "seagrass"); plot(predictors)
pred.vars <- names(predictors); pred.vars <- pred.vars[pred.vars != "aspect"] # remove aspect, because it has no biological relevance for fish

raster_df <- as.data.frame(predictors, xy = TRUE, na.rm = TRUE) %>%
  glimpse()

dat <- terra::extract(predictors, pres[, c("longitude", "latitude")], cells = TRUE) # Environmental values for where fish was found

dat <- as.data.frame(dat) %>% 
  mutate(ID = row_number()) %>% 
  glimpse()

if (study_site == "waatu") {
  dat <- pres %>%
    left_join(dat, by = "ID") %>%
    filter(as.numeric(depth_m) < 60) # remove deep drops here, otherwise the models are not happy.
} else {
  dat <- pres %>%
    left_join(dat, by = "ID")
}

head(dat)
saveRDS(dat, paste0("outputs/distribution_modelling_outputs/D03_", study_site, "_gam_model_data.rds")) # for rmarkdown


# Explore data quickly
ggplot(dat, aes(x = longitude, y = latitude)) +
  geom_point(aes(size = maxn), col = 'navy') + # all species
  scale_size_continuous(range = c(1, 7.5)) +
  facet_wrap(~fullspp)


## Model selection and prediction ---------------------------------------------

# Modify the dataframe to your species
dat_filtered <- dat[dat$fullspp %in% model_spp,] %>%
  na.omit() %>% # This removes any incomplete data
  glimpse()

ggplot(dat_filtered, aes(x = longitude, y = latitude)) +
  geom_point(aes(size = maxn), col = 'navy')

outdir <- ("outputs/distribution_modelling_outputs/D03_distribution_model_diagnostic_plots/")

use.dat <- na.omit(dat_filtered)
use.dat <- as.data.frame(use.dat)

## Testing which model family is best -----------------------------------------

Model1 <- gam(maxn ~
                s(seagrass, bs = 'cr'),
              family = poisson(link = "log"),  data = use.dat)

# Step 1: Generate candidate models using `generate.model.set`
model.set <- generate.model.set(
  use.dat = use.dat,
  test.fit = Model1,
  pred.vars.cont = pred.vars,
  k = 5,
  cov.cutoff = 0.7,
  max.predictors = 3
)

# Step 2: Fit all candidate models using Poisson and store AIC values
out.list <- fit.model.set(model.set,
                          max.models = 600,
                          parallel = T)


# Step 3: Select the best Poisson model based on AIC
out.list$failed.models # no failed models, all good
mod.table           <- out.list$mod.data.out
mod.table           <- mod.table[order(mod.table$AICc), ] # AICc value, the lower, the better the performance of the model
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc) # wi.AICc is the model weight, which compares the proportion of evidence for each model. The higher the better.
out.i               <- mod.table[which(mod.table$delta.AICc <= 2), ] # This compares the models that have comparable AICc values (here, within 2 AICc values)
out.i
best.model.name     <- as.character(out.i$modname)[1]; print(best.model.name)

# Step 4: Extract the formula of the best model
best_model <- out.list$success.models[[best.model.name]]
best_formula <- formula(best_model)
print(best_formula)

# Step 5: Define families for comparison
families <- list(
  poisson = poisson(link = "log"),
  NB = nb(link = "log"),
  tweedie = tw()
)

# Initialize an empty data frame to store the results
model_comparison <- data.frame(
  Model = character(),
  AIC = numeric(),
  Dispersion = numeric(),
  Deviance = numeric(),
  stringsAsFactors = FALSE
)

# Step 6: Fit the best formula to Poisson, Tweedie, and NB
for (fam_name in names(families)) {
  model <- gam(best_formula, family = families[[fam_name]], data = use.dat)
  gam.check(model)
  
  # Store model information
  model_info <- data.frame(
    Model = fam_name,
    AIC = AIC(model),
    Dispersion = deviance(model) / df.residual(model),
    Deviance = model$deviance,
    stringsAsFactors = FALSE
  )
  
  # Add to the results table
  model_comparison <- rbind(model_comparison, model_info)
}

# Print the final model comparison table
print(model_comparison)

# Choosing the family of the model
#fam <- poisson(link = "log")
#fam <- nb(link = "log")
fam <- tw(link = "log")

# We'll use Tweedie for all species because of convergence issues with NB models

## Run the model selection process --------------------------------------------

# Create a reference model to test all other models to.
Model1 <- gam(maxn ~ s(seagrass, bs = 'cr', k = 5), 
              family = fam, 
              data = use.dat, 
              method = 'REML')

summary(Model1)

# Generate sets of 3 variable combinations
model.set <- generate.model.set(use.dat = use.dat,
                                test.fit = Model1,
                                pred.vars.cont = pred.vars,
                                k = 5,
                                cov.cutoff = 0.7,
                                max.predictors = 3)

# Fit the sets of variables to the data
out.list <- fit.model.set(model.set,
                          max.models = 600,
                          parallel = T)
names(out.list)

# From all the combinations of predictors, we'll now look at how they compare to each other
out.list$failed.models # look for failed models
out.list$success.models # look for successful models

mod.table           <- out.list$mod.data.out
mod.table           <- mod.table[order(mod.table$AICc), ] # AICc value, the lower, the better the performance of the model

mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc) # wi.AICc is the model weight, which compares the proportion of evidence for each model. The higher the better.
out.i               <- mod.table[which(mod.table$delta.AICc <= 2), ] # This compares the models that have comparable AICc values (here, within 2 AICc values)
out.i

# here, models for some of the species/study-site combos have to be manually chosen because the top model does not reflect actual distribution
best.model.name     <- ifelse(study_site == "waatern" & model_spp == "Glaucosomatidae Glaucosoma hebraicum", "bathy+reef", as.character(out.i$modname)[1]); print(best.model.name)
best.model.name     <- ifelse(study_site == "waatu" & model_spp == "Labridae Ophthalmolepis lineolatus", "roughness+seagrass", best.model.name); print(best.model.name)

#saveRDS(best.model.name, file = sub("_\\.rds$", ".rds", paste("data/rmd/gam_best_model", dat_extent, paste(gsub(" ", "_", spp), collapse = "_"), ".rds", sep = "_")))

# We'll plot to see the spread of data in the model. Does it make biological sense? etc. very few observations where you expect few etc.
png(paste0("outputs/distribution_modelling_outputs/D03_distribution_model_diagnostic_plots/gam_fit_plot_", study_site, "_", model_spp, ".png"), width = 900, height = 900)  # or use jpeg() or pdf()
par(mfrow = c(3, 1))
best.model = out.list$success.models[[best.model.name]]
plot(best.model, all.terms = TRUE, pages = 1, residuals = TRUE, pch = 16)
title(main = paste("Model Diagnostics", model_spp), line = -6, cex.main = 2)
dev.off()

png(paste0("outputs/distribution_modelling_outputs/D03_gam_diagnostic_plot_", study_site, "_", model_spp, ".png"), width = 900, height = 900)  # or use jpeg() or pdf()
par(mfrow = c(2, 2))
best.model = out.list$success.models[[best.model.name]]
gam.check(best.model)
title(main = paste("Model Diagnostics", model_spp), line = -6, cex.main = 2)
dev.off()

# Using the best models above, fitting a GAM to fish
best.model.vars <- strsplit(best.model.name, "\\+")[[1]]; best.model.vars

# Format variables correctly for the gam to interpret
make_vars <- function(vars) {
  result <- sapply(vars, function(var) {
    # Check if 'aspect' is part of the variable name
    bs_type <- ifelse(grepl("aspect", var), "cc", "cr")
    
    # Create the s() format for each variable with the correct bs value
    paste0("s(", var, ", k = 3, bs = '", bs_type, "')")
  })
  
  # Join the results into a single string, separated by '+'
  return(paste(result, collapse = " + "))
}
best_vars <- make_vars(best.model.vars)

# Use a function to automatically fit the model, filtered to the model_spp data, and to the optimal variables
make_optimal_model <- function(vars, family) {
  # Create the model formula as a string
  formula_str <- paste("maxn ~", vars)
  
  # Fit the model using the formula string
  model <- gam(as.formula(formula_str),
               data = use.dat, method = "REML", family = family)
  return(model)
}
m_fish <- make_optimal_model(best_vars, fam)

summary(m_fish)
plot(m_fish, pages = 1, residuals = T, cex = 5)


## Remove variable values that are out of range -------------------------------

sel_rast <- predictors[[names(m_fish$model[2:length(names(m_fish$model))])]]; plot(sel_rast) # Subset the environmental raster to only the variables in the model

# Cropping the abundance prediction to observed env. variables only
extracted_values <- terra::extract(sel_rast, dat[, c("longitude", "latitude")]) %>%
  as_tibble() %>% 
  dplyr::select(!ID) %>% 
  glimpse() # Values of environmental variables observed on BRUVs (extracted from prediction rasters, not from direct annotations)

mask_raster <- rast(nrows = nrow(sel_rast), # Empty raster in which to store masks
                    ncols = ncol(sel_rast),
                    ext = ext(sel_rast),
                    crs = crs(sel_rast),
                    nlyrs = nlyr(sel_rast))


for (i in 1:nlyr(sel_rast)) {
  # Get the range for the current layer
  layer_min <- min(extracted_values[, i], na.rm = TRUE)
  layer_max <- max(extracted_values[, i], na.rm = TRUE)
  
  # Print the layer min and max
  cat("Layer", i, ", min:", layer_min, ", max:", layer_max, "\n")
  
  # Create a temporary vector for the current layer of sel_rast
  current_layer_values <- values(sel_rast[[i]])
  
  # Create a temporary vector for the mask layer, set to NA where out of range
  #mask_layer <- ifelse(current_layer_values < layer_min | current_layer_values > layer_max, NA, current_layer_values)
  mask_layer <- ifelse(is.na(current_layer_values), NA, current_layer_values)
  
  # Assign the mask_layer as a SpatRaster layer
  mask_raster[[i]] <- rast(nrows = nrow(sel_rast), ncols = ncol(sel_rast), vals = mask_layer, ext = ext(sel_rast), crs = crs(sel_rast))
}
plot(mask_raster)

# Keep only the cells that exist in all layers - this is to avoid extrapolating
valid_mask <- Reduce(`&`, lapply(1:length(mask_raster), function(i) !is.na(mask_raster[[i]])))
plot(valid_mask)
masked_raster <- mask(mask_raster, valid_mask)
masked_raster[!valid_mask] <- NA; plot(masked_raster)

ras_crop <- as.data.frame(masked_raster, xy = TRUE, na.rm = TRUE) %>% glimpse()


## Predicting mature distribution ---------------------------------------------

predicted_abundance <- cbind(ras_crop, "p_fish" = mgcv::predict.gam(m_fish, ras_crop, type = "response", se.fit = T)) %>%
  glimpse()

p_fish <- rast(predicted_abundance)

par(mfrow=c(1,1))
plot(p_fish$p_fish.fit, pch = 19,
     main = paste(model_spp, best.model.name),
     xlab = "Longitude", ylab = "Latitude")


## MEGAPLOT -------------------------------------------------------------------

# prepare for  saving the plot
png(paste0("outputs/distribution_modelling_outputs/D03_distribution_prediction_plots/D03_", study_site, "_distribution_prediction_plot_", paste(gsub(" ", "_", model_spp)), ".png"),
    width = 1200, height = 1200, res = 150)  # Open PNG device

# Set up the layout for a 3x3 grid
layout(matrix(c(1, 1, 3, 1, 1, 4, 2, 2, 5), nrow = 3, ncol = 3, byrow = TRUE))

# Adjust margins for better spacing
par(mar = c(10, 3, 3, 3))  # Minimal margins for individual plots
par(oma = c(1, 1, 1, 1))  # Minimal outer margins for overall layout
# Plot 1: p_mature$p_mature.fit (scatter plot) - Plot occupies positions 1, 1, and 1
plot(p_fish$p_fish.fit, pch = 25,
     main = paste(gsub(" ", "_", model_spp), collapse = "_"),
     xlab = "Longitude", ylab = "Latitude",
     xlim = ext_plot[2:1], ylim = ext_plot[4:3])

# Plot 2: Model summary (m_mature) - Plot occupies positions 2, 2
par(mar = c(0, 0, 0, 0))  # Remove margins
plot.new()  # Create a blank space for text
summary_info <- capture.output(summary(m_fish))  # Capture the summary output
text(0.1, 0.8, paste(summary_info, collapse = "\n"), cex = 1, adj = 0, family = "mono")

# Plot 3: final covariates
raster_names <- names(mask_raster)
for (i in 1:min(length(raster_names), 3)) {
  plot(mask_raster[[raster_names[i]]], 
       main = raster_names[i])
}
dev.off() # save for future reference


## Save predictions as rasters ------------------------------------------------

writeRaster(c(p_fish$p_fish.fit, p_fish$p_fish.se.fit), 
            filename = paste0("outputs/distribution_modelling_outputs/D03_distribution_rasters/D03_", study_site, "_predicted_distribution_raster_", paste(model_spp, collapse = "_"),
                              ".tif"), 
            overwrite = TRUE)


### END ###