# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare abundance detected)
# Data:    2024 MEGlab Geographe Bay BRUV habitat data.
# Task:    Predict habitats over Geographe Bay using BRUV observations and GAMs.
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    February 2025

# -----------------------------------------------------------------------------

# Load libraries
library(tidyverse) # for data manipulation
library(FSSgam) # for selecting the models
library(mgcv) # for making the models
library(predicts) # for making the prediction
library(terra) # for extracting bathy values at opcode locations
library(sf) # for dealing with polygons


rm(list = ls())

study_site <- "waatu" # waatu or waatern

ext_plot <- if(study_site == "waatern") {
  ext <- c(115.035, 115.68012207, -33.69743897, -33.35) 
} else { 
  ext <- c(114.7, 115, -34.25, -33.95)
} # Extent to crop layers to


## Files used in this script --------------------------------------------------

file_obs_hab            <- paste0("data/tidy/D01_", study_site, "_tidy_habitat.rds") # Observed habitat values, cleaned in 01_prep
file_bathy_predictors   <- paste0("data/tidy/D01_", study_site, "_bathy_predictors.rds") # made in 01_prep


# 250m habitat predictions ----------------------------------------------------

# -- Let's start with the 250m predictions. 

## Load data ------------------------------------------------------------------
habi <- readRDS(file_obs_hab) %>% glimpse()

# Extract bathymetry values to the location of the opcodes
bathy <- readRDS(file_bathy_predictors); plot(bathy)
bathy_values <- terra::extract(bathy, habi)

hab <- habi %>%
  dplyr::bind_cols(bathy_values) %>%   # Combine the opcode with the bathy_values dataframe
  group_by(opcode) %>%
  ungroup() %>%
  #na.omit() %>% # Removing drops that are beyond the lidar extent
  filter(!is.na(aspect_250m)) %>% # removing one NA value
  pivot_longer(cols = c("sand", "reef", "seagrass"), # specify the columns to pivot
               values_to = "number",
               names_to = "taxa") %>% 
  glimpse()
hab <- as.data.frame(hab) %>% glimpse()

names(hab)
pred.vars <- c("bathy_250m", #"aspect_250m", 
               "roughness_250m", "detrended_250m")

colSums(is.na(hab)) # There should be no NAs. We've removed them just above, so that only LiDAR extent data exists.

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)
round(cor(hab[, pred.vars]), 2) # All good, no high correlations

# Check to make sure Response vector has not more than 80% zeros
unique.vars = unique(as.character(hab$taxa))
head(hab)
unique.vars.use = character()
for(i in 1:length(unique.vars)){
  temp.dat = hab[which(hab$taxa == unique.vars[i]),]
  if(length(which(temp.dat$taxa == 0))/nrow(temp.dat)<0.9){
    unique.vars.use = c(unique.vars.use,unique.vars[i])}
}

unique.vars.use # All remain


## Run the full subset model selection ----------------------------------------
outdir    <- ("outputs/distribution_modelling_outputs/")
resp.vars <- unique.vars.use
out.all   <- list()
var.imp   <- list()


## Loop through the FSS function for each Abiotic taxa ------------------------
hab <- hab %>% na.omit() # removing lidar columns otherwise it wont run
colSums(is.na(hab))
for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  use.dat <- hab[hab$taxa == resp.vars[i],]
  use.dat   <- as.data.frame(use.dat)
  
  # Basic model to compare all the other combinations with
  Model1  <- gam(cbind(number, (total_pts - number)) ~ # Success and failure counts for each habitat
                   s(bathy_250m, bs = 'cr', k = 3),
                 family = binomial(link = "logit"),  data = use.dat)
  
  # Generate variable combinations to test in models
  model.set <- generate.model.set(use.dat = use.dat,
                                  test.fit = Model1,
                                  pred.vars.cont = pred.vars,
                                  cyclic.vars = c("aspect_250m"),
                                  k = 5,
                                  cov.cutoff = 0.7,
                                  max.predictors = 3
  )
  
  # Fit models to all variable combinations
  out.list <- fit.model.set(model.set,
                            max.models = 600,
                            parallel = T)
  names(out.list)
  
  # Examine outcomes
  out.list$failed.models # examine the list of failed models
  mod.table <- out.list$mod.data.out  # look at the model selection table
  mod.table <- mod.table[order(mod.table$AICc), ]
  mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
  out.i     <- mod.table[which(mod.table$delta.AICc <= 2), ] # Top models
  out.all   <- c(out.all, list(out.i))
  var.imp   <- c(var.imp, list(out.list$variable.importance$aic$variable.weights.raw))
  
  # Plot the top models
  for(m in 1:nrow(out.i)){
    best.model.name <- as.character(out.i$modname[m])
    
    png(file = paste(outdir, "D02", study_site, resp.vars[i], "250m_mod_fits.png", sep = "_"))
    if(best.model.name != "null"){
      par(mfrow = c(3, 1), mar = c(9, 4, 3, 1))
      best.model = out.list$success.models[[best.model.name]]
      plot(best.model, all.terms = T, pages = 1, residuals = T, pch = 16)
      mtext(side = 2, text = resp.vars[i], outer = F)}
    dev.off()
  }
}

## Model fits and importance --------------------------------------------------

names(out.all) <- resp.vars
names(var.imp) <- resp.vars
all.mod.fits <- list_rbind(out.all, names_to = "response")
all.var.imp  <- do.call("rbind", var.imp)
out.all
# Sand : aspect + bathy + roughness
# Reef : aspect + bathy + detrended
# Seagrass: aspect + bathy + roughness


## reload & format the data ---------------------------------------------------

hab_test <- hab %>% 
  pivot_wider(names_from = taxa, values_from = number, values_fill = list(number = 0)) %>% 
  glimpse()

hab_model <- hab_test %>%
  group_by(opcode) %>%
  summarize(
    across(c(#bathy_lidar, roughness_lidar, aspect_lidar, 
      bathy_250m, aspect_250m, roughness_250m, detrended_250m), 
      ~ first(.x, order_by = .x),  # Keep the first non-NA value
      .names = "{.col}"),  # This makes sure column names are kept
    reef = sum(reef, na.rm = F),
    sand = sum(sand, na.rm = F),
    seagrass = sum(seagrass, na.rm = F),
    total_pts = first(total_pts, order_by = total_pts),  # Keeps the first non-NA value
    .groups = "drop"  # Remove grouping after summarization
  ) %>% 
  glimpse()

hab_model$diff <- (hab_model$reef + hab_model$seagrass + hab_model$sand) - hab_model$total_pts
hab_model$diff # There should only be zeroes.

preds <- readRDS(file_bathy_predictors) %>% dplyr::select(-aspect_250m); plot(preds)
preddf <- as.data.frame(preds, xy = TRUE, na.rm = TRUE) %>% glimpse()


## Make models for each habitat and predict extent ----------------------------

fam <- binomial("logit")

make_vars <- function(vars) {
  sapply(vars, function(var) {
    bs_type <- ifelse(grepl("aspect", var), "cc", "cr")
    paste0("s(", var, ", k = 5, bs = '", bs_type, "')")
  }) |> paste(collapse = " + ")
}

# Function to fit the optimal model
make_optimal_model <- function(response_var, predictor_terms, data) {
  formula_str <- paste0(
    "cbind(", response_var, ", total_pts - ", response_var, ") ~ ", predictor_terms
  )
  
  gam(
    as.formula(formula_str),
    data = data,
    method = "REML",
    family = binomial
  )
}

# Container for storing the fitted models
fitted_models <- list()

for (habitat in names(out.all)) {
  
  # Step 1: Extract best model variable names
  best_model_name <- out.all[[habitat]]$modname[1]
  best_model_vars <- strsplit(best_model_name, "\\+")[[1]]
  
  # Step 2: Convert variables into GAM smooth terms
  best_vars <- make_vars(best_model_vars)
  
  # Step 3: Fit the model using cbind(successes, failures)
  fitted_model <- make_optimal_model(
    response_var = habitat,
    predictor_terms = best_vars,
    data = hab_model
  )
  
  # Step 4: Store model
  fitted_models[[habitat]] <- fitted_model
  
  # Optional: print quick summary
  message("Fitted binomial model for habitat: ", habitat)
  print(summary(fitted_model))
}

plot(fitted_models$reef, page = 1, residuals = T, cex = 5)
plot(fitted_models$sand, page = 1, residuals = T, cex = 5)
plot(fitted_models$seagrass, page = 1, residuals = T, cex = 5)




# 
# # Save each model (based on top models above)
# 
# # Reef
# out.all$reef[[1]]
# m_reef_250m <- gam(cbind(reef, total_pts - reef) ~
#                      s(aspect_250m, k = 5, bs = "cc")  +
#                      s(detrended_250m, k = 5, bs = "cr")  +
#                      s(bathy_250m, k = 5, bs = "cr"),
#                      s(roughness_250m, k = 5, bs = "cr"),
#                    data = hab_model, method = "REML", family = binomial("logit"))
# summary(m_reef_250m)
# plot(m_reef_250m, pages = 1, residuals = T, cex = 5)
# 
# # Seagrass
# out.all$seagrass[[1]]
# m_seagrass_250m <- gam(cbind(seagrass, total_pts - seagrass) ~
#                          s(aspect_250m, k = 5, bs = "cc")  +
#                          #s(detrended_250m, k = 5, bs = "cr")  +
#                          s(bathy_250m, k = 5, bs = "cr") +
#                          s(roughness_250m, k = 5, bs = "cr"),
#                        data = hab_model, method = "REML", family = binomial("logit"))
# summary(m_seagrass_250m)
# plot(m_seagrass_250m, pages = 1, residuals = T, cex = 5)
# 
# # Sand
# out.all$sand[[1]]
# m_sand_250m <- gam(cbind(sand, total_pts - sand) ~
#                      s(aspect_250m, k = 5, bs = "cc")  +
#                      #s(detrended_250m, k = 5, bs = "cr")  +
#                      s(bathy_250m, k = 5, bs = "cr") +
#                      s(roughness_250m, k = 5, bs = "cr"),
#                    data = hab_model, method = "REML", family = binomial("logit"))
# summary(m_sand_250m)
# plot(m_sand_250m, pages = 1, residuals = T, cex = 5)


## predict, rasterise and plot ------------------------------------------------

preddf_250m <- cbind(preddf,
                     "preef_250"     = predict(fitted_models$reef, preddf, type = "response", se.fit = T),
                     "psand_250"     = predict(fitted_models$sand, preddf, type = "response", se.fit = T),
                     "pseagrass_250" = predict(fitted_models$seagrass, preddf, type = "response", se.fit = T)) %>%
  glimpse()

prasts_250m <- rast(preddf_250m %>% dplyr::select(x, y, preef_250.fit, psand_250.fit, pseagrass_250.fit),
                    crs = crs(preds)) %>%
  crop(ext); plot(prasts_250m, range = c(0, 1))


## Removing parts of the prediction that are not observed ---------------------

# using MESS (Elith et al. 2010), we'll remove parts of the prediction that have
# combinations of predictors that were not observed. This makes sure we're not
# predicting overly confidently.

habi_sf <- st_as_sf(habi, wkt = "geometry")
xy <- as.data.frame(habi_sf)
xy$geometry <- st_coordinates(habi_sf)

xy <- xy %>% 
  mutate(
    x = xy$geometry[, 1],
    y = xy$geometry[, 2]
  ) %>% 
  dplyr::select(x, y) %>% 
  glimpse()


### 250m prediction predicts::MESS --------------------------------------------

resp.vars <- c("preef_250", "psand_250", "pseagrass_250")
models <- list(fitted_models$reef, fitted_models$sand, fitted_models$seagrass); summary(models)

predhab <- preddf_250m

# Extract predictor variables from the model's terms
model_terms <- terms(models[[1]])
model_vars <- attr(model_terms, "term.labels")

# Check the model variables
print(model_vars)

# loop for each 250m habitat model
for(i in 1:length(resp.vars)) {
  # select a model from the list of models above
  print(resp.vars[i])
  mod <- models[[i]]
  
  # select the prediction raster
  temppred <- predhab %>% 
    dplyr::select(x, y, paste0(resp.vars[i], '.fit'),
                  paste0(resp.vars[i], '.se.fit')) %>% 
    rast(crs = "epsg:4326")
  
  # select the predictors of the model within the whole set of predictors
  model_vars <- attr(model_terms, "term.labels")
  dat <- terra::extract(subset(bathy, model_vars), xy, ID = F)
  
  # run the MESS function
  messrast <- predicts::mess(subset(bathy, model_vars), dat) %>% 
    terra::clamp(lower = -0.01, values = F)
  
  # remove unobserved areas
  messrast <- terra::crop(messrast, temppred)
  temppred_m <- terra::mask(temppred, messrast)
  
  # then add the final rasters together
  if (i == 1) {
    preddf_m_250m <- temppred_m
  }
  else {
    preddf_m_250m <- rast(list(preddf_m_250m, temppred_m))
  }
}
plot(preddf_m_250m) # check where things were removed


## Make a big file stack with the cleaned predictions and the bathy layers ----

## Unclean versions of the predictions
#saveRDS(preddf_250m, paste0("outputs/Length_comparison/2024_geographe_predicted-habitat_250m.rds"))
#writeRaster(prasts_250m, filename = "outputs/Length_comparison/2024_geographe_predicted-habitat_250m.tif", overwrite = TRUE)
#saveRDS(preddf_lidar, paste0("outputs/Length_comparison/2024_geographe_predicted-habitat_LiDAR.rds")) # This takes ~20 seconds
#writeRaster(prasts_lidar, filename = "outputs/Length_comparison/2024_geographe_predicted-habitat_LiDAR.tif", overwrite = TRUE)

bathy <- readRDS(file_bathy_predictors); plot(bathy)

predictors <- list(aspect_250 = bathy[[2]], depth_250 = bathy[[1]], roughness_250 = bathy[[3]], detrended_250 = bathy[[4]],
                   reef_250m = preddf_m_250m[[1]], sand_250m = preddf_m_250m[[3]], seagrass_250m = preddf_m_250m[[5]]
)

for (i in 1:length(predictors)) { # Unifying extents and resampling so that they can be stacked
  reference_raster <- predictors[[1]] # Use the raster with the highest resolution as a reference. Otherwise the re-sampling will apply lower resolutions to everything.
  predictors[[i]] <- terra::resample(predictors[[i]], reference_raster, method = "bilinear")
  cat(paste("New extent of", i, ":", ext(predictors[[i]]), "\n"))
}

list2env(predictors, envir = .GlobalEnv) # Put the elements of the list in the environment

predictors <- c(
  depth_250, aspect_250, roughness_250, detrended_250, # 250m data
  reef_250m, sand_250m, seagrass_250m                 # Habitat predictions using 250m bathymetry
)

plot(predictors)

#predictors <- raster::stack(predictors) # convert to the right format before saving
#predictors <-  rast(predictors) # convert to rast
saveRDS(predictors, file = paste0("data/tidy/D02_", study_site, "_all_predictors.rds")) # Takes a while if running with LiDAR
# ^ needed for 03_gam


### END ###
