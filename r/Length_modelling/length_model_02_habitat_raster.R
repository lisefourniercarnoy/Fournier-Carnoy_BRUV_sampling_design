# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. length model)
# Data:    2024 MEGlab Geographe Bay BRUV habitat data.
# Task:    Predict habitats over Geographe Bay using BRUV observations and GAMs.
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    February 2025

# -----------------------------------------------------------------------------

# -- To predict fish abundance over Geographe Bay, we need raster layers of
# -- habitat. We fit GAMs of habitat to bathymetry covariates to predict
# -- habitat extent over our area.

# -----------------------------------------------------------------------------

# Load libraries
library(tidyverse) # for data manipulation
library(FSSgam) # for selecting the models
library(mgcv) # for making the models
library(predicts) # for making the prediction
library(terra) # for extracting bathy values at opcode locations
library(sf) # for dealing with polygons


rm(list = ls())

ext <- c(115.035, 115.68012207, -33.69743897, -33.35) # Extent to crop layers to


## Files used in this script --------------------------------------------------

file_obs_hab            <- "data/tidy/2024_geographe_tidy_habitat.rds" # Observed habitat values, cleaned in 01_prep
file_250m_bathy_deriv   <- "data/spatial/rasters/2014_geographe_bathymetry-derivatives.RDS" # made in 01_prep
file_lidar_bathy        <- "data/spatial/rasters/LiDAR_geo_compressed.tif" # Too big for git, new file in the workflow

file_bathy_predictors   <- "data/spatial/rasters/geographe_bathy_predictors.rds" # made in 01_prep


# 250m habitat predictions ----------------------------------------------------

# -- Let's start with the 250m predictions. 

## Load data ------------------------------------------------------------------

habi <- readRDS(file_obs_hab) %>% # Observed habitat
  
  # filter out unscorables
  filter(!(level_2 %in% c("Fishes", "Echinoderms"))) %>% 
  
  # make the correct categories
  mutate(taxa = case_when(level_2 %in% "Macroalgae" ~ "macro",
                          level_3 %in% "Unconsolidated (soft)" ~ "sand",
                          level_3 %in% "Consolidated (hard)" ~ "rock",
                          level_2 %in% "Seagrasses" ~ "seagrass",
                          level_2 %in% c("Sessile invertebrates", "Sponges", "Bryozoa", "Cnidaria") ~ "sessile_inverts")) %>%
  pivot_wider(names_from = taxa, values_from = number, values_fill = list(number = 0)) %>% 
  
  # make reef and remove the individual reef habitats
  mutate(reef = macro + rock + sessile_inverts) %>%
  dplyr::select(-c(sessile_inverts, macro, rock)) %>% 
  
  # sum rows to have a single row per opcode
  group_by(opcode, latitude_dd, longitude_dd) %>%  # Group by the `opcode`
  summarize(
    sand = sum(sand, na.rm = TRUE),
    seagrass = sum(seagrass, na.rm = TRUE),
    reef = sum(reef, na.rm = TRUE)
  ) %>%
  
  # add total points - these should add up to all the reef, sand, seagrass
  mutate(total_pts = sum(c(reef, sand, seagrass), na.rm = TRUE)) %>% 
  
  # Transform to a spatial object to extract bathy values
  st_as_sf(coords = c("longitude_dd", "latitude_dd"), crs = 4326) %>%  # Replace with the correct column names
  glimpse()


# habi <- readRDS(file_obs_hab) %>% # Observed habitat
#   
#   # Shuffling groups a bit to clear up and end up with fewer categories
#   mutate(level_2_3 = paste(level_2, level_3, sep = "_")) %>%
#   pivot_wider(names_from = level_2_3, values_from = number, values_fill = list(number = 0)) %>% 
#   mutate(unscorable = Fishes_NA + Echinoderms_NA,
#          sessile_invert = `Sessile invertebrates_NA` + Sponges_NA + Bryozoa_NA + Cnidaria_Corals + Cnidaria_Hydrocorals,
#          macro = `Macroalgae_Filamentous / filiform` + 
#            `Macroalgae_NA` + 
#            `Macroalgae_Erect fine branching` + 
#            `Macroalgae_Erect coarse branching` +
#            `Macroalgae_Encrusting` + 
#            Macroalgae_Laminate +
#            `Macroalgae_Large canopy-forming`,
#          sand = `Substrate_Unconsolidated (soft)`,
#          rock = `Substrate_Consolidated (hard)`,
#          seagrass = Seagrasses_NA + `Seagrasses_Strap-like leaves`
#   ) %>% 
#   mutate(reef = macro + rock + sessile_invert) %>% 
#   
#   # Sum rows to have a single row per opcode
#   group_by(opcode, longitude_dd, latitude_dd) %>%  # Group by the `opcode`
#   summarize(
#     unscorable = sum(unscorable, na.rm = TRUE),
#     sessile_invert = sum(sessile_invert, na.rm = TRUE),
#     macro = sum(macro, na.rm = TRUE),
#     sand = sum(sand, na.rm = TRUE),
#     rock = sum(rock, na.rm = TRUE),
#     seagrass = sum(seagrass, na.rm = TRUE),
#     reef = sum(reef, na.rm = TRUE)
#   ) %>%  
#   mutate(total_pts = sum(c(sessile_invert, macro, sand, rock, seagrass, na.rm = TRUE))) %>% 
#   
#   # Pivot it longer again
#   pivot_longer(cols = c("sand", "reef", "seagrass"), # specify the columns to pivot
#                values_to = "number", 
#                names_to = "taxa") %>%
#   ungroup() %>%   # Ungroup after summarizing to return to a regular dataframe
#   
#   dplyr::select("opcode", "longitude_dd", "latitude_dd", # Select final columns
#                 "taxa", "number", "total_pts") %>%
#   
#   # Transform to a spatial object to extract bathy values
#   st_as_sf(coords = c("longitude_dd", "latitude_dd"), crs = 4326) %>%  # Replace with the correct column names
#   glimpse()

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
pred.vars <- c("bathy_250m", "aspect_250m", "roughness_250m", "detrended_250m")

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
outdir    <- ("outputs/Length_comparison/")
resp.vars <- unique.vars.use
out.all   <- list()
var.imp   <- list()


## Loop through the FSS function for each Abiotic taxa ------------------------
hab <- hab[-c(5, 6, 7)] %>% na.omit() # removing lidar columns otherwise it wont run
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

    png(file = paste(outdir, m, resp.vars[i], "250m_mod_fits.png", sep = ""))
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
# Reef : aspect + bathy + roughness
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

preds <- readRDS(file_bathy_predictors)[[4:7]]; plot(preds)
preddf <- as.data.frame(preds, xy = TRUE, na.rm = TRUE) %>% glimpse()


## Make models for each habitat and predict extent ----------------------------

# Save each model (based on top models above)

# Reef
out.all$reef[[1]]
m_reef_250m <- gam(cbind(reef, total_pts - reef) ~
                     s(aspect_250m, k = 5, bs = "cc")  +
                     #s(detrended_250m, k = 5, bs = "cr")  +
                     s(bathy_250m, k = 5, bs = "cr") +
                     s(roughness_250m, k = 5, bs = "cr"),
                   data = hab_model, method = "REML", family = binomial("logit"))
summary(m_reef_250m)
plot(m_reef_250m, pages = 1, residuals = T, cex = 5)

# Seagrass
out.all$seagrass[[1]]
m_seagrass_250m <- gam(cbind(seagrass, total_pts - seagrass) ~
                         s(aspect_250m, k = 5, bs = "cc")  +
                         #s(detrended_250m, k = 5, bs = "cr")  +
                         s(bathy_250m, k = 5, bs = "cr") +
                         s(roughness_250m, k = 5, bs = "cr"),
                       data = hab_model, method = "REML", family = binomial("logit"))
summary(m_seagrass_250m)
plot(m_seagrass_250m, pages = 1, residuals = T, cex = 5)

# Sand
out.all$sand[[1]]
m_sand_250m <- gam(cbind(sand, total_pts - sand) ~
                     s(aspect_250m, k = 5, bs = "cc")  +
                     #s(detrended_250m, k = 5, bs = "cr")  +
                     s(bathy_250m, k = 5, bs = "cr") +
                     s(roughness_250m, k = 5, bs = "cr"),
                   data = hab_model, method = "REML", family = binomial("logit"))
summary(m_sand_250m)
plot(m_sand_250m, pages = 1, residuals = T, cex = 5)


## predict, rasterise and plot ------------------------------------------------
  
preddf_250m <- cbind(preddf,
                     "preef_250"     = predict(m_reef_250m, preddf, type = "response", se.fit = T),
                     "psand_250"     = predict(m_sand_250m, preddf, type = "response", se.fit = T),
                     "pseagrass_250" = predict(m_seagrass_250m, preddf, type = "response", se.fit = T)) %>%
  glimpse()

prasts_250m <- rast(preddf_250m %>% dplyr::select(x, y, preef_250.fit, psand_250.fit, pseagrass_250.fit),
                     crs = crs(preds)) %>%
  crop(ext); plot(prasts_250m, range = c(0, 1))

# LiDAR habitat predictions ---------------------------------------------------

# -- We need to create another predicted habitat map using the finer-scale LiDAR
# -- bathymetry data. We'll use the exact same method as above.

## Load data ------------------------------------------------------------------

habi <- readRDS(file_obs_hab) %>% # Observed habitat
  
  # filter out unscorables
  filter(!(level_2 %in% c("Fishes", "Echinoderms"))) %>% 
  
  # make the correct categories
  mutate(taxa = case_when(level_2 %in% "Macroalgae" ~ "macro",
                          level_3 %in% "Unconsolidated (soft)" ~ "sand",
                          level_3 %in% "Consolidated (hard)" ~ "rock",
                          level_2 %in% "Seagrasses" ~ "seagrass",
                          level_2 %in% c("Sessile invertebrates", "Sponges", "Bryozoa", "Cnidaria") ~ "sessile_inverts")) %>%
  pivot_wider(names_from = taxa, values_from = number, values_fill = list(number = 0)) %>% 
  
  # make reef and remove the individual reef habitats
  mutate(reef = macro + rock + sessile_inverts) %>%
  dplyr::select(-c(sessile_inverts, macro, rock)) %>% 
  
  # sum rows to have a single row per opcode
  group_by(opcode, latitude_dd, longitude_dd) %>%  # Group by the `opcode`
  summarize(
    sand = sum(sand, na.rm = TRUE),
    seagrass = sum(seagrass, na.rm = TRUE),
    reef = sum(reef, na.rm = TRUE)
  ) %>%
  
  # add total points - these should add up to all the reef, sand, seagrass
  mutate(total_pts = sum(c(reef, sand, seagrass), na.rm = TRUE)) %>% 
  
  # Transform to a spatial object to extract bathy values
  st_as_sf(coords = c("longitude_dd", "latitude_dd"), crs = 4326) %>%  # Replace with the correct column names
  glimpse()

# Check the result
head(habi)

# Extract bathymetry values to the location of the opcodes
bathy <- readRDS(file_bathy_predictors); plot(bathy)
bathy_values <- terra::extract(bathy, habi)

hab <- habi %>%
  dplyr::bind_cols(bathy_values) %>%   # Combine the opcode with the bathy_values dataframe
  group_by(opcode) %>%
  ungroup() %>%
  na.omit() %>% # Removing drops that are beyond the lidar extent
  pivot_longer(cols = c("sand", "reef", "seagrass"), # specify the columns to pivot
               values_to = "number",
               names_to = "taxa") %>% 
  glimpse()

hab <- as.data.frame(hab) %>% glimpse()

colSums(is.na(hab)) # There should be no NAs. We've removed them just above, so that only LiDAR extent data exists.

names(hab)
pred.vars <- c("bathy_lidar", "aspect_lidar", "roughness_lidar")

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)
round(cor(hab[, pred.vars]), 2) # All good, no high correlations

# Check to make sure Response vector has not more than 80% zeros
unique.vars = unique(as.character(hab$taxa))

unique.vars.use = character()
for(i in 1:length(unique.vars)){
  temp.dat = hab[which(hab$taxa == unique.vars[i]),]
  if(length(which(temp.dat$taxa == 0))/nrow(temp.dat)<0.9){
    unique.vars.use = c(unique.vars.use,unique.vars[i])}
}

unique.vars.use # All remain


## Run the full subset model selection ----
outdir    <- ("outputs/Length_comparison/")
resp.vars <- unique.vars.use
out.all   <- list()
var.imp   <- list()


## Loop through the FSS function for each Abiotic taxa ------------------------

for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  use.dat <- hab[hab$taxa == resp.vars[i],]
  use.dat   <- as.data.frame(use.dat)

  # Basic model to compare all the other combinations with
  Model1  <- gam(cbind(number, (total_pts - number)) ~ # Success and failure counts for each habitat
                   s(bathy_lidar, bs = 'cr', k = 3),
                 family = binomial(link = "logit"),  data = use.dat)
  # Generate variable combinations to test in models
  model.set <- generate.model.set(use.dat = use.dat,
                                  test.fit = Model1,
                                  pred.vars.cont = pred.vars,
                                  cyclic.vars = c("aspect"),
                                  k = 5,
                                  cov.cutoff = 0.7,
                                  max.predictors = 2 # reduce max number of preds. because I only have 3
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

    png(file = paste(outdir, m, resp.vars[i], "lidar_mod_fits.png", sep = ""))
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
# Sand :        aspect + bathy
# Reef :        aspect + roughness
# Seagrasses:   bathy + roughness


## reload & format the data ---------------------------------------------------

hab_test <- hab %>% 
  pivot_wider(names_from = taxa, values_from = number, values_fill = list(number = 0)) %>% 
  glimpse()

hab_model <- hab_test %>%
  group_by(opcode) %>%
  summarize(
    across(c(bathy_lidar, roughness_lidar, aspect_lidar, 
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
hab_model$diff # there should be only zeroes

preds <- readRDS(file_bathy_predictors)[[1:3]]; plot(preds)
preddf <- as.data.frame(preds, xy = TRUE, na.rm = TRUE) %>% glimpse()


## Make models for each habitat and predict extent ----------------------------

# Save each model (based on top models above)

# Reef
out.all$reef[[1]]
m_reef_lidar <- gam(cbind(reef, total_pts - reef) ~
                      s(aspect_lidar, k = 5, bs = "cc")  +
                      s(roughness_lidar, k = 5, bs = "cr"),
                      #s(bathy_lidar, k = 5, bs = "cr"),
                    data = hab_model, method = "REML", family = binomial("logit"))
summary(m_reef_lidar)
plot(m_reef_lidar, pages = 1, residuals = T, cex = 5)

# Sand
out.all$sand[[1]]
m_sand_lidar <- gam(cbind(sand, total_pts - sand) ~
                      s(aspect_lidar, k = 5, bs = "cc")  +
                      #s(roughness_lidar, k = 5, bs = "cr") +
                      s(bathy_lidar, k = 5, bs = "cr"),
                      data = hab_model, method = "REML", family = binomial("logit"))
summary(m_sand_lidar)
plot(m_sand_lidar, pages = 1, residuals = T, cex = 5)

# Seagrass
out.all$seagrass[[1]]
m_seagrass_lidar <- gam(cbind(seagrass, total_pts - seagrass) ~
                          #s(aspect_lidar, k = 5, bs = "cc")  +
                          s(roughness_lidar, k = 5, bs = "cr") +
                          s(bathy_lidar, k = 5, bs = "cr"),
                          data = hab_model, method = "REML", family = binomial("logit"))
summary(m_seagrass_lidar)
plot(m_seagrass_lidar, pages = 1, residuals = T, cex = 5)


## predict, rasterise and plot ------------------------------------------------

preddf_lidar <- cbind(preddf,
                      "preef_lidar"     = predict(m_reef_lidar, preddf, type = "response", se.fit = T),
                      "psand_lidar"     = predict(m_sand_lidar, preddf, type = "response", se.fit = T),
                      "pseagrass_lidar" = predict(m_seagrass_lidar, preddf, type = "response", se.fit = T)) %>%
  glimpse()

prasts_lidar <- rast(preddf_lidar %>% dplyr::select(x, y, preef_lidar.fit, psand_lidar.fit, pseagrass_lidar.fit),
               crs = crs(preds)) %>%
  crop(ext); plot(prasts_lidar, range = c(0, 1))


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
models <- list(m_reef_250m, m_sand_250m, m_seagrass_250m); summary(models)

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


### lidar prediction predicts::MESS -------------------------------------------


resp.vars <- c("preef_lidar", "psand_lidar", "pseagrass_lidar")
models <- list(m_reef_lidar, m_sand_lidar, m_seagrass_lidar); summary(models)

predhab <- preddf_lidar

# Extract predictor variables from the model's terms
model_terms <- terms(models[[1]])
model_vars <- attr(model_terms, "term.labels")

# Check the model variables
print(model_vars)

ext(temppred_m) == ext(rast(predhab))

# loop for each lidar habitat model
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
    preddf_m_lidar <- temppred_m
  }
  else {
    preddf_m_lidar <- rast(list(preddf_m_lidar, temppred_m))
  }
}
plot(preddf_m_lidar) # check where things were removed


## Make a big file stack with the cleaned predictions and the bathy layers ----

## Unclean versions of the predictions
#saveRDS(preddf_250m, paste0("outputs/Length_comparison/2024_geographe_predicted-habitat_250m.rds"))
#writeRaster(prasts_250m, filename = "outputs/Length_comparison/2024_geographe_predicted-habitat_250m.tif", overwrite = TRUE)
#saveRDS(preddf_lidar, paste0("outputs/Length_comparison/2024_geographe_predicted-habitat_LiDAR.rds")) # This takes ~20 seconds
#writeRaster(prasts_lidar, filename = "outputs/Length_comparison/2024_geographe_predicted-habitat_LiDAR.tif", overwrite = TRUE)

bathy <- readRDS(file_bathy_predictors); plot(bathy)

predictors <- list(bathy_lidar = bathy[[1]], roughness_lidar = bathy[[2]], aspect_lidar = bathy[[3]], # if including LiDAR data
                   reef_lidar = preddf_m_lidar[[1]], sand_lidar = preddf_m_lidar[[3]], seagrass_lidar = preddf_m_lidar[[5]],
                   aspect_250 = bathy[[4]], depth_250 = bathy[[5]], roughness_250 = bathy[[6]], detrended_250 = bathy[[7]],
                   reef_250m = preddf_m_250m[[1]], sand_250m = preddf_m_250m[[3]], seagrass_250m = preddf_m_250m[[5]]
)

#predictors <- list( # if not including LiDAR data
#  aspect_250 = bathy[[2]], depth_250 = bathy[[1]], roughness_250 = bathy[[3]], detrended_250 = bathy[[4]],
#  reef_250m = prasts_250m[[1]], sand_250m = prasts_250m[[2]], seagrass_250m = prasts_250m[[3]]
#)

for (i in 1:length(predictors)) { # Unifying extents and resampling so that they can be stacked
  reference_raster <- predictors[[1]] # Use the raster with the highest resolution as a reference. Otherwise the re-sampling will apply lower resolutions to everything.
  predictors[[i]] <- terra::resample(predictors[[i]], reference_raster, method = "bilinear")
  cat(paste("New extent of", i, ":", ext(predictors[[i]]), "\n"))
}

list2env(predictors, envir = .GlobalEnv) # Put the elements of the list in the environment

predictors <- c(
  bathy_lidar, roughness_lidar, aspect_lidar,          # LiDAR data
  depth_250, aspect_250, roughness_250, detrended_250, # 250m data
  reef_lidar, sand_lidar, seagrass_lidar,              # Habitat predictions using LiDAR bathymetry
  reef_250m, sand_250m, seagrass_250m                 # Habitat predictions using 250m bathymetry
)

# names(pred_rast) <- c("bathy_lidar", "roughness_lidar", "aspect_lidar",
#                       "aspect_250m", "bathy_250m",      "roughness_250m", "detrended_250m",
#                       "reef_lidar",  "sand_lidar",      "seagrass_lidar",
#                       "reef_250m",   "sand_250m",       "seagrass_250m")

plot(predictors)
predictors <- readRDS("data/spatial/rasters/geographe_all_predictors.rds")
class(predictors)

#predictors <- raster::stack(predictors) # convert to the right format before saving
#predictors <-  rast(predictors) # convert to rast
saveRDS(predictors, file = "data/spatial/rasters/geographe_all_predictors.rds") # Takes a while if running with LiDAR
# ^ needed for 03_gam


### END ###
