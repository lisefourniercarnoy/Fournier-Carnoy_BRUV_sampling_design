# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. length model)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Predict habitats over Geographe Bay using BRUV observations and GAMs.
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    September 2024

# -----------------------------------------------------------------------------

# -- For the spatial distribution model, we need habitat variables which will
# -- help to predict our fish's distribution. We don't have habitat data for the
# -- whole of the area (we only have point observation from BRUVs), so we need
# -- to predict it using bathymetry data. We will use different bathymetry
# -- resolutions to make a prediction: 250m resolution and LiDAR bathymetry
# -- and see what the model prefers.

# -----------------------------------------------------------------------------

# Load libraries
library(reshape2)
library(mgcv)
library(ggplot2)
library(terra)
library(raster)
library(predicts)
library(tidyverse)
library(sf)
library(nlraa)
library(patchwork)
library(FSSgam)

rm(list = ls())

ext <- c(115.035, 115.68012207, -33.69743897, -33.35) # Extent to crop layers to


## Files used in this script --------------------------------------------------

file_obs_hab            <- "data/tidy/2024_geographe_spabal_tidy_habitat.rds" # Observed habitat values
file_250m_bathy_deriv   <- "data/spatial/rasters/2014_geographe_bathymetry-derivatives.RDS"
#file_lidar_bathy        <- "data/spatial/rasters/LiDAR_geo_compressed.tif" # Too big for git

file_predictors         <- "data/spatial/rasters/geographe_bathy_predictors.rds"


# 250m habitat predictions ----------------------------------------------------

# -- Let's start with the 250m predictions.

## Load data ------------------------------------------------------------------

habi <- readRDS(file_obs_hab) %>% # Observed habitat
  pivot_longer(cols = c("Macroalgae", "Unconsolidated", "Consolidated", "Seagrasses"),
               values_to = "number", names_to = "taxa") %>%
  dplyr::mutate(depth = as.double(depth) * -1) %>%
  filter(longitude > 115.3 & OpCode != '12RC9') %>%
  na.omit() %>%
  mutate(taxa = case_when(
    taxa %in% c('Macroalgae', 'Consolidated') ~ 'Reef',
    TRUE ~ taxa
  )) %>%
  glimpse()


names(habi)
pred.vars <- c("gadepth", "aspect", "roughness", "detrended")

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)
round(cor(habi[ , pred.vars]), 2) # All good, no high correlations

# Check to make sure Response vector has not more than 80% zeros
unique.vars = unique(as.character(habi$taxa))

unique.vars.use = character()
for(i in 1:length(unique.vars)){
  temp.dat = habi[which(habi$taxa == unique.vars[i]),]
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

for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  use.dat <- habi[habi$taxa == resp.vars[i],]
  use.dat   <- as.data.frame(use.dat)

  # Basic model to compare all the other combinations with
  Model1  <- gam(cbind(number, (total_pts - number)) ~ # Success and failure counts for each habitat
                   s(gadepth, bs = 'cr', k = 3),
                 family = binomial(link = "logit"),  data = use.dat)

  # Generate variable combinations to test in models
  model.set <- generate.model.set(use.dat = use.dat,
                                  test.fit = Model1,
                                  pred.vars.cont = pred.vars,
                                  cyclic.vars = c("aspect"),
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
# Reef : detrended + gadepth + roughness
# Unconsolidated : detrended + gadepth + roughness
# Seagrass: detrended + gadepth + roughness


## reload & format the data ---------------------------------------------------

habi <- readRDS(file_obs_hab) %>%
  filter(OpCode != "12RC9") %>% # Also removing one OpCode with point count issues.
  mutate(Reef = Macroalgae + Consolidated) %>%  # Create a new column by summing Macroalgae and Rock
  dplyr::select(-Macroalgae, -Consolidated) %>% # Optionally remove the old columns if no longer needed
  glimpse()


preds <- readRDS(file_250m_bathy_deriv); plot(preds)

preddf <- as.data.frame(preds, xy = TRUE, na.rm = TRUE)


## Make models for each habitat and predict extent ----------------------------

# Save each model (based on top models above)

# Reef
m_reef_250m <- gam(cbind(Reef, total_pts - Reef) ~
                     #s(aspect, k = 5, bs = "cr")  +
                     s(detrended, k = 5, bs = "cr")  +
                     s(gadepth, k = 5, bs = "cr") +
                     s(roughness, k = 5, bs = "cr"),
                   data = habi, method = "REML", family = binomial("logit"))
summary(m_reef_250m)
plot(m_reef_250m, pages = 1, residuals = T, cex = 5)

# Seagrass
m_seagrass_250m <- gam(cbind(Seagrasses, total_pts - Seagrasses) ~
                         #s(aspect, k = 5, bs = "cr")  +
                         s(detrended, k = 5, bs = "cr")  +
                         s(gadepth, k = 5, bs = "cr") +
                         s(roughness, k = 5, bs = "cr"),
                       data = habi, method = "REML", family = binomial("logit"))
summary(m_seagrass_250m)
plot(m_seagrass_250m, pages = 1, residuals = T, cex = 5)

# Sand
m_sand_250m <- gam(cbind(Unconsolidated, total_pts - Unconsolidated) ~
                     #s(aspect, k = 5, bs = "cr")  +
                     s(detrended, k = 5, bs = "cr")  +
                     s(gadepth, k = 5, bs = "cr") +
                     s(roughness, k = 5, bs = "cr"),
                   data = habi, method = "REML", family = binomial("logit"))
summary(m_sand_250m)
plot(m_sand_250m, pages = 1, residuals = T, cex = 5)


## predict, rasterise and plot ------------------------------------------------

preddf <- cbind(preddf,
                "preef_250m"     = predict(m_reef_250m, preddf, type = "response", se.fit = T),
                "psand_250m"     = predict(m_sand_250m, preddf, type = "response", se.fit = T),
                "pseagrass_250m" = predict(m_seagrass_250m, preddf, type = "response", se.fit = T)) %>%
  glimpse()

prasts_250m <- rast(preddf %>% dplyr::select(x, y, preef_250m.fit, psand_250m.fit, pseagrass_250m.fit),
               crs = crs(preds)) %>%
  crop(ext)
plot(prasts_250m)

# Save the predictions for use in the length model
saveRDS(preddf, paste0("outputs/Length_comparison/2024_geographe_predicted-habitat_250m.rds"))
writeRaster(prasts_250m, filename = "outputs/Length_comparison/2024_geographe_predicted-habitat_250m.tif", overwrite = TRUE)



# LiDAR habitat predictions ---------------------------------------------------

# -- We need to create another predicted habitat map using the finer-scale LiDAR
# -- bathymetry data. We'll use the exact same method as above.

## Load data ------------------------------------------------------------------

habi <- readRDS(file_obs_hab) %>% # Observed habitat
  pivot_longer(cols = c("Macroalgae", "Unconsolidated", "Consolidated", "Seagrasses"),
               values_to = "number", names_to = "taxa") %>%
  dplyr::mutate(depth = as.double(depth) * -1) %>%
  filter(OpCode != '12RC9') %>%
  na.omit() %>%
  mutate(taxa = case_when(
    taxa %in% c('Macroalgae', 'Consolidated') ~ 'Reef',
    TRUE ~ taxa
  )) %>%
  ungroup() %>%
  mutate(ID = row_number()) %>%
  glimpse()

habi_sp <- st_as_sf(habi, coords = c("longitude", "latitude"), crs = 4326); plot(habi_sp)

lidar <- rast(file_predictors)[[1:3]]; names(lidar) <- c("bathy_lidar", "roughness_lidar", "aspect_lidar"); plot(lidar)

# Extract LiDAR values and include the ID
lidar_values <- extract(lidar, habi_sp, data.frame = TRUE) %>%
  mutate(ID = row_number()) %>%
  na.omit() %>%  # Optionally, remove NA values
  glimpse

# Now merge using the ID
habi <- habi %>%
  left_join(lidar_values, by = "ID") %>%
  na.omit() %>%
  glimpse()

names(habi)
pred.vars <- c("bathy_lidar", "aspect_lidar", "roughness_lidar")

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)
round(cor(habi[ , pred.vars]), 2) # All good, no high correlations

# Check to make sure Response vector has not more than 80% zeros
unique.vars = unique(as.character(habi$taxa))

unique.vars.use = character()
for(i in 1:length(unique.vars)){
  temp.dat = habi[which(habi$taxa == unique.vars[i]),]
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
  use.dat <- habi[habi$taxa == resp.vars[i],]
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
# Reef :              bathy + roughness
# Unconsolidated :    bathy + roughness
# Seagrasses:         bathy + roughness


## reload & format the data ---------------------------------------------------

habi <- readRDS(file_obs_hab) %>% # Observed habitat
  mutate(Reef = Macroalgae + Consolidated) %>%  # Create a new column by summing Macroalgae and Rock
  dplyr::select(-Macroalgae, -Consolidated) %>% # Optionally remove the old columns if no longer needed
  dplyr::mutate(depth = as.double(depth) * -1) %>%
  filter(OpCode != '12RC9') %>%
  na.omit() %>%
  ungroup() %>%
  mutate(ID = row_number()) %>%
  glimpse()

habi_sp <- st_as_sf(habi, coords = c("longitude", "latitude"), crs = 4326); plot(habi_sp)

lidar <- rast(file_predictors)[[1:3]]; names(lidar) <- c("bathy_lidar", "roughness_lidar", "aspect_lidar"); plot(lidar)

# Extract LiDAR values and include the ID
lidar_values <- extract(lidar, habi_sp, data.frame = TRUE) %>%
  mutate(ID = row_number()) %>%
  na.omit() %>%  # Optionally, remove NA values
  glimpse

# Now merge using the ID
habi <- habi %>%
  left_join(lidar_values, by = "ID") %>%
  na.omit() %>%
  glimpse()

preds <- rast(file_predictors); plot(preds)

preddf <- as.data.frame(preds, xy = TRUE, na.rm = TRUE)


## Make models for each habitat and predict extent ----------------------------

# Save each model (based on top models above)

# Reef
m_reef_lidar <- gam(cbind(Reef, total_pts - Reef) ~
                      #s(aspect_lidar, k = 5, bs = "cr")  +
                      s(roughness_lidar, k = 5, bs = "cr") +
                      s(bathy_lidar, k = 5, bs = "cr"),
                    data = habi, method = "REML", family = binomial("logit"))
summary(m_reef_lidar)
plot(m_reef_lidar, pages = 1, residuals = T, cex = 5)

# Sand
m_sand_lidar <- gam(cbind(Unconsolidated, total_pts - Unconsolidated) ~
                      #s(aspect_lidar, k = 5, bs = "cr")  +
                      s(roughness_lidar, k = 5, bs = "cr") +
                      s(bathy_lidar, k = 5, bs = "cr"),
                      data = habi, method = "REML", family = binomial("logit"))
summary(m_sand_lidar)
plot(m_sand_lidar, pages = 1, residuals = T, cex = 5)

# Seagrass
m_seagrass_lidar <- gam(cbind(Seagrasses, total_pts - Seagrasses) ~
                          #s(aspect_lidar, k = 5, bs = "cr")  +
                          s(roughness_lidar, k = 5, bs = "cr") +
                          s(bathy_lidar, k = 5, bs = "cr"),
                          data = habi, method = "REML", family = binomial("logit"))
summary(m_seagrass_lidar)
plot(m_seagrass_lidar, pages = 1, residuals = T, cex = 5)


## predict, rasterise and plot ------------------------------------------------

#preddf <- cbind(preddf,
#                "preef_lidar"     = predict(m_reef_lidar, preddf, type = "response", se.fit = T),
#                "psand_lidar"     = predict(m_sand_lidar, preddf, type = "response", se.fit = T),
#                "pseagrass_lidar" = predict(m_seagrass_lidar, preddf, type = "response", se.fit = T)) %>%
#  glimpse()

#prasts_lidar <- rast(preddf %>% dplyr::select(x, y, preef_lidar.fit, psand_lidar.fit, pseagrass_lidar.fit),
#               crs = crs(preds)) %>%
#  crop(ext); plot(prasts_lidar)

# Save the predictions for use in the length model
#saveRDS(preddf, paste0("outputs/Length_comparison/2024_geographe_predicted-habitat_LiDAR.rds"))
#writeRaster(prasts_lidar, filename = "outputs/Length_comparison/2024_geographe_predicted-habitat_LiDAR.tif", overwrite = TRUE)


## Making a file stack for all environmental data -----------------------------

bathy <- readRDS(file_predictors); plot(bathy) # aspect and bathy are switched but we'll fix that in a sec

#predictors <- list(bathy_lidar = bathy[[1]], roughness_lidar = bathy[[2]], aspect_lidar = bathy[[3]], # if including LiDAR data
#                   reef_lidar = prasts_lidar[[1]], sand_lidar = prasts_lidar[[2]], seagrass_lidar = prasts_lidar[[3]],
#                   aspect_250 = bathy[[4]], depth_250 = bathy[[5]], roughness_250 = bathy[[6]], detrended_250 = bathy[[7]],
#                   reef_250m = prasts_250m[[1]], sand_250m = prasts_250m[[2]], seagrass_250m = prasts_250m[[3]]
#                   )

predictors <- list( # if not including LiDAR data
  aspect_250 = bathy[[2]], depth_250 = bathy[[1]], roughness_250 = bathy[[3]], detrended_250 = bathy[[4]],
  reef_250m = prasts_250m[[1]], sand_250m = prasts_250m[[2]], seagrass_250m = prasts_250m[[3]]
)

for (i in 1:length(predictors)) { # Unifying extents and resampling so that they can be stacked
  reference_raster <- predictors[[1]] # Use the raster with the highest resolution as a reference. Otherwise the re-sampling will apply lower resolutions to everything.
  predictors[[i]] <- terra::resample(predictors[[i]], reference_raster, method = "bilinear")
  cat(paste("New extent of", i, ":", ext(predictors[[i]]), "\n"))
}

list2env(predictors, envir = .GlobalEnv) # Put the elements of the list in the environment

predictors <- c(
  #bathy_lidar, roughness_lidar, aspect_lidar,          # LiDAR data
  depth_250, aspect_250, roughness_250, detrended_250, # 250m data
  #reef_lidar, sand_lidar, seagrass_lidar,              # Habitat predictions using LiDAR bathymetry
  reef_250m, sand_250m, seagrass_250m)                 # Habitat predictions using 250m bathymetry

names(predictors) <- c(#"bathy_lidar", "roughness_lidar", "aspect_lidar", # Add clear names
                       "bathy_250m", "aspect_250m", "roughness_250m", "detrended_250m",
                       #"reef_lidar", "sand_lidar", "seagrass_lidar",
                       "reef_250m", "sand_250m", "seagrass_250m")
plot(predictors)

saveRDS(predictors, file = "data/spatial/rasters/geographe_all_predictors.rds") # May take a while if running with LiDAR
# ^ needed for 03_gam

### END ###
