### ARCHIVE FOR SAMPLING DESIGN COMPARISON ###

# from 03_gam -----------------------------------------------------------------
# non-streamlined model construction
## PINK SNAPPER ---------------------------------------------------------------

# Select the species you want to model
model_spp <- c(#"Choerodon rubescens",
  "Chrysophrys auratus"
  #"Glaucosoma hebraicum",
  #"Nemadactylus valenciennesi",
  #"Pseudocaranx spp"
)

dat_PS <- dat[dat$full_spp == model_spp,] %>%
  group_by(latitude, longitude, ID,
           bathy_250m, aspect_250m, roughness_250m, detrended_250m, reef_250m, sand_250m, seagrass_250m,
           bathy_lidar, roughness_lidar, aspect_lidar, reef_lidar, sand_lidar, seagrass_lidar
  ) %>%
  summarise(
    count_mature = sum(count_mature, na.rm = TRUE),
    count_immature = sum(count_immature, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  glimpse()

length(unique(dat$opcode))

### Fitting GAM, mature fish --------------------------------------------------

outdir <- ("outputs/Length_comparison/mature/")

use.dat <- na.omit(dat_PS)
use.dat <- as.data.frame(use.dat)
Model1 <- gam(count_mature ~ # Count of Mature individuals
                s(seagrass_250m, bs = 'cr'),
              family = poisson(link = "log"),  data = use.dat)

# Generate sets of 3 variable combinations
model.set <- generate.model.set(use.dat = use.dat,
                                test.fit = Model1,
                                pred.vars.cont = pred.vars,
                                cyclic.vars = c('aspect_250m', 'aspect_lidar'),
                                k = 5,
                                cov.cutoff = 0.7,
                                max.predictors = 3)

# Fit the sets of variables to the data
out.list <- fit.model.set(model.set,
                          max.models = 600,
                          parallel = T)
names(out.list)

# From all the combinations of predictors, we'll now look at how they compare to each other
out.list$failed.models # no failed models, all good
mod.table           <- out.list$mod.data.out
mod.table           <- mod.table[order(mod.table$AICc), ] # AICc value, the lower, the better the performance of the model
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc) # wi.AICc is the model weight, which compares the proportion of evidence for each model. The higher the better.
out.i               <- mod.table[which(mod.table$delta.AICc <= 2), ] # This compares the models that have comparable AICc values (here, the 2 best models)
out.i

best.model.name     <- as.character(out.i$modname); print(best.model.name)
saveRDS(best.model.name, file = "data/rmd/gam_best_model_PS.rds") # for Rmarkdown

png(file = paste(outdir, "mature_mod_fits.png", sep = ""))

# We'll plot to see the spread of data in the model. Does it make biological sense? etc. very few observations where you expect few etc.
for(m in 1:nrow(out.i)){
  best.model.name <- as.character(out.i$modname[m])
  
  png(file = paste(outdir, m, "PinkSnapper_mature", "mod_fits.png", sep = ""))
  if(best.model.name != "null"){
    par(mfrow = c(3, 1), mar = c(9, 4, 3, 1))
    best.model = out.list$success.models[[best.model.name]]
    plot(best.model, all.terms = T, pages = 1, residuals = T, pch = 16)
    mtext(side = 2, text = "mature", outer = F)}
  dev.off()
}

# Using the best models above, fitting a GAM to mature fish
out.i[[1]]

# TESTING POISSON
m_mature <- gam(count_mature ~
                  #s(bathy_lidar, k = 3, bs = "cr") +
                  #s(aspect_lidar, k = 3, bs = "cc") +
                  #s(roughness_lidar, k = 3, bs = "cr") +
                  
                  #s(bathy_250m, k = 3, bs = "cr") +
                  #s(aspect_250m, k = 3, bs = "cc") +
                  #s(roughness_250m, k = 3, bs = "cr") +
                  #s(detrended_250m, k = 3, bs = "cr") +
                  
                  s(reef_lidar, k = 3, bs = "cr") +
                  s(sand_lidar, k = 3, bs = "cr") +
                  #s(seagrass_lidar, k = 3, bs = "cr") +
                  
                  #s(reef_250m, k = 3, bs = "cr"),
                  s(seagrass_250m, k = 3, bs = "cr"),
                #s(sand_250m, k = 3, bs = "cr"),
                
                data = dat, method = "REML", family = poisson(link = "log"))
summary(m_mature)
plot(m_mature, pages = 1, residuals = T, cex = 5)

# TESTING TWEEDIE
m_mature <- gam(count_mature ~
                  #s(bathy_lidar, k = 3, bs = "cr") +
                  #s(aspect_lidar, k = 3, bs = "cc") +
                  #s(roughness_lidar, k = 3, bs = "cr") +
                  
                  #s(bathy_250m, k = 3, bs = "cr") +
                  #s(aspect_250m, k = 3, bs = "cc") +
                  #s(roughness_250m, k = 3, bs = "cr") +
                  #s(detrended_250m, k = 3, bs = "cr") +
                  
                  s(reef_lidar, k = 3, bs = "cr") +
                  s(sand_lidar, k = 3, bs = "cr") +
                  #s(seagrass_lidar, k = 3, bs = "cr") +
                  
                  #s(reef_250m, k = 3, bs = "cr"),
                  s(seagrass_250m, k = 3, bs = "cr"),
                #s(sand_250m, k = 3, bs = "cr"),
                
                data = dat, method = "REML", family = tw())
summary(m_mature)
plot(m_mature)
### Remove variable values that are out of range ------------------------------

sel_rast <- pred_rast[[names(m_mature$model[2:length(names(m_mature$model))])]]; plot(sel_rast) # Subset the environmental raster to only the variables in the model


# Cropping the abundance prediction to observed env. variables only
extracted_values <- terra::extract(sel_rast, dat[, c("longitude", "latitude")]) %>%
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
  mask_layer <- ifelse(current_layer_values < layer_min | current_layer_values > layer_max, NA, current_layer_values)
  
  # Assign the mask_layer as a SpatRaster layer
  mask_raster[[i]] <- rast(nrows = nrow(sel_rast), ncols = ncol(sel_rast), vals = mask_layer, ext = ext(sel_rast), crs = crs(sel_rast))
}
plot(mask_raster)

# Keep only the cells that exist in all layers - this is to avoid extrapolating
valid_mask <- !is.na(mask_raster[[1]]) & !is.na(mask_raster[[2]]) & !is.na(mask_raster[[3]]); plot(valid_mask)
masked_raster <- mask(mask_raster, valid_mask)
masked_raster[!valid_mask] <- NA; plot(masked_raster)

ras_crop <- as.data.frame(masked_raster, xy = TRUE, na.rm = TRUE) %>%
  glimpse()


### Predicting mature distribution --------------------------------------------

predicted_abundance <- cbind(ras_crop, "p_mature" = mgcv::predict.gam(m_mature, ras_crop, type = "response", se.fit = T)) %>%
  glimpse()

p_mature <- rast(predicted_abundance)

par(mfrow=c(1,1))
plot(p_mature$p_mature.fit, pch = 19,
     main = paste('Pink Snapper -', out.i$modname[1]),
     xlab = "Longitude", ylab = "Latitude",
     xlim = c(115.32, 115.5426), ylim = c(-33.60944, -33.40019))

# Save predictions as rasters
writeRaster(c(p_mature$p_mature.fit, p_mature$p_mature.se.fit), filename = "outputs/Length_comparison/predicted_abundance_rasters/mature_predicted_abundance_PS.tif", overwrite = TRUE)




## OTHER SPP ------------------------------------------------------------------

# Select the species you want to model
model_spp <- c("Choerodon rubescens"
               #"Chrysophrys auratus",
               #"Glaucosoma hebraicum"
               #"Nemadactylus valenciennesi"
               #"Pseudocaranx spp"
)

dat_BG <- dat[dat$full_spp == model_spp,] %>%
  group_by(latitude, longitude, ID,
           bathy_250m, aspect_250m, roughness_250m, detrended_250m, reef_250m, sand_250m, seagrass_250m,
           bathy_lidar, roughness_lidar, aspect_lidar, reef_lidar, sand_lidar, seagrass_lidar
  ) %>%
  summarise(
    count_mature = sum(count_mature, na.rm = TRUE),
    count_immature = sum(count_immature, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  ungroup() %>%
  glimpse()

ggplot(dat_BG, aes(x=longitude, y=latitude)) +
  geom_point(aes(size=count_mature), col = 'navy')


### Fitting GAM, mature fish --------------------------------------------------

outdir <- ("outputs/Length_comparison/mature/")

use.dat <- na.omit(dat_BG)
use.dat <- as.data.frame(use.dat)
Model1 <- gam(count_mature ~ # Count of Mature individuals
                s(seagrass_250m, bs = 'cr'),
              family = poisson(link = "log"),  data = use.dat)

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
out.list$failed.models # no failed models, all good
mod.table           <- out.list$mod.data.out
mod.table           <- mod.table[order(mod.table$AICc), ] # AICc value, the lower, the better the performance of the model
mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc) # wi.AICc is the model weight, which compares the proportion of evidence for each model. The higher the better.
out.i               <- mod.table[which(mod.table$delta.AICc <= 2), ] # This compares the models that have comparable AICc values (here, the 2 best models)
out.i
best.model.name     <- as.character(out.i$modname); print(best.model.name)
saveRDS(best.model.name, file = "data/rmd/gam_best_model_spp_complex.rds")

png(file = paste(outdir, "mature_mod_fits.png", sep = ""))

# We'll plot to see the spread of data in the model. Does it make biological sense? etc. very few observations where you expect few etc.
for(m in 1:nrow(out.i)){
  best.model.name <- as.character(out.i$modname[m])
  
  png(file = paste(outdir, m, "mature", "mod_fits.png", sep = ""))
  if(best.model.name != "null"){
    par(mfrow = c(3, 1), mar = c(9, 4, 3, 1))
    best.model = out.list$success.models[[best.model.name]]
    plot(best.model, all.terms = T, pages = 1, residuals = T, pch = 16)
    mtext(side = 2, text = "mature", outer = F)}
  dev.off()
}

# Using the best models above, fitting a GAM to mature fish
out.i[[1]]
m_mature <- gam(count_mature ~
                  #s(bathy_lidar, k = 3, bs = "cr") +
                  #s(aspect_lidar, k = 3, bs = "cr") +
                  #s(roughness_lidar, k = 3, bs = "cr") +
                  
                  #s(bathy_250m, k = 3, bs = "cr") +
                  #s(aspect_250m, k = 3, bs = "cr") +
                  s(roughness_250m, k = 3, bs = "cr") +
                  #s(detrended_250m, k = 3, bs = "cr") +
                  
                  s(reef_lidar, k = 3, bs = "cr") +
                  s(sand_lidar, k = 3, bs = "cr"),
                #s(seagrass_lidar, k = 3, bs = "cr"),
                
                #s(reef_250m, k = 3, bs = "cr") +
                #s(seagrass_250m, k = 3, bs = "cr"),
                #s(sand_250m, k = 3, bs = "cr"),
                
                data = dat, method = "REML", family = poisson(link = "log"))
summary(m_mature)
plot(m_mature, pages = 1, residuals = T, cex = 5)


## TWEEDIE

out.i[[1]]
m_mature <- gam(count_mature ~
                  #s(bathy_lidar, k = 3, bs = "cr") +
                  #s(aspect_lidar, k = 3, bs = "cr") +
                  #s(roughness_lidar, k = 3, bs = "cr") +
                  
                  #s(bathy_250m, k = 3, bs = "cr") +
                  #s(aspect_250m, k = 3, bs = "cr") +
                  s(roughness_250m, k = 3, bs = "cr") +
                  #s(detrended_250m, k = 3, bs = "cr") +
                  
                  s(reef_lidar, k = 3, bs = "cr") +
                  s(sand_lidar, k = 3, bs = "cr"),
                #s(seagrass_lidar, k = 3, bs = "cr"),
                
                #s(reef_250m, k = 3, bs = "cr") +
                #s(seagrass_250m, k = 3, bs = "cr"),
                #s(sand_250m, k = 3, bs = "cr"),
                
                data = dat, method = "REML", family = tw())
summary(m_mature)
plot(m_mature, pages = 1, residuals = T, cex = 5)

### Remove variable values that are out of range ------------------------------

sel_rast <- pred_rast[[names(m_mature$model[2:length(names(m_mature$model))])]]; plot(sel_rast) # Subset the environmental raster to only the variables in the model


# Cropping the abundance prediction to observed env. variables only
extracted_values <- terra::extract(sel_rast, dat[, c("longitude", "latitude")]) %>%
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
  mask_layer <- ifelse(current_layer_values < layer_min | current_layer_values > layer_max, NA, current_layer_values)
  
  # Assign the mask_layer as a SpatRaster layer
  mask_raster[[i]] <- rast(nrows = nrow(sel_rast), ncols = ncol(sel_rast), vals = mask_layer, ext = ext(sel_rast), crs = crs(sel_rast))
}
plot(mask_raster)

# Keep only the cells that exist in all layers - this is to avoid extrapolating
valid_mask <- !is.na(mask_raster[[1]]) & 
  !is.na(mask_raster[[2]]) #& !is.na(mask_raster[[3]])
plot(valid_mask)
masked_raster <- mask(mask_raster, valid_mask)
masked_raster[!valid_mask] <- NA; plot(masked_raster)

ras_crop <- as.data.frame(masked_raster, xy = TRUE, na.rm = TRUE) %>%
  glimpse()


### Predicting mature distribution --------------------------------------------

predicted_abundance <- cbind(ras_crop, "p_mature" = mgcv::predict.gam(m_mature, ras_crop, type = "response", se.fit = T)) %>%
  glimpse()

p_mature <- rast(predicted_abundance)

par(mfrow=c(1,1))
plot(p_mature$p_mature.fit, pch = 19,
     main = paste(model_spp, out.i$modname[1]),
     xlab = "Longitude", ylab = "Latitude",
     xlim = c(115.32, 115.5426), ylim = c(-33.60944, -33.40019))

# Save predictions as rasters
writeRaster(c(p_mature$p_mature.fit, p_mature$p_mature.se.fit), filename = "outputs/Length_comparison/predicted_abundance_rasters/mature_predicted_abundance_other.tif", overwrite = TRUE)


