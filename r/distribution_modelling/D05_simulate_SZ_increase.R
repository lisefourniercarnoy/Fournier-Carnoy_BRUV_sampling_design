# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare abundance detected)
# Data:    Modelled data 
# Task:    Simulate a Sanctuary Zone within which there is a __% increase in population
# Author:  Lise Fournier-Carnoy
# Date:    September 2024

# -----------------------------------------------------------------------------

library(tidyverse) # for manipulating data
library(sf) # for manipulating shapefiles
library(terra) # for manipulating rasters
library(tidyterra) # for manipulating rasters
library(stars)


# Clear memory
rm(list=ls())

study_site <- "waatern"

## Files used in this script --------------------------------------------------

file_land           <- "QGIS layers/clean/wadandi_land_highres.shp"
file_sim_SZ         <- paste0("QGIS layers/clean/", study_site, "_simulated_SZ.shp")
file_samp_area      <- paste0("QGIS layers/produced_from_code/D04_", study_site, "_sampling_area_deeper_than_7m.rds")
file_detrend_val    <- paste0("data/tidy/D01_", study_site, "_bathy_predictors.rds")
file_strata         <- paste0("outputs/distribution_modelling_outputs/D04_", study_site, "_sample_area_detrended_strata.rds")

# all the sampling designs
file_simp_spabal_sd       <- paste0("outputs/distribution_modelling_outputs/D04_sampling_designs/D04_", study_site, "_simple_spatially_balanced_designs.rds")

file_str_spabal_sd_2525   <- paste0("outputs/distribution_modelling_outputs/D04_sampling_designs/D04_", study_site, "_stratified_spatial_balance_designs_25_25.rds")
file_str_spabal_sd_2030   <- paste0("outputs/distribution_modelling_outputs/D04_sampling_designs/D04_", study_site, "_stratified_spatial_balance_designs_20_30.rds")

file_pref_sd_2525         <- paste0("outputs/distribution_modelling_outputs/D04_sampling_designs/D04_", study_site, "_preferential_designs_25_25.rds")
file_pref_sd_2030         <- paste0("outputs/distribution_modelling_outputs/D04_sampling_designs/D04_", study_site, "_preferential_designs_20_30.rds")

file_clump_sd             <- paste0("outputs/distribution_modelling_outputs/D04_sampling_designs/D04_", study_site, "_clustered_designs.rds")


crs <- 7850
crs_rast <- "epsg:7850"

choose_from <- list.files(
  path = "outputs/distribution_modelling_outputs/D03_distribution_rasters/", 
  pattern = paste0(study_site, "_predicted_distribution_raster_"), 
  full.names = TRUE
)
choose_from

file_dist_pred <- choose_from[[3]] # choose the fish prediction file

## Load files -----------------------------------------------------------------

ext <- c(115.035, 115.68012207, -33.69743897, -33.20243897)

aus <- st_read(file_land) %>% # Land
  st_transform(crs) %>%
  glimpse()

sim_sz <- st_read(file_sim_SZ) %>% # Simulated SZ
  st_make_valid() %>%
  st_transform(crs) %>%
  st_difference(st_union(aus)) %>%
  glimpse()

spp <- gsub(paste0("D03_", study_site, "_predicted_distribution_raster_|\\.tif"), "", basename(file_dist_pred))
samp_area <- readRDS(file_samp_area)
strata <- readRDS(file_strata) %>% 
  mutate(strata = case_when(
    detrended == '0' ~ "low",
    detrended == '1' ~ "mid",
    detrended == '2' ~ "high",
  )) %>% 
  glimpse()
detrend <- readRDS(file_detrend_val); detrend <- detrend$detrended_250m; plot(detrend)

fish <- rast(file_dist_pred) %>% project(crs_rast); plot(fish)

## Simulate an increase in the zone -------------------------------------------
plot(sim_sz$geometry)
inner_ntz <- st_buffer(sim_sz, -1000) # 1km inner buffer
ntz_edge_effect <- st_difference(sim_sz, inner_ntz)
plot(ntz_edge_effect$geometry, col = "blue"); plot(inner_ntz$geometry, col = "red", add =T)

masked_raster_inner <- mask(fish, inner_ntz); plot(masked_raster_inner$p_fish.fit)
masked_raster_edge <- mask(fish, ntz_edge_effect); plot(masked_raster_edge, add = T)


# Increase the values in the masked raster by 80%
inner_increased_raster_fit <- masked_raster_inner$p_fish.fit * 1.8; plot(inner_increased_raster_fit) # Fitted values
inner_increased_raster_se <- masked_raster_inner$p_fish.se.fit * 1.8; plot(inner_increased_raster_se) # Standard errors

edge_increased_raster_fit <- masked_raster_edge$p_fish.fit * 1.32; plot(edge_increased_raster_fit) # Fitted values
edge_increased_raster_se <- masked_raster_edge$p_fish.se.fit * 1.32; plot(edge_increased_raster_se) # Standard errors


# Combine the increased SZ and the rest of the raster
final_raster_fit <- fish$p_fish.fit  # Start with the original raster
final_raster_se <- fish$p_fish.se.fit  # Start with the original raster

# Replace the SZ values of the original raster with the increased values
final_raster_fit[!is.na(inner_increased_raster_fit)] <- inner_increased_raster_fit[!is.na(inner_increased_raster_fit)]
final_raster_se[!is.na(inner_increased_raster_se)] <- inner_increased_raster_se[!is.na(inner_increased_raster_se)]

final_raster_fit[!is.na(edge_increased_raster_fit)] <- edge_increased_raster_fit[!is.na(edge_increased_raster_fit)]
final_raster_se[!is.na(edge_increased_raster_se)] <- edge_increased_raster_se[!is.na(edge_increased_raster_se)]


plot(final_raster_fit)
plot(final_raster_se)


# smooth abundance edge effect test -------------------------------------------

# want to make a smooth increase in abundance from out to in SZ

sz_abundance_fit <- mask(fish$p_fish.fit, sim_sz); plot(sz_abundance_fit)
sz_abundance_se <- mask(fish$p_fish.se.fit, sim_sz); plot(sz_abundance_se)

sz_border <- st_boundary(sim_sz); plot(sz_border)

sz_centroid <- as.points(sz_abundance_fit) %>% st_as_sf()
plot(sz_centroid)

# Compute smallest distance from each point to the border
min_dist <- as.numeric(st_distance(sz_centroid, sz_border))
hist(min_dist)

# what the abundance increase looks like
test <- data.frame(dist = 0:max(min_dist))
max <- 0.8 # max abundance increase
mid <- 800 # midpoint of the curve
scale <- 500 # smoothness of curve
test$abundance_increase <- max * (tanh((test$dist - mid)/scale)/2 + 0.5)
plot(test)
min(test$abundance_increase)

# assign these distances to the raster
plot(sz_abundance_fit)
sz_abundance_fit$dist <- NA
sz_abundance_fit$dist[!is.na(sz_abundance_fit$p_fish.fit), ] <- min_dist
plot(sz_abundance_fit)

plot(sz_abundance_se)
sz_abundance_se$dist <- NA
sz_abundance_se$dist[!is.na(sz_abundance_se$p_fish.se.fit), ] <- min_dist
plot(sz_abundance_se)

# make the smooth abundance increase with logistic curve
dist_val <- values(sz_abundance_fit$dist)
fish_val <- values(sz_abundance_fit$p_fish.fit); plot(fish_val)
sz_val_fit <- fish_val * (1+ max * ((tanh((dist_val - mid)/scale)/2 + 0.5))); plot(sz_val_fit)
sz_abundance_fit$new_ab <- NA
sz_abundance_fit$new_ab <- sz_val_fit
plot(sz_abundance_fit)

dist_val <- values(sz_abundance_se$dist)
fish_val <- values(sz_abundance_se$p_fish.se.fit); plot(fish_val)
sz_val_se <- fish_val * (1+max * ((tanh((dist_val - mid)/scale)/2 + 0.5))); plot(sz_val_se)
sz_abundance_se$new_ab <- NA
sz_abundance_se$new_ab <- sz_val_se
plot(sz_abundance_se)

# join the 'outside-NTZ' raster back in
final_raster <- c(fish$p_fish.fit, fish$p_fish.se.fit); plot(final_raster)

final_raster$p_fish.fit[!is.na(sz_abundance_fit$new_ab)] <- sz_abundance_fit$new_ab[!is.na(sz_abundance_fit$new_ab)]
final_raster$p_fish.se.fit[!is.na(sz_abundance_se$new_ab)] <- sz_abundance_se$new_ab[!is.na(sz_abundance_se$new_ab)]

plot(final_raster)

# end test



# Full raster with both fitted values and standard errors
# final_raster <- c(final_raster_fit, final_raster_se)

saveRDS(final_raster, paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_smooth_edge_effect_", gsub(" ", "_", spp), ".rds"))

max_ab <- max(na.omit(values(crop(final_raster$p_fish.fit, ext(samp_area)))))
ggplot() +
  labs(title = "Location of the simulated SZ, \nwith 80% increase in abundance within") +
  geom_spatraster(data = final_raster$p_fish.fit) + # Predicted abundance
  geom_sf(data = aus) + # Land
  geom_sf(data = sim_sz, colour = "green", fill = NA, linewidth = 1) + # Simulated SZ
  geom_sf(data = samp_area, colour = "red", fill = NA, linewidth = 1) + # Sample area
  coord_sf(crs = crs_rast, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_fill_gradient(limits = c(0, max_ab), low = "white", high = "darkblue") + # Adjust raster color scale
  theme_minimal()


### Extract abundance values from sampling design locations -------------------

final_names <- c("design_id", "in_SZ", "strata", "fit_values", "se_values", "detrend_val","geometry", "SD")



#### Simple stratified sampling -----------------------------------------------

simp_spabal_sd <- readRDS(file_simp_spabal_sd)
names(simp_spabal_sd) <- paste0("simple_spabal_", 1:1000) # name each list element (simulated design)

extracted_values <- list()

for (i in seq_along(simp_spabal_sd)) {
  print(i)
  # design_name <- names(simp_spabal_sd)[[i]]
  # points <- simp_spabal_sd[[i]]  # Keep as sf
  # 
  # # Convert to SpatVector just for extraction
  # points_vect <- vect(points)
  # 
  # # Extract values and take only the second column (value), not the ID
  # fit_values <- suppressWarnings(terra::extract(final_raster$p_fish.fit, points_vect)[, 2])
  # se_values  <- suppressWarnings(terra::extract(final_raster$p_fish.se.fit, points_vect)[, 2])
  # detrend_val <- suppressWarnings(terra::extract(detrend$detrended, points_vect)[, 2]) # add detrended bathymetry value
  # points <- st_join(points, strata[, c("strata")], left = TRUE) # add strata type
  # 
  # # Add extracted values directly to the sf object
  # points$fit_values <- fit_values
  # points$se_values <- se_values
  # points$detrend_val <- detrend_val
  # 
  # points$design_id <- design_name # to identify which simulation
  # 
  # # Store the result with geometry intact
  # extracted_values[[design_name]] <- points
  design_name <- names(simp_spabal_sd)[i]
  points <- st_transform(simp_spabal_sd[[i]], crs)  # keep as sf
  points_vect <- vect(points)
  
  # Extract all raster values with ID to maintain alignment
  fit <- terra::extract(final_raster[["p_fish.fit"]], points_vect)
  se <- terra::extract(final_raster[["p_fish.se.fit"]], points_vect)
  det <- terra::extract(detrend[["detrended_250m"]], points_vect)
  
  # Initialize with NA and assign by ID
  points$fit_values    <- NA
  points$se_values     <- NA
  points$detrend_val   <- NA
  
  points$fit_values[fit$ID]  <- fit[[2]]
  points$se_values[se$ID]    <- se[[2]]
  points$detrend_val[det$ID] <- det[[2]]
  
  points <- st_join(points, strata["strata"], left = TRUE)
  points$design_id <- design_name
  
  extracted_values[[design_name]] <- points
}
summary(extracted_values)

combined_data <- bind_rows(extracted_values) %>% glimpse()
combined_data$SD <- "simple_spabal"

# select the clean columns
combined_data <- combined_data %>% dplyr::select(all_of(final_names))
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_simple_spatially_balanced_", gsub(" ", "_", spp), "_extracted_values_smooth_edge_effect.rds"))
readRDS(paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_simple_spatially_balanced_", gsub(" ", "_", spp), "_extracted_values_smooth_edge_effect.rds"))


#### Stratified Spatial Balance 25-25 -----------------------------------------

str_spabal_sd_2525 <- readRDS(file_str_spabal_sd_2525); str_spabal_sd_2525 <- lapply(str_spabal_sd_2525, function(x) dplyr::select(x, -strata)) # remove numbered strata column
names(str_spabal_sd_2525) <- paste0("stratified_spabal_2525_", 1:1000) # name each list element

# Extract abundance values to each SD point
extracted_values <- list()

for (i in seq_along(str_spabal_sd_2525)) {
  print(i)
  # design_name <- names(str_spabal_sd_2525)[i]
  # points <- str_spabal_sd_2525[[i]]  # keep as sf
  # 
  # # Convert to SpatVector just for extraction
  # points_vect <- vect(points)
  # 
  # # Extract values and take only the second column (value), not the ID
  # fit_values <- terra::extract(final_raster$p_fish.fit, points_vect)[, 2]
  # se_values  <- terra::extract(final_raster$p_fish.se.fit, points_vect)[, 2]
  # detrend_val <- terra::extract(detrend$detrended, points_vect)[, 2] # add detrended bathymetry value
  # #points <- st_join(points, strata[, c("strata")], left = TRUE) # add strata type
  # 
  # # Add extracted values directly to the sf object
  # points$fit_values <- fit_values
  # points$se_values <- se_values
  # points$detrend_val <- detrend_val
  # 
  # points$design_id <- design_name # to identify which simulation
  # 
  # # Store the result with geometry intact
  # extracted_values[[design_name]] <- points
  design_name <- names(str_spabal_sd_2525)[i]
  points <- st_transform(str_spabal_sd_2525[[i]], crs)  # keep as sf
  points_vect <- vect(points)
  
  # Extract all raster values with ID to maintain alignment
  fit <- terra::extract(final_raster[["p_fish.fit"]], points_vect)
  se <- terra::extract(final_raster[["p_fish.se.fit"]], points_vect)
  det <- terra::extract(detrend[["detrended_250m"]], points_vect)
  
  # Initialize with NA and assign by ID
  points$fit_values    <- NA
  points$se_values     <- NA
  points$detrend_val   <- NA
  
  points$fit_values[fit$ID]  <- fit[[2]]
  points$se_values[se$ID]    <- se[[2]]
  points$detrend_val[det$ID] <- det[[2]]
  
  points <- st_join(points, strata["strata"], left = TRUE)
  points$design_id <- design_name
  
  extracted_values[[design_name]] <- points
}
summary(extracted_values)

combined_data <- bind_rows(extracted_values) %>% glimpse()
combined_data$SD <- "stratified_spabal_2525"

# select the clean columns
combined_data <- combined_data %>% dplyr::select(all_of(final_names))
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_stratified_spatially_balanced_2525_", gsub(" ", "_", spp), "_extracted_values_smooth_edge_effect.rds"))


#### Stratified Spatial Balance 20-30 -------------------------------------------------

str_spabal_sd_2030 <- readRDS(file_str_spabal_sd_2030); str_spabal_sd_2030 <- lapply(str_spabal_sd_2030, function(x) dplyr::select(x, -strata)) # remove numbered strata column
names(str_spabal_sd_2030) <- paste0("stratified_spabal_2030_", 1:1000) # name each list element

# Extract abundance values to each SD point
extracted_values <- list()

for (i in seq_along(str_spabal_sd_2030)) {
  print(i)
  # design_name <- names(str_spabal_sd_2030)[i]
  # points <- str_spabal_sd_2030[[i]]  # keep as sf
  # 
  # # Convert to SpatVector just for extraction
  # points_vect <- vect(points)
  # 
  # # Extract values and take only the second column (value), not the ID
  # fit_values <- terra::extract(final_raster$p_fish.fit, points_vect)[, 2]
  # se_values  <- terra::extract(final_raster$p_fish.se.fit, points_vect)[, 2]
  # detrend_val <- terra::extract(detrend$detrended, points_vect)[, 2] # add detrended bathymetry value
  # points <- st_join(points, strata[, c("strata")], left = TRUE) # add strata type
  # 
  # # Add extracted values directly to the sf object
  # points$fit_values <- fit_values
  # points$se_values <- se_values
  # points$detrend_val <- detrend_val
  # 
  # points$design_id <- design_name # to identify which simulation
  # 
  # # Store the result with geometry intact
  # extracted_values[[design_name]] <- points
  design_name <- names(str_spabal_sd_2030)[i]
  points <- st_transform(str_spabal_sd_2030[[i]], crs)  # keep as sf
  points_vect <- vect(points)
  
  # Extract all raster values with ID to maintain alignment
  fit <- terra::extract(final_raster[["p_fish.fit"]], points_vect)
  se <- terra::extract(final_raster[["p_fish.se.fit"]], points_vect)
  det <- terra::extract(detrend[["detrended_250m"]], points_vect)
  
  # Initialize with NA and assign by ID
  points$fit_values    <- NA
  points$se_values     <- NA
  points$detrend_val   <- NA
  
  points$fit_values[fit$ID]  <- fit[[2]]
  points$se_values[se$ID]    <- se[[2]]
  points$detrend_val[det$ID] <- det[[2]]
  
  points <- st_join(points, strata["strata"], left = TRUE)
  points$design_id <- design_name
  
  extracted_values[[design_name]] <- points
}
summary(extracted_values)

combined_data <- bind_rows(extracted_values) %>% glimpse()
combined_data$SD <- "stratified_spabal_2030"

# select the clean columns
combined_data <- combined_data %>% dplyr::select(all_of(final_names))
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_stratified_spatially_balanced_2030_", gsub(" ", "_", spp), "_extracted_values_smooth_edge_effect.rds"))


#### Preferential 25-25 -------------------------------------------------------

pref_sd_2525 <- readRDS(file_pref_sd_2525)
names(pref_sd_2525) <- paste0("pref_2525_", 1:1000) # name each list element

# Extract abundance values to each SD point
extracted_values <- list()

for (i in seq_along(pref_sd_2525)) {
  print(i)
  # design_name <- names(pref_sd_2525)[i]
  # points <- st_transform(pref_sd_2525[[i]], crs)  # keep as sf
  # 
  # # Convert to SpatVector just for extraction
  # points_vect <- vect(points)
  # 
  # # Extract values and take only the second column (value), not the ID
  # fit_values <- terra::extract(final_raster$p_fish.fit, points_vect)[, 2]
  # se_values  <- terra::extract(final_raster$p_fish.se.fit, points_vect)[, 2]
  # detrend_val <- terra::extract(detrend$detrended, points_vect)[, 2] # add detrended bathymetry value
  # points <- st_join(points, strata[, c("strata")], left = TRUE) # add strata type
  # 
  # # Add extracted values directly to the sf object
  # points$fit_values <- fit_values
  # points$se_values <- se_values
  # points$detrend_val <- detrend_val
  # 
  # points$in_SZ <- as.logical(st_within(points$geometry, SZ, sparse = FALSE)[,1])
  # 
  # points$design_id <- design_name # to identify which simulation
  # 
  # # Store the result with geometry intact
  # extracted_values[[design_name]] <- points
  
  design_name <- names(pref_sd_2525)[i]
  points <- st_transform(pref_sd_2525[[i]], crs)  # keep as sf
  points_vect <- vect(points)
  
  # Extract all raster values with ID to maintain alignment
  fit <- terra::extract(final_raster[["p_fish.fit"]], points_vect)
  se <- terra::extract(final_raster[["p_fish.se.fit"]], points_vect)
  det <- terra::extract(detrend[["detrended_250m"]], points_vect)
  
  # Initialize with NA and assign by ID
  points$fit_values    <- NA
  points$se_values     <- NA
  points$detrend_val   <- NA
  
  points$fit_values[fit$ID]  <- fit[[2]]
  points$se_values[se$ID]    <- se[[2]]
  points$detrend_val[det$ID] <- det[[2]]
  
  #points <- st_join(points, strata["strata"], left = TRUE)
  points$design_id <- design_name
  
  extracted_values[[design_name]] <- points
}
summary(extracted_values)

combined_data <- bind_rows(extracted_values) %>% glimpse()
combined_data$SD <- "pref_2525"

# select the clean columns
combined_data <- combined_data %>% dplyr::select(all_of(final_names))
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_preferential_2525_", gsub(" ", "_", spp), "_extracted_values_smooth_edge_effect.rds"))


#### Preferential 20-30 -------------------------------------------------------

pref_sd_2030 <- readRDS(file_pref_sd_2030)
names(pref_sd_2030) <- paste0("pref_2030_", 1:1000) # name each list element

# Extract abundance values to each SD point
extracted_values <- list()

for (i in seq_along(pref_sd_2030)) {
  print(i)
  # Extract all raster values with ID to maintain alignment
  # fit <- terra::extract(final_raster$p_fish.fit, points_vect)
  # se <- terra::extract(final_raster$p_fish.se.fit, points_vect)
  # det <- terra::extract(detrend$detrended, points_vect)
  # 
  # # Initialize with NA and assign by ID
  # points$fit_values    <- NA
  # points$se_values     <- NA
  # points$detrend_val   <- NA
  # 
  # points$fit_values[fit$ID]  <- fit[[2]]
  # points$se_values[se$ID]    <- se[[2]]
  # points$detrend_val[det$ID] <- det[[2]]
  # 
  # points <- st_join(points, strata["strata"], left = TRUE)
  # points$design_id <- design_name
  # 
  
  design_name <- names(pref_sd_2030)[i]
  points <- st_transform(pref_sd_2030[[i]], crs)  # keep as sf
  points_vect <- vect(points)
  
  # Extract all raster values with ID to maintain alignment
  fit <- terra::extract(final_raster[["p_fish.fit"]], points_vect)
  se <- terra::extract(final_raster[["p_fish.se.fit"]], points_vect)
  det <- terra::extract(detrend[["detrended_250m"]], points_vect)
  
  # Initialize with NA and assign by ID
  points$fit_values    <- NA
  points$se_values     <- NA
  points$detrend_val   <- NA
  
  points$fit_values[fit$ID]  <- fit[[2]]
  points$se_values[se$ID]    <- se[[2]]
  points$detrend_val[det$ID] <- det[[2]]
  
  #points <- st_join(points, strata["strata"], left = TRUE)
  points$design_id <- design_name
  
  extracted_values[[design_name]] <- points
}
summary(extracted_values)

combined_data <- bind_rows(extracted_values) %>% glimpse()
combined_data$SD <- "pref_2030"

# select the clean columns
combined_data <- combined_data %>% dplyr::select(all_of(final_names))
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_preferential_2030_", gsub(" ", "_", spp), "_extracted_values_smooth_edge_effect.rds"))


#### Clustered ----------------------------------------------------------------

clump_sd <- readRDS(file_clump_sd)
names(clump_sd) <- paste0("clustered_", 1:1000)

# Extract abundance values to each SD point
extracted_values <- list()

# for (i in seq_along(clump_SD)) { # THIS TAKES A WHILE
#   # Keep track of which point and which design each value is from
#   design_name <- names(clump_SD)[i]
#   points <- clump_SD[[i]]
# 
#   # Extract fitted and SE values at each point (after SZ increase)
#   fit_values_after <- raster::extract(final_raster$p_fish.fit, points)
#   se_values_after <- raster::extract(final_raster$p_fish.se.fit, points)
# 
#   # Select the relevant columns
#   in_SZ <- data.frame(in_SZ = st_drop_geometry(points)$in_SZ)
# 
#   # combine
#   result <- cbind(in_SZ,
#                   cluster_id = st_drop_geometry(points)$cluster_id,
#                   fit_value_after = fit_values_after$p_fish.fit,
#                   se_value_after = se_values_after$p_fish.se.fit)
#   extracted_values[[design_name]] <- result
# }

for (i in seq_along(clump_sd)) {
  print(i)
  # design_name <- names(simp_spabal_sd)[[i]]
  # points <- simp_spabal_sd[[i]]  # Keep as sf
  # 
  # # Convert to SpatVector just for extraction
  # points_vect <- vect(points)
  # 
  # # Extract values and take only the second column (value), not the ID
  # fit_values <- suppressWarnings(terra::extract(final_raster$p_fish.fit, points_vect)[, 2])
  # se_values  <- suppressWarnings(terra::extract(final_raster$p_fish.se.fit, points_vect)[, 2])
  # detrend_val <- suppressWarnings(terra::extract(detrend$detrended, points_vect)[, 2]) # add detrended bathymetry value
  # points <- st_join(points, strata[, c("strata")], left = TRUE) # add strata type
  # 
  # # Add extracted values directly to the sf object
  # points$fit_values <- fit_values
  # points$se_values <- se_values
  # points$detrend_val <- detrend_val
  # 
  # points$design_id <- design_name # to identify which simulation
  # 
  # # Store the result with geometry intact
  # extracted_values[[design_name]] <- points
  # 
  
  design_name <- names(clump_sd)[i]
  points <- st_transform(clump_sd[[i]], crs)  # keep as sf
  points_vect <- vect(points)
  
  # Extract all raster values with ID to maintain alignment
  fit <- terra::extract(final_raster[["p_fish.fit"]], points_vect)
  se <- terra::extract(final_raster[["p_fish.se.fit"]], points_vect)
  det <- terra::extract(detrend[["detrended_250m"]], points_vect)
  
  # Initialize with NA and assign by ID
  points$fit_values    <- NA
  points$se_values     <- NA
  points$detrend_val   <- NA
  
  points$fit_values[fit$ID]  <- fit[[2]]
  points$se_values[se$ID]    <- se[[2]]
  points$detrend_val[det$ID] <- det[[2]]
  
  points <- st_join(points, strata["strata"], left = TRUE)
  points$design_id <- design_name
  
  extracted_values[[design_name]] <- points
}

combined_data <- bind_rows(extracted_values) %>% glimpse()
combined_data$SD <- "clustered"
combined_data <- combined_data %>% rename(in_SZ = zone) %>% mutate(in_SZ = ifelse(in_SZ == "SZ", TRUE, FALSE)) %>% glimpse()

# select the clean columns
combined_data <- combined_data %>% dplyr::select(all_of(final_names), "cluster_id")
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, paste0("outputs/distribution_modelling_outputs/D05_extracted_values/D05_", study_site, "_clustered_", gsub(" ", "_", spp), "_extracted_values_smooth_edge_effect.rds"))

### END ###
