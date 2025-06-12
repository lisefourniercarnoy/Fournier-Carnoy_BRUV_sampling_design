# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare Lengths detected)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Simulate a Sanctuary Zone within which there is a __% increase in population
# Author:  Lise Fournier-Carnoy
# Date:    September 2024

# -----------------------------------------------------------------------------

# Status:  Trying to extract abundance values to all sampling points of both SD

# -----------------------------------------------------------------------------

library(tidyverse) # for manipulating data
library(sf) # for manipulating shapefiles
library(terra) # for manipulating rasters
library(tidyterra) # for manipulating rasters
library(stars)


# Clear memory
rm(list=ls())


## Files used in this script --------------------------------------------------

file_land           <- "data/spatial/shapefiles/wadandi_land.shp"
file_sim_SZ         <- "data/spatial/shapefiles/simulated_SZ.shp"
file_samp_area      <- "data/spatial/shapefiles/L04_sampling_area_deeper_than_7m.rds"
file_detrend_val    <- "data/spatial/rasters/2024_geographe_bathymetry-derivatives.rds"
file_strata         <- "data/rmd/L04_samp_area_detrended.rds"

# all the sampling designs
file_simp_spabal_sd       <- "outputs/Length_comparison/sampling_designs/L04_simple_spatially_balanced_designs.rds"

file_str_spabal_sd_2525   <- "outputs/Length_comparison/sampling_designs/L04_stratified_spatially_balanced_designs_25_25.rds"
file_str_spabal_sd_2030   <- "outputs/Length_comparison/sampling_designs/L04_stratified_spatially_balanced_designs_20_30.rds"

file_pref_sd_2525         <- "outputs/Length_comparison/sampling_designs/L04_preferential_designs_25_25.rds"
file_pref_sd_2030         <- "outputs/Length_comparison/sampling_designs/L04_preferential_designs_20_30.rds"

file_clump_sd             <- "outputs/Length_comparison/sampling_designs/L04_clustered_designs.rds"


dat_extent <-  "full_data"
crs <- 7850
crs_rast <- "epsg:7850"

choose_from <- list.files(
  path = "outputs/Length_comparison/predicted_abundance_rasters/", 
  pattern = "mature_predicted_abundance_", 
  full.names = TRUE
)
choose_from

file_mature_pred <- choose_from[[1]] # choose the fish prediction file

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

spp <- gsub("mature_predicted_abundance_|\\.tif", "", basename(file_mature_pred))
samp_area <- readRDS(file_samp_area)
strata <- readRDS(file_strata) %>% 
  mutate(strata = case_when(
    detrended == '0' ~ "low",
    detrended == '1' ~ "mid",
    detrended == '2' ~ "high",
  )) %>% 
  glimpse()
detrend <- readRDS(file_detrend_val); detrend <- detrend$detrended; plot(detrend)

mature <- rast(file_mature_pred) %>% project(crs_rast); plot(mature)

## Simulate an increase in the zone -------------------------------------------

masked_raster <- mask(mature, sim_sz); plot(masked_raster)

# Increase the values in the masked raster by 80%
increased_raster_fit <- masked_raster$p_mature.fit * 1.8; plot(increased_raster_fit) # Fitted values
increased_raster_se <- masked_raster$p_mature.se.fit * 1.8; plot(increased_raster_se) # Standard errors


# Combine the increased SZ and the rest of the raster
final_raster_fit <- mature$p_mature.fit  # Start with the original raster
final_raster_se <- mature$p_mature.se.fit  # Start with the original raster

# Replace the SZ values of the original raster with the increased values
final_raster_fit[!is.na(increased_raster_fit)] <- increased_raster_fit[!is.na(increased_raster_fit)]
final_raster_se[!is.na(increased_raster_se)] <- increased_raster_se[!is.na(increased_raster_se)]

plot(final_raster_fit)
plot(final_raster_se)


# Full raster with both fitted values and standard errors
final_raster <- c(final_raster_fit, final_raster_se)

saveRDS(final_raster, paste0("outputs/Length_comparison/predicted_abundance_rasters/SZ_increase/L05_mature_predicted_abundance_", spp, "SZ_increase.rds"))

ggplot() +
  labs(title = "Location of the simulated SZ, \nwith 80% increase in abundance within") +
  geom_spatraster(data = final_raster$p_mature.fit) + # Predicted abundance
  geom_sf(data = aus) + # Land
  geom_sf(data = sim_sz, colour = "green", fill = NA, linewidth = 1) + # Simulated SZ
  geom_sf(data = samp_area, colour = "red", fill = NA, linewidth = 1) + # Sample area
  coord_sf(crs = crs_rast, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_fill_gradient(limits = c(0, 5), low = "white", high = "darkblue") + # Adjust raster color scale
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
  # fit_values <- suppressWarnings(terra::extract(final_raster$p_mature.fit, points_vect)[, 2])
  # se_values  <- suppressWarnings(terra::extract(final_raster$p_mature.se.fit, points_vect)[, 2])
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
  fit <- terra::extract(final_raster$p_mature.fit, points_vect)
  se <- terra::extract(final_raster$p_mature.se.fit, points_vect)
  det <- terra::extract(detrend$detrended, points_vect)
  
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
saveRDS(combined_data, paste0("outputs/Length_comparison/sampling_designs/L05_simple_spatially_balanced_", spp, "_extracted_values.rds"))
readRDS(paste0("outputs/Length_comparison/sampling_designs/L05_simple_spatially_balanced_", spp, "_extracted_values.rds"))


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
  # fit_values <- terra::extract(final_raster$p_mature.fit, points_vect)[, 2]
  # se_values  <- terra::extract(final_raster$p_mature.se.fit, points_vect)[, 2]
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
  fit <- terra::extract(final_raster$p_mature.fit, points_vect)
  se <- terra::extract(final_raster$p_mature.se.fit, points_vect)
  det <- terra::extract(detrend$detrended, points_vect)
  
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
saveRDS(combined_data, paste0("outputs/Length_comparison/sampling_designs/L05_stratified_spatially_balanced_2525_", spp, "_extracted_values.rds"))


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
  # fit_values <- terra::extract(final_raster$p_mature.fit, points_vect)[, 2]
  # se_values  <- terra::extract(final_raster$p_mature.se.fit, points_vect)[, 2]
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
  fit <- terra::extract(final_raster$p_mature.fit, points_vect)
  se <- terra::extract(final_raster$p_mature.se.fit, points_vect)
  det <- terra::extract(detrend$detrended, points_vect)
  
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
saveRDS(combined_data, paste0("outputs/Length_comparison/sampling_designs/L05_stratified_spatially_balanced_2030_", spp, "_extracted_values.rds"))


#### Preferential 25-25 -------------------------------------------------------

pref_sd_2525 <- readRDS(file_pref_sd_2525)
names(pref_sd_2525) <- paste0("pref_2525_", 1:1000) # name each list element
try <- pref_sd_2525[[1]]
crs(try)

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
  # fit_values <- terra::extract(final_raster$p_mature.fit, points_vect)[, 2]
  # se_values  <- terra::extract(final_raster$p_mature.se.fit, points_vect)[, 2]
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
  fit <- terra::extract(final_raster$p_mature.fit, points_vect)
  se <- terra::extract(final_raster$p_mature.se.fit, points_vect)
  det <- terra::extract(detrend$detrended, points_vect)
  
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
combined_data$SD <- "pref_2525"

# select the clean columns
combined_data <- combined_data %>% dplyr::select(all_of(final_names))
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, paste0("outputs/Length_comparison/sampling_designs/L05_preferential_2525_", spp, "_extracted_values.rds"))


#### Preferential 20-30 -------------------------------------------------------

pref_sd_2030 <- readRDS(file_pref_sd_2030)
names(pref_sd_2030) <- paste0("pref_2030_", 1:1000) # name each list element

# Extract abundance values to each SD point
extracted_values <- list()

for (i in seq_along(pref_sd_2030)) {
  print(i)
  # Extract all raster values with ID to maintain alignment
  # fit <- terra::extract(final_raster$p_mature.fit, points_vect)
  # se <- terra::extract(final_raster$p_mature.se.fit, points_vect)
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
  fit <- terra::extract(final_raster$p_mature.fit, points_vect)
  se <- terra::extract(final_raster$p_mature.se.fit, points_vect)
  det <- terra::extract(detrend$detrended, points_vect)
  
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
combined_data$SD <- "pref_2030"

# select the clean columns
combined_data <- combined_data %>% dplyr::select(all_of(final_names))
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, paste0("outputs/Length_comparison/sampling_designs/L05_preferential_2030_", spp, "_extracted_values.rds"))


#### Clustered ----------------------------------------------------------------

clump_sd <- readRDS(file_clump_sd)
clump_sd[[1]]
# Extract abundance values to each SD point
extracted_values <- list()

# for (i in seq_along(clump_SD)) { # THIS TAKES A WHILE
#   # Keep track of which point and which design each value is from
#   design_name <- names(clump_SD)[i]
#   points <- clump_SD[[i]]
# 
#   # Extract fitted and SE values at each point (after SZ increase)
#   fit_values_after <- raster::extract(final_raster$p_mature.fit, points)
#   se_values_after <- raster::extract(final_raster$p_mature.se.fit, points)
# 
#   # Select the relevant columns
#   in_SZ <- data.frame(in_SZ = st_drop_geometry(points)$in_SZ)
# 
#   # combine
#   result <- cbind(in_SZ,
#                   cluster_id = st_drop_geometry(points)$cluster_id,
#                   fit_value_after = fit_values_after$p_mature.fit,
#                   se_value_after = se_values_after$p_mature.se.fit)
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
  # fit_values <- suppressWarnings(terra::extract(final_raster$p_mature.fit, points_vect)[, 2])
  # se_values  <- suppressWarnings(terra::extract(final_raster$p_mature.se.fit, points_vect)[, 2])
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
  fit <- terra::extract(final_raster$p_mature.fit, points_vect)
  se <- terra::extract(final_raster$p_mature.se.fit, points_vect)
  det <- terra::extract(detrend$detrended, points_vect)
  
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

# select the clean columns
combined_data <- combined_data %>% dplyr::select(all_of(final_names), "cluster_id")
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, paste0("outputs/Length_comparison/sampling_designs/L05_clustered_", spp, "_extracted_values.rds"))

### END ###
