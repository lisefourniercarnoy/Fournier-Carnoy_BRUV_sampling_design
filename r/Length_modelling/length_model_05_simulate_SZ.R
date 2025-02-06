# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare Lengths detected)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Simulate a Sanctuary Zone within which there is a __% increase in population
# Author:  Lise Fournier-Carnoy
# Date:    September 2024

# -----------------------------------------------------------------------------

# Status:  Trying to extract abundance values to all sampling points of both SD

# -----------------------------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(stars)


# Clear memory
rm(list=ls())


## Files used in this script --------------------------------------------------

file_land <- "data/spatial/shapefiles/wadandi_land.shp"
file_sim_SZ <-"data/spatial/shapefiles/simulated_SZ.shp"
file_samp_area    <- "data/spatial/shapefiles/sampling_area_deeper_than_7m.rds"

file_spabal_SD <- "outputs/Length_comparison/sampling_designs/spatially_balanced_designs.rds"
file_pref_SD <- "outputs/Length_comparison/sampling_designs/preferential_designs.rds"

file_mature_pred <- "outputs/Length_comparison/predicted_abundance_rasters/mature_predicted_abundance.tif"


## Load files -----------------------------------------------------------------

ext <- c(115.035, 115.68012207, -33.69743897, -33.20243897)


aus <- st_read(file_land) %>% # Land
  st_transform(4326) %>%
  glimpse()

sim_sz <- st_read(file_sim_SZ) %>% # Simulated SZ
  st_make_valid() %>%
  st_transform(4326) %>%
  st_difference(st_union(aus)) %>%
  glimpse()

mature <- rast(file_mature_pred) %>%
  glimpse() %>%
  crop(ext) %>%
  project(st_crs(sim_sz)$proj4string)

samp_area <- readRDS(file_samp_area)


ggplot() +
  labs(title = "Location of the simulated SZ") +
  geom_spatraster(data = mature, aes(fill = p_mature.fit)) + # Predicted abundance
  geom_sf(data = aus) + # Land
  geom_sf(data = sim_sz, colour = "green", fill = NA, linewidth = 1) + # Simulated SZ
  coord_sf(crs = 4326, xlim = c(115.0, 115.65), ylim = c(-33.2, -33.7)) +
  theme_minimal()



## Mature fish ----------------------------------------------------------------

### Simulate an increase in the zone ------------------------------------------

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

saveRDS(final_raster, "outputs/Length_comparison/predicted_abundance_rasters/mature_predicted_abundance_SZ_increase.rds")

ggplot() +
  labs(title = "Location of the simulated SZ, \nwith 80% increase in abundance within") +
  geom_spatraster(data = final_raster$p_mature.fit) + # Predicted abundance
  geom_sf(data = aus) + # Land
  geom_sf(data = sim_sz, colour = "green", fill = NA, linewidth = 1) + # Simulated SZ
  geom_sf(data = samp_area, colour = "red", fill = NA, linewidth = 1) + # Sample area
  coord_sf(crs = 4326, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_fill_gradient(limits = c(0, 40), low = "lightblue", high = "darkblue") + # Adjust raster color scale
  theme_minimal()


### Extract abundance values from sampling design locations -------------------

#### Spatially balanced -------------------------------------------------------

spabal_SD <- readRDS(file_spabal_SD)

spabal_SD <- lapply(spabal_SD, function(sf_design) { # Categorise each point by whether it's in a zone or not
    sf_design$in_SZ <- st_within(sf_design, sim_sz, sparse = FALSE) %>%
      rowSums() > 0  # TRUE if point is in the SZ, otherwise FALSE
    return(sf_design)
  })
summary(spabal_SD$spabal_design_1$in_SZ)
plot(spabal_SD$spabal_design_1)

# Extract abundance values to each SD point
extracted_values <- list()

for (i in seq_along(spabal_SD)) { # THIS TAKES A WHILE
  # Keep track of which point and which design each value is from
  design_name <- names(spabal_SD)[i]
  points <- spabal_SD[[i]]

  # Extract fitted and SE values at each point (before SZ increase)
  fit_values_before <- raster::extract(mature$p_mature.fit, points)
  se_values_before <- raster::extract(mature$p_mature.se.fit, points)

  # Extract fitted and SE values at each point (after SZ increase)
  fit_values_after <- raster::extract(final_raster$p_mature.fit, points)
  se_values_after <- raster::extract(final_raster$p_mature.se.fit, points)

  # Add values in respective columns
  result <- cbind(st_drop_geometry(points), fit_values_before, se_values_before$p_mature.se.fit, fit_values_after, se_values_after$p_mature.se.fit)
  colnames(result) <- c("in_SZ", "ID_1", "fit_value_before", "se_value_before", "ID_2", "fit_value_after", "se_value_after")
  extracted_values[[design_name]] <- result
}

glimpse(extracted_values)

# Combine data into one

combined_data <- bind_rows(
  lapply(seq_along(extracted_values), function(i) {
    df <- extracted_values[[i]]
    df$design_id <- paste("Design", i)  # Add a column to identify design
    return(df)
  })
)
combined_data$SD <- "spabal"

glimpse(combined_data)

# Save for next script
saveRDS(combined_data, "outputs/Length_comparison/sampling_designs/spatially_balanced_extracted_values.rds")


#### Preferential -------------------------------------------------------------

pref_SD <- readRDS(file_pref_SD)

pref_SD <- lapply(pref_SD, function(sf_design) { # Categorise each point by whether it's in a zone or not
  sf_design$in_SZ <- st_within(sf_design, sim_sz, sparse = FALSE) %>%
    rowSums() > 0  # TRUE if point is in the SZ, otherwise FALSE
  return(sf_design)
})
summary(pref_SD$spabal_design_1$in_SZ)

# Extract abundance values to each SD point
extracted_values <- list()

for (i in seq_along(pref_SD)) {
  # Keep track of which point and which design each value is from
  design_name <- names(pref_SD)[i]
  points <- pref_SD[[i]]

  # Extract fitted and SE values at each point (before SZ increase)
  fit_values_before <- terra::extract(mature$p_mature.fit, points)
  se_values_before <- terra::extract(mature$p_mature.se.fit, points)

  # Extract fitted and SE values at each point (after SZ increase)
  fit_values_after <- terra::extract(final_raster$p_mature.fit, points)
  se_values_after <- terra::extract(final_raster$p_mature.se.fit, points)

  # Add values in respective columns
  result <- cbind(st_drop_geometry(points), fit_values_before, se_values_before$p_mature.se.fit, fit_values_after, se_values_after$p_mature.se.fit)
  colnames(result) <- c("in_SZ", "ID_1", "fit_value_before", "se_value_before", "ID_2", "fit_value_after", "se_value_after")
  extracted_values[[design_name]] <- result
}

glimpse(extracted_values)


# Combine data into one
combined_data <- bind_rows(
  lapply(seq_along(extracted_values), function(i) {
    df <- extracted_values[[i]]
    df$design_id <- paste("Design", i)  # Add a column to identify design
    return(df)
  })
)
combined_data$SD <- "pref"
glimpse(combined_data)

# Save for next script
saveRDS(combined_data, "outputs/Length_comparison/sampling_designs/preferential_extracted_values.rds")


### END ###
