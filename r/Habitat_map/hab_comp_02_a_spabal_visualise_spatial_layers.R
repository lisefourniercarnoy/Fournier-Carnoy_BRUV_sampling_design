# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (C. Habitat map comp.)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Create bathymetry derivatives and finalise habitat dataframe
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  Done for now. To be changed with 2024 data once it's processed.

# -----------------------------------------------------------------------------


# Clear memory
rm(list=ls())

# Load libraries
library(sf)
library(terra)
library(stars)
library(starsExtra)
library(tidyverse)
library(tidyterra)
library(patchwork)

# Set the study name.
name <- '2024_geographe'
sampling_design <- "_spabal"


#### Bathymetry ----

# Load the bathymetry data (GA 250m resolution).
e <- ext(114.9, 115.7, -33.7, -33.2)

bathy <- rast("data/spatial/rasters/GB-SW_250mBathy.tif") %>%
  crop(e) %>%
  clamp(upper = 0, lower = -250, values = F) %>%
  trim()
plot(bathy)


#### Bathymetry derivatives ----

# Create terrain metrics (bathymetry derivatives).
preds <- terrain(bathy, neighbors = 8,
                 v = c("aspect", "roughness"),
                 unit = "degrees")

# Create detrended bathymetry.
zstar <- st_as_stars(bathy)
detre <- detrend(zstar, parallel = 8)
detre <- as(object = detre, Class = "SpatRaster")
names(detre) <- c("detrended", "lineartrend")

# Join depth, terrain metrics and detrended bathymetry.
preds <- rast(list(bathy, preds, detre[[1]]))
names(preds)[1] <- "gadepth"

# Save the bathymetry derivatives.
saveRDS(preds,
        file = paste0("data/spatial/rasters/",
                      name, sampling_design, "_bathymetry-derivatives.rds"))


#### Sample points ----

# Read in the metadata.
metadata <- read.csv(paste0("data/tidy/", name, "_spabal_tidy_metadata.csv")) %>%
  dplyr::mutate(longitude = as.numeric(longitude),
                latitude = as.numeric(latitude)) %>%
  glimpse()

# Convert metadata to a spatial file and check alignment with bathymetry.
metadata.vect <- vect(metadata, geom = c("longitude", "latitude"), crs = "epsg:4326")

plot(preds[[1]])
plot(metadata.vect, add = T)

# Extract bathymetry derivatives at each of the samples.
metadata.bathy.derivatives   <- cbind(metadata,
                                      terra::extract(preds, metadata.vect)) %>%
  filter_at(vars(gadepth, aspect, roughness, detrended), all_vars(!is.na(.))) %>% # Removes samples missing bathymetry derivatives
  glimpse()

# Save the metadata bathymetry derivatives.
saveRDS(metadata.bathy.derivatives, paste0("data/tidy/", name, sampling_design, "_metadata-bathymetry-derivatives.rds"))

aus <- st_read("data/spatial/shapefiles/WAcoastline.shp") #data/spatial/shp/cstauscd_r.mif")                            # geodata 100k coastline available: https://data.gov.au/dataset/ds-ga-a05f7892-eae3-7506-e044-00144fdd4fa6/

# Load marine park boundaries, state and commonwealth.
marine.parks <- st_read("data/spatial/shapefiles/CAPAD_SW.shp") %>%
  dplyr::filter(NAME %in% c("South-west Corner", "Eastern Recherche")) %>%
  dplyr::mutate(ZONE_TYPE = str_replace_all(ZONE_TYPE, " \\s*\\([^\\)]+\\)", "")) %>%
  glimpse()

# Plot bathymetry derivatives.
plot(preds)

depth <- ggplot() +
  geom_spatraster(data = preds, aes(fill = gadepth)) +
  scale_fill_viridis_c(option = "A", na.value = "transparent",
                       name = "Depth") +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  geom_sf(data = marine.parks %>% dplyr::filter(RES_NUMBER %in% c("swswcnpz15", "swearnpz02")),
          fill = NA, colour = "#7bbc63", linewidth = 0.6) +
  theme_minimal() +
  coord_sf(crs = 4326, xlim = c(ext(preds)[1], ext(preds)[2]),
           ylim = c(ext(preds)[4], ext(preds)[3]))

roughness <- ggplot() +
  geom_spatraster(data = preds, aes(fill = roughness)) +
  scale_fill_viridis_c(option = "F", na.value = "transparent",
                       name = "Roughness", direction = -1) +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  geom_sf(data = marine.parks %>% dplyr::filter(RES_NUMBER %in% c("swswcnpz15", "swearnpz02")),
          fill = NA, colour = "#7bbc63", linewidth = 0.6) +
  theme_minimal() +
  coord_sf(crs = 4326, xlim = c(ext(preds)[1], ext(preds)[2]),
           ylim = c(ext(preds)[4], ext(preds)[3]))

aspect <- ggplot() +
  geom_spatraster(data = preds, aes(fill = aspect)) +
  scale_fill_viridis_c(option = "C", na.value = "transparent",
                       name = "Aspect") +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  geom_sf(data = marine.parks %>% dplyr::filter(RES_NUMBER %in% c("swswcnpz15", "swearnpz02")),
          fill = NA, colour = "#7bbc63", linewidth = 0.6) +
  theme_minimal() +
  coord_sf(crs = 4326, xlim = c(ext(preds)[1], ext(preds)[2]),
           ylim = c(ext(preds)[4], ext(preds)[3]))

detrended <- ggplot() +
  geom_spatraster(data = preds, aes(fill = detrended)) +
  scale_fill_viridis_c(option = "D", na.value = "transparent",
                       name = "Detrended") +
  geom_sf(data = aus, fill = "seashell2", colour = "grey80", size = 0.1) +
  geom_sf(data = marine.parks %>% dplyr::filter(RES_NUMBER %in% c("swswcnpz15", "swearnpz02")),
          fill = NA, colour = "#7bbc63", linewidth = 0.6) +
  theme_minimal() +
  coord_sf(crs = 4326, xlim = c(ext(preds)[1], ext(preds)[2]),
           ylim = c(ext(preds)[4], ext(preds)[3]))

depth / roughness / aspect / detrended + plot_layout(guides = "collect")
dev.off()


## Habitat data ----

# Load and format the tidy habitat data into a format suitable for modelling.
# For this we need predictors for every datapoint, and habitat classes in tidy "broad" classes
final.habitat <- readRDS(paste0("data/tidy/", name, sampling_design, "_tidy_habitat.rds")) %>%
  dplyr::mutate(longitude = as.double(longitude),
                latitude = as.double(latitude)) %>%
  left_join(metadata.bathy.derivatives) %>%
  dplyr::filter(!is.na(gadepth)) %>%
  dplyr::select(-ID) %>%
  glimpse()

# Save the output for use in the modelling scripts
saveRDS(final.habitat, file = paste0("data/tidy/", name, sampling_design, "_tidy_habitat.rds"))


