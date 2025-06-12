# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare Lengths detected)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Simulate sampling designs
# Author:  Lise Fournier-Carnoy, heavily inspired by Claude Spencer
# Date:    September 2024

# -----------------------------------------------------------------------------

# Status:

# -----------------------------------------------------------------------------

library(spsurvey) # for making random sampling designs
library(tidyverse) # for tidy data manipulation and plotting and everything
library(sf) # for dealing with shapefiles
library(terra) # for dealing with shapefiles
library(stars) # for dealing with shapefiles
library(starsExtra) # for dealing with shapefiles
library(tidyterra) # for dealing with shapefiles
library(ggnewscale)
library(ggspatial) # for scalebars on plots
library(nngeo)
library(cubelyr)
library(units) # for cluster design, making sure cluster points aren't too close.
library(spatstat.geom) # for random design with min distance
library(spatstat.random) # for random design with min distance

# Clear memory
rm(list=ls())

set.seed(12345)

## Load data ------------------------------------------------------------------

# Files used in this script

file_bathy_deriv        <- "data/spatial/rasters/2024_geographe_bathymetry-derivatives.rds"

file_land               <- "data/spatial/shapefiles/wadandi_land.shp"
file_sim_SZ             <- "data/spatial/shapefiles/simulated_SZ.shp"
file_pref_samp_area     <- "QGIS layers/polys/simulated_SZ_preferential_sampling_area.shp"

file_abund_crop         <- "outputs/Length_comparison/predicted_abundance_rasters/mature_predicted_abundance.tif"


# Load bathymetry and SZs

# Extent to crop layers
ext <- c(115.035, 115.68012207, -33.69743897, -33.20243897)
crs <- 7850 # projected in m

preds <- readRDS(file_bathy_deriv) %>%
  crop(ext); plot(preds)

# Load shapefiles
aus <- st_read(file_land) %>% # Land
  st_transform(crs) %>%
  glimpse(); plot(aus)

SZ <- st_read(file_sim_SZ) %>% # Simulated SZ
  st_make_valid() %>%
  st_transform(crs) %>%
  st_difference(st_union(aus)) %>%
  glimpse(); plot(SZ$geometry)


# Cut out 7m depth out of the sampling zone
samp_area <- st_read(file_pref_samp_area) %>% # Sampling area
  st_transform(crs) %>%
  st_difference(st_union(aus)) %>%
  st_make_valid()

depth_mask <- preds$gadepth < -7 # Select the depth to cut out
plot(depth_mask)
depth_polygons <- as.polygons(depth_mask, na.rm = TRUE) # Convert to polygon
depth_sf <- st_as_sf(depth_polygons) %>% # Convert again for intersection
  st_transform(crs) %>%
  st_make_valid()

samp_area <- st_intersection(samp_area, depth_sf) %>% # Select only deeper-than-5m
  slice(2); plot(samp_area)

saveRDS(samp_area, "data/spatial/shapefiles/L04_sampling_area_deeper_than_7m.rds")

## SIMPLE SPATIAL BALANCED ----------------------------------------------------

### Run the sampling design ---------------------------------------------------
samp_area$temp <- 'temp' # temporary column
samp_temp <- tibble(temp = 50) # for GRTS to know how many samples to put

sample.design <- grts(st_transform(samp_area, 7850),
                      n_base = samp_temp,
                      n_over = 10,
                      stratum_var = "temp",
                      DesignID = "LFC-GEO",
                      mindis = 500,
                      maxtry = 20)
sample.design <- st_sf(sample.design$sites_base)
sample.design$in_SZ <- as.logical(st_within(sample.design$geometry, SZ, sparse = FALSE)[,1])

tempdat <- st_nn(sample.design$sites_base, sample.design$sites_base, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = 500) %>% # measure only the ones closer than 500m
  glimpse()

samples_sf_simsb <- sample.design$sites_base %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 1)) %>%
  # If everything works well (and the grts mindis call above is set to 500), there should be only zeroes in dists.
  # If not, remove the ones that have dist between 0-500m (0 excluded)
  glimpse()
summary(samples_sf_simsb$dists) # all zero. perfect.


### Plot example design -------------------------------------------------------

ggplot() +
  geom_sf(data = sf_all, aes(fill = as.factor(detrended)), colour = NA) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#E0E0E0", "2" = "#9E9E9E"),
                    na.value = NA,
                    labels = c("0" = "Strata 1 - low",
                               "1" = "Strata 2 - medium",
                               "2" = "Strata 3 - high")) +
  labs(fill = "Detrended bathymetry",
       title = "a. Simple spatially balanced sampling design") +
  new_scale_fill() +
  geom_sf(data = aus) +
  geom_sf(data = SZ, colour = "#7f6000", fill = NA, linewidth = 1) +
  geom_sf(data = samples_sf_simsb, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples_sf_simsb, 'data/rmd/L04_example_sampling_design_simple_spatial_balance.rds') # for use in Rmarkdown


### Loop for 1000 designs -----------------------------------------------------

num_designs <- 1000 # Number of SB sampling designs to create
sampling_designs <- vector("list", num_designs)

for(i in 1:num_designs) {
  # Create a sampling design
  sampling_design <- grts(st_transform(samp_area, 7850),
                        n_base = samp_temp,
                        n_over = 10,
                        stratum_var = "temp",
                        DesignID = "LFC-GEO",
                        mindis = 500,
                        maxtry = 20)
  sampling_design <- st_sf(sampling_design$sites_base)
  sampling_design$in_SZ <- as.logical(st_within(sampling_design$geometry, SZ, sparse = FALSE)[,1])
  
  sampling_designs[[i]] <- sampling_design
}
head(sampling_designs[[1]])

# Save the list of sampling designs
saveRDS(sampling_designs, file = "outputs/Length_comparison/sampling_designs/L04_simple_spatially_balanced_designs.rds")



## Make strata ----------------------------------------------------------------

# Using detrended bathymetry - same as 2024 SB design
hist(preds$detrended)
detrended_qs <- c(0, 0.5, 0.8, 1)
detrended_cuts <- global(preds$detrended, probs = detrended_qs, fun = quantile, na.rm = T)
cat_detrended <- classify(preds$detrended, rcl = as.numeric(detrended_cuts[1,]))
par(mfrow = c(1, 1)); plot(cat_detrended)

# Look at roughness just for shits n gigs
hist(preds$roughness)
roughness_qs <- c(0, 0.5, 0.8, 1)
roughness_cuts <- global(preds$roughness, probs = roughness_qs, fun = quantile, na.rm = T)
cat_roughness <- classify(preds$roughness, rcl = as.numeric(roughness_cuts[1,]))
par(mfrow = c(1, 1)); plot(cat_roughness)


sf_detrended <- as.factor(cat_detrended) %>% # Detrended bathymetry strata
  as.polygons() %>% 
  st_as_sf() %>% 
  st_transform(crs); plot(sf_detrended)

sf_SZ <- st_intersection(sf_detrended, SZ) %>% # Strata within SZ
  mutate(in_SZ = TRUE); plot(sf_SZ)

sf_detrended <- as.factor(cat_detrended) %>% # Strata without SZ
  as.polygons() %>%
  st_as_sf() %>%
  st_transform(crs) %>% 
  st_difference(st_union(SZ)) %>%
  mutate(in_SZ = FALSE); plot(sf_detrended)

sf_all <- bind_rows(sf_SZ, sf_detrended) %>%
  st_make_valid; plot(sf_all) # Strata stitched back together

inp_stars <- st_as_stars(sf_all); plot(inp_stars) # Convert into a 'stars' object, for later

sf_all <- st_intersection(sf_all, samp_area) %>% # Crop to the preferential sampling area
  st_make_valid() %>%
  st_collection_extract("POLYGON"); plot(sf_all)

saveRDS(sf_all, 'data/rmd/L04_samp_area_detrended.rds')

inp_stars <- st_as_stars(sf_all); plot(inp_stars)


# Set the number of samples in each strata
inp_sf <- st_as_sf(inp_stars) %>%
  st_make_valid() %>%
  mutate(strata = as.integer(detrended) + 1) %>%
  mutate( # Group polygons by strata and SZ status
    zone_strata = case_when(
      in_SZ == TRUE & strata == 1 ~ "SZ_strata_1",
      in_SZ == TRUE & strata == 2 ~ "SZ_strata_2",
      in_SZ == TRUE & strata == 3 ~ "SZ_strata_3",
      in_SZ == FALSE & strata == 1 ~ "out_strata_1",
      in_SZ == FALSE & strata == 2 ~ "out_strata_2",
      in_SZ == FALSE & strata == 3 ~ "out_strata_3"
    ),
    nsamps_sb_25_25 = case_when(
      zone_strata == "SZ_strata_1" ~ 10,
      zone_strata == "SZ_strata_2" ~ 10,
      zone_strata == "SZ_strata_3" ~ 5,
      zone_strata == "out_strata_1" ~ 10,
      zone_strata == "out_strata_2" ~ 10,
      zone_strata == "out_strata_3" ~ 5
    ),
    nsamps_sb_20_30 = case_when(
      zone_strata == "SZ_strata_1" ~ 8,
      zone_strata == "SZ_strata_2" ~ 8,
      zone_strata == "SZ_strata_3" ~ 4,
      zone_strata == "out_strata_1" ~ 12,
      zone_strata == "out_strata_2" ~ 12,
      zone_strata == "out_strata_3" ~ 6
    ),
    nsamps_p_25_25 = case_when(
      zone_strata == "SZ_strata_3" ~ 25,
      zone_strata == "out_strata_3" ~ 25
    ),
    nsamps_p_20_30 = case_when(
      zone_strata == "SZ_strata_3" ~ 20,
      zone_strata == "out_strata_3" ~ 30
    )
  ) %>%
  # group_by(zone_strata) %>%
  # summarise(
  #   geometry = st_union(geometry),
  #   nsamps = first(nsamps),  # Retrieve nsamps for each zone_strata
  #   .groups = 'drop'
  # ) %>%
  dplyr::mutate(area = st_area(.),
                mindis = 500) %>%
  st_transform(crs = crs) %>%
  glimpse()
unique(inp_sf)
plot(inp_sf) # Check that number of samples in each group is all good

saveRDS(inp_sf, "data/spatial/shapefiles/L04_sampling_design_zone_strata_all.rds") # Save to use in other scripts


## SPATIALLY BALANCED SAMPLING DESIGN, 25-25 ----------------------------------

# GRTS needs the number of samples in this horrible wide format for some reason
base_samps <- data.frame(nsamps = inp_sf$nsamps_sb_25_25,
                         zone_strata = inp_sf$zone_strata) %>%
  pivot_wider(names_from = zone_strata,
              values_from = nsamps) %>%
  glimpse()
base_samps <- tibble(
  "SZ_strata_1" = 10,
  "SZ_strata_2" = 10,
  "SZ_strata_3" = 5,
  "out_strata_1" = 10,
  "out_strata_2" = 10,
  "out_strata_3" = 5
)

### Run the sampling design ---------------------------------------------------

sample.design <- grts(inp_sf,
                      n_base = base_samps,
                      n_over = 10,
                      stratum_var = "zone_strata",
                      DesignID = "LFC-GEO",
                      mindis = 500,
                      maxtry = 20)

plot(sample.design)


### Filter out points that are too close to each other ------------------------

tempdat <- st_nn(sample.design$sites_base, sample.design$sites_base, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = 500) %>% # measure only the ones closer than 500m
  glimpse()


samples_sf_sb <- sample.design$sites_base %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 1)) %>%
  # If everything works well (and the grts mindis call above is set to 500), there should be only zeroes in dists.
  # If not, remove the ones that have dist between 0-500m (0 excluded)
  glimpse()
summary(samples_sf_sb$dists) # all zero. perfect.


### Plot the sampling design --------------------------------------------------

ggplot() +
  geom_sf(data = sf_all, aes(fill = as.factor(detrended)), colour = NA) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#E0E0E0", "2" = "#9E9E9E"),
                    na.value = NA,
                    labels = c("0" = "Strata 1 - low",
                               "1" = "Strata 2 - medium",
                               "2" = "Strata 3 - high")) +
  labs(fill = "Detrended bathymetry",
       title = "a. Spatially balanced sampling design  25 in SZ, 25 outside the SZ") +
  new_scale_fill() +
  geom_sf(data = aus) +
  geom_sf(data = SZ, colour = "#7f6000", fill = NA, linewidth = 1) +
  geom_sf(data = samples_sf_sb, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples_sf_sb, 'data/rmd/L04_example_sampling_design_stratified_spatial_balance_25_25.rds') # for use in Rmarkdown


### Loop for 1000 designs -----------------------------------------------------

num_designs <- 1000 # Number of SB sampling designs to create
sampling_designs <- vector("list", num_designs)

for (i in 1:num_designs) {
  print(i)
  # Create a sampling design
  sampling_design <- grts(inp_sf,
                          n_base = base_samps,
                          n_over = 10,
                          stratum_var = "zone_strata",
                          DesignID = paste("LFC-GEO", i, sep = "-"),
                          mindis = 100,
                          maxtry = 20)

  # Add strata and zone_strata information to sampled points
  sampled_points <- sampling_design$sites_base %>%
    dplyr::mutate(
      zone_strata = st_join(., inp_sf)$zone_strata,
    )

  sampling_designs[[i]] <- sampled_points
}
head(sampling_designs[[1]])
# Save the list of sampling designs
saveRDS(sampling_designs, file = "outputs/Length_comparison/sampling_designs/L04_stratified_spatial_balanced_designs_25_25.rds")

## SPATIALLY BALANCED SAMPLING DESIGN, 20-30 ----------------------------------

# GRTS needs the number of samples in this horrible wide format for some reason
base_samps <- tibble(
  "SZ_strata_1" = 8,
  "SZ_strata_2" = 8,
  "SZ_strata_3" = 4,
  "out_strata_1" = 12,
  "out_strata_2" = 12,
  "out_strata_3" = 6
)

### Run the sampling design ---------------------------------------------------
sample.design <- grts(inp_sf,
                      n_base = base_samps,
                      n_over = 10,
                      stratum_var = "zone_strata",
                      DesignID = "LFC-GEO",
                      mindis = 500,
                      maxtry = 20)

plot(sample.design)


### Filter out points that are too close to each other ------------------------

tempdat <- st_nn(sample.design$sites_base, sample.design$sites_base, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = 500) %>% # measure only the ones closer than 500m
  glimpse()


samples_sf_sb <- sample.design$sites_base %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 2)) %>%
  # If everything works well (and the grts mindis call above is set to 500), there should be only zeroes in dists.
  # If not, remove the ones that have dist between 0-500m (0 excluded)
  glimpse()
summary(samples_sf_sb$dists) # all zero. perfect.


### Plot the sampling design --------------------------------------------------

ggplot() +
  geom_sf(data = sf_all, aes(fill = as.factor(detrended)), colour = NA) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#E0E0E0", "2" = "#9E9E9E"),
                    na.value = NA,
                    labels = c("0" = "Strata 1 - low",
                               "1" = "Strata 2 - medium",
                               "2" = "Strata 3 - high")) +
  labs(fill = "Detrended bathymetry",
       title = "a. Spatially balanced sampling design") +
  new_scale_fill() +
  geom_sf(data = aus) +
  geom_sf(data = SZ, colour = "#7f6000", fill = NA, linewidth = 1) +
  geom_sf(data = samples_sf_sb, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples_sf_sb, 'data/rmd/L04_example_sampling_design_stratified_spatial_balance_20_30.rds') # for use in Rmarkdown


### Loop for 1000 designs -----------------------------------------------------

num_designs <- 1000 # Number of SB sampling designs to create
sampling_designs <- vector("list", num_designs)

for (i in 1:num_designs) {
  print(i)
  # Create a sampling design
  sampling_design <- grts(inp_sf,
                          n_base = base_samps,
                          n_over = 10,
                          stratum_var = "zone_strata",
                          DesignID = paste("LFC-GEO", i, sep = "-"),
                          mindis = 100,
                          maxtry = 20)
  
  # Add strata and zone_strata information to sampled points
  sampled_points <- sampling_design$sites_base %>%
    dplyr::mutate(
      zone_strata = st_join(., inp_sf)$zone_strata,
    )
  
  sampling_designs[[i]] <- sampled_points
}
head(sampling_designs[[1]])
# Save the list of sampling designs
saveRDS(sampling_designs, file = "outputs/Length_comparison/sampling_designs/L04_stratified_spatial_balanced_designs_20_30.rds.rds")


## PREFERENTIAL SAMPLING DESIGN 25-25 -----------------------------------------

# GRTS needs the number of samples in this horrible wide format for some reason
base_samps <- tibble(
  "SZ_strata_3" = 25,
  "out_strata_3" = 25
)

### Run the sampling design ----

sample.design <- grts(inp_sf,
                      n_base = base_samps,
                      n_over = 10,
                      stratum_var = "zone_strata",
                      DesignID = "LFC-GEO",
                      mindis = 10,
                      maxtry = 20)

plot(sample.design)


### Filter out points that are too close to each other ------------------------

tempdat <- st_nn(sample.design$sites_base, sample.design$sites_base, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = 500) %>% # measure only the ones closer than 500m
  glimpse()


samples_sf_pref <- sample.design$sites_base %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 1)) %>%
  # If everything works well (and the grts mindis call above is set to 500), there should be only zeroes in dists.
  # If not, remove the ones that have dist between 0-500m (0 excluded)
  glimpse()


### Plot the sampling design --------------------------------------------------

tempdat <- st_nn(samples_sf_pref, samples_sf_pref,
                 returnDist = T, progress = F, k = 5, maxdist = 500)


ggplot() +
  geom_sf(data = sf_all, aes(fill = as.factor(detrended)), colour = NA) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#E0E0E0", "2" = "#9E9E9E"),
                    na.value = NA,
                    labels = c("0" = "Strata 1 - low",
                               "1" = "Strata 2 - medium",
                               "2" = "Strata 3 - high")) +
  labs(fill = "Detrended bathymetry",
       title = "b. Preferential sampling design") +
  new_scale_fill() +
  geom_sf(data = aus) +
  geom_sf(data = SZ, colour = "#7f6000", fill = NA, linewidth = 1) +
  geom_sf(data = samples_sf_pref, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples_sf_pref, 'data/rmd/L04_example_sampling_design_preferential_25_25.rds') # for use in Rmarkdown



### Loop for 1000 designs ------------------------------------------------------

num_designs <- 1000 # Number of preferential sampling designs to create
sampling_designs <- vector("list", num_designs) # Create an empty list in which to store designs

for (i in 1:num_designs) { # Create SB sampling designs
  print(i)
  # Create a sampling design
  sampling_design <- grts(inp_sf,
                          n_base = base_samps,
                          n_over = 10,
                          stratum_var = "zone_strata",
                          DesignID = paste("LFC-GEO", i, sep = "-"),
                          mindis = 100,
                          maxtry = 20)
  
  # Add strata and zone_strata information to sampled points
  sampled_points <- sampling_design$sites_base %>%
    dplyr::mutate(
      zone_strata = st_join(., inp_sf)$zone_strata,
    )
  
  sampling_designs[[i]] <- sampled_points
}


# Check the sampling designs, looking to see that they're all different
plot(sampling_designs[[1]])
plot(sampling_designs[[2]])

# Convert the list into usable sf objects
sampling_design_coord <- lapply(sampling_designs, function(x) {
  data.frame(
    lat_WGS84 = x$sites_base$lat_WGS84,
    lon_WGS84 = x$sites_base$lon_WGS84
  )
})
names(sampling_design_coord) <- paste0("pref_2525_design_", seq_along(sampling_design_coord)) # Rename each design

plot(sf_list[[2]])

# Save the list
saveRDS(sf_list, file = "outputs/Length_comparison/sampling_designs/L04_preferential_designs_25_25.rds")



## PREFERENTIAL SAMPLING DESIGN 20-30 -----------------------------------------

# GRTS needs the number of samples in this horrible wide format for some reason
base_samps <- tibble(
  "SZ_strata_3" = 20,
  "out_strata_3" = 30
)


### Run the sampling design ----
sample.design <- grts(inp_sf,
                      n_base = base_samps,
                      n_over = 10,
                      stratum_var = "zone_strata",
                      DesignID = "LFC-GEO",
                      mindis = 10,
                      maxtry = 20)

plot(sample.design)


### Filter out points that are too close to each other ------------------------

tempdat <- st_nn(sample.design$sites_base, sample.design$sites_base, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = 20) %>% # measure only the ones closer than 500m
  glimpse()


samples_sf_pref <- sample.design$sites_base %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 1)) %>%
  # If everything works well (and the grts mindis call above is set to 500), there should be only zeroes in dists.
  # If not, remove the ones that have dist between 0-500m (0 excluded)
  glimpse()


### Plot the sampling design --------------------------------------------------
tempdat <- st_nn(samples_sf_pref, samples_sf_pref,
                 returnDist = T, progress = F, k = 5, maxdist = 500)


ggplot() +
  geom_sf(data = sf_all, aes(fill = as.factor(detrended)), colour = NA) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#E0E0E0", "2" = "#9E9E9E"),
                    na.value = NA,
                    labels = c("0" = "Strata 1 - low",
                               "1" = "Strata 2 - medium",
                               "2" = "Strata 3 - high")) +
  labs(fill = "Detrended bathymetry",
       title = "b. Preferential sampling design") +
  new_scale_fill() +
  geom_sf(data = aus) +
  geom_sf(data = SZ, colour = "#7f6000", fill = NA, linewidth = 1) +
  geom_sf(data = samples_sf_pref, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples_sf_pref, 'data/rmd/L04_example_sampling_design_preferential_20_30.rds') # for use in Rmarkdown


### Loop for 1000 designs ------------------------------------------------------

num_designs <- 1000 # Number of preferential sampling designs to create
sampling_designs <- vector("list", num_designs) # Create an empty list in which to store designs

for (i in 1:num_designs) { # Create SB sampling designs
  print(i)
  # Create a sampling design
  sampling_design <- grts(inp_sf,
                          n_base = base_samps,
                          n_over = 10,
                          stratum_var = "zone_strata",
                          DesignID = paste("LFC-GEO", i, sep = "-"),
                          mindis = 100,
                          maxtry = 20)
  
  # Add strata and zone_strata information to sampled points
  sampled_points <- sampling_design$sites_base %>%
    dplyr::mutate(
      zone_strata = st_join(., inp_sf)$zone_strata,
    )
  
  sampling_designs[[i]] <- sampled_points
}


# Check the sampling designs, looking to see that they're all different
plot(sampling_designs[[1]])
plot(sampling_designs[[2]])

# Convert the list into usable sf objects
sampling_design_coord <- lapply(sampling_designs, function(x) {
  data.frame(
    lat_WGS84 = x$sites_base$lat_WGS84,
    lon_WGS84 = x$sites_base$lon_WGS84
  )
})
names(sampling_design_coord) <- paste0("pref_2030_design_", seq_along(sampling_design_coord)) # Rename each design

plot(sf_list[[2]])

# Save the list
saveRDS(sf_list, file = "outputs/Length_comparison/sampling_designs/L04_preferential_designs_20_30.rds")


## CLUMPED SAMPLING DESIGN ----------------------------------------------------

# Reorganise layers to suit clumped sampling
sf_all_temp <- st_transform(sf_all, 4326)
samp_area_temp <- st_transform(samp_area, 4326)
sf::sf_use_s2(FALSE)

sf_all_clump <- st_intersection(sf_all_temp, samp_area_temp) %>% # Crop to the preferential sampling area
  st_make_valid() %>%
  st_collection_extract("POLYGON") %>%
  dplyr::filter(detrended == 2) %>%  # For clumped design, we will choose reef (high detrended bathy)
  # Now we need to separate the fished zones into South and North (so that we can have 2 clumps either side of the SZ)
  st_cast("POLYGON") %>%                                    # Split the features
  mutate(centroid = st_centroid(geometry),                  # Create centroids for filtering
         latitude = st_coordinates(centroid)[,2])           # Extract latitude

plot(sf_all_clump$geometry)

# Split into fished north, fished south and SZ, then union and put back together
sf_north  <- sf_all_clump %>% filter(in_SZ == FALSE & latitude > -33.51) %>% dplyr::select(-centroid, -latitude) %>% st_union() %>% st_sf(in_SZ = "FALSE_north", geometry = .); plot(sf_north$geometry)
sf_south  <- sf_all_clump %>% filter(in_SZ == FALSE & latitude <= -33.51) %>% dplyr::select(-centroid, -latitude) %>% st_union() %>% st_sf(in_SZ = "FALSE_south", geometry = .); plot(sf_south$geometry)
sf_SZ     <- sf_all_clump %>% filter(in_SZ == TRUE) %>% dplyr::select(-centroid, -latitude) %>% st_union() %>% st_sf(in_SZ = "TRUE", geometry = .); plot(sf_SZ$geometry)

# Combine them back together
sf_clump <- bind_rows(sf_north, sf_south, sf_SZ); plot(sf_clump)


### Create cluster centres ----------------------------------------------------

clump_area <- sf_clump

# Function to generate random points within a polygon
generate_points_in_polygon <- function(polygon, num_points, min_distance = 2000) {
  polygon <- st_transform(polygon, crs = crs)  # Ensure CRS consistency
  points <- data.frame(x = numeric(num_points), y = numeric(num_points))
  cluster_centers <- list()  # To store the cluster centers
  
  # Convert min_distance to a "units" object in meters
  min_distance <- set_units(min_distance, "m")
  
  i <- 1
  while (i <= num_points) {
    # Generate random point
    point <- st_sfc(st_point(c(runif(1, st_bbox(polygon)[1], st_bbox(polygon)[3]),
                               runif(1, st_bbox(polygon)[2], st_bbox(polygon)[4]))),
                    crs = crs)
    
    # Check if the point is inside the polygon
    if (length(st_within(point, polygon)[[1]]) > 0 && st_within(point, polygon)[[1]]) {
      # Check the distance to existing points
      valid_point <- TRUE
      for (existing_point in cluster_centers) {
        distance <- st_distance(point, existing_point)
        
        # Ensure distance comparison works by converting the distance to the same unit (meters)
        if (distance < min_distance) {
          valid_point <- FALSE
          break  # Stop checking if it's too close
        }
      }
      
      # If the point is valid, add it to the cluster centers
      if (valid_point) {
        points[i, ] <- c(st_coordinates(point)[1], st_coordinates(point)[2])
        cluster_centers[[i]] <- point  # Store the valid point as a cluster center
        i <- i + 1
      }
    }
  }
  return(points)
}

## Create cluster centres, around which points will be.
preferential_area_SZ    <- clump_area[clump_area$in_SZ == "TRUE", ]
preferential_area_north <- clump_area[clump_area$in_SZ == "FALSE_north", ]
preferential_area_south <- clump_area[clump_area$in_SZ == "FALSE_south", ]


# Generate 2 clusters within the SZ area, and 4 outside the SZ
num_clusters <- 2
num_clusters_tot <- num_clusters * 3 # x in the SZ, x south, x north = x*3

cluster_centers_true <- generate_points_in_polygon(preferential_area_SZ, num_clusters)
cluster_centers_north <- generate_points_in_polygon(preferential_area_north, num_clusters)
cluster_centers_south <- generate_points_in_polygon(preferential_area_south, num_clusters)

cluster_centers <- st_as_sf(rbind(cluster_centers_true, cluster_centers_north, cluster_centers_south), coords = c("x", "y"), crs = st_crs(clump_area)); plot(cluster_centers)


### Create points around the cluster centres ----------------------------------

# Parameters for the clusters
cluster_radius          <- 400  # The radius of each clump (in meters)
points_per_cluster_out  <- 8  # Number of points per cluster outside SZ
points_per_cluster_SZ   <- 9  # Number of points per cluster inside SZ (this is to add up to 50 total points)
min_distance            <- 200  # Minimum distance between cluster points in meters

# Convert cluster centers to UTM (meters)
cluster_centers_sf <- st_as_sf(cluster_centers, coords = c("X", "Y"), crs = crs)

# Prepare layers for the function
cluster_centers_utm       <- as.data.frame(st_coordinates(cluster_centers_sf))
samp_area_utm             <- st_transform(samp_area, crs)
preferential_area_SZ_utm  <- st_transform(preferential_area_SZ, crs)

generate_clumped_points <- function(cluster_centers_utm, samp_area_utm, points_per_cluster, cluster_radius, preferential_area_SZ, min_distance) {
  clumped_points <- list()
  
  # Loop through each cluster center
  for (i in 1:nrow(cluster_centers_utm)) {
    cluster_points <- data.frame(
      x = numeric(0),  
      y = numeric(0),  
      cluster_id = integer(0),
      in_SZ = logical(0)  # New column to indicate if the point is in SZ
    )
    
    # Convert cluster center to sf point
    cluster_sf <- st_sfc(st_point(c(cluster_centers_utm$X[i], cluster_centers_utm$Y[i])), crs = st_crs(preferential_area_SZ))
    
    # Check if cluster center is inside SZ
    is_in_SZ <- lengths(st_intersects(cluster_sf, preferential_area_SZ)) > 0
    points_per_cluster <- ifelse(is_in_SZ, points_per_cluster_SZ, points_per_cluster_out)
    
    # Generate points while maintaining minimum distance constraint
    while (nrow(cluster_points) < points_per_cluster) {
      # Generate random point within cluster radius
      x_random <- rnorm(1, mean = cluster_centers_utm$X[i], sd = cluster_radius)
      y_random <- rnorm(1, mean = cluster_centers_utm$Y[i], sd = cluster_radius)
      
      # Create sf point
      point <- st_sfc(st_point(c(x_random, y_random)), crs = st_crs(preferential_area_SZ))
      
      # Ensure the point is inside the sampling area
      if (lengths(st_intersects(point, samp_area_utm)) > 0) {
        
        # Check if the point is inside the SZ
        point_in_SZ <- lengths(st_intersects(point, preferential_area_SZ)) > 0
        
        # Ensure the point is inside or outside SZ based on cluster category
        if (is_in_SZ) {
          if (!point_in_SZ) next  # Skip if not inside SZ
        } else {
          if (point_in_SZ) next  # Skip if inside SZ
        }
        
        # Convert existing points in the cluster to sf object
        if (nrow(cluster_points) > 0) {
          existing_points_sf <- st_as_sf(cluster_points, coords = c("x", "y"), crs = st_crs(preferential_area_SZ))
          
          # Compute distances to existing points
          distances <- st_distance(point, existing_points_sf)
          
          # Ensure new point is at least min_distance from all existing points
          if (any(distances < units::set_units(min_distance, "m"))) next # Skip this point if too close
        }
        
        # Add point to cluster
        cluster_points <- rbind(cluster_points, data.frame(x = x_random, y = y_random, cluster_id = i, in_SZ = point_in_SZ))
      }
    }
    
    # Store generated cluster points
    clumped_points[[i]] <- cluster_points
  }
  
  # Combine all clusters into a single data frame
  clumped_points_df <- do.call(rbind, clumped_points)
  return(clumped_points_df)
}


# Generate the clumped points within the sampling area
clumped_points_sf <- generate_clumped_points(cluster_centers_utm, samp_area_utm, points_per_cluster, cluster_radius, preferential_area_SZ_utm, min_distance)
clumped_points_sf <- st_as_sf(clumped_points_sf, coords = c("x", "y"), crs = crs)

# Convert back to spatial features in WGS84 CRS
samples_sf_cluster <- st_transform(clumped_points_sf, 4326)

# View results
print(samples_sf_cluster)

### Plot the sampling design --------------------------------------------------

ggplot() +
  geom_sf(data = sf_all, aes(fill = as.factor(detrended)), colour = NA) +
  scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#E0E0E0", "2" = "#9E9E9E"),
                    na.value = NA,
                    labels = c("0" = "Strata 1 - low",
                               "1" = "Strata 2 - medium",
                               "2" = "Strata 3 - high")) +
  labs(fill = "Detrended bathymetry",
       title = "Clustered sampling design") +
  new_scale_fill() +
  geom_sf(data = aus) +
  geom_sf(data = SZ, colour = "#7f6000", fill = NA, linewidth = 1) +
  geom_sf(data = samples_sf_cluster, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
saveRDS(samples_sf_cluster, 'data/rmd/L04_example_sampling_design_clustered.rds') # for use in Rmarkdown


### Loop for 1000 designs -----------------------------------------------------

# Number of designs you want to generate
num_designs <- 1000
sampling_designs <- vector("list", num_designs)  # Create an empty list to store the designs

# Loop to generate 1000 designs
for (i in 1:num_designs) {
  # Set a random seed for each iteration to ensure different cluster centers and points
  set.seed(i)  # This ensures a different seed for each design, which leads to different random points
  print(i)
  
  # Generate cluster centers for each design
  cluster_centers_true <- generate_points_in_polygon(preferential_area_SZ, num_clusters)
  cluster_centers_north <- generate_points_in_polygon(preferential_area_north, num_clusters)
  cluster_centers_south <- generate_points_in_polygon(preferential_area_south, num_clusters)

  # Combine cluster centers from all areas
  cluster_centers <- st_as_sf(rbind(cluster_centers_true, cluster_centers_north, cluster_centers_south),
                              coords = c("x", "y"), crs = st_crs(clump_area))

  # Generate clumped points around these cluster centers
  clumped_points_sf <- generate_clumped_points(cluster_centers_utm, samp_area_utm, points_per_cluster, cluster_radius, preferential_area_SZ_utm, min_distance)
  clumped_points_sf <- st_as_sf(clumped_points_sf, coords = c("x", "y"), crs = crs)

  # Convert clumped points back to WGS84 CRS
  clumped_points_sf_wgs84 <- st_transform(clumped_points_sf, 4326)

  # Save the current design in the list
  sampling_designs[[i]] <- clumped_points_sf_wgs84
}

# Check the sampling designs
plot(sampling_designs[[1]])
plot(sampling_designs[[2]])


names(sampling_designs) <- paste0("cluster_design_", seq_along(sampling_designs))

# If they're already sf objects, no need to convert. But ensure they have the correct CRS.
sf_list <- lapply(sampling_designs, function(x) {
  st_set_crs(x, 4326)  # Set to WGS 84 if not already set
})

# Now you can access and plot directly with attributes preserved
plot(sf_list[[5]])
plot(sf_list[[15]])
plot(sf_list[[2]])

# Save the list of designs
saveRDS(sf_list, file = "outputs/Length_comparison/sampling_designs/L04_clustered_designs.rds")

### END ###