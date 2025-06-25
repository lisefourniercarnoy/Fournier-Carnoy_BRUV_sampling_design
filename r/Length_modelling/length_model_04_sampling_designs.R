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

# Clear memory
rm(list=ls())

study_site <- "waatu"

min_dis <- 250 # minimum distance between points for all sampling designs

set.seed(12345)

## Load data ------------------------------------------------------------------

# Files used in this script

file_bathy_deriv        <- paste0("data/tidy/L01_", study_site, "_bathy_predictors.rds")

file_land               <- "QGIS layers/clean/wadandi_land_highres.shp"
file_sim_SZ             <- paste0("QGIS layers/clean/", study_site, "_simulated_SZ.shp")
file_pref_samp_area     <- paste0("QGIS layers/clean/", study_site, "_simulated_SZ_sampling_area.shp")


# Load bathymetry and SZs

# Extent to crop layers
ext <- if(study_site == "waatern") {
  ext <- c(115.035, 115.68012207, -33.69743897, -33.35) 
} else { 
  ext <- c(114.7, 115, -34.25, -33.95)
} # Extent to crop layers to
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

depth_mask <- preds$bathy_250m < -7 # Select the depth to cut out
plot(depth_mask)
depth_polygons <- as.polygons(depth_mask, na.rm = TRUE) # Convert to polygon
depth_sf <- st_as_sf(depth_polygons) %>% # Convert again for intersection
  st_transform(crs) %>%
  st_make_valid()

samp_area <- st_intersection(samp_area, depth_sf) %>% # Select only deeper-than-5m
  slice(2); plot(samp_area)

saveRDS(samp_area, paste0("QGIS layers/produced_from_code/L04_", study_site, "_sampling_area_deeper_than_7m.rds"))


## SIMPLE SPATIAL BALANCED ----------------------------------------------------

### Run the sampling design ---------------------------------------------------
samp_area$zone_strata <- 'all' # temporary column
base_samps <- tibble("all" = 50) # for GRTS to know how many samples to put

sample.design <- grts(st_transform(samp_area, 7850),
                      n_base = base_samps,
                      n_over = 10,
                      stratum_var = "zone_strata",
                      DesignID = "LFC-",
                      mindis = min_dis,
                      maxtry = 20)

sample.design <- st_sf(sample.design$sites_base)
sample.design$in_SZ <- as.logical(st_within(sample.design$geometry, SZ, sparse = FALSE)[,1])

tempdat <- st_nn(sample.design$geometry, sample.design$geometry, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = min_dis) %>% # measure only the ones closer than 500m
  glimpse()

samples <- sample.design %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 2)) %>%
  # If everything works well (and the grts mindis call above is set to 500), there should be only zeroes in dists.
  # If not, remove the ones that have dist between 0-500m (0 excluded)
  glimpse()
summary(samples$dists) # all zero. perfect.


# remove points too close to each other (if any)
samples <- samples %>% mutate(row_id = row_number())
nn_list <- tempdat$nn
dist_list <- tempdat$dist
remove_points <- c()

# Check if there are any neighbors beyond self (length > 1 in any nn_list element)
has_neighbors <- any(sapply(nn_list, function(x) length(x) > 1))

if(has_neighbors) {
  for(i in seq_along(nn_list)) {
    neighbours <- nn_list[[i]]
    distances <- dist_list[[i]]
    
    if(length(neighbours) < 2) next # skip if no neighbors except self
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(i, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  if(length(remove_points) > 0) {
    samples_cleaned <- samples[-remove_points, ]
  } else {
    samples_cleaned <- samples
  }
  
} else {
  # No neighbors beyond self, so no removal
  samples_cleaned <- samples
} # this identifies one point from the pair of points within min_dis of each other to remove

if(length(remove_points) > 0) {
  samples_cleaned <- samples[-remove_points, ] # Remove points
  
  # Report
  cat("Removed points:", length(remove_points), "\n")
  cat("Remaining points:", nrow(samples_cleaned), "\n")
} else { 
  # No points to remove, keep original data
  samples_cleaned <- samples
  
  cat("No points removed.\n")
  cat("Total points:", nrow(samples_cleaned), "\n")
} # this removes the points identified above, if there are any.


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
  geom_sf(data = samples_cleaned, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples_cleaned, paste0('data/rmd/L04_', study_site, '_example_sampling_design_simple_spatial_balance.rds')) # for use in Rmarkdown


### Loop for 1000 designs -----------------------------------------------------

num_designs <- 1000 # Number of SB sampling designs to create
sampling_designs <- vector("list", num_designs)
for (i in 1:num_designs) {
  cat("Generating simple spatially balanced design", i, "\n")
  
  # Step 1: Generate GRTS design
  raw_design <- grts(
    st_transform(samp_area, 7850),
    n_base = base_samps,
    n_over = 10,
    stratum_var = "zone_strata",
    DesignID = paste0("LFC-", i),
    mindis = min_dis,
    maxtry = 20
  )
  
  # Step 2: Extract sites_base and convert to sf
  samples <- st_sf(raw_design$sites_base)
  
  # Step 3: Nearest neighbor check
  tempdat <- st_nn(
    samples,
    samples,
    returnDist = TRUE,
    progress = FALSE,
    k = 5,
    maxdist = min_dis
  )
  
  samples <- samples %>% mutate(row_id = row_number())
  nn_list <- tempdat$nn
  dist_list <- tempdat$dist
  
  remove_points <- c()
  has_neighbors <- any(sapply(nn_list, function(x) length(x) > 1))
  
  if (has_neighbors) {
    for (j in seq_along(nn_list)) {
      neighbors <- nn_list[[j]]
      distances <- dist_list[[j]]
      
      if (length(neighbors) < 2) next
      
      close_neighbors <- neighbors[distances > 0 & distances < min_dis]
      for (nbr in close_neighbors) {
        pair <- sort(c(j, nbr))
        remove_points <- c(remove_points, pair[1])
      }
    }
    remove_points <- unique(remove_points)
    if (length(remove_points) > 0) {
      samples <- samples[-remove_points, ]
    }
  }
  
  # Step 4: Add spatial indicator
  samples$in_SZ <- as.logical(st_within(samples$geometry, SZ, sparse = FALSE)[, 1])
  
  # Step 5: Store final design
  sampling_designs[[i]] <- samples
}
plot(sampling_designs[[1]]$geometry)


# Save the list of sampling designs
saveRDS(sampling_designs, file = paste0("outputs/Length_comparison/sampling_designs/L04_", study_area, "simple_spatially_balanced_designs.rds"))



## Make strata ----------------------------------------------------------------

# Using detrended bathymetry - same as 2024 SB design
hist(preds$detrended_250m)
detrended_qs <- c(0, 0.5, 0.9, 1)
detrended_cuts <- global(preds$detrended_250m, probs = detrended_qs, fun = quantile, na.rm = T)
cat_detrended <- classify(preds$detrended_250m, rcl = as.numeric(detrended_cuts[1,]))
par(mfrow = c(1, 1)); plot(cat_detrended)

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
  st_collection_extract("POLYGON") %>% 
  rename(detrended = detrended_250m); plot(sf_all)

saveRDS(sf_all, paste0("data/rmd/L04_", study_site, "_sample_area_detrended_strata.rds"))

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
  dplyr::mutate(area = st_area(.),
                mindis = 500) %>%
  st_transform(crs = crs) %>%
  glimpse()
unique(inp_sf)
plot(inp_sf) # Check that number of samples in each group is all good

saveRDS(inp_sf, "outputs/L04_sampling_design_strata_sample_numbers.rds") # Save to use in other scripts


## SPATIALLY BALANCED SAMPLING DESIGN, 25-25 ----------------------------------

# GRTS needs the number of samples in this horrible wide format for some reason
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
                      DesignID = "LFC-",
                      mindis = min_dis,
                      maxtry = 20)

plot(sample.design)


### Filter out points that are too close to each other ------------------------

tempdat <- st_nn(sample.design$sites_base, sample.design$sites_base, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = min_dis) %>% # measure only the ones closer than min_dis
  glimpse()

samples <- sample.design$sites_base %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 2)) %>%
  # If everything works well (and the grts mindis call above is set above), there should be only zeroes in dists.
  # If not, remove the ones that have dist within the mindis (0 excluded)
  glimpse()
summary(samples$dists) # there may be some pairs of points too close, which we will remove below.

# remove points too close to each other
samples <- samples %>% mutate(row_id = row_number())
nn_list <- tempdat$nn
dist_list <- tempdat$dist
remove_points <- c()

# Check if there are any neighbors beyond self (length > 1 in any nn_list element)
has_neighbors <- any(sapply(nn_list, function(x) length(x) > 1))

if(has_neighbors) {
  for(i in seq_along(nn_list)) {
    neighbours <- nn_list[[i]]
    distances <- dist_list[[i]]
    
    if(length(neighbours) < 2) next # skip if no neighbors except self
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(i, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  if(length(remove_points) > 0) {
    samples_cleaned <- samples[-remove_points, ]
  } else {
    samples_cleaned <- samples
  }
  
} else {
  # No neighbors beyond self, so no removal
  samples_cleaned <- samples
} # this identifies one point from the pair of points within min_dis of each other to remove

if(length(remove_points) > 0) {
  samples_cleaned <- samples[-remove_points, ] # Remove points
  
  # Report
  cat("Removed points:", length(remove_points), "\n")
  cat("Remaining points:", nrow(samples_cleaned), "\n")
} else { 
  # No points to remove, keep original data
  samples_cleaned <- samples
  
  cat("No points removed.\n")
  cat("Total points:", nrow(samples_cleaned), "\n")
} # this removes the points identified above, if there are any.

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
  geom_sf(data = samples_cleaned, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples_cleaned, paste0('data/rmd/L04_', study_site, '_example_sampling_design_stratified_spatial_balance_25_25.rds')) # for use in Rmarkdown


### Loop for 1000 designs -----------------------------------------------------

num_designs <- 1000
sampling_designs <- vector("list", num_designs)

for (i in 1:num_designs) {
  cat("Generating stratified spatially balanced (25/25) design", i, "\n")
  
  # Generate sampling design
  sampling_design <- grts(inp_sf,
                          n_base = base_samps,
                          n_over = 10,
                          stratum_var = "zone_strata",
                          DesignID = paste("LFC-", i, sep = "-"),
                          mindis = min_dis,
                          maxtry = 20)
  
  # Add zone_strata info by spatial join
  sampled_points <- sampling_design$sites_base %>%
    dplyr::mutate(zone_strata = st_join(., inp_sf)$zone_strata) %>%
    dplyr::mutate(row_id = row_number())  # Add row ID for filtering
  
  # Calculate nearest neighbors and distances within min_dis
  tempdat <- st_nn(sampled_points, sampled_points, 
                   returnDist = TRUE, 
                   progress = FALSE, 
                   k = 5, 
                   maxdist = min_dis)
  
  nn_list <- tempdat$nn
  dist_list <- tempdat$dist
  
  remove_points <- integer()  # Initialize empty vector
  
  # Find points to remove (first point in each close pair)
  for(j in seq_along(nn_list)) {
    neighbours <- nn_list[[j]]
    distances <- dist_list[[j]]
    
    if(length(neighbours) < 2) next
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(j, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  # Remove those points from sampled_points
  if(length(remove_points) > 0) {
    samples_cleaned <- sampled_points[-remove_points, ]
    cat("  → Removed", length(remove_points), "points. Remaining:", nrow(samples_cleaned), "\n")
  } else {
    samples_cleaned <- sampled_points
    cat("  → No points removed. All", nrow(samples_cleaned), "retained.\n")
  }
  
  # Store the cleaned sample in the list
  sampling_designs[[i]] <- samples_cleaned
}
summary(sapply(sampling_designs, nrow)) # maximum 3 removed points, should be alright
head(sampling_designs[[1]])


# Save the list of sampling designs
saveRDS(sampling_designs, file = paste0("outputs/Length_comparison/sampling_designs/L04_", study_site, "_stratified_spatial_balanced_designs_25_25.rds"))

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
                      DesignID = "LFC-",
                      mindis = min_dis,
                      maxtry = 20)

plot(sample.design)


### Filter out points that are too close to each other ------------------------

tempdat <- st_nn(sample.design$sites_base, sample.design$sites_base, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = min_dis) %>% # measure only the ones closer than min_dis
  glimpse()

samples <- sample.design$sites_base %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 2)) %>%
  # If everything works well (and the grts mindis call above is set above), there should be only zeroes in dists.
  # If not, remove the ones that have dist within the mindis (0 excluded)
  glimpse()
summary(samples$dists) # there may be some pairs of points too close, which we will remove below.

# remove points too close to each other
samples <- samples %>% mutate(row_id = row_number())
nn_list <- tempdat$nn
dist_list <- tempdat$dist
remove_points <- c()

# Check if there are any neighbors beyond self (length > 1 in any nn_list element)
has_neighbors <- any(sapply(nn_list, function(x) length(x) > 1))

if(has_neighbors) {
  for(i in seq_along(nn_list)) {
    neighbours <- nn_list[[i]]
    distances <- dist_list[[i]]
    
    if(length(neighbours) < 2) next # skip if no neighbors except self
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(i, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  if(length(remove_points) > 0) {
    samples_cleaned <- samples[-remove_points, ]
  } else {
    samples_cleaned <- samples
  }
  
} else {
  # No neighbors beyond self, so no removal
  samples_cleaned <- samples
} # this identifies one point from the pair of points within min_dis of each other to remove

if(length(remove_points) > 0) {
  samples_cleaned <- samples[-remove_points, ] # Remove points
  
  # Report
  cat("Removed points:", length(remove_points), "\n")
  cat("Remaining points:", nrow(samples_cleaned), "\n")
} else { 
  # No points to remove, keep original data
  samples_cleaned <- samples
  
  cat("No points removed.\n")
  cat("Total points:", nrow(samples_cleaned), "\n")
} # this removes the points identified above, if there are any.


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
  geom_sf(data = samples_cleaned, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples, paste0('data/rmd/L04_', study_site, '_example_sampling_design_stratified_spatial_balance_20_30.rds')) # for use in Rmarkdown


### Loop for 1000 designs -----------------------------------------------------

num_designs <- 1000
sampling_designs <- vector("list", num_designs)

for (i in 1:num_designs) {
  cat("Generating stratified spatially balanced (20/30) design", i, "\n")
  
  # Generate sampling design
  sampling_design <- grts(inp_sf,
                          n_base = base_samps,
                          n_over = 10,
                          stratum_var = "zone_strata",
                          DesignID = paste("LFC-", i, sep = "-"),
                          mindis = min_dis,
                          maxtry = 20)
  
  # Add zone_strata info by spatial join
  sampled_points <- sampling_design$sites_base %>%
    dplyr::mutate(zone_strata = st_join(., inp_sf)$zone_strata) %>%
    dplyr::mutate(row_id = row_number())  # Add row ID for filtering
  
  # Calculate nearest neighbors and distances within min_dis
  tempdat <- st_nn(sampled_points, sampled_points, 
                   returnDist = TRUE, 
                   progress = FALSE, 
                   k = 5, 
                   maxdist = min_dis)
  
  nn_list <- tempdat$nn
  dist_list <- tempdat$dist
  
  remove_points <- integer()  # Initialize empty vector
  
  # Find points to remove (first point in each close pair)
  for(j in seq_along(nn_list)) {
    neighbours <- nn_list[[j]]
    distances <- dist_list[[j]]
    
    if(length(neighbours) < 2) next
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(j, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  # Remove those points from sampled_points
  if(length(remove_points) > 0) {
    samples_cleaned <- sampled_points[-remove_points, ]
    cat("  → Removed", length(remove_points), "points. Remaining:", nrow(samples_cleaned), "\n")
  } else {
    samples_cleaned <- sampled_points
    cat("  → No points removed. All", nrow(samples_cleaned), "retained.\n")
  }
  
  # Store the cleaned sample in the list
  sampling_designs[[i]] <- samples_cleaned
}
summary(sapply(sampling_designs, nrow)) # maximum 3 removed points, should be alright
head(sampling_designs[[1]])

# Save the list of sampling designs
saveRDS(sampling_designs, file = paste0("outputs/Length_comparison/sampling_designs/L04_", study_site, "_stratified_spatial_balanced_designs_20_30.rds.rds"))


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
                      mindis = min_dis,
                      maxtry = 20)

plot(sample.design)


### Filter out points that are too close to each other ------------------------


tempdat <- st_nn(sample.design$sites_base, sample.design$sites_base, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = min_dis) %>% # measure only the ones closer than min_dis
  glimpse()

samples <- sample.design$sites_base %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 2)) %>%
  # If everything works well (and the grts mindis call above is set above), there should be only zeroes in dists.
  # If not, remove the ones that have dist within the mindis (0 excluded)
  glimpse()
summary(samples$dists) # there may be some pairs of points too close, which we will remove below.

# remove points too close to each other
samples <- samples %>% mutate(row_id = row_number())
nn_list <- tempdat$nn
dist_list <- tempdat$dist
remove_points <- c()

# Check if there are any neighbors beyond self (length > 1 in any nn_list element)
has_neighbors <- any(sapply(nn_list, function(x) length(x) > 1))

if(has_neighbors) {
  for(i in seq_along(nn_list)) {
    neighbours <- nn_list[[i]]
    distances <- dist_list[[i]]
    
    if(length(neighbours) < 2) next # skip if no neighbors except self
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(i, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  if(length(remove_points) > 0) {
    samples_cleaned <- samples[-remove_points, ]
  } else {
    samples_cleaned <- samples
  }
  
} else {
  # No neighbors beyond self, so no removal
  samples_cleaned <- samples
} # this identifies one point from the pair of points within min_dis of each other to remove

if(length(remove_points) > 0) {
  samples_cleaned <- samples[-remove_points, ] # Remove points
  
  # Report
  cat("Removed points:", length(remove_points), "\n")
  cat("Remaining points:", nrow(samples_cleaned), "\n")
} else { 
  # No points to remove, keep original data
  samples_cleaned <- samples
  
  cat("No points removed.\n")
  cat("Total points:", nrow(samples_cleaned), "\n")
} # this removes the points identified above, if there are any.


### Plot the sampling design --------------------------------------------------

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
  geom_sf(data = samples, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples_sf_pref, paste0('data/rmd/L04_', study_site, '_example_sampling_design_preferential_25_25.rds')) # for use in Rmarkdown



### Loop for 1000 designs ------------------------------------------------------

num_designs <- 1000
sampling_designs <- vector("list", num_designs)

for (i in 1:num_designs) {
  cat("Generating preferential (25/25) design", i, "\n")
  
  # Generate sampling design
  sampling_design <- grts(inp_sf,
                          n_base = base_samps,
                          n_over = 10,
                          stratum_var = "zone_strata",
                          DesignID = paste("LFC-", i, sep = "-"),
                          mindis = min_dis,
                          maxtry = 20)
  
  # Add zone_strata info by spatial join
  sampled_points <- sampling_design$sites_base %>%
    dplyr::mutate(zone_strata = st_join(., inp_sf)$zone_strata) %>%
    dplyr::mutate(row_id = row_number())  # Add row ID for filtering
  
  # Calculate nearest neighbors and distances within min_dis
  tempdat <- st_nn(sampled_points, sampled_points, 
                   returnDist = TRUE, 
                   progress = FALSE, 
                   k = 5, 
                   maxdist = min_dis)
  
  nn_list <- tempdat$nn
  dist_list <- tempdat$dist
  
  remove_points <- integer()  # Initialize empty vector
  
  # Find points to remove (first point in each close pair)
  for(j in seq_along(nn_list)) {
    neighbours <- nn_list[[j]]
    distances <- dist_list[[j]]
    
    if(length(neighbours) < 2) next
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(j, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  # Remove those points from sampled_points
  if(length(remove_points) > 0) {
    samples_cleaned <- sampled_points[-remove_points, ]
    cat("  → Removed", length(remove_points), "points. Remaining:", nrow(samples_cleaned), "\n")
  } else {
    samples_cleaned <- sampled_points
    cat("  → No points removed. All", nrow(samples_cleaned), "retained.\n")
  }
  
  # Store the cleaned sample in the list
  sampling_designs[[i]] <- samples_cleaned
}
summary(sapply(sampling_designs, nrow)) # maximum 3 removed points, should be alright
head(sampling_designs[[1]])

# Save the list
saveRDS(sf_list, file = paste0("outputs/Length_comparison/sampling_designs/L04_", study_site, "_preferential_designs_25_25.rds"))



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
                      mindis = min_dis,
                      maxtry = 20)

plot(sample.design)


### Filter out points that are too close to each other ------------------------


tempdat <- st_nn(sample.design$sites_base, sample.design$sites_base, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = min_dis) %>% # measure only the ones closer than min_dis
  glimpse()

samples <- sample.design$sites_base %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 2)) %>%
  # If everything works well (and the grts mindis call above is set above), there should be only zeroes in dists.
  # If not, remove the ones that have dist within the mindis (0 excluded)
  glimpse()
summary(samples$dists) # there may be some pairs of points too close, which we will remove below.

# remove points too close to each other
samples <- samples %>% mutate(row_id = row_number())
nn_list <- tempdat$nn
dist_list <- tempdat$dist
remove_points <- c()

# Check if there are any neighbors beyond self (length > 1 in any nn_list element)
has_neighbors <- any(sapply(nn_list, function(x) length(x) > 1))

if(has_neighbors) {
  for(i in seq_along(nn_list)) {
    neighbours <- nn_list[[i]]
    distances <- dist_list[[i]]
    
    if(length(neighbours) < 2) next # skip if no neighbors except self
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(i, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  if(length(remove_points) > 0) {
    samples_cleaned <- samples[-remove_points, ]
  } else {
    samples_cleaned <- samples
  }
  
} else {
  # No neighbors beyond self, so no removal
  samples_cleaned <- samples
} # this identifies one point from the pair of points within min_dis of each other to remove

if(length(remove_points) > 0) {
  samples_cleaned <- samples[-remove_points, ] # Remove points
  
  # Report
  cat("Removed points:", length(remove_points), "\n")
  cat("Remaining points:", nrow(samples_cleaned), "\n")
} else { 
  # No points to remove, keep original data
  samples_cleaned <- samples
  
  cat("No points removed.\n")
  cat("Total points:", nrow(samples_cleaned), "\n")
} # this removes the points identified above, if there are any.


### Plot the sampling design --------------------------------------------------

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
  geom_sf(data = samples, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples, paste0('data/rmd/L04_', study_site,'_example_sampling_design_preferential_20_30.rds')) # for use in Rmarkdown



### Loop for 1000 designs ------------------------------------------------------

num_designs <- 1000
sampling_designs <- vector("list", num_designs)

for (i in 1:num_designs) {
  cat("Generating preferential (20/30) design", i, "\n")
  
  # Generate sampling design
  sampling_design <- grts(inp_sf,
                          n_base = base_samps,
                          n_over = 10,
                          stratum_var = "zone_strata",
                          DesignID = paste("LFC-", i, sep = "-"),
                          mindis = min_dis,
                          maxtry = 20)
  
  # Add zone_strata info by spatial join
  sampled_points <- sampling_design$sites_base %>%
    dplyr::mutate(zone_strata = st_join(., inp_sf)$zone_strata) %>%
    dplyr::mutate(row_id = row_number())  # Add row ID for filtering
  
  # Calculate nearest neighbors and distances within min_dis
  tempdat <- st_nn(sampled_points, sampled_points, 
                   returnDist = TRUE, 
                   progress = FALSE, 
                   k = 5, 
                   maxdist = min_dis)
  
  nn_list <- tempdat$nn
  dist_list <- tempdat$dist
  
  remove_points <- integer()  # Initialize empty vector
  
  # Find points to remove (first point in each close pair)
  for(j in seq_along(nn_list)) {
    neighbours <- nn_list[[j]]
    distances <- dist_list[[j]]
    
    if(length(neighbours) < 2) next
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(j, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  # Remove those points from sampled_points
  if(length(remove_points) > 0) {
    samples_cleaned <- sampled_points[-remove_points, ]
    cat("  → Removed", length(remove_points), "points. Remaining:", nrow(samples_cleaned), "\n")
  } else {
    samples_cleaned <- sampled_points
    cat("  → No points removed. All", nrow(samples_cleaned), "retained.\n")
  }
  
  # Store the cleaned sample in the list
  sampling_designs[[i]] <- samples_cleaned
}
summary(sapply(sampling_designs, nrow)) # maximum 3 removed points, should be alright
head(sampling_designs[[1]])

# Save the list
saveRDS(sampling_designs, file = paste0("outputs/Length_comparison/sampling_designs/L04_", study_site, "_preferential_designs_20_30.rds"))


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
cluster_radius          <- 300  # The radius of each clump (in meters)
points_per_cluster_out  <- 8  # Number of points per cluster outside SZ
points_per_cluster_SZ   <- 9  # Number of points per cluster inside SZ (this is to add up to 50 total points)
min_distance            <- min_dis  # Minimum distance between cluster points in meters

# Convert cluster centers to UTM (meters)
cluster_centers_sf <- st_as_sf(cluster_centers, coords = c("X", "Y"), crs = crs)

# Prepare layers for the function
cluster_centers_utm       <- as.data.frame(st_coordinates(cluster_centers_sf))
samp_area_utm             <- st_transform(samp_area, crs)
preferential_area_SZ_utm  <- st_transform(preferential_area_SZ, crs)

generate_clumped_points <- function(cluster_centers_utm, samp_area_utm,
                                    points_per_cluster, cluster_radius,
                                    preferential_area_SZ,
                                    preferential_area_north,
                                    preferential_area_south,
                                    SZ,
                                    min_distance,
                                    points_per_cluster_SZ = 9,
                                    points_per_cluster_out = 8) {
  clumped_points <- list()
  
  for (i in 1:nrow(cluster_centers_utm)) {
    cluster_points <- data.frame(x = numeric(0), y = numeric(0), cluster_id = integer(0), zone = character(0))
    
    # Create sf point for the cluster center
    cluster_sf <- st_sfc(st_point(c(cluster_centers_utm$X[i], cluster_centers_utm$Y[i])), crs = st_crs(preferential_area_SZ))
    
    # Determine which zone the cluster center falls in
    if (lengths(st_intersects(cluster_sf, preferential_area_SZ)) > 0) {
      current_zone <- "SZ"
      zone_polygon <- SZ  # Use strict SZ polygon
      points_needed <- points_per_cluster_SZ
    } else if (lengths(st_intersects(cluster_sf, preferential_area_north)) > 0) {
      current_zone <- "north"
      zone_polygon <- st_difference(samp_area_utm, st_make_valid(SZ))
      points_needed <- points_per_cluster_out
    } else if (lengths(st_intersects(cluster_sf, preferential_area_south)) > 0) {
      current_zone <- "south"
      zone_polygon <- st_difference(samp_area_utm, st_make_valid(SZ))
      points_needed <- points_per_cluster_out
    } else {
      warning(paste("Cluster center", i, "does not fall into any known zone. Skipping."))
      next
    }
    
    # Generate cluster points
    while (nrow(cluster_points) < points_needed) {
      x_random <- rnorm(1, mean = cluster_centers_utm$X[i], sd = cluster_radius)
      y_random <- rnorm(1, mean = cluster_centers_utm$Y[i], sd = cluster_radius)
      point <- st_sfc(st_point(c(x_random, y_random)), crs = st_crs(zone_polygon))
      
      # Check if point is in the correct zone and in sampling area
      if (lengths(st_intersects(point, zone_polygon)) > 0 &&
          lengths(st_intersects(point, samp_area_utm)) > 0) {
        
        # Enforce minimum distance between points in the same cluster
        if (nrow(cluster_points) > 0) {
          existing_points_sf <- st_as_sf(cluster_points, coords = c("x", "y"), crs = st_crs(zone_polygon))
          distances <- st_distance(point, existing_points_sf)
          
          if (any(distances < units::set_units(min_distance, "m"))) next
        }
        
        # Accept point
        cluster_points <- rbind(cluster_points, data.frame(
          x = x_random,
          y = y_random,
          cluster_id = i,
          zone = current_zone
        ))
      }
    }
    
    clumped_points[[i]] <- cluster_points
  }
  
  clumped_points_df <- do.call(rbind, clumped_points)
  return(clumped_points_df)
}

clumped_points_sf <- generate_clumped_points(
  cluster_centers_utm,
  samp_area_utm,
  points_per_cluster = NULL,
  cluster_radius,
  preferential_area_SZ_utm,
  preferential_area_north = st_transform(preferential_area_north, crs),
  preferential_area_south = st_transform(preferential_area_south, crs),
  SZ = SZ,
  min_distance
)

plot(clumped_points_sf$x, clumped_points_sf$y)

# Generate the clumped points within the sampling area
clumped_points_sf <- st_as_sf(clumped_points_sf, coords = c("x", "y"), crs = crs)

# Convert back to spatial features in WGS84 CRS
samples <- st_transform(clumped_points_sf, 4326)

# View results
print(samples)
plot(samples)

### Filter out points that are too close to each other ------------------------

tempdat <- st_nn(samples$geometry, samples$geometry, # Measure the distance between points and their nearest neighbour
                 returnDist = T, progress = F, k = 5, maxdist = min_dis) %>% # measure only the ones closer than min_dis
  glimpse()

samples <- samples %>%
  dplyr::mutate(nn = sapply(tempdat[[1]], "[", 1),
                dists = sapply(tempdat[[2]], "[", 2)) %>%
  # If everything works well (and the grts mindis call above is set above), there should be only zeroes in dists.
  # If not, remove the ones that have dist within the mindis (0 excluded)
  glimpse()
summary(samples$dists) # there may be some pairs of points too close, which we will remove below.

# remove points too close to each other
samples <- samples %>% mutate(row_id = row_number())
nn_list <- tempdat$nn
dist_list <- tempdat$dist
remove_points <- c()

# Check if there are any neighbors beyond self (length > 1 in any nn_list element)
has_neighbors <- any(sapply(nn_list, function(x) length(x) > 1))

if(has_neighbors) {
  for(i in seq_along(nn_list)) {
    neighbours <- nn_list[[i]]
    distances <- dist_list[[i]]
    
    if(length(neighbours) < 2) next # skip if no neighbors except self
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    
    for(nbr in close_neighbours) {
      pair <- sort(c(i, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  if(length(remove_points) > 0) {
    samples_cleaned <- samples[-remove_points, ]
  } else {
    samples_cleaned <- samples
  }
  
} else {
  # No neighbors beyond self, so no removal
  samples_cleaned <- samples
} # this identifies one point from the pair of points within min_dis of each other to remove

if(length(remove_points) > 0) {
  samples_cleaned <- samples[-remove_points, ] # Remove points
  
  # Report
  cat("Removed points:", length(remove_points), "\n")
  cat("Remaining points:", nrow(samples_cleaned), "\n")
} else { 
  # No points to remove, keep original data
  samples_cleaned <- samples
  
  cat("No points removed.\n")
  cat("Total points:", nrow(samples_cleaned), "\n")
} # this removes the points identified above, if there are any.

plot(samples_cleaned)


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
  geom_sf(data = samples_cleaned, aes(colour = "Simulated sample points")) +
  coord_sf(crs = crs, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  scale_color_manual(name = '', values = c('Simulated sample points' = 'red')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

saveRDS(samples_sf_cluster, paste0('data/rmd/L04_', study_site, '_example_sampling_design_clustered.rds')) # for use in Rmarkdown


### Loop for 1000 designs -----------------------------------------------------

num_designs <- 1000
sampling_designs <- vector("list", length = num_designs)
par(mfrow = c(3, 3)) # the loop has a plot check for each simulation, setting up a 3x3 layout

for (i in 1:num_designs) {
  cat("Iteration:", i, "\n")
  
  # Generate cluster centers
  
  preferential_area_SZ    <- clump_area[clump_area$in_SZ == "TRUE", ]
  preferential_area_north <- clump_area[clump_area$in_SZ == "FALSE_north", ]
  preferential_area_south <- clump_area[clump_area$in_SZ == "FALSE_south", ]
  
  num_clusters <- 2
  cluster_centers_true  <- generate_points_in_polygon(preferential_area_SZ, num_clusters)
  cluster_centers_north <- generate_points_in_polygon(preferential_area_north, num_clusters)
  cluster_centers_south <- generate_points_in_polygon(preferential_area_south, num_clusters)
  
  cluster_centers <- st_as_sf(rbind(cluster_centers_true, cluster_centers_north, cluster_centers_south), coords = c("x", "y"), crs = crs)
  cluster_centers_sf <- st_transform(cluster_centers, crs) # for plotting within the loop only
  cluster_centers_utm <- as.data.frame(st_coordinates(cluster_centers)) 
  
  
  # Generate cluster points 
  
  samp_area_utm <- st_transform(samp_area, crs)
  preferential_area_SZ_utm <- st_transform(preferential_area_SZ, crs)
  SZ <- st_transform(SZ, crs)
  
  clumped_points_sf <- generate_clumped_points(
    cluster_centers_utm,
    samp_area_utm,
    points_per_cluster = NULL,
    cluster_radius,
    preferential_area_SZ_utm,
    preferential_area_north = st_transform(preferential_area_north, crs),
    preferential_area_south = st_transform(preferential_area_south, crs),
    SZ,
    min_distance
  )

  clumped_points_sf <- st_as_sf(clumped_points_sf, coords = c("x", "y"), crs = crs)
  samples <- st_transform(clumped_points_sf, 4326)

  
  # Remove points too close to each other

  tempdat <- st_nn(samples$geometry, samples$geometry, returnDist = TRUE, progress = FALSE, k = 5, maxdist = min_dis)
  
  samples <- samples %>%
    mutate(nn = sapply(tempdat[[1]], `[`, 1),
           dists = sapply(tempdat[[2]], `[`, 2),
           row_id = row_number())
  
  remove_points <- c()
  nn_list <- tempdat$nn
  dist_list <- tempdat$dist
  
  for(j in seq_along(nn_list)) {
    neighbours <- nn_list[[j]]
    distances <- dist_list[[j]]
    
    if (length(neighbours) < 2) next
    
    close_neighbours <- neighbours[distances > 0 & distances < min_dis]
    for (nbr in close_neighbours) {
      pair <- sort(c(j, nbr))
      remove_points <- c(remove_points, pair[1])
    }
  }
  
  remove_points <- unique(remove_points)
  
  samples_cleaned <- if (length(remove_points) > 0) {
    samples[-remove_points, ]
  } else {
    samples
  }
  
  cat("Removed:", length(remove_points), "points\n")
  cat("Remaining:", nrow(samples_cleaned), "points\n\n")
  
  
  # Store the result
  
  sampling_designs[[i]] <- samples_cleaned

  
  # Plot check 
  
  cluster_centers_wgs <- st_transform(cluster_centers, crs = 4326)
  SZ_wgs <- st_transform(SZ, crs = 4326)
  
  plot(cluster_centers_wgs$geometry,
       col = "#A3320B",
       pch = 4,
       cex = 1.2,
       main = paste("Simulation number", i))
  plot(samples$geometry,
       col = "#EFA00B",
       pch = 20,
       add = TRUE)
  plot(SZ_wgs$geometry,
       border = "#3E6990",
       lwd = 2,
       add = TRUE)
}

saveRDS(sampling_designs, file = paste0("outputs/Length_comparison/sampling_designs/L04_", study_site, "_clustered_designs.rds"))

### END ###