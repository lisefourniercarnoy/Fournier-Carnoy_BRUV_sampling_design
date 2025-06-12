# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. length model)
# Data:    2024 MEGlab Geographe Bay BRUV habitat & abundance data.
# Task:    Make tidy files for future scripts
# Author:  Lise Fournier-Carnoy
# Date:    February 2025

# -----------------------------------------------------------------------------

# Status:  Final data, all clean.

# -----------------------------------------------------------------------------

library(dismo)
library(tidyverse)
library(sf)
library(terra)
library(raster)
library(CheckEM)

# Clear memory
rm(list=ls())


ext <- c(115.035, 115.68012207, -33.69743897, -33.35) # Extent to crop layers to

## Files used in this script --------------------------------------------------

file_fish_length        <- "data/tidy/2024-04_Geographe_stereo-BRUVs_complete-length.rds" # obtained from 00_checkem
file_metadata           <- "data/raw/2024_commonwealth/2024-04_Geographe_stereo-BRUVs_Metadata.rds" # obtained from MEGlab labsheets

file_obs_hab            <- "data/raw/2024_commonwealth/2024-04_Geographe_stereo-BRUVs_forwards_Dot Point Measurements.txt" # output form TransectMeasure

file_250m_bathy         <- "data/spatial/rasters/GB-SW_250mBathy.tif" # from Geoscience Australia
file_250m_bathy_deriv   <- "data/spatial/rasters/2014_geographe_bathymetry-derivatives.RDS" # From Claude's code
file_lidar_bathy        <- "data/spatial/rasters/LiDAR_geo_compressed.tif" # Too big for git


## Presence/Absence data ------------------------------------------------------

L50 <- c('Chrysophrys auratus' = 0, #566, # L50 from Wakefield et al. 2015 - males 585, females 566
         #'Pseudocaranx spp' = 279, # L50 from male P. dentex, Farmer et al. 2005
         'Glaucosoma hebraicum' = 0, #301, # L50 from Hesp 2002 - males 320, females 301 - but using all fish because papers find NTZ makes dhuies more abundant not bigger
         #'Choerodon rubescens' = 264, # L50 from females, Nardi 2006
         #'Nemadactylus valenciennesi' = 400,# L50 from females, Coulson et al. 2010
         #"Epinephelides armatus" = 285, # L50 from females west coast, Buckland 2022
         "Ophthalmolepis lineolatus" = 0 #184 # l50 from females, Morton 2008
         ) 

metadata <- readRDS(file_metadata) %>% 
  rename(opcode = sample) %>% 
  dplyr::filter(successful_count == "Yes" & successful_length == "Yes") %>% 
  dplyr::select("opcode", "latitude_dd", "longitude_dd", "status", "depth_m") %>% 
  glimpse()

length_data <- readRDS(file_fish_length) %>%
  rename(opcode = sample) %>% 
  mutate(full_spp = paste(genus, species)) %>% 
  glimpse()

dat <- length_data %>% 
  # Select species present in L50
  filter(full_spp %in% names(L50)) %>% 
  
  # Join the metadata with the length data
  left_join(metadata %>% distinct(opcode, longitude_dd, latitude_dd), by = "opcode") %>% 
  
  # Separate into size categories, treating missing lengths as immature
  group_by(longitude_dd, latitude_dd, full_spp, opcode) %>% 
  summarise(
    count_mature = sum(coalesce(length_mm, 0) > L50[full_spp], na.rm = TRUE),
    count_immature = sum(coalesce(length_mm, 0) <= L50[full_spp], na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  
  # Ensure that all full_spp are present even if no data is available for the species
  complete(opcode, full_spp, fill = list(count_mature = 0, count_immature = 0)) %>% 
  complete(opcode = unique(metadata$opcode), full_spp = full_spp, fill = list(count_mature = 0, count_immature = 0)) %>% 
  
  # Fill in missing longitude and latitude
  left_join(metadata %>% distinct(opcode, longitude_dd, latitude_dd), by = "opcode") %>% 
  mutate(
    longitude = coalesce(longitude_dd.x, longitude_dd.y),
    latitude = coalesce(latitude_dd.x, latitude_dd.y)
  ) %>%
  dplyr::select(-longitude_dd.x, -latitude_dd.x, -longitude_dd.y, -latitude_dd.y) %>% 
  
  # Create an ID column
  mutate(ID = row_number()) %>% 
  glimpse()

length(unique(dat$opcode)) == length(unique(metadata$opcode)) # should be true. if not there are opcodes missing.

ggplot(dat, aes(x = longitude, y = latitude, size = count_mature)) +
  facet_wrap(~full_spp) +
  geom_point()

saveRDS(dat, "data/tidy/mature_presence_latlong_all_lengths.rds") # Needed for 03_gam

summary(as.factor(dat$full_spp[dat$count_mature>0]))


## Environmental data ---------------------------------------------------------

### Observed habitat formatting -----------------------------------------------

# workflow from C. Spencer, https://globalarchivemanual.github.io/CheckEM/articles/r-workflows/check-habitat.html
metadata <- read_metadata(here::here("data/raw/2024_commonwealth/")) %>%
  dplyr::select(campaignid, sample, longitude_dd, latitude_dd, date_time, location, site, depth_m, successful_count, successful_length, successful_habitat_forwards, successful_habitat_backwards) %>%
  glimpse()

points <- read_TM(here::here("data/raw/2024_commonwealth/"),
                  sample = "opcode")

habitat <- points %>%
  dplyr::filter(relief_annotated %in% "No") %>%
  dplyr::select(campaignid, sample, starts_with("level"), scientific) %>%
  glimpse()

relief <- points %>%
  dplyr::filter(relief_annotated %in% "Yes") %>%
  dplyr::select(campaignid, sample, starts_with("level"), scientific) %>%
  glimpse()

num.points <- 20

# Check nothing's missing, that everything is where it should be
wrong.points.habitat <- habitat %>%
  group_by(campaignid, sample) %>%
  summarise(points.annotated = n()) %>%
  left_join(metadata) %>%
  dplyr::mutate(expected = case_when(successful_habitat_forwards %in% "Yes" & successful_habitat_backwards %in% "Yes" ~ num.points * 2, 
                                     successful_habitat_forwards %in% "Yes" & successful_habitat_backwards %in% "No" ~ num.points * 1, 
                                     successful_habitat_forwards %in% "No" & successful_habitat_backwards %in% "Yes" ~ num.points * 1, 
                                     successful_habitat_forwards %in% "No" & successful_habitat_backwards %in% "No" ~ num.points * 0)) %>%
  dplyr::filter(!points.annotated == expected) %>%
  glimpse() # if there are no rows, all good.

wrong.points.relief <- relief %>%
  group_by(campaignid, sample) %>%
  summarise(points.annotated = n()) %>%
  left_join(metadata) %>%
  dplyr::mutate(expected = case_when(successful_habitat_forwards %in% "Yes" & successful_habitat_backwards %in% "Yes" ~ num.points * 2, 
                                     successful_habitat_forwards %in% "Yes" & successful_habitat_backwards %in% "No" ~ num.points * 1, 
                                     successful_habitat_forwards %in% "No" & successful_habitat_backwards %in% "Yes" ~ num.points * 1, 
                                     successful_habitat_forwards %in% "No" & successful_habitat_backwards %in% "No" ~ num.points * 0)) %>%
  dplyr::filter(!points.annotated == expected) %>%
  glimpse() # if there are no rows, all good.

habitat.missing.metadata <- anti_join(habitat, metadata, by = c("sample")) %>%
  glimpse() # if there are no rows, all good.

metadata.missing.habitat <- anti_join(metadata, habitat, by = c("sample")) %>%
  glimpse() # here there's one drop that doesn't have habitat, normal because the files are overall missing. it's the only one in Geographe. Normal. All chill.



tidy.habitat <- habitat %>%
  dplyr::mutate(number = 1) %>%                                     
  left_join(catami) %>%
  dplyr::select(campaignid, sample, number, starts_with("level"), family, genus, species) %>%
  dplyr::filter(!level_2 %in% c("","Unscorable", NA)) %>%  
  group_by(campaignid, sample, across(starts_with("level")), family, genus, species) %>%
  dplyr::tally(number, name = "number") %>%
  ungroup() %>%                                                     
  dplyr::select(campaignid, sample, level_1, everything()) %>%
  left_join(metadata %>% dplyr::select(sample, longitude_dd, latitude_dd), by = "sample") %>%
  rename(opcode = sample) %>% 
  glimpse()

# Save for 02_habitat_raster
saveRDS(tidy.habitat, file = here::here(paste0("data/tidy/2024_geographe_tidy_habitat.rds")))


### Bathymetry ----------------------------------------------------------------

# 250m bathymetry
bathy_250m <- rast(file_250m_bathy) %>%
  raster::crop(as(extent(ext), 'SpatialPolygons')); plot(bathy_250m)

# 250m Bathymetry derivatives
bathy_250m_d <- readRDS(file_250m_bathy_deriv) %>%
  raster::crop(ext); plot(bathy_250m_d)

depth_250       <- bathy_250m_d[[1]]; plot(depth_250)#;      writeRaster(bathy_250m_d[[1]], "data/spatial/rasters/geographe_250m_depth.tif", overwrite = T)
aspect_250      <- bathy_250m_d[[2]]; plot(aspect_250)#;     writeRaster(bathy_250m_d[[2]], "data/spatial/rasters/geographe_250m_aspect.tif", overwrite = T)
roughness_250   <- bathy_250m_d[[3]]; plot(roughness_250)#;  writeRaster(bathy_250m_d[[3]], "data/spatial/rasters/geographe_250m_roughness.tif", overwrite = T)
detrended_250   <- bathy_250m_d[[4]]; plot(detrended_250)#;  writeRaster(bathy_250m_d[[4]], "data/spatial/rasters/geographe_250m_detrended_bathy.tif", overwrite = T)

# LiDAR bathymetry
bathy_lidar <- rast(file_lidar_bathy); bathy_lidar <- terra::project(bathy_lidar, crs("+proj=longlat +datum=WGS84 +no_defs"))
bathy_lidar <- crop(bathy_lidar, ext); plot(bathy_lidar)

# LiDAR bathymetry derivatives (can't make detrended bathy so skipping that)
bathy_lidar_d <- raster::terrain(bathy_lidar[[1]], neighbors = 8,
                                 v = c("aspect", "roughness"),
                                 unit = "degrees"); plot(bathy_lidar_d)

aspect_lidar      <- bathy_lidar_d[[1]]; plot(aspect_lidar)
roughness_lidar   <- bathy_lidar_d[[2]]; plot(roughness_lidar)


## Making a file stack for all bathymetry data --------------------------------

predictors <- list(bathy_lidar = bathy_lidar[[1]], roughness_lidar = roughness_lidar, aspect_lidar = aspect_lidar,
                   aspect_250 = aspect_250, depth_250 = depth_250, roughness_250 = roughness_250, detrended_250 = detrended_250)

for (i in 1:length(predictors)) { # Unifying extents and resampling so that they can be stacked
  reference_raster <- predictors[[1]] # Use the raster with the highest resolution as a reference. Otherwise the re-sampling will apply lower resolutions to everything.
  predictors[[i]] <- terra::resample(predictors[[i]], reference_raster, method = "bilinear")
  cat(paste("New extent of", i, ":", ext(predictors[[i]]), "\n"))
}

list2env(predictors, envir = .GlobalEnv) # Put the elements of the list in the environment

predictors <- c(
  bathy_lidar, roughness_lidar, aspect_lidar,          # LiDAR data
  depth_250, aspect_250, roughness_250, detrended_250) # 250m data

names(predictors) <- c("bathy_lidar", "roughness_lidar", "aspect_lidar", # Add clear names
                       "bathy_250m", "aspect_250m", "roughness_250m", "detrended_250m")
plot(predictors)

saveRDS(predictors, file = "data/spatial/rasters/geographe_bathy_predictors.rds") # This takes a while


### END ###
