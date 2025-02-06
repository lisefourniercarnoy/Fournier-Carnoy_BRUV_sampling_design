# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (B. Compare communities detected)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Make tidy files for future scripts
# Author:  Lise Fournier-Carnoy
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  LIDAR ELEMENTS ARE COMMENTED OUT BECAUSE FILES TOO BIG FOR GIT

# -----------------------------------------------------------------------------

library(dismo)
library(tidyverse)
library(sf)
library(terra)
library(raster)

# Clear memory
rm(list=ls())


ext <- c(115.035, 115.68012207, -33.69743897, -33.35) # Extent to crop layers to

## Files used in this script --------------------------------------------------

file_fish_length        <- "data/raw/2024_commonwealth/2024_geographe_lengths.txt"
file_metadata           <- "data/raw/2024_commonwealth/2024_geographe_metadata.csv"

file_hab_ras            <- "outputs/Length_comparison/2024_geographe_predicted-habitat.tif"

file_250m_bathy         <- "data/spatial/rasters/GB-SW_250mBathy.tif"
file_250m_bathy_deriv   <- "data/spatial/rasters/2014_geographe_bathymetry-derivatives.RDS"
#file_lidar_bathy        <- "data/spatial/rasters/LiDAR_geo_compressed.tif" # Too big for git


## Presence/Absence data ------------------------------------------------------

L50 <- c(#'Chrysophrys auratus' = 375, # L50 from Wakefield et al. 2015
         #'Pseudocaranx spp' = 279, # L50 from male P. dentex, Farmer et al. 2005
         'Glaucosoma hebraicum' = 298, # L50 from Evans-Powell 2022
         'Achoerodus gouldii' = 653, # L50 for females, Coulson et al 2007
         'Choerodon rubescens' = 264, # L50 from females, Nardi 2006
         'Epinephelides armatus' = 250.6, # L50 from males, Moore et al 2007
         'Nemadactylus valenciennesi' = 400) # L50 from females, Coulson et al. 2010


data <- read.table(file_metadata, header = TRUE, sep = ",", stringsAsFactors = FALSE)

length_data <- read.table(file_fish_length, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% glimpse()

data <- left_join(length_data, data, by = c("OpCode" = "opcode")) %>%
  # Select species
  mutate(full_spp = paste(Genus, Species)) %>%
  filter(full_spp %in% names(L50) & !is.na(Length)) %>%
  group_by(longitude_dd, latitude_dd) %>%

  # Separate into size categories
  summarise(
    count_mature = sum(Length > L50[full_spp], na.rm = TRUE),
    count_immature = sum(Length < L50[full_spp], na.rm = TRUE),
    .groups = 'drop'
  ) %>%

  # Remove duplicate rows
  distinct(longitude_dd, .keep_all = TRUE) %>%

  # Add an ID column
  mutate(ID = row_number()) %>%
  rename(longitude = longitude_dd,
         latitude = latitude_dd) %>%
  glimpse()

# Quickly plot where they are
ggplot(data, aes(x = longitude, y = latitude, size = count_mature)) +
  geom_point()


## Environmental data ---------------------------------------------------------

# Habitat - script 02_habitat_raster must be run before this
hab <- stack(file_hab_ras) %>%
  raster::crop(ext) %>%
  glimpse(); plot(hab)

reef      <- rast(hab[[1]])%>%raster::crop(ext); plot(reef)#;      writeRaster(hab[[1]], "data/spatial/rasters/geographe_reef.tif", overwrite = T)
sand      <- rast(hab[[2]])%>%raster::crop(ext); plot(sand)#;      writeRaster(hab[[6]], "data/spatial/rasters/geographe_sand.tif", overwrite = T)
seagrass  <- rast(hab[[3]])%>%raster::crop(ext); plot(seagrass)#;  writeRaster(hab[[3]], "data/spatial/rasters/geographe_seagrass.tif", overwrite = T)

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
#bathy_lidar <- rast(file_lidar_bathy); bathy_lidar <- terra::project(bathy_lidar, crs("+proj=longlat +datum=WGS84 +no_defs"))
#bathy_lidar <- crop(bathy_lidar, ext); plot(bathy_lidar)

# LiDAR bathymetry derivatives (can't make detrended bathy so skipping that)
#bathy_lidar_d <- raster::terrain(bathy_lidar[[1]], neighbors = 8,
#                                 v = c("aspect", "roughness"),
#                                 unit = "degrees"); plot(bathy_lidar_d)

#aspect_lidar      <- bathy_lidar_d[[1]]; plot(aspect_lidar)
#roughness_lidar   <- bathy_lidar_d[[2]]; plot(roughness_lidar)


## Making a file stack for all environmental data -----------------------------

predictors <- list(#bathy_lidar = bathy_lidar[[1]], roughness_lidar = roughness_lidar, aspect_lidar = aspect_lidar,
                   aspect_250 = aspect_250, depth_250 = depth_250, roughness_250 = roughness_250, detrended_250 = detrended_250)

for (i in 1:length(predictors)) { # Unifying extents and resampling so that they can be stacked
  reference_raster <- predictors[[1]] # Use the raster with the highest resolution as a reference. Otherwise the re-sampling will apply lower resolutions to everything.
  predictors[[i]] <- terra::resample(predictors[[i]], reference_raster, method = "bilinear")
  cat(paste("New extent of", i, ":", ext(predictors[[i]]), "\n"))
}

list2env(predictors, envir = .GlobalEnv) # Put the elements of the list in the environment

predictors <- c(
  #bathy_lidar, roughness_lidar, aspect_lidar,          # LiDAR data
  depth_250, aspect_250, roughness_250, detrended_250) # 250m data

names(predictors) <- c(#"bathy_lidar", "roughness_lidar", "aspect_lidar", # Add clear names
                       "aspect_250m", "bathy_250m", "roughness_250m", "detrended_250m")
plot(predictors)

saveRDS(predictors, file = "data/spatial/rasters/geographe_bathy_predictors.rds") # This takes a while

saveRDS(data, "data/tidy/mature_presence_latlong.rds") # Needed for 03_gam

### END ###
