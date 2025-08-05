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
library(stars) # for detrended bathymetry
library(starsExtra) # for detrended bathymetry

# Clear memory
rm(list=ls())

study_site <- "waatu" # waatu or waatern

ext <- if(study_site == "waatern") {
  ext <- c(115.035, 115.68012207, -33.69743897, -33.35) 
} else { 
  ext <- c(114.7, 115.05, -34.25, -33.95)
} # Extent to crop layers to


## Files used in this script --------------------------------------------------

file_metadata <- if (study_site == "waatern") {
  "data/raw/2024_waatern_commonwealth/2024-04_Geographe_stereo-BRUVs_Metadata.rds"
} else {
  "data/raw/2024_waatu_commonwealth/2024-10_SwC_stereo-BRUVs_Metadata.csv"
}

file_maxn <- ifelse(study_site == "waatern", 
                    "data/tidy/2024_waatern_all_tidy_maxn.csv",
                    "data/tidy/2024_waatu_all_tidy_maxn.csv")

file_obs_hab            <- ifelse(study_site == "waatern", 
                                  "data/raw/2024_waatern_commonwealth/2024-04_Geographe_stereo-BRUVs_forwards_Dot Point Measurements.txt",
                                  "data/raw/2024_waatu_state/benthos_summarise.RDS")

file_250m_bathy         <- "QGIS layers/clean/wadandi_250m_bathy.tif" # from Geoscience Australia


## Presence/Absence data ------------------------------------------------------

metadata <- if (study_site == "waatern") {
  metadata <- readRDS(file_metadata)
} else {
  metadata <- read.csv(file_metadata) %>% dplyr::filter(!is.na(ip))
} 

if (study_site == "waatern") {
  metadata <- metadata %>% rename(opcode = sample)
}

metadata <- metadata %>% # reading the right format depending on the study site
  rename(depth = depth_m) %>% 
  dplyr::filter(str_detect(opcode, "BV")) %>%  # take only the spatially balanced drops. this is for the distribution model
  dplyr::filter(successful_count == "Yes") %>% 
  dplyr::select("opcode", "latitude_dd", "longitude_dd", "status", "depth") %>% 
  glimpse()

dat <- read.csv(file_maxn) %>% 
  dplyr::filter(fullspp %in% c('Sparidae Chrysophrys auratus', 
                               'Labridae Ophthalmolepis lineolatus', 
                               'Glaucosomatidae Glaucosoma hebraicum'),
                opcode %in% metadata$opcode) %>% 

  # Join the metadata with the length data
  left_join(metadata %>% distinct(opcode, longitude_dd, latitude_dd), by = "opcode") %>% 
  
  # Ensure that all fullspp are present even if no data is available for the species
  complete(opcode = unique(metadata$opcode), fullspp = fullspp, fill = list(MaxN = 0)) %>% 
  
  # Fill in missing longitude and latitude
  left_join(metadata %>% distinct(opcode, longitude_dd, latitude_dd, depth, status), by = "opcode") %>%
  mutate(
    longitude = coalesce(longitude_dd.x, longitude_dd.y),
    latitude = coalesce(latitude_dd.x, latitude_dd.y),
    depth   = coalesce(depth.x, as.numeric(depth.y)),  # convert to numeric if needed
    status  = coalesce(status.x, status.y),
  ) %>%
  dplyr::select(-longitude_dd.x, -latitude_dd.x, -longitude_dd.y, -latitude_dd.y, 
                -depth.x, -depth.y, -status.x, -status.y, -sd, -X) %>% 
  
  # Create an ID column
  mutate(ID = row_number()) %>% 
  glimpse()

length(unique(dat$opcode)) == length(unique(metadata$opcode)) # should be true. if not there are opcodes missing.

ggplot(dat, aes(x = longitude, y = latitude, size = MaxN)) +
  facet_wrap(~fullspp) +
  geom_point()

saveRDS(dat, paste0("data/tidy/L01_", study_site, "_presence_latlong.rds")) # Needed for 03_gam


## Environmental data ---------------------------------------------------------

### Observed habitat formatting -----------------------------------------------

# waatu-specific workflow
if (study_site == "waatu") {
  library(SQAPI)
  token <- readRDS("secrets/api_token.rds") %>%
    writeClipboard()
  api <- SQAPI$new()
  
  filter <- query_filter(name = "annotation_set_id", op = "eq", val = "17486") # this requests to download data from annotation_set_id that equal 17486, i.e the 2024-10-SWC
  r <- export(
    api = api,
    endpoint = "api/annotation/export",
    query_filters = filter,
    template = "data.csv"
  )
  waatu_hab <- parse_api(r, filetype = "csv")
} # this statement retrieves habitat data for Waatu/capes from Squidle. This only needs to be done once.

if (study_site == "waatu") {
  test_hab <- waatu_hab %>% 
    dplyr::select(label.name, point.media.deployment.name) %>% 
    rename(opcode = point.media.deployment.name,
           habitat = label.name) %>% 
    glimpse()
  
  tidy.habitat <- merge(test_hab, metadata, by = "opcode") %>% 
    glimpse()
  
  # Save for 02_habitat_raster
  saveRDS(tidy.habitat, file = here::here(paste0("data/tidy/L01_", study_site, "_tidy_habitat.rds")))
} # this statement cleans up the habitat data in a format we can use, and saves it for L02.

# waatern-specific workflow
if (study_site == "waatern") {
  # workflow from C. Spencer, https://globalarchivemanual.github.io/CheckEM/articles/r-workflows/check-habitat.html
  metadata <- read_metadata(here::here(paste0("data/raw/2024_", study_site, "_commonwealth/"))) %>%
    dplyr::select(campaignid, sample, longitude_dd, latitude_dd, date_time, location, site, depth_m, successful_count, successful_length, successful_habitat_forwards, successful_habitat_backwards) %>%
    glimpse()
  
  points <- read_TM(here::here(paste0("data/raw/2024_", study_site, "_commonwealth/")),
                    sample = "opcode"
  )
  
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
  saveRDS(tidy.habitat, file = here::here(paste0("data/tidy/L01_", study_site, "_tidy_habitat.rds")))
} # this statement reads waatern habitat data and checks it.

### Bathymetry ----------------------------------------------------------------

# 250m bathymetry
bathy_250 <- rast(file_250m_bathy) %>%
  raster::crop(as(extent(ext), 'SpatialPolygons')); plot(bathy_250)

# 250m Bathymetry derivatives
aspect_250 <- terra::terrain(bathy_250, "aspect", unit = "degrees", neighbors = 8); plot(aspect_250)
roughness_250 <- terrain(bathy_250, "roughness", neighbors = 8); plot(roughness_250)

zstar <- st_as_stars(bathy_250)
detrended_250 <- detrend(zstar, parallel = 8) %>%
  rast()
names(detrended_250) <- c("detrended", "lineartrend")
detrended_250 <- detrended_250[["detrended"]]; plot(detrended_250)


## Making a file stack for all bathymetry data --------------------------------

predictors <- list(aspect_250 = aspect_250, depth_250 = bathy_250, roughness_250 = roughness_250, detrended_250 = detrended_250)

for (i in 1:length(predictors)) { # Unifying extents and resampling so that they can be stacked
  reference_raster <- predictors[[1]] # Use the raster with the highest resolution as a reference. Otherwise the re-sampling will apply lower resolutions to everything.
  predictors[[i]] <- terra::resample(predictors[[i]], reference_raster, method = "bilinear")
  cat(paste("New extent of", i, ":", ext(predictors[[i]]), "\n"))
}

list2env(predictors, envir = .GlobalEnv) # Put the elements of the list in the environment

predictors <- c(depth_250, aspect_250, roughness_250, detrended_250) # 250m data

names(predictors) <- c("bathy_250m", "aspect_250m", "roughness_250m", "detrended_250m")
plot(predictors)

saveRDS(predictors, file = paste0("data/tidy/L01_", study_site, "_bathy_predictors.rds"))


### END ###
