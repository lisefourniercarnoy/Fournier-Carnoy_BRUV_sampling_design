# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (B. Compare communities detected)
# Data:    2024 MEGlab BRUV data.
# Task:    Tidy-up dataframes.
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    August 2024

# -----------------------------------------------------------------------------

library(tidyverse)
library(sf)

# Clear memory
rm(list=ls())


# Waatern ---------------------------------------------------------------------

file_spabal_metadata  <- "data/tidy/2024_waatern_commonwealth/2024_waatern_BRUVs_metadata.rds"
file_spabal_maxn      <- "data/tidy/00_2024_waatern_BRUVs_tidy_counts.rds"
file_spabal_habitat <-  "data/tidy/L01_waatern_tidy_habitat.rds"

file_pref_metadata    <- "data/raw/2024_waatern_state/2024_02_NgariCapes.MP.Monitoring_stereoBRUVs_Metadata.csv"
file_pref_maxn        <- "data/raw/2024_waatern_state/2024_02_NgariCapes.MP.Monitoring_stereoBRUVs_Points 1.txt"
file_pref_habitat <- "data/raw/2024_waatu_state/2024-02_NgariCapes.MP.Monitoring_stereoBRUVs_benthos.csv"

## Commonwealth data setup ----------------------------------------------------

# -- Using MEGLAB BRUV data, I'm combining the metadata information to the MaxNs

# Metadata
spabal_metadata <- readRDS(file_spabal_metadata) %>%
  mutate(longitude_dd = as.numeric(longitude_dd),
         latitude_dd = as.numeric(latitude_dd),
         depth_m = as.numeric(depth_m)) %>% 
  rename(longitude = longitude_dd,
         latitude = latitude_dd,
         depth = depth_m,
         opcode = sample) %>%  
  dplyr::select(opcode, longitude, latitude, depth, status) %>%
  glimpse()


# MaxN
maxn <- readRDS(file_spabal_maxn) %>%
  rename(opcode = sample,
         fullspp = scientific) %>%
  group_by(opcode, fullspp) %>%
  summarise(maxn = sum(count, na.rm = TRUE)) %>%  # Get MaxN
  filter(!(fullspp == " " | fullspp == " spp") & maxn > 0) %>%  # Remove empty/incomplete spp names.
  glimpse()

# Habitat
hab <- readRDS(file_spabal_habitat) %>% 
  st_drop_geometry() %>% 
  glimpse

# Tidy-up
cw_tidy_maxn <- maxn %>%
  left_join(spabal_metadata) %>%
  left_join(hab) %>% 
  group_by(opcode) %>%
  dplyr::select(opcode, fullspp, maxn, longitude, latitude, depth, status, names(hab)) %>%
  mutate(sd = ifelse(grepl('P', opcode) & longitude>115.3, "preferential", "spatially balanced")) %>%
  glimpse()


## State data setup -----------------------------------------------------------

# -- Using DBCA BRUV data, I'm combining the metadata information to the MaxNs

# Metadata
pref_metadata <- read.csv(file_pref_metadata) %>%
  dplyr::select(opcode, longitude_dd, latitude_dd, depth_m, status) %>%
  mutate(depth_m = as.numeric(depth_m)) %>%
  rename(longitude = longitude_dd,
         latitude = latitude_dd,
         depth = depth_m) %>%
  distinct(opcode, .keep_all = TRUE) %>%  # Keep the first occurrence of each 'opcode', some of them are duplicated fsr
  glimpse()

# MaxN
maxn <- read.table(file_pref_maxn, sep = "\t", header = T) %>%
  mutate(fullspp = paste(Family, Genus, Species)) %>%
  mutate(fullspp = case_when(
    str_detect(fullspp, "spilomelanurus|vittiger") ~ "Monacanthidae Acanthaluteres sp1",  # change these two into sp1, for consistency with cw data
    TRUE ~ fullspp
  )) %>%
  filter(!(fullspp == " " | fullspp == " spp")) %>% # Remove empty/incomplete spp names.
  group_by(Frame, fullspp, OpCode) %>%
  summarise(maxn = n()) %>%
  ungroup() %>% # Counting each frame's number of each species ...
  group_by(fullspp, OpCode) %>%
  slice(which.max(maxn))%>% # ... then selecting the maxn per opcode
  ungroup() %>%
  dplyr::select(OpCode, fullspp, maxn) %>%
  rename(opcode = OpCode) %>%
  filter(!(fullspp == "  " | fullspp == " spp") & maxn > 0) %>%  # Remove empty/incomplete spp names.
  glimpse()

# Habitat
hab <- read.csv(file_pref_habitat, skip = 4, header = T) %>% 
  rename(hab1 = X, 
         hab2 = X.1) %>% 
  dplyr::select(Filename, hab1, hab2) %>% 
  filter(grepl("EGB", Filename)) %>% 
  filter(!(hab1 %in% c("Unscorable", "Fishes"))) %>% 
  mutate(
    opcode = sub("NCMP_2024_(EGB[0-9]+-[0-9]+)_*\\.jpg", "\\1", Filename),
    taxa = case_when(hab1 %in% c("Macroalgae") ~ "macro",
                     hab2 %in% c("Encrusting", "Erect fine branching", "Large canopy-forming",
                                 "Macroalgae", "Erect coarse branching", "Ecklonia radiata", 
                                 "Scytothalia dorycarpa", "Filamentous / filiform") ~ "macro",
                     hab2 %in% "Unconsolidated (soft)" ~ "sand",
                     hab2 %in% "Consolidated (bare rock)" ~ "rock",
                     hab1 %in% "Seagrasses" ~ "seagrass",
                     hab2 %in% c("Cup-likes", "Hydroids", "Massive forms", "Sponges", "Crusts", "Hydrocorals",
                                 "Black & Octocorals", "Mixed sessile invertebrates (small unidentifiable)",
                                 "Bryozoa", "Ascidians", "Erect forms") ~ "sessile_inverts"),
    number = 1
  ) %>%
  pivot_wider(names_from = taxa, values_from = number, values_fn = sum, values_fill = 0) %>% 
  
  # make reef and remove the individual reef habitats
  mutate(reef = macro) %>% #  + rock + sessile_inverts - but there are none in the state data
  dplyr::select(-c(macro, hab1, hab2, Filename)) %>% 
  
  # sum rows to have a single row per opcode
  group_by(opcode) %>%  # Group by opcode
  summarize(
    seagrass = sum(seagrass, na.rm = TRUE),
    sand = sum(sand, na.rm = TRUE),
    reef = sum(reef, na.rm = TRUE)
  ) %>%
  
  # add total points - these should add up to all the reef, sand, seagrass
  mutate(total_pts = reef + sand + seagrass) %>% 
  glimpse()

# Tidy-up
state_tidy_maxn <- maxn %>%
  left_join(pref_metadata) %>%
  left_join(hab) %>% 
  group_by(longitude, latitude) %>%
  dplyr::select(opcode, fullspp, maxn, longitude, latitude, depth, status, names(hab)) %>%
  dplyr::filter(grepl('EGB', opcode, fixed = T)) %>%
  mutate(sd = "preferential") %>%
  glimpse()


## Combining into one ---------------------------------------------------------

names(state_tidy_maxn); names(cw_tidy_maxn) # should be the exact same
dat <- bind_rows(state_tidy_maxn, cw_tidy_maxn)
dat <- st_as_sf(dat, coords = c("longitude", "latitude"))
saveRDS(dat, file = paste0("data/tidy/C01_2024_waatern_all_tidy_maxn.rds"))
st_write(dat, dsn = "data/tidy", layer = "C01_2024_waatern_all_tidy_maxn", driver = "ESRI Shapefile", delete_layer = T)

# Waatu -----------------------------------------------------------------------

file_spabal_metadata  <- "data/raw/2024_waatu_commonwealth/2024-10_SwC_stereo-BRUVs_Metadata.csv"
file_spabal_maxn      <- "data/raw/2024_waatu_commonwealth/2024-10_SwC_stereo-BRUVs_points.txt"
file_spabal_habitat   <- "data/tidy/D01_waatu_tidy_habitat.RDS"


file_pref_metadata    <- "data/raw/2024_waatu_state/metadata.RDS"
file_pref_maxn        <- "data/raw/2024_waatu_state/count.RDS"
file_pref_habitat     <- "data/raw/2024_waatu_state/2024-02_NgariCapes.MP.Monitoring_stereoBRUVs_benthos.csv"

## Commonwealth data setup ----------------------------------------------------

# -- Using MEGLAB BRUV data, I'm combining the metadata information to the MaxNs

# Metadata
spabal_metadata <- read.csv(file_spabal_metadata) %>%
  dplyr::select(opcode, longitude_dd, latitude_dd, depth_m, status) %>%
  rename(longitude = longitude_dd,
         latitude = latitude_dd,
         depth = depth_m) %>%
  glimpse()

# MaxN
maxn <- read_tsv(file_spabal_maxn, quote = "\"", show_col_types = FALSE) %>%
  rename(opcode = OpCode) %>%
  group_by(opcode, Family, Genus, Species) %>%
  summarise(maxn = sum(Number, na.rm = TRUE)) %>%  # Get MaxN
  mutate(fullspp = paste(Family, Genus, Species),
         opcode = gsub("_", "-", opcode)) %>%
  filter(!(fullspp == " " | fullspp == " spp")
         & maxn > 0) %>%  # Remove empty/incomplete spp names.
  filter(!opcode %in% c('SWC-BV-133', "SWC-BV-022")) %>% # remove two drops which can't be maxn'ed
  filter(!grepl("IO", opcode)) %>% # these are preferential drops that i'm removing
  filter(!grepl("AP", opcode)) %>% # these are preferential drops that i'm removing
  glimpse()

# Habitat
hab <- readRDS(file_spabal_habitat) %>%
  st_drop_geometry() %>% 
  glimpse()

# Tidy-up
cw_tidy_maxn <- maxn %>%
  left_join(spabal_metadata) %>%
  left_join(hab) %>% 
  group_by(opcode) %>%
  dplyr::select(opcode, fullspp, maxn, longitude, latitude, depth, status, names(hab)) %>%
  mutate(sd = ifelse(grepl('P', opcode) & longitude>115.3, "preferential", "spatially balanced")) %>%
  glimpse()


## State data setup -----------------------------------------------------------

# -- Using DBCA BRUV data, I'm combining the metadata information to the MaxNs

# Metadata
pref_metadata <- readRDS(file_pref_metadata) %>%
  filter(date_time >= as.Date("2024-01-01") & date_time < as.Date("2025-01-01")) %>% # select only 2024 data
  dplyr::select(sample, longitude_dd, latitude_dd, depth_m, status, sample_url) %>%
  mutate(depth_m = as.numeric(depth_m)) %>%
  rename(opcode = sample,
         longitude = longitude_dd,
         latitude = latitude_dd,
         depth = depth_m) %>%
    filter(str_detect(opcode, "CF")) %>% # select cape freycinet points
  glimpse()

# MaxN
maxn <- readRDS(file_pref_maxn) %>%
  filter(sample_url %in% pref_metadata$sample_url) %>% # select only the cafe freycinet drops. original data has Geographe as well
  left_join(pref_metadata, by = "sample_url") %>%  
  rename(
    Family = family,
    Genus = genus,
    Species = species
  ) %>% 
  mutate(fullspp = paste(Family, Genus, Species)) %>%
  mutate(fullspp = case_when(
    str_detect(fullspp, "spilomelanurus|vittiger") ~ "Monacanthidae Acanthaluteres sp1",  # change these two into sp1, for consistency with cw data
    TRUE ~ fullspp
  )) %>%
  filter(!(fullspp == " " | fullspp == " spp")) %>% # Remove empty/incomplete spp names.
  dplyr::select(opcode, fullspp, count) %>%
  rename(maxn = count) %>% 
  mutate(opcode = str_replace_all(opcode, "_", "-")) %>% # format opcode so it's the same as metadata
  glimpse()

# Habitat
hab <- read.csv(file_pref_habitat, skip = 4, header = T) %>% 
  rename(hab1 = X, 
         hab2 = X.1) %>% 
  dplyr::select(Filename, hab1, hab2) %>% 
  filter(grepl("CF", Filename)) %>% 
  filter(!(hab1 %in% c("Unscorable", "Fishes"))) %>% 
  mutate(
    Filename = gsub(" ", "", Filename),  # removes all spaces
    opcode = sub("NCMP_2024_(CF[0-9]+-[0-9]+)\\.jpg", "\\1", Filename),
    taxa = case_when(hab1 %in% c("Macroalgae") ~ "macro",
                     hab2 %in% c("Encrusting", "Erect fine branching", "Large canopy-forming",
                                 "Macroalgae", "Erect coarse branching", "Ecklonia radiata", 
                                 "Scytothalia dorycarpa", "Filamentous / filiform") ~ "macro",
                     hab2 %in% "Unconsolidated (soft)" ~ "sand",
                     hab2 %in% "Consolidated (bare rock)" ~ "rock",
                     hab2 %in% "Thalassodendron spp" ~ "seagrass",
                     hab2 %in% c("Cup-likes", "Hydroids", "Massive forms", "Sponges", "Crusts", "Hydrocorals",
                                 "Black & Octocorals", "Mixed sessile invertebrates (small unidentifiable)",
                                 "Bryozoa", "Ascidians", "Erect forms") ~ "sessile_inverts"),
    number = 1
  ) %>%
  pivot_wider(names_from = taxa, values_from = number, values_fn = sum, values_fill = 0) %>% 
  
  # make reef and remove the individual reef habitats
  mutate(reef = macro) %>% #  + rock + sessile_inverts - but there are none in the state data
  dplyr::select(-c(macro, hab1, hab2, Filename)) %>% 
  
  # sum rows to have a single row per opcode
  group_by(opcode) %>%  # Group by opcode
  summarize(
    sand = sum(sand, na.rm = TRUE),
    reef = sum(reef, na.rm = TRUE)
  ) %>%
  
  # add total points - these should add up to all the reef, sand, seagrass
  mutate(total_pts = reef + sand) %>% 
  glimpse()
hab$seagrass <- 0; glimpse(hab) # adding this so that we can join state and commonwealth

# Tidy-up
state_tidy_maxn <- maxn %>%
  left_join(pref_metadata) %>%
  left_join(hab) %>% 
  group_by(longitude, latitude) %>%
  dplyr::select(opcode, fullspp, maxn, longitude, latitude, depth, status, names(hab)) %>%
  dplyr::filter(grepl('CF', opcode, fixed = TRUE) & !grepl('NCMP', opcode)) %>% # keep the 2024 CF points, not the 2022 ones
  mutate(sd = "preferential") %>%
  glimpse()


## Combining into one ---------------------------------------------------------

names(state_tidy_maxn); names(cw_tidy_maxn) # should be the exact same
dat <- bind_rows(state_tidy_maxn, cw_tidy_maxn) %>% mutate(status = ifelse(status == "Fished", "Fished", "No-take"))
saveRDS(dat, file = paste0("data/tidy/C01_2024_waatu_all_tidy_maxn.rds"))


names(state_tidy_maxn); names(cw_tidy_maxn) # should be the exact same
dat <- bind_rows(state_tidy_maxn, cw_tidy_maxn)
dat <- st_as_sf(dat, coords = c("longitude", "latitude"))
saveRDS(dat, file = paste0("data/tidy/C01_2024_waatu_all_tidy_maxn.rds"))
st_write(dat, dsn = "data/tidy", layer = "C01_2024_waatu_all_tidy_maxn", driver = "ESRI Shapefile", delete_layer = T)


### END ###