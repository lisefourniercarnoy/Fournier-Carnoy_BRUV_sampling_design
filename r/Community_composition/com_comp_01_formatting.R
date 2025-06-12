# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (B. Compare communities detected)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Tidy-up dataframes.
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    August 2024

# -----------------------------------------------------------------------------

# -- The aim here is to format the MaxNs and combine into one, so that it can be
# -- used in NMDS. Some things need to be done, including adding metadata and
# -- selecting only EGB drops from the state data.

# -----------------------------------------------------------------------------

library(tidyverse)


# Clear memory
rm(list=ls())


# Waatern ---------------------------------------------------------------------

file_spabal_metadata <- "data/raw/2024_waatern_commonwealth/2024-04_Geographe_stereo-BRUVs_Metadata.csv"
file_spabal_maxn <- "data/raw/2024_waatern_commonwealth/2024_geographe_points.txt"

file_pref_metadata <- "data/raw/2024_waatern_state/2024_02_NgariCapes.MP.Monitoring_stereoBRUVs_Metadata.csv"
file_pref_maxn <- "data/raw/2024_waatern_state/2024_02_NgariCapes.MP.Monitoring_stereoBRUVs_Points 1.txt"


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
maxn <- read.table(file_spabal_maxn, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  rename(opcode = OpCode) %>%
  group_by(opcode, Family, Genus, Species) %>%
  summarise(MaxN = sum(Number, na.rm = TRUE)) %>%  # Get MaxN
  mutate(fullspp = paste(Family,Genus, Species)) %>%
  filter(!(fullspp == " " | fullspp == " spp")
         & MaxN > 0) %>%  # Remove empty/incomplete spp names.
  glimpse()

# Tidy-up
cw_tidy_maxn <- maxn %>%
  left_join(spabal_metadata) %>%
  group_by(opcode) %>%
  dplyr::select(opcode, fullspp, MaxN, longitude, latitude, depth, status) %>%
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
  summarise(MaxN = n()) %>%
  ungroup() %>% # Counting each frame's number of each species ...
  group_by(fullspp, OpCode) %>%
  slice(which.max(MaxN))%>% # ... then selecting the maxn per opcode
  ungroup() %>%
  dplyr::select(OpCode, fullspp, MaxN) %>%
  rename(opcode = OpCode) %>%
  glimpse()


# Tidy-up
state_tidy_maxn <- maxn %>%
  left_join(pref_metadata) %>%
  group_by(longitude, latitude) %>%
  dplyr::select(opcode, fullspp, MaxN, longitude, latitude, depth, status) %>%
  dplyr::filter(grepl('EGB', opcode, fixed = T)) %>%
  mutate(sd = "preferential") %>%
  glimpse()


## Combining into one ---------------------------------------------------------

names(state_tidy_maxn); names(cw_tidy_maxn) # should be the exact same
dat <- bind_rows(state_tidy_maxn, cw_tidy_maxn)
write.csv(dat, file = paste0("data/tidy/2024_waatern_all_tidy_maxn.csv"))

# Waatu -----------------------------------------------------------------------
file_spabal_metadata <- "data/raw/2024_waatu_commonwealth/2024-10_SwC_stereo-BRUVs_Metadata.csv"
file_spabal_maxn <- "data/raw/2024_waatu_commonwealth/2024-10_SwC_stereo-BRUVs_points.txt"

file_pref_metadata <- "data/raw/2024_waatu_state/metadata.RDS"
file_pref_maxn <- "data/raw/2024_waatu_state/count.RDS"

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
maxn <- read.table(file_spabal_maxn, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  rename(opcode = OpCode) %>%
  group_by(opcode, Family, Genus, Species) %>%
  summarise(MaxN = sum(Number, na.rm = TRUE)) %>%  # Get MaxN
  mutate(fullspp = paste(Family,Genus, Species)) %>%
  filter(!(fullspp == " " | fullspp == " spp")
         & MaxN > 0) %>%  # Remove empty/incomplete spp names.
  glimpse()

# Tidy-up
cw_tidy_maxn <- maxn %>%
  left_join(spabal_metadata) %>%
  group_by(opcode) %>%
  dplyr::select(opcode, fullspp, MaxN, longitude, latitude, depth, status) %>%
  mutate(sd = ifelse(grepl('P', opcode) & longitude>115.3, "preferential", "spatially balanced")) %>%
  glimpse()


## State data setup -----------------------------------------------------------

# -- Using DBCA BRUV data, I'm combining the metadata information to the MaxNs

# Metadata
test <- pref_metadata[pref_metadata$date_time >= as.Date("2024-01-01") & 
                pref_metadata$date_time < as.Date("2025-01-01"), ]

pref_metadata <- readRDS(file_pref_metadata) %>%
  filter(date_time >= as.Date("2024-01-01") & date_time < as.Date("2025-01-01")) %>% # select only 2024 data
  dplyr::select(sample, longitude_dd, latitude_dd, depth_m, status) %>%
  mutate(depth_m = as.numeric(depth_m)) %>%
  rename(opcode = sample,
         longitude = longitude_dd,
         latitude = latitude_dd,
         depth = depth_m) %>%
    filter(str_detect(opcode, "CF")) %>% # select cape freycinet points
  
  distinct(opcode, .keep_all = TRUE) %>%  # Keep the first occurrence of each 'opcode', some of them are duplicated fsr
  glimpse()


# MaxN
maxn <- readRDS(file_pref_maxn) %>%
  left_join(readRDS(file_pref_metadata) %>% dplyr::select(c(sample, sample_url)) %>% rename(opcode = sample), by = "sample_url") %>% 
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
  rename(MaxN = count) %>% 
  mutate(opcode = str_replace_all(opcode, "_", "-")) %>% # format opcode so it's the same as metadata
  glimpse()


# Tidy-up
state_tidy_maxn <- maxn %>%
  left_join(pref_metadata) %>%
  group_by(longitude, latitude) %>%
  dplyr::select(opcode, fullspp, MaxN, longitude, latitude, depth, status) %>%
  dplyr::filter(grepl('CF', opcode, fixed = T)) %>%
  mutate(sd = "preferential") %>%
  glimpse()


## Combining into one ---------------------------------------------------------

names(state_tidy_maxn); names(cw_tidy_maxn) # should be the exact same
dat <- bind_rows(state_tidy_maxn, cw_tidy_maxn)
write.csv(dat, file = paste0("data/tidy/2024_waatu_all_tidy_maxn.csv"))



### END ###