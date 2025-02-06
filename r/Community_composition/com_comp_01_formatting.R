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

library(dplyr)


# Clear memory
rm(list=ls())


# Files used in this script

file_spabal_metadata <- "data/raw/2024_geographe_metadata.csv"
file_spabal_maxn <- "data/raw/2024_geographe_points.txt"

file_pref_metadata <- "data/raw/2024_state/2024_02_NgariCapes.MP.Monitoring_stereoBRUVs_Metadata.csv" # Final dataset
file_pref_maxn <- "data/raw/2024_state/2024_02_NgariCapes.MP.Monitoring_stereoBRUVs_Points 1.txt" # Final dataset


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
write.csv(dat, file = paste0("data/tidy/2024_geographe_all_tidy_maxn.csv"))



### END ###
