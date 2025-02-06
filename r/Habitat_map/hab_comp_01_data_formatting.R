# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (C. Habitat map comp.)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Tidy-up dataframes.
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  Done for now. Just need Commonwealth data once it's processed

# -----------------------------------------------------------------------------


# Clear memory
rm(list=ls())

# Load libraries
library(CheckEM)
library(tidyverse)
library(ggplot2)
library(arrow)

# Set the study name
name <- "2024_geographe"


#### SPATIALLY BALANCED needs actual data ----

# Load and format metadata
spabal_metadata <- read.csv("data/raw/2014-Geographe-stereo-BRUVs.checked.metadata.csv") %>%
  dplyr::select(sample, longitude, latitude, depth, status) %>%
  rename(OpCode = sample) %>%
  glimpse()

write.csv(spabal_metadata, file = paste0("data/tidy/", name, "_spabal_tidy_metadata.csv"))

spabal_habitat <- read.table("data/raw/2014-12_Geographe.Bay_Habitat.point.score.txt", sep = "\t", header = T) %>%
  dplyr::select(OpCode, Relief0, Relief1, Relief2, Relief3, Relief4,
                fieldofview.Limited, fieldofview.Open,
                Consolidated, Macroalgae, Seagrasses, Stony.corals, Unconsolidated) %>%
  rename(total_pts = fieldofview.Open) # Not sure if this is right, but 2014 data doesn't have total_pts
plot(spabal_metadata$longitude, spabal_metadata$latitude)

# Tidy-up

names(spabal_habitat); names(spabal_metadata) # 'opcode' should match in order to join the two

spabal_tidy_habitat <- spabal_habitat %>%
  left_join(spabal_metadata) %>%
  group_by(OpCode, longitude, latitude, depth) %>%
  glimpse()
saveRDS(spabal_tidy_habitat, file = paste0("data/tidy/", name, "_spabal_tidy_habitat.rds"))
plot(spabal_tidy_habitat$longitude, spabal_tidy_habitat$latitude)


#### PREFERENTIAL DATA ----

# Load and format metadata

pref_metadata <- read.csv("data/raw/2024_state/2024_02_NgariCapes.MP.Monitoring_stereoBRUVs_Metadata.csv") %>%
  dplyr::select(opcode, longitude_dd, latitude_dd, depth_m, status) %>%
  rename(longitude = longitude_dd,
         latitude = latitude_dd,
         depth = depth_m) %>%
  glimpse()
write.csv(pref_metadata, file = paste0("data/tidy/", name, "_pref_tidy_metadata.csv"))

cols <- c("Filename", "remove", "Image row", "Image col", "CampaignID", "remove2", "opcode", "Transect", "level1", "level2", "Radius %")

pref_habitat <- read.table("data/raw/2024_state/NCMP_2024_Dot Point Measurements.txt", sep = "\t", header = F, skip = 5, fill = T) %>%
  select(where(~ !(all(is.na(.)) | all(. == ""))))  %>%
  setNames(cols) %>%
  mutate(habitat = case_when( # Categorising habitat to match Spatially Balanced dataset
    level1 == "Macroalgae" ~ "Macroalgae",
    level1 == "Substrate" & level2 %in% c("Consolidated (hard)") ~ "Consolidated",
    level1 == "Substrate" & level2 %in% c("Unconsolidated (soft)") ~ "Unconsolidated",
    level1 == "Cnidaria" ~ "Stony.corals",
    level1 == "Fishes" ~ "Unscorable",
    level1 == "Seagrasses" ~ "Seagrasses",
    level1 == "Sessile invertebrates" ~ "Invertebrate",
    level1 == "Unscorable" ~ "Unscorable",
    level1 == "Sponges" ~ "Invertebrate",
    TRUE ~ "Other"
  )) %>%
  select(opcode, habitat) %>%
  group_by(opcode, habitat) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ungroup() %>%
  filter(habitat != "Unscorable")

total_pts_per_opcode <- pref_habitat %>%
  group_by(opcode) %>%
  summarise(total_pts = sum(count), .groups = 'drop')
pref_habitat <- pref_habitat %>%
  left_join(total_pts_per_opcode, by = "opcode") %>%
  arrange(opcode, habitat) %>%
  pivot_wider(
    names_from = habitat,   # Column names for the new columns
    values_from = count,    # Values for the new columns
    values_fill = list(count = 0)  # Fill missing values with 0
  ) %>%
  glimpse()


# Tidy-up

names(pref_habitat); names(pref_metadata) # 'opcode' should match in order to join the two

pref_tidy_habitat <- pref_habitat %>%
  left_join(pref_metadata) %>%
  group_by(opcode, longitude, latitude, depth) %>%
  glimpse()
saveRDS(pref_tidy_habitat, file = paste0("data/tidy/", name, "_pref_tidy_habitat.rds"))
