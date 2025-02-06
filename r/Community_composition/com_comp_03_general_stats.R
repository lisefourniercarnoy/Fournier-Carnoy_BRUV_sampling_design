# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (B. Compare communities detected)
# Data:    2024 BRUV MaxN data. (MEGLAB + DBCA)
# Task:    Get general numbers for results section
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    December 2024

# -----------------------------------------------------------------------------

# Status:  NMDS done, first try at full-community manyglm()

# -----------------------------------------------------------------------------


# Clear memory
rm(list=ls())


#### Load libraries ----

library(vegan)
library(ggplot2)
library(tidyverse)
library(mvabund)
library(sf)


dat <- read.csv("data/tidy/2024_geographe_all_tidy_maxn.csv") %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  glimpse()

# Total number of spp

tot <- length(unique(dat$fullspp)); tot # total number of species
unk <- length(grep("spp|sp1|sp2|unknown|sp", unique(dat$fullspp), ignore.case = TRUE)); unk # number of uncertain species (ones that are spp, sp1, sp2 or unknown)
tot-unk # number of certain species

# Number of spp unique to one or the other sampling

preferential_species <- dat %>%
  filter(sd == "preferential") %>%
  pull(fullspp) %>%
  unique()

# Get species in "spatially balanced"
spatially_balanced_species <- dat %>%
  filter(sd == "spatially balanced") %>%
  pull(fullspp) %>%
  unique()

# Find species that are only in "preferential"
only_preferential <- setdiff(preferential_species, spatially_balanced_species)

# Find species that are only in "spatially balanced"
only_spatially_balanced <- setdiff(spatially_balanced_species, preferential_species)

# Get the count of unique species in each group
length(only_preferential)
length(only_spatially_balanced)


# Number of species in the community composition box
crs <- st_crs(4326)
box <- st_read("QGIS layers/polys/com_comp_analysis_limits.shp") %>%
  st_make_valid() %>%
  st_transform(crs)
plot(box)
dat2 <- st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326); plot(dat2)

dat2 <- st_intersection(dat2, box); plot(dat2)




# Get unique species for all sd groups
all_species <- dat2 %>%
  pull(fullspp) %>%
  unique()

# Get unique species for preferential sd
preferential_species <- dat2 %>%
  filter(sd == "preferential") %>%
  pull(fullspp) %>%
  unique()

# Get unique species for spatially balanced sd
spatially_balanced_species <- dat2 %>%
  filter(sd == "spatially balanced") %>%
  pull(fullspp) %>%
  unique()

# Find species that are only in preferential
only_preferential <- setdiff(preferential_species, spatially_balanced_species)

# Find species that are only in spatially balanced
only_spatially_balanced <- setdiff(spatially_balanced_species, preferential_species)

# Create the summary table
summary_table <- tibble(
  Total_Species_All_SD = length(all_species),
  Total_Species_Preferential = length(preferential_species),
  Total_Species_Spatially_Balanced = length(spatially_balanced_species),
  Species_Only_Preferential = length(only_preferential),
  Species_Only_Spatially_Balanced = length(only_spatially_balanced)
)

# Print the summary table
print(summary_table)

