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

colour_palette <- readLines("chapter_colours.txt"); colour_palette <- eval(parse(text = colour_palette))
source("custom_theme.R")

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


## Habitat sampled differences -------------------------------------------------

hab_pref <- readRDS("data/tidy/2024_geographe_pref_tidy_habitat.rds") %>%
  dplyr::filter(opcode %in% dat2$opcode) %>%
  mutate(Reef = Macroalgae + Consolidated + Invertebrate + Stony.corals) %>%
  dplyr::select(c(opcode, total_pts, Unconsolidated, Seagrasses, Reef)) %>%
  glimpse()


hab_spabal <- readRDS("data/tidy/2024_geographe_tidy_habitat.rds") %>%
  # filter out unscorables
  filter(!(level_2 %in% c("Fishes", "Echinoderms"))) %>% 
  
  # make the correct categories
  mutate(taxa = case_when(level_2 %in% "Macroalgae" ~ "macro",
                          level_3 %in% "Unconsolidated (soft)" ~ "sand",
                          level_3 %in% "Consolidated (hard)" ~ "rock",
                          level_2 %in% "Seagrasses" ~ "seagrass",
                          level_2 %in% c("Sessile invertebrates", "Sponges", "Bryozoa", "Cnidaria") ~ "sessile_inverts")) %>%
  pivot_wider(names_from = taxa, values_from = number, values_fill = list(number = 0)) %>% 
  
  # make reef and remove the individual reef habitats
  mutate(reef = macro + rock + sessile_inverts) %>%
  dplyr::select(-c(sessile_inverts, macro, rock)) %>% 
  
  # sum rows to have a single row per opcode
  group_by(opcode, latitude_dd, longitude_dd) %>%  # Group by the `opcode`
  summarize(
    sand = sum(sand, na.rm = TRUE),
    seagrass = sum(seagrass, na.rm = TRUE),
    reef = sum(reef, na.rm = TRUE)
  ) %>%
  
  # add total points - these should add up to all the reef, sand, seagrass
  mutate(total_pts = sum(c(reef, sand, seagrass), na.rm = TRUE)) %>% 
  
  # Transform to a spatial object to extract bathy values
  st_as_sf(coords = c("longitude_dd", "latitude_dd"), crs = 4326) %>%  # Replace with the correct column names
  dplyr::filter(opcode %in% dat2$opcode) %>%
  glimpse()

# Calculate the total proportions for `hab_pref` across the entire dataset
hab_pref_proportions_total <- hab_pref %>%
  ungroup() %>%  # Remove any existing grouping
  summarise(
    total_sand = sum(Unconsolidated, na.rm = TRUE),
    total_seagrass = sum(Seagrasses, na.rm = TRUE),
    total_reef = sum(Reef, na.rm = TRUE),
    total_pts = sum(total_pts, na.rm = TRUE)
  ) %>%
  mutate(
    sand_prop = total_sand / total_pts,
    seagrass_prop = total_seagrass / total_pts,
    reef_prop = total_reef / total_pts
  )

# Calculate the total proportions for `hab_spabal` across the entire dataset
hab_spabal_proportions_total <- hab_spabal %>%
  ungroup() %>%  # Remove any existing grouping
  summarise(
    total_sand = sum(sand, na.rm = TRUE),
    total_seagrass = sum(seagrass, na.rm = TRUE),
    total_reef = sum(reef, na.rm = TRUE),
    total_pts = sum(total_pts, na.rm = TRUE)
  ) %>%
  mutate(
    sand_prop = total_sand / total_pts,
    seagrass_prop = total_seagrass / total_pts,
    reef_prop = total_reef / total_pts
  )

# View the results
hab_pref_proportions_total
hab_spabal_proportions_total

# Tidy the hab_spabal proportions by explicitly selecting only the numeric columns
hab_spabal_tidy <- hab_spabal_proportions_total %>%
  st_drop_geometry() %>%  # This will remove the geometry column
  select(sand_prop, seagrass_prop, reef_prop) %>%  # Explicitly select only the numeric columns
  #pivot_longer(cols = everything(), names_to = "habitat", values_to = "proportion") %>% 
  mutate(dataset = "Spatially Balanced")  # Add a column for the dataset name

# Tidy the hab_pref proportions (No changes needed here)
hab_pref_tidy <- hab_pref_proportions_total %>% 
  select(sand_prop, seagrass_prop, reef_prop) %>% 
  #pivot_longer(cols = everything(), names_to = "habitat", values_to = "proportion") %>% 
  mutate(dataset = "Clustered")  # Add a column for the dataset name

# Combine both datasets
tidy_proportions <- bind_rows(hab_pref_tidy, hab_spabal_tidy) %>%
  rename(Sand = "sand_prop",
         Seagrass = "seagrass_prop",
         Reef = "reef_prop",
         `Sampling Design` = "dataset")

# View the tidied data
tidy_proportions

saveRDS(tidy_proportions, "data/rmd/general_stats_habitats_sampled_proportion.rds")



### END ###