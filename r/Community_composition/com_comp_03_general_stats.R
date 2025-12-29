# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison
# Data:    2024 BRUV MaxN data. (MEGLAB + DBCA)
# Task:    Get general numbers for results section
# Author:  Lise Fournier-Carnoy
# Date:    November 2025

# -----------------------------------------------------------------------------

rm(list=ls())

library(tidyverse) # general dat manipulation
library(sf) # deal with spatial objects
library(gt) # pretty tables for article

colour_palette <- readLines("chapter_colours.txt"); colour_palette <- eval(parse(text = colour_palette))
source("custom_theme.R")

# WAATERN ---------------------------------------------------------------------

study_site <- "waatern"
dat <- readRDS(paste0("data/tidy/C01_2024_", study_site, "_all_tidy_maxn.rds")) %>%
  glimpse()

# percentage of habitat observed across each design
crs <- st_crs(4326)
box <- st_read(paste0("QGIS layers/clean/", study_site, "_community_composition_analysis_limits.shp")) %>%
  st_make_valid() %>%
  st_transform(crs)
plot(box$geometry)
dat2 <- st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326) %>% st_set_crs(4326); plot(dat2$geometry)
st_crs(dat2)
dat2 <- st_intersection(dat2, box); plot(dat2$geometry)

hab_waatern <- readRDS("data/tidy/D01_waatern_tidy_habitat.rds") %>%
  dplyr::filter(opcode %in% dat2$opcode) %>%
  dplyr::select(c(opcode, total_pts, sand, seagrass, reef)) %>%
  mutate(sd = ifelse(grepl("P", opcode), "Targeted", "Stratified spatially balanced (assemblage analysis dataset)")) %>% 
  glimpse()


hab_waatern_tot <- hab_waatern %>% 
  group_by(sd) %>%
  summarise(
    total_pts = sum(total_pts, na.rm = TRUE),
    sand = sum(sand, na.rm = TRUE),
    seagrass = sum(seagrass, na.rm = TRUE),
    reef = sum(reef, na.rm = TRUE)
  ) %>%
  mutate(
    sand_pct = round(100 * sand / total_pts, 1),
    seagrass_pct = round(100 * seagrass / total_pts, 1),
    reef_pct = round(100 * reef / total_pts, 1)
  )

print(hab_waatern_tot)


# WAATU -----------------------------------------------------------------------
study_site <- "waatu"

dat <- readRDS(paste0("data/tidy/C01_2024_", study_site, "_all_tidy_maxn.rds")) %>%
  glimpse()

# percentage of habitat observed across each design
crs <- st_crs(4326)
box <- st_read(paste0("QGIS layers/clean/", study_site, "_community_composition_analysis_limits.shp")) %>%
  st_make_valid() %>%
  st_transform(crs)
plot(box$geometry)
dat2 <- st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326) %>% st_set_crs(4326); plot(dat2$geometry)
dat2 <- st_intersection(dat2, box); plot(dat2$geometry)

hab_waatu <- readRDS("data/tidy/D01_waatu_tidy_habitat.rds") %>%
  dplyr::filter(opcode %in% dat2$opcode) %>%
  dplyr::select(c(opcode, total_pts, sand, seagrass, reef)) %>%
  mutate(sd = ifelse(grepl("P", opcode), "Targeted", "Stratified spatially balanced (assemblage analysis dataset)")) %>% 
  st_drop_geometry() %>% 
  glimpse()

hab_waatu2 <- read.csv("data/raw/2024_waatu_state/2024-02_NgariCapes.MP.Monitoring_stereoBRUVs_benthos.csv", skip = 4) %>% 
  dplyr::select(opcode = Filename, hab_1 = X, hab_2 = X.1) %>% 
  dplyr::filter(grepl("CF", opcode)) %>% 
  mutate(hab = case_when(
    hab_2 == "Unconsolidated (soft)" ~ "sand",
    hab_2 == "Unconsolidated (hard)" ~ "reef",
    hab_1 == "Seagrasses" ~ "seagrass",
    hab_1 == "Unscorable" ~ NA,
    hab_1 == "Fishes" ~ NA,
    is.na(hab_1) ~ NA,
    TRUE ~ 'reef'
  )) %>% 
  count(opcode, hab, name = "n_hab") %>%
  pivot_wider(
    names_from = hab,
    values_from = n_hab,
    names_prefix = "",
    values_fill = 0
  ) %>% 
  mutate(sd = "Targeted",
         total_pts = sand + reef) %>% 
  glimpse()



hab_waatu_tot <- rbind(hab_waatu, hab_waatu2) %>% 
  group_by(sd) %>%
  summarise(
    total_pts = sum(total_pts, na.rm = TRUE),
    sand = sum(sand, na.rm = TRUE),
    seagrass = sum(seagrass, na.rm = TRUE),
    reef = sum(reef, na.rm = TRUE)
  ) %>%
  mutate(
    sand_pct = round(100 * sand / total_pts, 1),
    seagrass_pct = round(100 * seagrass / total_pts, 1),
    reef_pct = round(100 * reef / total_pts, 1)
  )

print(hab_waatu_tot)


## rbind together -------------------------------------------------------------

table <- rbind(hab_waatu_tot %>% dplyr::select(sd, sand_pct, seagrass_pct, reef_pct) %>% mutate(site = "Waatu/Capes Region"),
      st_drop_geometry(hab_waatern_tot) %>% dplyr::select(sd, sand_pct, seagrass_pct, reef_pct) %>% mutate(site = "Waatern/Geographe Bay"))

table %>%
  gt(groupname_col = "site") %>%
  cols_label(
    sand_pct = "Sand (%)",
    sd = "",
    seagrass_pct = "Seagrass (%)",
    reef_pct = "Reef (%)"
  ) %>%
  fmt_number(columns = where(is.integer), decimals = 0) %>%
  tab_style(
    style = cell_text(weight = "bold",
                      color = colour_palette[2]),
    locations = cells_row_groups()
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold",
                      color = colour_palette[2]),
    locations = cells_column_labels()
  ) %>% 
  tab_source_note(
    source_note = "% refers to the benthos cover across all samples.\n
    The 'assemblage analysis dataset' refers to data defined by the red polygon in Figure 1D-E, used in the assemblage analysis."
  ) 


## Simulated areas ------------------------------------------------------------

# calculating areas,  Waatern
land <- st_read(paste0("QGIS layers/clean/wadandi_land_highres.shp")) %>% st_transform(4326)
plot(land$geometry)
st_crs(box)
study_site <- "waatern"
box <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ_sampling_area.shp"))
box <- rmapshaper::ms_erase(box, land)
plot(box$geometry)
st_area(st_union(box))/1000000 # area of the simulated sampling area

SZ <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ.shp"))
SZ <- rmapshaper::ms_erase(SZ, land)

plot(SZ$geometry, add = T)
st_area(SZ)/1000000 # area of the simulated SZ


# calculating areas,  Waatu
land <- st_read(paste0("QGIS layers/clean/wadandi_land_highres.shp")) %>% st_transform(4326)
plot(land$geometry)
st_crs(box)
study_site <- "waatu"
box <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ_sampling_area.shp"))
box <- rmapshaper::ms_erase(box, land)
plot(box$geometry)
st_area(st_union(box))/1000000 # area of the simulated sampling area

SZ <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ.shp"))
SZ <- rmapshaper::ms_erase(SZ, land)

plot(SZ$geometry, add = T)
st_area(SZ)/1000000 # area of the simulated SZ

### END ###
