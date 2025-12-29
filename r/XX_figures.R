# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare abundance detected)
# Data:    everything made in previous scripts
# Task:    Make figures
# Author:  Lise Fournier-Carnoy, heavily inspired by Claude Spencer
# Date:    November 2025

# -----------------------------------------------------------------------------

# Status: the entirety of the script (for 1 site) should yield 15 files 
#         (all in outputs/distribution_modelling_outputs/D04_sampling_designs/)

# -----------------------------------------------------------------------------

rm(list=ls())

colour_palette <- readLines("chapter_colours.txt"); colour_palette <- eval(parse(text = colour_palette))
source("custom_theme.R")

library(tidyverse)
library(patchwork)
library(sf)
library(terra)
library(ggpubr)
library(cowplot)
library(imager)
library(grid) # for images
library(magick) # for images
library(gt)

## Figure 1. the study area ---------------------------------------------------


# load data from both study sites
land <- st_read("QGIS layers/clean/wadandi_land_highres.shp"); plot(land$geometry)
sz_state <- st_read("QGIS layers/clean/dbca_wadandi_sz.shp"); plot(sz_state$geometry)
sz_cw <- st_read("QGIS layers/clean/wadandi_commonwealth_marine_parks.shp"); sz_cw <- sz_cw[sz_cw$ZoneName == "National Park Zone",]; plot(sz_cw$geometry)
waatern_com_zone <- st_read("QGIS layers/clean/waatern_community_composition_analysis_limits.shp"); plot(waatern_com_zone$geometry)
waatu_com_zone <- st_read("QGIS layers/clean/waatu_community_composition_analysis_limits.shp"); plot(waatu_com_zone$geometry)

waatern_clustered_maxn <- st_read("QGIS layers/clean/waatern_2024_state_BRUVs.shp")
waatu_clustered_maxn <- st_read("QGIS layers/clean/waatu_2024_state_BRUVs.shp")

waatern_spabal_maxn <- st_read("QGIS layers/clean/waatern_2024_commonwealth_BRUVs.shp")
waatu_spabal_maxn <- st_read("QGIS layers/clean/waatu_2024_commonwealth_BRUVs.shp")

waatern_bbox <- st_bbox(c(xmin = 115.05, xmax = 115.54, ymax = -33.67, ymin = -33.35), crs = st_crs(4326))
waatern <- ggplot() +
  geom_sf(data = waatern_com_zone, color = colour_palette[6], alpha = 0) +
  geom_sf(data = land, fill = "grey90", color = "grey", size = 0.2) +
  geom_sf(data = sz_state, color = colour_palette[1], fill = colour_palette[2], alpha = 0.3) +
  geom_sf(data = sz_cw, color = colour_palette[1], fill = colour_palette[2], alpha = 0.3) +
  #  geom_sf(data = waatern_clustered_maxn, color =  colour_palette[4]) +
  geom_sf(data = waatern_spabal_maxn, size = 0.5, color = ifelse(waatern_spabal_maxn$SD == "spabal", colour_palette[6], colour_palette[4])) +
  theme_minimal() +
  coord_sf(
    xlim = waatern_bbox[c(1, 3)], ylim = waatern_bbox[c(2, 4)],   # Set based on your desired latitude range
    expand = FALSE
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top"
  ) + 
  ggspatial::annotation_scale(location = 'br') +
  labs(x = "Waatern/Geographe Bay")
waatern


waatu_clustered_maxn$type <- "Targeted"
waatu_spabal_maxn$type <- "Stratified spatially balanced"
waatu_bbox <- st_bbox(c(xmin = 114.7, xmax = 115.05, ymax = -34.24, ymin = -33.95), crs = st_crs(4326))
waatu <- ggplot() +
  geom_sf(data = waatu_com_zone, color = colour_palette[6], alpha = 0) +
  geom_sf(data = land, fill = "grey90", color = "grey", size = 0.2) +
  geom_sf(data = sz_state, color = colour_palette[1], fill = colour_palette[2], alpha = 0.3) +
  geom_sf(data = sz_cw, color = colour_palette[1], fill = colour_palette[2], alpha = 0.3) +
  geom_sf(data = waatu_clustered_maxn, size = 0.5, aes(color = type)) +
  geom_sf(data = waatu_spabal_maxn, size = 0.5, aes(color = type)) +
  coord_sf(
    xlim = waatu_bbox[c(1, 3)],
    ylim = waatu_bbox[c(2, 4)],
    expand = FALSE
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "Targeted" = colour_palette[4],
      "Stratified spatially balanced" = colour_palette[6]
    )
  ) +
  
  # Clean theme and legend placement
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top"
  ) + 
  ggspatial::annotation_scale(location = 'br') +
  labs(x = "Waatu/Capes Region")
waatu

aus <- st_read("QGIS layers/clean/AUS_2021_AUST_GDA2020.shp")

aus_plot <- ggplot() +
  geom_sf(data = aus,
          fill = "grey90", color = "grey") +
  geom_sf(data = st_as_sfc(st_bbox(c(xmin = 114.8, xmax = 115.6, ymax = -34.4, ymin = -33.4), crs = st_crs(4326))),
          fill = NA, colour = colour_palette[5]) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top"
  ) +
  coord_sf(
    xlim = c(110, 155),
    ylim = c(-45, -10),
    expand = FALSE
  )
aus_plot

overview <- ggplot() +
  geom_sf(data = waatu_com_zone,
          aes(color = 'Assemblage analysis limits'),
          fill = NA, size = 0.5) +
  geom_sf(data = waatern_com_zone,
          aes(color = 'Assemblage analysis limits'),
          fill = NA, size = 0.5) +
  geom_sf(data = land,
          fill = "grey90", color = "grey", size = 0.2) +
  geom_sf(data = sz_state,
          aes(fill = 'No-take zones'),
          color = colour_palette[2], alpha = 0.1) +
  
  geom_sf(data = st_as_sfc(waatern_bbox),
          color = colour_palette[5], fill = NA) +
  geom_sf(data = st_as_sfc(waatu_bbox),
          color = colour_palette[5], fill = NA) +
  
  geom_sf(data = sz_cw,
          aes(fill = 'No-take zones'),
          color = colour_palette[2], alpha = 0.1) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(
    xlim = c(114.8, 115.6),
    ylim = c(-34.4, -33.4),
    expand = FALSE
  ) +
  scale_fill_manual(name = "", values = c('No-take zones' = colour_palette[2])) +
  scale_color_manual(name = "", values = c('Assemblage analysis limits' = colour_palette[6]))
overview


image1 <- rasterGrob(load.image("photos/wadandi_ranger_BRUV.JPG"), interpolate = TRUE)
image2 <- rasterGrob(load.image("photos/BRUV_G.hebraicum_C.auratus_Waatu.PNG"), interpolate = TRUE)
image_plot1 <- ggplot() + annotation_custom(image1) + theme_void()
image_plot2 <- ggplot() + annotation_custom(image2) + theme_void() 

images <- ggarrange(
  ggarrange(image_plot1, 
            image_plot2, 
            ncol = 1),
  widths = c(1,1),
  nrow = 1, 
  heights = c(1,1)
)
images

final_plot <- (images | (overview | aus_plot) / (waatu |  waatern)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
final_plot

ggsave("plots/fig1.png", plot = final_plot, width = 8, height = 5, dpi = 1500)



## Table 1. summary of deployment numbers -------------------------------------

waatern_clustered_maxn <- st_read("QGIS layers/clean/waatern_2024_state_BRUVs.shp")
waatu_clustered_maxn <- st_read("QGIS layers/clean/waatu_2024_state_BRUVs.shp")

waatern_spabal_maxn <- st_read("QGIS layers/clean/waatern_2024_commonwealth_BRUVs.shp")
waatu_spabal_maxn <- st_read("QGIS layers/clean/waatu_2024_commonwealth_BRUVs.shp")

waatern_com_zone <- st_read("QGIS layers/clean/waatern_community_composition_analysis_limits.shp"); plot(waatern_com_zone$geometry)
waatu_com_zone <- st_read("QGIS layers/clean/waatu_community_composition_analysis_limits.shp"); plot(waatu_com_zone$geometry)
length(st_intersection(waatu_spabal_maxn, waatu_com_zone)$opcode)

sample_summary <- data.frame(
  Region = c('A. Waatu/Capes Region', 'A. Waatu/Capes Region', 'A. Waatu/Capes Region',
             'B. Waatern/Geographe Bay', 'B. Waatern/Geographe Bay', 'B. Waatern/Geographe Bay'),
  Design = c('Stratified spatially balanced sampling (distribution model dataset)', 'Stratified spatially balanced sampling (assemblage analysis dataset)', 'Targeted', 
             'Stratified spatially balanced sampling (distribution model dataset)', 'Stratified spatially balanced sampling (assemblage analysis dataset)', 'Targeted'),
  Total_Samples = as.integer(c(
    length(waatu_spabal_maxn$opcode),
    length(st_intersection(waatu_spabal_maxn, waatu_com_zone)$opcode), #subset
    length(waatu_clustered_maxn$opcode),
    
    length(waatern_spabal_maxn$opcode),
    length(st_intersection(waatern_spabal_maxn, waatern_com_zone)$opcode), #subset
    length(waatern_clustered_maxn$opcode)
  )),
  NTZ_Samples = as.integer(c(
    sum(waatu_spabal_maxn$status == "No-take", na.rm = TRUE),
    sum((st_intersection(waatu_spabal_maxn, waatu_com_zone))$status == "No-take", na.rm = TRUE),
    sum(waatu_clustered_maxn$status == "No-take", na.rm = TRUE),
    
    sum(waatern_spabal_maxn$status == "No-take", na.rm = TRUE),
    sum((st_intersection(waatern_spabal_maxn, waatern_com_zone))$status == "No-take", na.rm = TRUE),
    sum(waatern_clustered_maxn$status == "No-take", na.rm = TRUE)
  )),
  Fished_Samples = as.integer(c(
    sum(waatu_spabal_maxn$status == "Fished", na.rm = TRUE),
    sum((st_intersection(waatu_spabal_maxn, waatu_com_zone))$status == "Fished", na.rm = TRUE),
    sum(waatu_clustered_maxn$status == "Fished", na.rm = TRUE),
    
    sum(waatern_spabal_maxn$status == "Fished", na.rm = TRUE),
    sum((st_intersection(waatern_spabal_maxn, waatern_com_zone))$status == "Fished", na.rm = TRUE),
    sum(waatern_clustered_maxn$status == "Fished", na.rm = TRUE)
  ))
)
sample_summary <- sample_summary[, c("Region", "Design", "NTZ_Samples", "Fished_Samples", "Total_Samples")]

# View table
print(sample_summary)
library(gt)

sample_summary %>%
  gt(groupname_col = "Region") %>%
  cols_label(
    Region = "Region",
    Design = "",
    Total_Samples = "Total",
    NTZ_Samples = "NTZ",
    Fished_Samples = "Fished"
  ) %>%
  fmt_number(columns = where(is.integer), decimals = 0) %>%
  tab_style(
    style = cell_text(weight = "bold", color = colour_palette[2]),
    locations = cells_row_groups()
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold", color = colour_palette[2]),
    locations = cells_column_labels()
  ) %>% 
  tab_source_note(
    source_note = "The 'assemblage analysis dataset' refers to stratified spatially balanced samples located in the red polygons in Figure 1, used for assemblage comparison analysis \n
    The 'distribution model dataset' refers to all stratified spatially balanced data in Figure 1, used for species distribution modelling"
  ) 


## Table 2. observed habitat -------------------------------------------------
# WAATERN

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


# WAATU
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


# rbind together
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


## Figure A2. simulated number of samples ---------------------------------------

waatern_forced <- tibble(
  ` `                      = c("Simple Spatially Balanced", "Stratified Spatially Balanced (25/25)", "Stratified Spatially Balanced (20/30)", "Preferential (25/25)", "Preferential (20/30)", "Clustered"),
  
  `in NTZ (strata 1)`      = c(NA, 7, 5, 25, 20, NA),
  `in NTZ (strata 2)`      = c(NA, 8, 7, 0,  0,  NA),
  `in NTZ (strata 3)`      = c(NA, 10, 8, 0,  0,  NA),
  
  `outside NTZ (strata 1)` = c(NA, 7, 8, 25, 30, NA),
  `outside NTZ (strata 2)` = c(NA, 8, 10, 0,  0,  NA),
  `outside NTZ (strata 3)` = c(NA, 10, 12, 0,  0,  NA),
  
  `Total samples`          = rep(50, 6),
  
  `Region`               = rep("A. Waatern/Geographe Bay", 6)
)

waatu_forced <- tibble(
  ` `                      = c("Simple Spatially Balanced", "Stratified Spatially Balanced (25/25)", "Stratified Spatially Balanced (20/30)", "Preferential (25/25)", "Preferential (20/30)", "Clustered"),
  `in NTZ (strata 1)`      = c(NA, 5, 3, 18, 14, NA),
  `in NTZ (strata 2)`      = c(NA, 6, 5, 0,  0,  NA),
  `in NTZ (strata 3)`      = c(NA, 7, 6, 0,  0,  NA),
  
  `outside NTZ (strata 1)` = c(NA, 5, 6, 18, 22, NA),
  `outside NTZ (strata 2)` = c(NA, 6, 7, 0,  0,  NA),
  `outside NTZ (strata 3)` = c(NA, 7, 9, 0,  0,  NA),
  
  `Total samples`          = rep(36, 6),
  
  `Region`               = rep("B. Waatu/Capes Region", 6)
)

sample_n_forced <- rbind(waatern_forced, waatu_forced)

sample_n_forced %>% 
  gt(groupname_col = "Region") %>% 
  cols_label(
    ` ` = "",
    `outside NTZ (strata 1)` = "strata 1",
    `outside NTZ (strata 2)` = "strata 2",
    `outside NTZ (strata 3)` = "strata 3",
    `in NTZ (strata 1)` = "strata 1",
    `in NTZ (strata 2)` = "strata 2",
    `in NTZ (strata 3)` = "strata 3",
    `Total samples` = "Total samples"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold", color = colour_palette[2]),
    locations = cells_row_groups()
  ) %>% 
  tab_spanner(
    label = "Samples forced in NTZ",
    columns = c(2:4)
  ) %>% 
  tab_spanner(
    label = "Samples forced outside NTZ",
    columns = c(5:7)
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold", 
                      size = px(20)),
    locations = cells_column_spanners()
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold", color = colour_palette[2]),
    locations = list(cells_column_labels(), cells_column_spanners())
  ) %>% 
  tab_style_body(
    fn = function(x) is.na(x),
    style = cell_text(color = "gray")
  ) %>% 
  tab_style_body(
    fn = function(x) x==0,
    style = cell_text(color = "gray")
  ) %>% 
  tab_source_note(
    source_note = "Numbers here are the target numbers, although the design generation process may remove a point or more if the minimum distance between points (250m) is violated."
  ) 


## Figure A3. realised abundance ratios ---------------------------------------

realised <- tibble(
  ` `                        = c(rep("Waatern/Geographe Bay", 3), rep("Waatu/Capes Region", 3)),
  
  `Species`                  = c(rep(c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"), 2)),
  `Realised abundance (inside NTZ)`   = c(3.29, 2.65, 0.14, 3.12, 16.1, 0.9),
  `Realised abundance (outside NTZ)`  = c(1.99, 1.53, 0.1,  1.7,  9.93, 0.47),
  `Realised abundance ratio`          = c(1.64, 1.73, 1.37, 1.83, 1.62, 1.91),
)

realised %>% 
  gt(groupname_col = "Species") %>% 
  tab_style(
    style = cell_text(weight = "bold", color = colour_palette[2]),
    locations = cells_row_groups()
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold", 
                      size = px(20)),
    locations = cells_column_spanners()
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold", color = colour_palette[2]),
    locations = list(cells_column_labels())
  ) %>% 
  tab_style_body(
    fn = function(x) is.na(x),
    style = cell_text(color = "gray")
  ) %>% 
  tab_style_body(
    fn = function(x) x==0,
    style = cell_text(color = "gray")
  ) %>% 
  tab_source_note(
    source_note = "The realised abundance is the mean predicted abundance of the entire area (inside or outside NTZ) within which sampling designs were simulated."
  ) 

## Figure 3. Example sampling designs -----------------------------------------

land <- st_read("QGIS layers/clean/wadandi_land_highres.shp"); plot(land$geometry)

create_sampling_design_plot <- function(sf_data, sampling_area, SZ, aus, samples_sf, colour_palette, legend = TRUE) {
  p <- ggplot() +
    geom_sf(data = sf_data, aes(fill = as.factor(detrended)), colour = NA) +
    scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#E0E0E0", "2" = "#9E9E9E"),
                      na.value = NA,
                      labels = c("0" = "Strata 1\n(low)",
                                 "1" = "Strata 2\n(medium)",
                                 "2" = "Strata 3\n(high)")) +
    labs(fill = "Detrended\nbathymetry") +
    geom_sf(data = SZ, aes(colour = "Simulated\nNTZ"), fill = NA, linewidth = 1) +
    geom_sf(data = land, color = "darkgray", fill = 'lightgray', size = 0.2) +
    geom_sf(data = samples_sf, aes(colour = "Simulated\nsample points"), size = 1) +
    coord_sf(crs = 4326, xlim = ext(sampling_area)[1:2], ylim = ext(sampling_area)[3:4]) +
    scale_color_manual(name = '', values = c('Simulated\nsample points' = colour_palette[6], 
                                             'Simulated\nNTZ' = colour_palette[2])) +
    custom_theme +
    theme(
      axis.text.x = element_blank(),  # Explicitly remove x-axis text
      axis.text.y = element_blank(),  # Explicitly remove y-axis text
      axis.ticks = element_blank(),   # Remove axis ticks
      axis.title.x = element_blank(), # Remove x-axis title
      axis.title.y = element_blank(),  # Remove y-axis title
      legend.position = "top",
      plot.margin = margin(0, 0, 0, 0, "cm")
    )
}


# WAATERN example designs

study_site <- "waatern"
target_crs <- 4326

SZ                  <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ.shp")) %>% st_transform(crs=target_crs) # simulated SZ shapefile
sampling_area_waatern <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ_sampling_area.shp")) %>% st_transform(crs=target_crs) # simulation area

sf_all              <- readRDS(paste0("outputs/distribution_modelling_outputs/D04_", study_site, "_sample_area_detrended_strata.rds")) %>% st_transform(crs=target_crs) # detrended bathymetry strata in the simulation area

samples_sf_rand     <- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site, '_example_sampling_design_simple_spatial_balance.rds')) %>% st_transform(crs=target_crs)# example points for simulated random design
samples_sf_pref_2525<- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site, '_example_sampling_design_preferential_25_25.rds')) %>% st_transform(crs=target_crs) # example points for simulated preferential design
samples_sf_pref_2030<- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site,'_example_sampling_design_preferential_20_30.rds')) %>% st_transform(crs=target_crs) # example points for simulated preferential design
samples_sf_sb_2525  <- readRDS(paste0("outputs/distribution_modelling_outputs/D04_sampling_designs/D04_", study_site, "_example_sampling_design_stratified_spatial_balance_25_25.rds")) %>% st_transform(crs=target_crs) # example points for simulated spatially balanced design
samples_sf_sb_2030  <- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site, '_example_sampling_design_stratified_spatial_balance_20_30.rds')) %>% st_transform(crs=target_crs) # example points for simulated spatially balanced design
samples_sf_clump    <- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site, '_example_sampling_design_clustered.rds')) %>% st_transform(crs=target_crs) # example points for simulated clustered design

# make the subplots
waatern_ex_rand <- create_sampling_design_plot(sf_all, sampling_area_waatern, SZ, aus, samples_sf_rand, colour_palette, legend = TRUE)
waatern_ex_spabal_2525 <- create_sampling_design_plot(sf_all, sampling_area_waatern, SZ, aus, samples_sf_sb_2525, colour_palette, legend = TRUE)
waatern_ex_spabal_2030 <- create_sampling_design_plot(sf_all, sampling_area_waatern, SZ, aus, samples_sf_sb_2030, colour_palette, legend = TRUE)
waatern_ex_pref_2525 <- create_sampling_design_plot(sf_all, sampling_area_waatern, SZ, aus, samples_sf_pref_2525, colour_palette, legend = TRUE)
waatern_ex_pref_2030 <- create_sampling_design_plot(sf_all, sampling_area_waatern, SZ, aus, samples_sf_pref_2030, colour_palette, legend = TRUE)
waatern_ex_clump <- create_sampling_design_plot(sf_all, sampling_area_waatern, SZ, aus, samples_sf_clump, colour_palette, legend = FALSE)


# WAATU example designs

# Load layers
study_site <- "waatu"
target_crs <- 4326

SZ                  <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ.shp")) %>% st_transform(crs=target_crs) # simulated SZ shapefile
sampling_area_waatu <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ_sampling_area.shp")) %>% st_transform(crs=target_crs) # simulation area

sf_all              <- readRDS(paste0("outputs/distribution_modelling_outputs/D04_", study_site, "_sample_area_detrended_strata.rds")) %>% st_transform(crs=target_crs) # detrended bathymetry strata in the simulation area

samples_sf_rand     <- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site, '_example_sampling_design_simple_spatial_balance.rds')) %>% st_transform(crs=target_crs)# example points for simulated random design
samples_sf_pref_2525<- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site, '_example_sampling_design_preferential_25_25.rds')) %>% st_transform(crs=target_crs) # example points for simulated preferential design
samples_sf_pref_2030<- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site,'_example_sampling_design_preferential_20_30.rds')) %>% st_transform(crs=target_crs) # example points for simulated preferential design
samples_sf_sb_2525  <- readRDS(paste0("outputs/distribution_modelling_outputs/D04_sampling_designs/D04_", study_site, "_example_sampling_design_stratified_spatial_balance_25_25.rds")) %>% st_transform(crs=target_crs) # example points for simulated spatially balanced design
samples_sf_sb_2030  <- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site, '_example_sampling_design_stratified_spatial_balance_20_30.rds')) %>% st_transform(crs=target_crs) # example points for simulated spatially balanced design
samples_sf_clump    <- readRDS(paste0('outputs/distribution_modelling_outputs/D04_sampling_designs/D04_', study_site, '_example_sampling_design_clustered.rds')) %>% st_transform(crs=target_crs) # example points for simulated clustered design

pred_PS      <- as.data.frame(readRDS(paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_Sparidae_Chrysophrys_auratus.rds")), xy=TRUE)
pred_SMW     <- as.data.frame(readRDS(paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_Labridae_Ophthalmolepis_lineolatus.rds")), xy=TRUE)
pred_WAD     <- as.data.frame(readRDS(paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_Glaucosomatidae_Glaucosoma_hebraicum.rds")), xy=TRUE)

# make the subplots
waatu_ex_rand <- create_sampling_design_plot(sf_all, sampling_area_waatu, SZ, aus, samples_sf_rand, colour_palette, legend = TRUE)
waatu_ex_spabal_2525 <- create_sampling_design_plot(sf_all, sampling_area_waatu, SZ, aus, samples_sf_sb_2525, colour_palette, legend = TRUE)
waatu_ex_spabal_2030 <- create_sampling_design_plot(sf_all, sampling_area_waatu, SZ, aus, samples_sf_sb_2030, colour_palette, legend = TRUE)
waatu_ex_pref_2525 <- create_sampling_design_plot(sf_all, sampling_area_waatu, SZ, aus, samples_sf_pref_2525, colour_palette, legend = TRUE)
waatu_ex_pref_2030 <- create_sampling_design_plot(sf_all, sampling_area_waatu, SZ, aus, samples_sf_pref_2030, colour_palette, legend = TRUE)
waatu_ex_clump <- create_sampling_design_plot(sf_all, sampling_area_waatu, SZ, aus, samples_sf_clump, colour_palette, legend = FALSE)


# combine all example SDs

ex_maps <- ggarrange(
  waatern_ex_rand + theme(legend.position = "top") + annotate("text", x = ext(sampling_area_waatern)[1],y = ext(sampling_area_waatern)[4], label = "simple sp. balanced", hjust = 0, vjust = 0.5, size = 3), 
  waatern_ex_spabal_2030 + theme(legend.position = "top") + annotate("text", x = ext(sampling_area_waatern)[1], y = ext(sampling_area_waatern)[4], label = "str. sp. balanced (20/30)", hjust = 0, vjust = 0.5, size = 3),  
  waatern_ex_spabal_2525 + theme(legend.position = "top") + annotate("text", x = ext(sampling_area_waatern)[1], y = ext(sampling_area_waatern)[4], label = "str. sp. balanced (25/25)", hjust = 0, vjust = 0.5, size = 3),  
  
  waatu_ex_rand + theme(legend.position = "top") + annotate("text", x = ext(sampling_area_waatu)[2], y = ext(sampling_area_waatu)[3], label = "simple sp. balanced", hjust = 0, vjust = -0.5, size = 3, angle = 90), 
  waatu_ex_spabal_2030 + theme(legend.position = "top") + annotate("text", x = ext(sampling_area_waatu)[2], y = ext(sampling_area_waatu)[3], label = "str. sp. balanced (20/30)", hjust = 0, vjust = -0.5, size = 3, angle = 90),  
  waatu_ex_spabal_2525 + theme(legend.position = "top") + annotate("text", x = ext(sampling_area_waatu)[2], y = ext(sampling_area_waatu)[3], label = "str. sp. balanced (25/25)", hjust = 0, vjust = -0.5,size = 3, angle = 90),  
  
  
  waatern_ex_pref_2030 + theme(legend.position = "none") + annotate("text", x = ext(sampling_area_waatern)[1], y = ext(sampling_area_waatern)[4], label = " preferential (20/30)", hjust = 0, vjust = 0.5, size = 3), 
  waatern_ex_pref_2525 + theme(legend.position = "none") + annotate("text", x = ext(sampling_area_waatern)[1], y = ext(sampling_area_waatern)[4],  label = " preferential (25/25)", hjust = 0, vjust = 0.5, size = 3), 
  waatern_ex_clump + theme(legend.position = "none") + annotate("text", x = ext(sampling_area_waatern)[1], y = ext(sampling_area_waatern)[4], label = "clustered", hjust = 0, vjust = 0.5, size = 3), 
  
  waatu_ex_pref_2030 + theme(legend.position = "none") + annotate("text", x = ext(sampling_area_waatu)[2], y = ext(sampling_area_waatu)[3], label = " preferential (20/30)", hjust = 0, vjust = -0.5, size = 3, angle = 90), 
  waatu_ex_pref_2525 + theme(legend.position = "none") + annotate("text", x = ext(sampling_area_waatu)[2], y = ext(sampling_area_waatu)[3], label = " preferential (25/25)", hjust = 0, vjust = -0.5, size = 3, angle = 90), 
  waatu_ex_clump + theme(legend.position = "none") + annotate("text", x = ext(sampling_area_waatu)[2], y = ext(sampling_area_waatu)[3], label = "clustered", hjust = 0, vjust = -0.5, size = 3, angle = 90), 
  
  nrow = 2, ncol = 6,
  widths = c(1, 1, 1),  # Adjust widths for spacing
  labels = c("A.", "B.", "C.", "D.", "E.", "F.",
             "G.", "H.", "I.", "J.", "K.", "L."),
  common.legend = TRUE
)
ex_maps

save_plot(paste0("plots/all_example_sampling_designs.png"), ex_maps, base_width = 10, base_height = 5)


# Define the function for creating the raster plot
create_raster_plot <- function(data, sampling_area, SZ, land, colour_palette) {
  
  bbox <- sf::st_bbox(sampling_area)
  
  # Filter raster data to only include points within the sampling area bounding box
  data_subset <- data[data$x >= bbox["xmin"] & data$x <= bbox["xmax"] &
                        data$y >= bbox["ymin"] & data$y <= bbox["ymax"], ]
  
  # Calculate min and max only from the subset
  fill_min <- min(data_subset$p_fish.fit, na.rm = TRUE)
  fill_max <- max(data_subset$p_fish.fit, na.rm = TRUE)
  
  
  ggplot() +
    geom_raster(data = data, aes(x = x, y = y, fill = p_fish.fit), interpolate = TRUE) +
    geom_sf(data = sampling_area, aes(color = 'Simulation area'),
            fill = 'transparent', linewidth = 1, show.legend = FALSE) +
    geom_sf(data = SZ, aes(color = 'Simulated NTZ'),
            fill = 'transparent', linewidth = 1, show.legend = FALSE) +
    scale_color_manual(name = '',
                       values = c('Simulation area' = "darkgray",
                                  'Simulated NTZ' = colour_palette[2])) +
    geom_sf(data = land, color = "darkgray", fill = 'lightgray', size = 0.2) +
    custom_theme +
    scale_fill_gradient(
      name = "",
      low = "white", high = colour_palette[4],
      limits = c(fill_min, fill_max)  # Use min/max from cropped area
    ) + 
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0, "cm")
    ) +
    coord_sf(crs = 4326, xlim = ext(sampling_area)[1:2], ylim = ext(sampling_area)[3:4])
}

# WAATERN

study_site <- "waatern"
SZ <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ.shp")) %>% st_transform(crs=target_crs) # simulated SZ shapefile

pred_PS      <- as.data.frame(terra::project(readRDS(paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_edge_effect_Sparidae_Chrysophrys_auratus.rds")), "EPSG:4326"), xy=TRUE)
pred_SMW     <- as.data.frame(terra::project(readRDS(paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_edge_effect_Labridae_Ophthalmolepis_lineolatus.rds")), "EPSG:4326"), xy=TRUE)
pred_WAD     <- as.data.frame(terra::project(readRDS(paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_edge_effect_Glaucosomatidae_Glaucosoma_hebraicum.rds")), "EPSG:4326"), xy=TRUE)

PS_raster_waatern <- create_raster_plot(pred_PS, sampling_area_waatern, SZ, land, colour_palette)
SMW_raster_waatern <- create_raster_plot(pred_SMW, sampling_area_waatern, SZ, land, colour_palette)
WAD_raster_waatern <- create_raster_plot(pred_WAD, sampling_area_waatern, SZ, land, colour_palette)


# WAATU
study_site <- "waatu"
SZ <- st_read(paste0("QGIS layers/clean/", study_site, "_simulated_SZ.shp")) %>% st_transform(crs=target_crs) # simulated SZ shapefile

pred_PS      <- as.data.frame(terra::project(readRDS(paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_edge_effect_Sparidae_Chrysophrys_auratus.rds")), "EPSG:4326"), xy=TRUE)
pred_SMW     <- as.data.frame(terra::project(readRDS(paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_edge_effect_Labridae_Ophthalmolepis_lineolatus.rds")), "EPSG:4326"), xy=TRUE)
pred_WAD     <- as.data.frame(terra::project(readRDS(paste0("outputs/distribution_modelling_outputs/D05_SZ_increase_distribution_rasters/D05_", study_site, "_SZ_increase_predicted_abundance_edge_effect_Glaucosomatidae_Glaucosoma_hebraicum.rds")), "EPSG:4326"), xy=TRUE)

PS_raster_waatu <- create_raster_plot(pred_PS, sampling_area_waatu, SZ, land, colour_palette)
SMW_raster_waatu <- create_raster_plot(pred_SMW, sampling_area_waatu, SZ, land, colour_palette)
WAD_raster_waatu <- create_raster_plot(pred_WAD, sampling_area_waatu, SZ, land, colour_palette)


# Arrange raster plots (Bottom Row) with individual legends
raster_plots <- ggarrange(
  PS_raster_waatern + annotate("text", x = ext(sampling_area_waatern)[1], y = ext(sampling_area_waatern)[4], label = "Yijarup/Pink \nSnapper", hjust = -0.1, vjust = 1.1, size = 3), 
  SMW_raster_waatern + annotate("text", x = ext(sampling_area_waatern)[1], y = ext(sampling_area_waatern)[4], label = "Southern \nMaori \nWrasse", hjust = -0.1, vjust = 1.1, size = 3), 
  WAD_raster_waatern + annotate("text", x = ext(sampling_area_waatern)[1], y = ext(sampling_area_waatern)[4], label = "Djubitj/West \nAustralian \nDhufish", hjust = -0.1, vjust = 1.1, size = 3), 
  
  PS_raster_waatu + annotate("text", x = ext(sampling_area_waatu)[2], y = ext(sampling_area_waatu)[3], label = "Yijarup/Pink Snapper", hjust = 0, vjust = -0.5, size = 3, angle = 90), 
  SMW_raster_waatu + annotate("text", x = ext(sampling_area_waatu)[2], y = ext(sampling_area_waatu)[3], label = "Southern Maori Wrasse", hjust = 0, vjust = -0.5, size = 3, angle = 90), 
  WAD_raster_waatu + annotate("text", x = ext(sampling_area_waatu)[2], y = ext(sampling_area_waatu)[3], label = "Djubitj/West Australian Dhufish", hjust = 0, vjust = -0.5, size = 3, angle = 90),
  nrow = 1, widths = c(1, 1, 1),
  labels = c("M.", "N.", "O.", "P.", "Q.", "R."),
  common.legend = FALSE
)
raster_plots
# make a big combined plot with example SDs + rasters

save_plot("plots/species_prediction_plot.png", raster_plots, base_width = 10, base_height = 2)


# full stitched plot

all_plots <- ggarrange(ex_maps, raster_plots, ncol =1, heights = c(2, 1), align = "hv")
all_plots
save_plot("plots/fig3.png", all_plots, base_width = 8, base_height = 7)

## Figure A1. detailed performance comparison ---------------------------------

species_lookup <- data.frame( # lookup table so that different dataframes recognise the species names
  species_formatted = c(
    "Yijarup/Pink Snapper",
    "Djubitj/West Australian Dhufish",
    "Southern Maori Wrasse"
  ),
  species = c(
    "Sparidae_Chrysophrys_auratus",
    "Glaucosomatidae_Glaucosoma_hebraicum",
    "Labridae_Ophthalmolepis_lineolatus"
  )
)

study_site <- "waatu"

raw <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_raw_data_smooth_edge_effect.rds"))%>% 
  left_join(species_lookup, by = "species") %>% 
  glimpse()
raw$species_formatted <- factor(raw$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

sig <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_significance_testing_model_results_smooth_edge_effect.rds")) %>% 
  left_join(species_lookup, by = "species") %>%
  dplyr::select(c(SD, design_id, p_value, estimate, species_formatted)) %>% 
  glimpse()
sig$species_formatted <- factor(sig$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

mb <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_mean_bias_smooth_edge_effect.rds")) %>% 
  left_join(species_lookup, by = "species") %>%
  glimpse()
mb$species_formatted <- factor(mb$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

rmse <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_RMSE_smooth_edge_effect.rds")) %>% 
  left_join(species_lookup, by = "species") %>%
  glimpse()
rmse$species_formatted <- factor(rmse$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

int_cov <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_interval_coverage_smooth_edge_effect.rds")) %>% 
  left_join(species_lookup, by = "species") %>%
  glimpse()
int_cov$species_formatted <- factor(int_cov$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

realised_ab <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_realised_abundances_smooth_edge_effect.rds"))
realised_long <- realised_ab %>%
  dplyr::select(species_raw, #species_formatted = species_formatted.x,
                SZ = mean_abundance_inside, fished = mean_abundance_outside) %>%
  pivot_longer(cols = c(SZ, fished),
               names_to = "status",
               values_to = "realised_abundance_mean") %>% 
  left_join(species_lookup, by = c("species_raw" = "species")) %>% 
  glimpse()
realised_long$species_formatted <- factor(realised_long$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))



gradient_colors <- colorRampPalette(c(colour_palette[4], colour_palette[6]))(6)


### raw data ------------------------------------------------------------------

p_raw <- ggplot(raw, aes(x = SD, y = obs_ab_mean, fill = status)) +
  geom_boxplot(width = 1, color = "black", alpha = 0.5, position = position_dodge(0.9), outlier.size = 0.5) +
  
  # Add horizontal dashed lines for realised abundance
  geom_hline(data = realised_long,
             aes(yintercept = realised_abundance_mean,
                 color = status,
                 group = interaction(species_formatted, status)),
             linetype = "dashed",
             linewidth = 0.8) +
  facet_wrap(~ species_formatted, scales = "free_y") +
  custom_theme +
  labs(
    fill = "Status",
    x = "Sampling Design",
    y = "Observed\nAbundance"
  ) +
  scale_fill_manual(values = colour_palette[c(6, 1)], labels = c("fished", "NTZ")) +
  scale_color_manual(values = colour_palette[c(6, 1)], guide = "none") +  # hide redundant legend
  theme(legend.position = "none",
        axis.text.x = element_blank(),  # Removes only x-axis text
        axis.title.x = element_blank()  # Removes x-axis title
  )

### mean bias -----------------------------------------------------------------

p_mb <- ggplot(mb, aes(x = SD, y = mean_bias, fill = SD)) +
  geom_boxplot(position = position_dodge(0.9), outlier.size = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = colour_palette[6], size = 1) +  
  facet_grid(cols = vars(species_formatted)) +
  labs(y = "Mean Bias") +
  custom_theme +
  theme(
    axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank(),  # Removes x-axis title
    strip.text.x = element_blank()  # Removes column facet titles (species names
  ) +
  scale_fill_manual(values = gradient_colors) +
  scale_color_manual(values = c("white" = colour_palette[6], "black" = "black")) +
  coord_cartesian(ylim = c(NA, 1))

## rmse -----------------------------------------------------------------------

p_rmse <- ggplot(rmse, aes(x = SD, y = rmse, fill = SD)) +
  geom_boxplot(position = position_dodge(0.9), outlier.size = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = colour_palette[6], size = 1) +  
  facet_grid(cols = vars(species_formatted)) +
  labs(y = "RMSE") +
  custom_theme +
  theme(
    axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank(),  # Removes x-axis title
    strip.text.x = element_blank()  # Removes column facet titles (species names
  ) +
  scale_fill_manual(values = gradient_colors) +
  scale_color_manual(values = c("white" = colour_palette[4], "black" = "black")) +
  scale_y_continuous(limits = c(0, NA)) +  # Cut off values below 0 - rmse cannot be negative, but the violin plot makes it look like there are neg values even though there aren't
  coord_cartesian(ylim = c(0, 0.75))

## interval coverage ----------------------------------------------------------

p_ic <- ggplot(int_cov, aes(x = SD, y = int_coverage, fill = SD)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~species_formatted) +
  labs(y = "Coverage") +
  (custom_theme %+replace% theme(
    axis.title.x = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x  = element_text(angle = 90, colour = "#08415C")
  )) +
  geom_hline(data = int_cov, 
             aes(yintercept = 0.95), linetype = "dashed", color = colour_palette[6], size = 1) +  
  scale_fill_manual(values = gradient_colors) +
  coord_cartesian(ylim = c(0.5, 1))



## p-value --------------------------------------------------------------------

species_lookup_mod <- species_lookup %>%
  mutate(short_species = str_extract(species, "[^_]+_[^_]+$"))

sig <- sig %>% 
  mutate(SD = recode(SD,
                     "simple_spabal" = "Si_SB",
                     "stratified_spabal_2030" = "St_SB_2030",
                     "stratified_spabal_2525" = "St_SB_2525",
                     "pref_2030" = "P_2030",
                     "pref_2525" = "P_2525",
                     "clustered" = "C"
  )) %>% 
  left_join(species_lookup_mod %>% select(species_formatted, species),
            by = "species_formatted")

summary_stats <- sig %>%
  group_by(species_formatted, SD) %>%
  summarise(
    total = n(),
    under_0.05 = sum(p_value < 0.05, na.rm = TRUE),
    perc_under_0.05 = round(100 * under_0.05 / total, 2)
  ) %>%
  arrange(desc(perc_under_0.05))

p_sig <- ggplot(sig, aes(x = SD, y = p_value)) +
  geom_boxplot(outlier.size = 0.5, varwidth = T, aes(fill = SD)) +
  facet_wrap(~ species_formatted) +
  labs(
    title = "",
    x = "",
    y = "P-value for\nNTZ effect"
  ) +
  custom_theme +
  theme(
    axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank(),  # Removes x-axis title
    strip.text.x = element_blank()  # Removes column facet titles (species names
  ) +
  scale_fill_manual(values = gradient_colors) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.01, ymax = 0.05, alpha = 0.25, fill = colour_palette[1])


# model estimate --------------------------------------------------------------

test <- realised_ab %>%
  mutate(ra = mean_abundance_outside/mean_abundance_inside)

test2 <- sig %>% 
  glimpse()

lookup <- tibble::tribble(
  ~species_formatted, ~species_raw,
  "Djubitj/West Australian Dhufish", "Glaucosomatidae_Glaucosoma_hebraicum",
  "Yijarup/Pink Snapper", "Sparidae_Chrysophrys_auratus",
  "Southern Maori Wrasse", "Labridae_Ophthalmolepis_lineolatus"
)

# Join to attach realised abundance
test2 <- test2 %>%
  left_join(lookup, by = "species_formatted") %>%
  left_join(
    realised_ab %>% dplyr::select(species_raw, realised_abundance_ratio),
    by = "species_raw"
  ) %>%
  rename(ra = realised_abundance_ratio)



# Plot
p_est <- ggplot(test2, aes(x = SD, y = estimate)) +
  geom_boxplot(outlier.shape = NA, varwidth = TRUE, aes(fill = SD)) +
  facet_wrap(~ species_formatted) +
  geom_hline(aes(yintercept = ra), linetype = "dashed", color = colour_palette[1], size = 1) +
  labs(
    x = "",
    y = "Estimate for\nNTZ effect"
  ) +
  custom_theme +
  scale_fill_manual(values = gradient_colors) +
  theme(
    axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank(),  # Removes x-axis title
    strip.text.x = element_blank()  # Removes column facet titles (species names
  ) +
  coord_cartesian(ylim = c(0, 2))


# all plots -------------------------------------------------------------------

if (study_site == "waatu") {
  waatu_detail <- (p_raw / p_sig / p_est / p_mb / p_rmse / p_ic) + 
    plot_layout(heights = rep(1, 6)) + 
    plot_annotation(
      title = '        Waatu/Capes Region',
      tag_levels = "A", tag_suffix = "."
    ) & 
    theme(plot.margin = margin(rep(0, 4)))  # reduce plot margins for each subplot
}
if (study_site == "waatern") {
  waatern_detail <- (p_raw / p_sig / p_est / p_mb / p_rmse / p_ic) + 
    plot_layout(heights = rep(1, 6)) + 
    plot_annotation(
      title = '        Waatern/Geographe Bay',
      tag_levels = "A", tag_suffix = "."
    ) & 
    theme(plot.margin = margin(rep(0, 4)))  # reduce plot margins for each subplot
}
detail <- waatu_detail | waatern_detail
ggsave("plots/figA1.png", plot = detail, width = 20, height = 10, dpi = 1500)


## Figure 4. overall ranking --------------------------------------------------

# I want to rank the sampling designs based on the metrics. 

# coverage
glimpse(int_cov)
int_cov_rank <- int_cov %>% 
  group_by(species) %>% 
  mutate(
    dist_to_ideal = abs(0.95 - int_coverage), # distance to the ideal 0.95 value
    int_cov_rank = scales::rescale(-dist_to_ideal, to = c(0, 1))) %>% 
  dplyr::select(SD, species, int_cov_rank)
int_cov_rank

# rmse
glimpse(rmse)
rmse_rank <- rmse %>% 
  group_by(species, SD) %>% 
  summarise(
    rmse_mean = mean(rmse),
    rmse_se = sd(rmse, na.rm = TRUE) / sqrt(sum(!is.na(rmse)))
  ) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  mutate(
    rmse_mean_rank = scales::rescale(-rmse_mean, to = c(0, 1)), # 'minus' because higher rank goes to lower values
    rmse_se_rank = scales::rescale(-rmse_se, to = c(0, 1)) # 'minus' because higher rank goes to lower values
  ) %>% 
  dplyr::select(SD, species, rmse_mean_rank, rmse_se_rank)
rmse_rank

# mean bias
glimpse(mb)
mb_rank <- mb %>% 
  group_by(species, SD) %>% 
  summarise(
    mb_mean = abs(mean(mean_bias)), # distance to zero
    mb_se = sd(mean_bias, na.rm = TRUE) / sqrt(sum(!is.na(mean_bias)))
  ) %>%
  ungroup() %>% 
  group_by(species) %>% 
  mutate(
    mb_mean_rank = scales::rescale(-abs(mb_mean), to = c(0, 1)),
    mb_se_rank = scales::rescale(-mb_se, to = c(0, 1))
  )
mb_rank  

# model estimate
glimpse(sig)

sig2 <- sig %>% # add realised abundance in the dataframe to compare to the observed abundance 
  left_join(lookup, by = "species_formatted") %>%
  left_join(
    realised_ab %>% dplyr::select(species_raw, realised_abundance_ratio),
    by = "species_raw"
  ) %>%
  dplyr::select(-species_raw) %>% 
  rename(ra = realised_abundance_ratio)
glimpse(sig2)

est_rank <- sig2 %>% 
  group_by(species_formatted, SD, ra) %>% 
  summarise(
    est_mean = mean(estimate),
    est_se = sd(estimate, na.rm = TRUE) / sqrt(sum(!is.na(estimate)))
  ) %>% 
  ungroup() %>% 
  group_by(species_formatted) %>% 
  mutate(
    est_mean_rank = scales::rescale((est_mean - ra), to = c(0, 1)),
    est_se_rank = scales::rescale(-(est_se - ra), to = c(0, 1))
  )
est_rank

# significance
glimpse(sig)
sig_rank <- sig %>% 
  group_by(species_formatted, SD) %>% 
  summarise(
    n_mod = n(),
    n_mod_noNA = sum(!is.na(p_value)),
    sig_tot = sum(p_value <= 0.05, na.rm = TRUE) / n_mod_noNA
  ) %>% 
  ungroup() %>% 
  group_by(species_formatted) %>% 
  mutate(
    sig_rank = scales::rescale(sig_tot, to = c(0, 1))
  )
sig_rank

# join all ranks together
rank1 <- merge(mb_rank, rmse_rank)
rank2 <- merge(sig_rank, est_rank)
rank3 <- merge(rank1, int_cov_rank) %>% 
  mutate(SD = fct_recode(SD, 
                    "pref_2030" = "P_2030",
                    "pref_2525" = "P_2525",
                    "stratified_spabal_2030" = "St_SB_2030",
                    "stratified_spabal_2525" = "St_SB_2525",
                    "simple_spabal" = "Si_SB",
                    "clustered" = "C"
  )) %>% 
  left_join(species_lookup, by = "species")
  

rank <- merge(rank2, rank3) %>% 
  dplyr::select(-c(mb_mean, mb_se, est_mean, est_se, ra, n_mod, sig_tot, n_mod_noNA, species)) %>% 
  mutate(
    tot_rank = scales::rescale((mb_mean_rank + mb_se_rank + rmse_mean_rank + rmse_se_rank + int_cov_rank + 
                                  1.5 * (est_mean_rank + est_se_rank + sig_rank)), to = c(0, 1))
  )
rank

rank_long <- pivot_longer(rank, cols = sig_rank:tot_rank, names_to = "metric")
rank_long$metric <- factor(rank_long$metric, levels = rev(c("tot_rank", 
                                                            "sig_rank", "est_mean_rank", "est_se_rank", 
                                                            "mb_mean_rank", "mb_se_rank",
                                                            "rmse_mean_rank", "rmse_se_rank",
                                                            "int_cov_rank")))

if (study_site == "waatern"){
  waatern_ranking_long <- rank_long %>% mutate(site = "Waatern\nGeographe Bay")
} else {
  waatu_ranking_long <- rank_long %>% mutate(site = "Waatu\nCapes Region")
}

# RUN FROM THE 'ranking' TITLE WITH BOTH SITES
all_ranking <- rbind(waatern_ranking_long, waatu_ranking_long)

all_ranking <- all_ranking %>% 
  left_join(species_lookup, by = "species_formatted") %>%
  mutate(
        species_formatted = factor(
      species_formatted,
      levels = c(
        "Yijarup/Pink Snapper",
        "Southern Maori Wrasse",
        "Djubitj/West Australian Dhufish"
      )
    ),
    species_formatted = fct_recode(species_formatted, 
                                   "Djubitj/West \nAustralian Dhufish" = "Djubitj/West Australian Dhufish")
  ) %>%
  mutate(
    metric = case_when(metric == "sig_rank" ~ "NTZ effect significance *",
                       metric == "est_mean_rank" ~ "Mean effect size *",
                       metric == "est_se_rank" ~ "Effect size SE *",
                       metric == "mb_mean_rank" ~ "Mean bias",
                       metric == "mb_se_rank" ~ "Mean bias SE",
                       metric == "rmse_mean_rank" ~ "Mean RMSE",
                       metric == "rmse_se_rank" ~ "RMSE SE",
                       metric == "int_cov_rank" ~ "Interval Coverage",
                       metric == "tot_rank" ~ "Overall rank",
                       .default = metric
    )
  ) %>% 
  mutate(
    metric = factor(
      metric,
      levels = rev(c("NTZ effect significance *", "Mean effect size *", "Effect size SE *", 
                     "Mean bias", "Mean bias SE", "Mean RMSE", "RMSE SE", "Interval Coverage", "Overall rank"))),
    SD = fct_recode(SD, 
                    "Preferential (20/30)" = "pref_2030",
                    "Preferential (25/25)" = "pref_2525",
                    "Str. sp. balanced (20/30)" = "stratified_spabal_2030",
                    "Str. sp. balanced (25/25)" = "stratified_spabal_2525",
                    "Simple sp. balanced" = "simple_spabal",
                    "Clustered" = "clustered"
    )
  )

ranking_plot <- ggplot() +
  # First: plot all non-total metrics
  geom_tile(
    data = all_ranking %>% filter(metric != "Overall rank"),
    aes(x = metric, y = SD, fill = value),
    color = "white"
  ) +
  scale_fill_gradient(high = colour_palette[2], low = "white", name = "Rank") +
  
  # New fill scale for total_rank
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = all_ranking %>% filter(metric == "Overall rank"),
    aes(x = metric, y = SD, fill = value),
    color = "white"
  ) +
  scale_fill_gradient(high = colour_palette[6], low = "white", name = "Overall rank") +
  facet_grid(site ~ species_formatted) +
  labs(
    x = "Metric",
    y = "Sampling Design"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  custom_theme +
  coord_flip()
ranking_plot
ggsave("plots/fig4a.png", plot = ranking_plot, width = 8, height = 6.5, dpi = 600)


sites_overall_rank <- all_ranking %>% 
  dplyr::filter(metric == "Overall rank") %>% 
  group_by(SD, site) %>% 
  summarise(rank_all = sum(value))

sites_rank_plot <- ggplot(sites_overall_rank, aes(x=SD, y=rank_all, fill = SD)) + 
  geom_bar(stat = "identity") +
  custom_theme +
  #  theme(axis.text.x = element_blank(),
  #        axis.text.y = element_blank(),
  #        strip.text.x = element_blank()) +
  #oord_flip() +
  labs(x = "", y = "relative\nperf.") +
  #facet_wrap(~site) +
  scale_fill_manual(values = gradient_colors)
sites_rank_plot

spp_overall_rank <- all_ranking %>% 
  dplyr::filter(metric == "Overall rank") %>% 
  group_by(SD, species) %>% 
  summarise(rank_all = sum(value))

spp_rank_plot <- ggplot(spp_overall_rank, aes(x=SD, y=rank_all, fill = SD)) + 
  geom_bar(stat = "identity") +
  custom_theme +
  #    theme(axis.text.x = element_blank(),
  #          axis.text.y = element_blank(),
  #          strip.text.x = element_blank()) +
  #    coord_flip() +
  labs(x = "", y = "relative\nperf.") +
  facet_wrap(~species) +
  scale_fill_manual(values = gradient_colors)
spp_rank_plot


(spp_rank_plot | sites_rank_plot) +
  plot_layout(heights = c(1, 2), widths = c(2, 1))


overall_rank <- all_ranking %>% 
  dplyr::filter(metric == "Overall rank") %>% 
  group_by(SD, species, site) %>% 
  summarise(rank_all = sum(value))
overall <- ggplot(overall_rank, aes(x=SD, y=rank_all, fill = SD)) + 
  geom_bar(stat = "identity") +
  custom_theme +
  #theme(axis.text.x = element_blank()) +
  #coord_flip() +
  labs(x = "", y = "relative\nperf.") +
  facet_grid(site~species) +
  scale_fill_manual(values = gradient_colors)
overall

ggsave("plot/fig4b.png", plot = overall, width = 3, height = 2, dpi = 600)






# below was wrong ---






# we'll rank each SD by how it performs for each spp. 
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# Summarise mean_bias by SD and species
mb_summary <- mb %>% group_by(SD, species) %>% summarise(mean_bias_mean = mean(mean_bias, na.rm = TRUE), mean_bias_se = se(mean_bias), .groups = "drop")

# Summarise rmse by SD and species
rmse_summary <- rmse %>% group_by(SD, species) %>% summarise(rmse_mean = mean(rmse, na.rm = TRUE), rmse_se = se(rmse), .groups = "drop")

# Summarise estimate and p-value by SD and species
sig_summary <- test2 %>% group_by(SD, species_raw) %>% 
  summarise(estimate_mean = mean(estimate, na.rm = TRUE), estimate_se = se(estimate), n_total = n(), n_signif = sum(p_value < 0.05, na.rm = TRUE), prop_signif = n_signif / n_total * 100, .groups = "drop")
ra_lookup <- test2 %>% dplyr::select(species_raw, ra) %>% distinct()
sig_summary <- sig_summary %>% left_join(ra_lookup, by = "species_raw") %>% mutate(est_diff = abs(estimate_mean - ra)) %>% group_by(species_raw) %>% arrange(species_raw) %>% rename(species = species_raw) %>% ungroup()

# Summarise interval coverage by SD and species
int_cov_summary <- int_cov %>% group_by(SD, species) %>% summarise(int_cov_mean = mean(int_coverage, na.rm = TRUE), .groups = "drop")

summary_by_SD_species <- mb_summary %>%
  full_join(rmse_summary, by = c("SD", "species")) %>%
  full_join(int_cov_summary, by = c("SD", "species")) %>% 
  full_join(sig_summary, by = c("SD", "species"))

realised_ratios <- readRDS(paste0("outputs/distribution_modelling_outputs/D06_sampling_design_performance/D06_", study_site, "_sampling_design_performance_realised_abundances.rds"))

# rank table
percentile_rank <- function(x) ecdf(x)(x)

ranking_table <- summary_by_SD_species %>%
  left_join(realised_ratios, by = c("species" = "species_raw")) %>%
  group_by(species) %>%
  mutate(
    estimate_error = abs(estimate_mean - realised_abundance_ratio),
    # if the goal is a small value, use "1 - ..."
    mean_bias_rank     = 1 - scales::rescale(abs(mean_bias_mean), to = c(0, 1)),
    mean_bias_se_rank  = 1 - scales::rescale(mean_bias_se, to = c(0, 1)),
    rmse_rank          = 1 - scales::rescale(rmse_mean, to = c(0, 1)),
    rmse_se_rank       = 1 - scales::rescale(rmse_se, to = c(0, 1)),
    int_cov_rank       = 1 - scales::rescale(abs(int_cov_mean - 0.95), to = c(0, 1)),
    estimate_mean_rank = scales::rescale(estimate_error, to = c(0, 1)),
    estimate_se_rank   = 1 - scales::rescale(estimate_se, to = c(0, 1)),
    sig_rank           = scales::rescale(prop_signif, to = c(0, 1)),
    
    total_rank = mean_bias_rank +
      mean_bias_se_rank +
      rmse_rank +
      rmse_se_rank +
      int_cov_rank +
      1.5 * estimate_mean_rank +
      1.5 * estimate_se_rank +
      1.5 * sig_rank
  ) %>%
  arrange(species, desc(total_rank)) %>%
  ungroup()

ranking_long <- ranking_table %>% group_by(species) %>% 
  ungroup() %>% 
  dplyr::select( # select only rank columns 
    species, SD, 
    mean_bias_rank, mean_bias_se_rank, 
    rmse_rank, rmse_se_rank, 
    int_cov_rank, 
    estimate_mean_rank, estimate_se_rank, 
    sig_rank, total_rank ) %>% 
  pivot_longer( cols = ends_with("_rank"), names_to = "metric", values_to = "rank" ) %>% 
  glimpse()

ranking_long <- ranking_long %>%
  mutate(metric = factor(metric, levels = c(
    "total_rank",    
    "int_cov_rank", 
    "rmse_se_rank", 
    "rmse_rank", 
    "mean_bias_se_rank", 
    "mean_bias_rank", 
    "estimate_se_rank", 
    "estimate_mean_rank", 
    "sig_rank"
  )), # total rank will always be at the top, but other than that it goes bottom up
  # Rename the labels for display in plot
  metric = fct_recode(metric,
                      "NTZ effect significance *" = "sig_rank",
                      "Effect size SE *" = "estimate_se_rank",
                      "Mean effect size *" = "estimate_mean_rank",
                      "Mean bias" = "mean_bias_rank",
                      "Mean bias SE" = "mean_bias_se_rank",
                      "Mean RMSE" = "rmse_rank",
                      "RMSE SE" = "rmse_se_rank",
                      "Interval Coverage" = "int_cov_rank",
                      "Total rank" = "total_rank"
  ),
  species = fct_recode(as.factor(species), 
                       "Yijarup\nPink Snapper" = "Sparidae_Chrysophrys_auratus", 
                       "Djubitj\nWest Australian\nDhufish" = "Glaucosomatidae_Glaucosoma_hebraicum", 
                       "Southern\nMaori Wrasse" = "Labridae_Ophthalmolepis_lineolatus"
  ),
  SD = fct_recode(SD, 
                  "Preferential (20/30)" = "P_2030",
                  "Preferential (25/25)" = "P_2525",
                  "Str. sp. balanced (20/30)" = "St_SB_2030",
                  "Str. sp. balanced (25/25)" = "St_SB_2525",
                  "Simple sp. balanced" = "Si_SB",
                  "Clustered" = "C"
  )
  )

if (study_site == "waatern"){
  waatern_ranking_long <- ranking_long %>% mutate(site = "Waatern\nGeographe Bay")
} else {
  waatu_ranking_long <- ranking_long %>% mutate(site = "Waatu\nCapes Region")
} 

# RUN FROM THE 'ranking' TITLE WITH BOTH SITES
all_ranking <- rbind(waatern_ranking_long, waatu_ranking_long) %>%
  group_by(site, species, metric) %>%
  mutate(
    rank_scaled = scales::rescale(rank, to = c(0, 1))
  ) %>%
  ungroup()

all_ranking$species <- factor(all_ranking$species, levels = c("Yijarup\nPink Snapper",  "Southern\nMaori Wrasse", "Djubitj\nWest Australian\nDhufish"))
ranking_plot <- ggplot() +
  # First: plot all non-total metrics
  geom_tile(
    data = all_ranking %>% filter(metric != "Total rank"),
    aes(x = metric, y = SD, fill = rank_scaled),
    color = "white"
  ) +
  scale_fill_gradient(high = colour_palette[2], low = "white", name = "Rank") +
  
  # New fill scale for total_rank
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = all_ranking %>% filter(metric == "Total rank"),
    aes(x = metric, y = SD, fill = rank_scaled),
    color = "white"
  ) +
  scale_fill_gradient(high = colour_palette[6], low = "white", name = "Total Rank") +
  
  facet_grid(site ~ species) +
  labs(
    x = "Metric",
    y = "Sampling Design"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  custom_theme +
  coord_flip()
ranking_plot
ggsave("plots/fig4a.png", plot = ranking_plot, width = 8, height = 6.5, dpi = 600)


overall_rank <- all_ranking %>% 
  dplyr::filter(metric == "Total rank") %>% 
  group_by(SD) %>% 
  summarise(rank_all = sum(rank))
overall <- ggplot(overall_rank, aes(x=SD, y=rank_all, fill = SD)) + 
  geom_bar(stat = "identity") +
  custom_theme +
  theme(axis.text.x = element_blank()
  ) +
  coord_flip() +
  labs(x = "", y = "relative\nperf.") +
  scale_fill_manual(values = gradient_colors)
overall
ggsave("plots/fig4b.png", plot = overall, width = 3, height = 2, dpi = 600)



## Figure 2. little plot of strata --------------------------------------------

strat <- readRDS(paste0("outputs/distribution_modelling_outputs/D04_waatu_sample_area_detrended_strata.rds")) %>% 
  mutate(detrended = case_when(detrended == 0 ~ "Strata 1 (low)",
                               detrended == 1 ~ "Strata 2 (medium)",
                               detrended == 2 ~ "Strata 3 (high)")) %>% 
  rename('detrended bathymetry strata' = detrended) %>% 
  glimpse()

strata_plot <- ggplot(strat) +
  geom_sf(aes(fill = `detrended bathymetry strata`), col = NA) +
  scale_fill_manual(values=c(colour_palette[3:1])) +
  coord_sf(ylim = c(6225877.51105516, 6232686.747557), xlim = c(308922.377479875, 313518.321854704)) +
  theme_minimal() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank())

save_plot(paste0("plots/fig2 (part).png"), strata_plot, base_width = 3, base_height = 3)

## END ##
