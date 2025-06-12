# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare Lengths detected)
# Data:    Geographe Bay Shapefiles and BRUV data.
# Task:    Plot figures for paper
# Author:  Lise Fournier-Carnoy
# Date:    December 2024

# -----------------------------------------------------------------------------

# Status:  Starting

# -----------------------------------------------------------------------------

library(sf)
library(raster)
library(ggplot2)
library(tmap)
library(terra)

colour_palette <- c("#6B0504", "#08415C", "#3E6990", "#fff2cc", "#f1c232", "#7f6000")
saveRDS(colour_palette, "chapter_colours.rds")

## Basic map of study area ----------------------------------------------------

aus                   <- st_read("data/spatial/shapefiles/wadandi_land.shp")
sim_sz                <- st_read("data/spatial/shapefiles/simulated_SZ.shp")
sampling_area         <- readRDS("data/spatial/shapefiles/sampling_area_deeper_than_7m.rds")
mature_pred           <- raster("outputs/Length_comparison/predicted_abundance_rasters/mature_predicted_abundance_SZ_increase.rds")


raster_df <- as.data.frame(mature_pred, xy = TRUE)

ggplot() +
  coord_sf(crs = 4326, xlim = ext(sampling_area)[1:2], ylim = ext(sampling_area)[3:4]) +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = p_mature.fit), interpolate = TRUE) +
  geom_sf(data = sampling_area, aes(color = 'Sampling Area'),
          fill = 'transparent', linewidth = 1) +
  geom_sf(data = sim_sz, aes(color = 'Simulated NTZ'),
          fill = 'transparent', linewidth = 1) +
  scale_color_manual(name = '',
                     values = c('Sampling Area' = '#f1c232', 'Simulated NTZ' = '#7f6000')) +
  geom_sf(data = aus,
          color = "darkgray", size = 0.2) +
  scale_fill_gradient(name = 'Predicted Abundance',
                      low = "#fff2cc", high = "#7f6000", limits = c(0, 2), na.value=NA) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())









## Example sampling design
ggplot() +
  geom_sf(data = sf_all, aes(fill = as.factor(detrended)), colour = NA) +
  scale_fill_manual(values = c("0" = "#fff2cc", "1" = "#ffd966", "2" = "#f1c232"),
                    na.value = NA,
                    labels = c("0" = "Strata 1 - low",
                               "1" = "Strata 2 - medium",
                               "2" = "Strata 3 - high")) +
  labs(fill = "Detrended bathymetry",
       title = "a. Example spatially balanced sampling design") +
  new_scale_fill() +
  geom_sf(data = aus) +
  geom_sf(data = SZ, colour = "#7f6000", fill = NA, linewidth = 1) +
  geom_sf(data = samples_sf, colour = "red") +
  coord_sf(crs = 4326, xlim = ext(samp_area)[1:2], ylim = ext(samp_area)[3:4]) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)  # Set the longitude (x-axis) labels to horizontal
) +
  theme(legend.position = "bottom",
         legend.box = "horizontal",
         legend.box.spacing = unit(0, "cm"),
         legend.key.size = unit(0.5, "cm"),
         legend.text = element_text(size = 10))


# Plot real BRUV location
samp_area <- st_read("QGIS layers/polys/com_comp_analysis_limits.shp")
AMP <- st_read("QGIS layers/polys/AMP_local.shp")
aus <- st_read("data/spatial/shapefiles/wadandi_land.shp")
pref <- st_read("QGIS layers/2024_pref_BRUVs.shp")
spabal <- st_read("QGIS layers/2024_spabal_BRUVs.shp")
MPA <- st_read("QGIS layers/polys/swc_sanctuaryzones.shp"); MPA <- st_zm(MPA); st_crs(MPA) <- 4326

library(terra)
ggplot() +
  geom_sf(data = aus, fill = "lightgrey", color = "darkgray", size = 0.2) +
  geom_sf(data = samp_area,
          aes(fill = 'Community \ncomposition \nanalysis limits'), color = "#f1c232", size = 5) +
  geom_sf(data = AMP[AMP$ZoneName == 'National Park Zone',],
          aes(fill = 'No-take \nzones'), color = "#ffd966", size = 0.5) +
  geom_sf(data = MPA,
          aes(fill = 'No-take \nzones'), color = "#ffd966", size = 0.2) +
  geom_sf(data = pref[pref$longitude > 115.3,],
          aes(color = 'Preferential \nsampling'), size = 1, shape = 16) +
  geom_sf(data = spabal[spabal$SD == 'spabal',],
          aes(color = 'Spatially \nbalanced \nsampling'), size = 1, shape = 16) +
  geom_sf(data = spabal[spabal$SD == 'pref',],
          aes(color = 'Preferential \nsampling'), size = 1, shape = 16) +
  coord_sf(crs = 4326, xlim = ext(spabal)[1:2], ylim = ext(spabal)[3:4]) +
  theme_minimal() +
  scale_fill_manual(values = c('No-take \nzones' = "#fff2cc", 'Community \ncomposition \nanalysis limits' = 'transparent')) +
  scale_color_manual(values = c('Preferential \nsampling' = "#f1c232", 'Spatially \nbalanced \nsampling' = "#7f6000")) +
  labs(fill = '', color = '') +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.spacing = unit(0, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 10))


# real bruv locations (correct colours)


aus             <- st_read("data/spatial/shapefiles/wadandi_land.shp") # land shapefile
com_samp_area   <- st_read("QGIS layers/polys/com_comp_analysis_limits.shp")
AMP             <- st_read("QGIS layers/polys/AMP_local.shp") # real SZs
MPA             <- st_read("QGIS layers/polys/swc_sanctuaryzones.shp"); MPA <- st_zm(MPA); st_crs(MPA) <- 4326 # real SZs

pref            <- st_read("QGIS layers/2024_pref_BRUVs.shp") # real BRUV state points
spabal          <- st_read("QGIS layers/2024_spabal_BRUVs.shp") # real BRUV commonwealth points
MaxN            <- read.csv("data/tidy/2024_geographe_all_tidy_maxn.csv") # all MaxNs


library(ggpubr) # for ggarrange
library(grid)
library(imager)
library(ggspatial) # for scale bars

palette <- c("#E5B25D", "#A3320B", "#6B0504", "#394F49", "#65743A")

# Load images as raster grobs
image1 <- rasterGrob(load.image("photos/wadandi_ranger_BRUV.JPG"), interpolate = TRUE)
image2 <- rasterGrob(load.image("photos/BRUV_photo_yijarup_skippy.PNG"), interpolate = TRUE)

# Arrange images into ggplot objects
image_plot1 <- ggplot() + annotation_custom(image1) + theme_void()
image_plot2 <- ggplot() + annotation_custom(image2) + theme_void() 

plot <- ggplot() +
  geom_sf(data = com_samp_area,
          aes(fill = 'Community \ncomposition \nanalysis limits'), color = "#A3320B", lwd = 1) +
  geom_sf(data = AMP[AMP$ZoneName == 'National Park Zone',],
          aes(fill = 'No-take \nzones'), alpha = 0.1, color = "#65743A", lwd = 1) +
  geom_sf(data = MPA,
          aes(fill = 'No-take \nzones'), alpha = 0.1, color = "#65743A", lwd = 1) +
  geom_sf(data = pref[pref$longitude > 115.3,],
          aes(color = 'Preferential \nsampling'), size = 1, shape = 16) +
  geom_sf(data = spabal[spabal$SD == 'spabal',],
          aes(color = 'Spatially \nbalanced \nsampling'), size = 1, shape = 16) +
  geom_sf(data = spabal[spabal$SD == 'pref',],
          aes(color = 'Preferential \nsampling'), size = 1, shape = 16) +
  geom_sf(data = aus, fill = "lightgrey", color = "darkgray", lwd = 1) +  # Move 'aus' layer to the top
  coord_sf(crs = 4326, xlim = ext(spabal)[1:2], ylim = ext(spabal)[3:4]) +
  theme_minimal() +
  scale_fill_manual(values = c('No-take \nzones' = "#65743A", 'Community \ncomposition \nanalysis limits' = 'transparent')) +
  scale_color_manual(values = c('Preferential \nsampling' = "#E5B25D", 'Spatially \nbalanced \nsampling' = "#6B0504")) +
  labs(fill = '', color = '') +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.spacing = unit(0, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 10)) +
  annotation_scale(location = "br", width_hint = 0.2)  # Adds a scale bar



plot_zoom <- ggplot() + # zoomed in
  geom_sf(data = com_samp_area,
          aes(fill = 'Community \ncomposition \nanalysis limits'), color = "#A3320B", lwd = 1) +
  geom_sf(data = AMP[AMP$ZoneName == 'National Park Zone',],
          aes(fill = 'No-take \nzones'), alpha = 0.1, color = "#65743A", lwd = 1) +
  geom_sf(data = MPA,
          aes(fill = 'No-take \nzones'), alpha = 0.1, color = "#65743A", lwd = 1) +
  geom_sf(data = pref[pref$longitude > 115.3,],
          aes(color = 'Preferential \nsampling'), size = 1, shape = 16) +
  geom_sf(data = spabal[spabal$SD == 'spabal',],
          aes(color = 'Spatially \nbalanced \nsampling'), size = 1, shape = 16) +
  geom_sf(data = spabal[spabal$SD == 'pref',],
          aes(color = 'Preferential \nsampling'), size = 1, shape = 16) +
  geom_sf(data = aus, fill = "lightgrey", color = "darkgray", lwd = 1) +
  coord_sf(crs = 4326, xlim = ext(com_samp_area)[1:2], ylim = ext(com_samp_area)[3:4]) +
  theme_minimal() +
  scale_fill_manual(values = c('No-take \nzones' = "#65743A", 'Community \ncomposition \nanalysis limits' = 'transparent')) +
  scale_color_manual(values = c('Preferential \nsampling' = "#E5B25D", 'Spatially \nbalanced \nsampling' = "#6B0504")) +
  labs(fill = '', color = '') +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.spacing = unit(0, "cm"),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        axis.text.x = element_blank(),  # Remove x-axis 
        axis.text.y = element_blank()   # Remove y-axis 
  ) +
  annotation_scale(location = "br", width_hint = 0.2)  # Adds a scale bar

# Combine the two plots with a shared legend
combined_plots <- ggarrange(
  plot, plot_zoom,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom",
  labels = c("C", "D"),
  widths = c(1.45, 1)  # Equal heights for both plots
)

# Arrange everything
final_layout <- ggarrange(
  ggarrange(image_plot1, 
            image_plot2, 
            ncol = 2, 
            labels = c("A", "B"),
            widths = c(1, 1)),
  combined_plots,
  nrow = 2, 
  heights = c(1.1, 1.6)
)

# Display the final layout
final_layout



### size distributions ----

## Files used in this script --------------------------------------------------

file_fish_length        <- "data/raw/2024_state/2024_02_NgariCapes.MP.Monitoring_stereoBRUVs_Lengths.txt" # obtained from 00_checkem
file_metadata           <- "data/raw/2024_state/2024_02_NgariCapes.MP.Monitoring_stereoBRUVs_Metadata.csv" # obtained from MEGlab labsheets

file_obs_hab            <- "data/raw/2024_commonwealth/2024-04_Geographe_stereo-BRUVs_forwards_Dot Point Measurements.txt" # output form TransectMeasure

file_250m_bathy         <- "data/spatial/rasters/GB-SW_250mBathy.tif" # from Geoscience Australia
file_250m_bathy_deriv   <- "data/spatial/rasters/2014_geographe_bathymetry-derivatives.RDS" # From Claude's code
file_lidar_bathy        <- "data/spatial/rasters/LiDAR_geo_compressed.tif" # Too big for git


## Presence/Absence data ------------------------------------------------------

L50 <- c('Chrysophrys auratus' = 566, # L50 from Wakefield et al. 2015 - males 585, females 566
         'Glaucosoma hebraicum' = 0, #301, # L50 from Hesp 2002 - males 320, females 301 - but using all fish because papers find NTZ makes dhuies more abundant not bigger
         "Ophthalmolepis lineolatus" = 184 # l50 from females, Morton 2008
) 

metadata <- read.csv(file_metadata) %>% 
  dplyr::filter(successful_count == "Yes" & successful_length == "Yes",
                grepl("EGB", opcode)) %>% 
  dplyr::select("opcode", "latitude_dd", "longitude_dd", "status", "depth_m") %>% 
  glimpse()

length_data <- read.table(file_fish_length, sep = "\t", header = T) %>%
  rename(opcode = OpCode) %>% 
  mutate(full_spp = paste(Genus, Species)) %>% 
  dplyr::filter(grepl("EGB", opcode)) %>% 
  glimpse()

names(length_data)
unique(length_data$full_spp)
dat <- length_data %>% 
  # Select species present in L50
  filter(full_spp %in% names(L50)) %>% 
  
  # Join the metadata with the length data
  left_join(metadata %>% distinct(opcode, longitude_dd, latitude_dd), by = "opcode") #%>% 
  
  # Separate into size categories, treating missing lengths as immature
  # group_by(longitude_dd, latitude_dd, full_spp, opcode) %>% 
  # summarise(
  #   count_mature = sum(coalesce(Length, 0) > L50[full_spp], na.rm = TRUE),
  #   count_immature = sum(coalesce(Length, 0) <= L50[full_spp], na.rm = TRUE),
  #   .groups = "drop"
  # ) %>% 
  # 
  # # Ensure that all full_spp are present even if no data is available for the species
  # complete(opcode, full_spp, fill = list(count_mature = 0, count_immature = 0)) %>% 
  # complete(opcode = unique(metadata$opcode), full_spp = full_spp, fill = list(count_mature = 0, count_immature = 0)) %>% 
  
  # Fill in missing longitude and latitude
  # left_join(metadata %>% distinct(opcode, longitude_dd, latitude_dd), by = "opcode") %>% 
  # mutate(
  #   longitude = coalesce(longitude_dd.x, longitude_dd.y),
  #   latitude = coalesce(latitude_dd.x, latitude_dd.y)
  # ) %>%
  # dplyr::select(-longitude_dd.x, -latitude_dd.x, -longitude_dd.y, -latitude_dd.y) %>% 
  # 
  # # Create an ID column
  # mutate(ID = row_number()) %>% 
  # glimpse()

length(unique(dat$opcode)) == length(unique(metadata$opcode)) # should be true. if not there are opcodes missing.

ggplot(dat, aes(x = longitude_dd, y = latitude_dd, size = Length)) +
  facet_wrap(~full_spp) +
  geom_point()

saveRDS(dat, "data/tidy/mature_presence_latlong.rds") # Needed for 03_gam

summary(as.factor(dat$full_spp[dat$count_mature>0]))


