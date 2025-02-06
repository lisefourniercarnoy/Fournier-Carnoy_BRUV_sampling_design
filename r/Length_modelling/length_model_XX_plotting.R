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

