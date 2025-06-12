# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare Lengths detected)
# Data:    Simulated data
# Task:    Plot and save final plots so that the RMarkdown knit goes faster
# Author:  Lise Fournier-Carnoy
# Date:    April 2025

# -----------------------------------------------------------------------------

# Status:  

# -----------------------------------------------------------------------------
rm(list = ls())  # Clear environment

# Example sampling designs

library(tidyverse) # for data manipulation
library(sf) # for dealing with shapefiles
library(raster) # for dealing with rasters
library(terra) # for dealing with rasters
library(ggpubr) # for plot details
library(ggtext) # for plot details
library(cowplot) # for plot saving
library(magick) # for final plot stitching


colour_palette <- eval(parse(text = readLines("chapter_colours.txt")))
source("custom_theme.R")

# Load layers
target_crs <- 4326

aus                 <- st_read("data/spatial/shapefiles/wadandi_land.shp") %>% st_transform(crs=target_crs) # land shapefile
SZ                  <- st_read("data/spatial/shapefiles/simulated_SZ.shp") %>% st_transform(crs=target_crs) # simulated SZ shapefile
sampling_area       <- readRDS("data/spatial/shapefiles/sampling_area_deeper_than_7m.rds") %>% st_transform(crs=target_crs) # simulation area

sf_all              <- readRDS('data/rmd/samp_area_detrended.rds') %>% st_transform(crs=target_crs) # detrended bathymetry strata in the simulation area

samples_sf_rand     <- readRDS('data/rmd/samples_sf_srs.rds') %>% st_transform(crs=target_crs) %>% dplyr::filter(design_id == 1) # example points for simulated random design
samples_sf_pref_2525<- readRDS('data/rmd/samples_sf_pref_25_25.rds') %>% st_transform(crs=target_crs) # example points for simulated preferential design
samples_sf_pref_2030<- readRDS('data/rmd/samples_sf_pref_20_30.rds') %>% st_transform(crs=target_crs) # example points for simulated preferential design
samples_sf_sb_2525  <- readRDS('data/rmd/samples_sf_sb_25in_25out.rds') %>% st_transform(crs=target_crs) # example points for simulated spatially balanced design
samples_sf_sb_2030  <- readRDS('data/rmd/samples_sf_sb_20in_30out.rds') %>% st_transform(crs=target_crs) # example points for simulated spatially balanced design
samples_sf_clump    <- readRDS('data/rmd/samples_sf_clump.rds') %>% st_transform(crs=target_crs) # example points for simulated clustered design

mature_pred_PS      <- as.data.frame(readRDS("outputs/Length_comparison/predicted_abundance_rasters/SZ_increase/mature_predicted_abundance_all_lengths_Chrysophrys_auratus_full_data_SZ_increase.rds"), xy=TRUE)
mature_pred_SMW     <- as.data.frame(readRDS("outputs/Length_comparison/predicted_abundance_rasters/SZ_increase/mature_predicted_abundance_all_lengths_Ophthalmolepis_lineolatus_full_data_SZ_increase.rds"), xy=TRUE)
mature_pred_WAD     <- as.data.frame(readRDS("outputs/Length_comparison/predicted_abundance_rasters/SZ_increase/mature_predicted_abundance_all_lengths_Glaucosoma_hebraicum_full_data_SZ_increase.rds"), xy=TRUE)


# Define a function for creating the sampling design plots --------------------

create_sampling_design_plot <- function(sf_data, sampling_area, SZ, aus, samples_sf, colour_palette, legend = TRUE) {
  p <- ggplot() +
    geom_sf(data = sf_data, aes(fill = as.factor(detrended)), colour = NA) +
    scale_fill_manual(values = c("0" = "#F5F5F5", "1" = "#E0E0E0", "2" = "#9E9E9E"),
                      na.value = NA,
                      labels = c("0" = "Strata 1 <br> (low)",
                                 "1" = "Strata 2 <br> (medium)",
                                 "2" = "Strata 3 <br> (high)")) +
    labs(fill = "Detrended <br> bathymetry") +
    geom_sf(data = SZ, aes(colour = "Simulated <br> NTZ"), fill = NA, linewidth = 1) +
    geom_sf(data = aus, color = "darkgray", fill = 'lightgray', size = 0.2) +
    geom_sf(data = samples_sf, aes(colour = "Simulated <br> sample points")) +
    coord_sf(crs = 4326, xlim = ext(sampling_area)[1:2], ylim = ext(sampling_area)[3:4]) +
    scale_color_manual(name = '', values = c('Simulated <br> sample points' = colour_palette[1], 
                                             'Simulated <br> NTZ' = colour_palette[3])) +
    custom_theme +
    theme(
      axis.text.x = element_blank(),  # Explicitly remove x-axis text
      axis.text.y = element_blank(),  # Explicitly remove y-axis text
      axis.ticks = element_blank(),   # Remove axis ticks
      axis.title.x = element_blank(), # Remove x-axis title
      axis.title.y = element_blank(),  # Remove y-axis title
    legend.position = "top"
    )
}

ex_rand <- create_sampling_design_plot(sf_all, sampling_area, SZ, aus, samples_sf_rand, colour_palette, legend = TRUE)
ex_spabal_2525 <- create_sampling_design_plot(sf_all, sampling_area, SZ, aus, samples_sf_sb_2525, colour_palette, legend = TRUE)
ex_spabal_2030 <- create_sampling_design_plot(sf_all, sampling_area, SZ, aus, samples_sf_sb_2030, colour_palette, legend = TRUE)
ex_pref_2525 <- create_sampling_design_plot(sf_all, sampling_area, SZ, aus, samples_sf_pref_2525, colour_palette, legend = TRUE)
ex_pref_2030 <- create_sampling_design_plot(sf_all, sampling_area, SZ, aus, samples_sf_pref_2030, colour_palette, legend = TRUE)
ex_clump <- create_sampling_design_plot(sf_all, sampling_area, SZ, aus, samples_sf_clump, colour_palette, legend = FALSE)


# Arrange the example sampling design plots (Top Row)
ex_maps <- ggarrange(
  ex_rand + theme(legend.position = "top") + annotate("text", 
                                                      x = ext(sampling_area)[1], 
                                                      y = ext(sampling_area)[4], 
                                                      label = "simple random", 
                                                      hjust = 0, vjust = 0.5,
                                                      size = 4), 
  ex_spabal_2030 + theme(legend.position = "top") + annotate("text", 
                                                             x = ext(sampling_area)[1], 
                                                             y = ext(sampling_area)[4], 
                                                             label = "sp. balanced (20/30)", 
                                                             hjust = 0, vjust = 0.5,
                                                             size = 4),  
  ex_spabal_2525 + theme(legend.position = "top") + annotate("text", 
                                                             x = ext(sampling_area)[1], 
                                                             y = ext(sampling_area)[4], 
                                                             label = "sp. balanced (25/25)", 
                                                             hjust = 0, vjust = 0.5,
                                                             size = 4),  
  ex_pref_2030 + theme(legend.position = "none") + annotate("text", 
                                                            x = ext(sampling_area)[1], 
                                                            y = ext(sampling_area)[4], 
                                                            label = "preferential (20/30)", 
                                                            hjust = 0, vjust = 0.5,
                                                            size = 4), 
  ex_pref_2525 + theme(legend.position = "none") + annotate("text", 
                                                            x = ext(sampling_area)[1], 
                                                            y = ext(sampling_area)[4], 
                                                            label = "preferential (25/25)", 
                                                            hjust = 0, vjust = 0.5,
                                                            size = 4), 
  ex_clump + theme(legend.position = "none") + annotate("text", 
                                                        x = ext(sampling_area)[1], 
                                                        y = ext(sampling_area)[4], 
                                                        label = "clustered", 
                                                        hjust = 0, vjust = 0.5,
                                                        size = 4), 
  nrow = 2, ncol = 3,
  widths = c(1, 1, 1),  # Adjust widths for spacing
  labels = c("A.", "B.", "C.", 
             "D.", "E.", "F."),
  common.legend = TRUE
)

save_plot("plots/example_sampling_designs.png", ex_maps, base_width = 8, base_height = 5)


# Define the function for creating the raster plot ----------------------------

create_raster_plot <- function(mature_data, sampling_area, SZ, aus, colour_palette) {
  ggplot() +
    geom_raster(data = mature_data, aes(x = x, y = y, fill = p_mature.fit), interpolate = TRUE) +
    geom_sf(data = sampling_area, aes(color = 'Simulation area'),
            fill = 'transparent', linewidth = 1, show.legend = F) +
    geom_sf(data = SZ, aes(color = 'Simulated NTZ'),
            fill = 'transparent', linewidth = 1, show.legend = F) +
    scale_color_manual(name = '',
                       values = c('Simulation area' = colour_palette[1], 'Simulated NTZ' = colour_palette[3])) +
    geom_sf(data = aus, color = "darkgray", fill = 'lightgray', size = 0.2) +
        custom_theme +
    scale_fill_gradient(
      name = "",
      low = "white", high = colour_palette[4]
      ) + 
    theme(
      axis.text.x = element_blank(),  # Explicitly remove x-axis text
      axis.text.y = element_blank(),  # Explicitly remove y-axis text
      axis.ticks = element_blank(),   # Remove axis ticks
      axis.title.x = element_blank(), # Remove x-axis title
      axis.title.y = element_blank(),  # Remove y-axis title
      legend.position = "none"
      ) +
    coord_sf(
      xlim = c(st_bbox(sampling_area)["xmin"], st_bbox(sampling_area)["xmax"]),
      ylim = c(st_bbox(sampling_area)["ymin"], st_bbox(sampling_area)["ymax"]),
      expand = FALSE
    )
}

PS_raster <- create_raster_plot(mature_pred_PS, sampling_area, SZ, aus, colour_palette)
SMW_raster <- create_raster_plot(mature_pred_SMW, sampling_area, SZ, aus, colour_palette)
WAD_raster <- create_raster_plot(mature_pred_WAD, sampling_area, SZ, aus, colour_palette)


# Arrange raster plots (Bottom Row) with individual legends
raster_plots <- ggarrange(
  PS_raster + annotate("text", 
                       x = ext(sampling_area)[1], 
                       y = ext(sampling_area)[4], 
                       label = "Yijarup/Pink \nSnapper", 
                       hjust = -0.1, vjust = 1.1,
                       size = 4), 
  SMW_raster + annotate("text", 
                        x = ext(sampling_area)[1], 
                        y = ext(sampling_area)[4], 
                        label = "Southern \nMaori \nWrasse", 
                        hjust = -0.1, vjust = 1.1,
                        size = 4), 
  WAD_raster + annotate("text", 
                        x = ext(sampling_area)[1], 
                        y = ext(sampling_area)[4], 
                        label = "Djubitj/West \nAustralian \nDhufish", 
                        hjust = -0.1, vjust = 1.1,
                        size = 4), 
  nrow = 1, widths = c(1, 1, 1), 
  labels = c("G.", "H.", "I.")
  )
save_plot("plots/species_prediction_plot.png", raster_plots, base_width = 8, base_height = 2.25)


## final plot -----------------------------------------------------------------

sampling_plot <- image_read("plots/example_sampling_designs.png")
raster_plot <- image_read("plots/species_prediction_plot.png")

combined_image <- image_append(c(sampling_plot, raster_plot), stack = TRUE)
image_write(combined_image, path = "plots/example_sampling_design_predicted_abundance.png", format = "png")


## single species performance -------------------------------------------------

species_lookup <- data.frame( # lookup table so that different dataframes recognise the species names
  species_formatted = c(
    "Yijarup/Pink Snapper",
    "Djubitj/West Australian Dhufish",
    "Southern Maori Wrasse"
  ),
  species = c(
    "Chrysophrys_auratus",
    "Glaucosoma_hebraicum",
    "Ophthalmolepis_lineolatus"
  )
)

raw <- readRDS("data/rmd/L06_sampling_design_performance_raw_data.rds") %>% 
  left_join(species_lookup, by = "species") %>% 
  glimpse()
raw$species_formatted <- factor(raw$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

sig <- readRDS("data/rmd/L06_significance_testing_model_results.rds") %>% 
  left_join(species_lookup, by = "species") %>%
  dplyr::select(c(SD, design_id, p_value, estimate, species_formatted)) %>% 
  glimpse()
sig$species_formatted <- factor(sig$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

mb <- readRDS("data/rmd/L06_sampling_design_performance_mean_bias.rds") %>% 
  left_join(species_lookup, by = "species") %>%
  glimpse()
mb$species_formatted <- factor(mb$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

rmse <- readRDS("data/rmd/L06_sampling_design_performance_rmse.rds") %>% 
  left_join(species_lookup, by = "species") %>%
  glimpse()
rmse$species_formatted <- factor(rmse$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

int_cov <- readRDS("data/rmd/L06_sampling_design_performance_interval_coverage.rds") %>% 
  left_join(species_lookup, by = "species") %>%
  glimpse()
int_cov$species_formatted <- factor(int_cov$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

realised_ab <- readRDS(paste0("data/rmd/sampling_design_performance_realised_abundances_all_lengths.rds"))
realised_long <- realised_abundance_df %>%
  dplyr::select(species, #species_formatted = species_formatted.x,
                SZ = mean_abundance_inside, fished = mean_abundance_outside) %>%
  pivot_longer(cols = c(SZ, fished),
               names_to = "status",
               values_to = "realised_abundance_mean") %>% 
  left_join(species_lookup, by = "species") %>% 
  glimpse()
realised_long$species_formatted <- factor(realised_long$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))






### raw data ------------------------------------------------------------------
p_raw <- ggplot(raw, aes(x = SD, y = obs_ab_mean, fill = status)) +
  geom_boxplot(width = 1, color = "black", alpha = 0.5, position = position_dodge(0.9), outlier.size = 0.5) +
  
  # Add horizontal dashed lines for realised abundance
  geom_hline(data = realised_long,
             aes(yintercept = realised_abundance_mean,
                 color = status,
                 group = interaction(species, status)),
             linetype = "dashed",
             size = 0.8) +
  facet_wrap(~ species_formatted, scales = "free_y") +
  labs(
    fill = "Status",
    x = "Sampling Design",
    y = "Observed<br>Abundance"
  ) +
  scale_fill_manual(values = colour_palette[1:2], labels = c("fished", "NTZ")) +
  scale_color_manual(values = colour_palette[1:2], guide = "none") +  # hide redundant legend
  custom_theme +
  theme(legend.position = "none",
        axis.text.x = element_blank(),  # Removes only x-axis text
        axis.title.x = element_blank()  # Removes x-axis title
  )


### mean bias -----------------------------------------------------------------

p_mb <- ggplot(mb, aes(x = SD, y = mean_bias, fill = SD)) +
  geom_boxplot(position = position_dodge(0.9), outlier.size = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = colour_palette[1], size = 1) +  
  facet_grid(cols = vars(species_formatted)) +
  labs(y = "Mean Bias") +
  custom_theme +
  theme(
    axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank(),  # Removes x-axis title
    strip.text.x = element_blank()  # Removes column facet titles (species names
  ) +
  scale_fill_manual(values = gradient_colors) +
  scale_color_manual(values = c("white" = colour_palette[4], "black" = "black")) +
  coord_cartesian(ylim = c(NA, 2))

## rmse -----------------------------------------------------------------------

p_rmse <- ggplot(rmse, aes(x = SD, y = rmse, fill = SD)) +
  geom_boxplot(position = position_dodge(0.9), outlier.size = 0.5) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = colour_palette[1], size = 1) +  
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
  coord_cartesian(ylim = c(0, 2))

## interval coverage ----------------------------------------------------------

p_ic <- ggplot(int_cov, aes(x = SD, y = int_coverage, fill = SD)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~species) +
  labs(y = "Coverage") +
  custom_theme +
  theme(
    #axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank(),  # Removes x-axis title
    strip.text.x = element_blank()  # Removes column facet titles (species names
  ) +
  geom_hline(data = interval_coverage, 
             aes(yintercept = 0.95), linetype = "dashed", color = colour_palette[1], size = 1) +  
  scale_fill_manual(values = gradient_colors) +
  coord_cartesian(ylim = c(0.5, 1))



## p-value --------------------------------------------------------------------

sig <- readRDS("data/rmd/L06_significance_testing_model_results.rds") %>% 
  mutate(SD = recode(SD,
                     "simple_spabal" = "Si_SB",
                     "stratified_spabal_2030" = "St_SB_2030",
                     "stratified_spabal_2525" = "St_SB_2525",
                     "pref_2030" = "P_2030",
                     "pref_2525" = "P_2525",
                     "clustered" = "C"
                     )
  ) %>% 
  left_join(species_lookup, by = "species") %>% 
  glimpse()
sig$species_formatted <- factor(sig$species_formatted, levels = c("Yijarup/Pink Snapper", "Southern Maori Wrasse", "Djubitj/West Australian Dhufish"))

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
  geom_text(
    data = summary_stats,
    aes(x = SD, y = 1.05,
        label = paste0("n=", total, "\nsig.=", round(perc_under_0.05, 1), "%")),
    size = 3,
    angle = 90,
    hjust = 1,
    inherit.aes = FALSE
  ) +
  labs(
    title = "",
    x = "",
    y = "P-value for<br>SZ effect"
  ) +
  custom_theme +
  theme(
    axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank(),  # Removes x-axis title
    strip.text.x = element_blank()  # Removes column facet titles (species names
  ) +
  scale_fill_manual(values = gradient_colors) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.01, ymax = 0.05), fill = NA, linetype = "dashed", color = colour_palette[1], size = 1)


# model estimate --------------------------------------------------------------

test <- realised_abundance_df %>%
  mutate(ra = mean_abundance_outside/mean_abundance_inside)

test <- sig %>%
  left_join(test, by = "species")
head(test)
# Plot
p_est <- ggplot(test, aes(x = SD, y = estimate)) +
  geom_boxplot(outlier.shape = NA, varwidth = TRUE, aes(fill = SD)) +
  facet_wrap(~ species_formatted) +
  geom_hline(aes(yintercept = ra), linetype = "dashed", color = colour_palette[1], size = 1) +
  labs(
    x = "",
    y = "Estimate for<br>SZ effect"
  ) +
  custom_theme +
  scale_fill_manual(values = gradient_colors) +
  theme(
    axis.text.x = element_blank(),  # Removes only x-axis text
    axis.title.x = element_blank(),  # Removes x-axis title
    strip.text.x = element_blank()  # Removes column facet titles (species names
  ) +
  coord_cartesian(ylim = c(0.3, 1.25))


# all plots -------------------------------------------------------------------

p_raw / p_sig / p_est / p_mb / p_rmse / p_ic + plot_layout(heights = c(1, 1, 1, 1, 1, 1))  




### END ###
