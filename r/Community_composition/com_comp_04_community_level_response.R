# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (B. Compare communities detected)
# Data:    2024 BRUV data (MEGLAB and DBCA)
# Task:    Calculate metrics of community-level response to NTZ protection
# Author:  Lise Fournier-Carnoy
# Date:    December 2024

# -----------------------------------------------------------------------------

# The aim here is to compare two sampling designs in whether they find the same effects of NTZ at the community-level.
# I'll be using recommended metrics from Soykan et al. 2015 to assess community 

# -----------------------------------------------------------------------------

rm(list=ls())

library(tidyverse)
library(sf)
library(vegan)
library(e1071) # for skewness
library(purrr)

colour_palette <- readLines("chapter_colours.txt"); colour_palette <- eval(parse(text = colour_palette))
source("custom_theme.R")

## Load data and calculate metrics --------------------------------------------

dat <- readRDS('data/rmd/community_composition_analysis_maxn.rds') %>% 
  st_drop_geometry() %>% 
  glimpse()

# prepare data
meta_cols <- c("opcode", "depth", "status", "sd")
species_data <- dat %>% dplyr::select(-c(all_of(meta_cols), "id"))
meta_data <- dat %>% dplyr::select(c(all_of(meta_cols), -depth)) %>% 
  mutate(sd = ifelse(sd == "preferential", "Clustered", "Spatially Balanced"))
as.data.frame(colnames(dat))
dat_long <- bind_cols(meta_data, species_data)

# compute metrics
compute_metrics <- function(df) {
  species <- df %>% dplyr::select(where(is.numeric))  # all numeric = species data
  
  `Total abundance` <- sum(species)
  
  species_vector <- colSums(species)
  species_nonzero <- species_vector[species_vector > 0]
  
  # Log-normal mu
  `log mu` <- mean(log1p(species_nonzero))
  
  # Log skew
  `log skew` <- skewness(log1p(species_nonzero))

  `Species richness` <- length(species_nonzero)
  
  tibble(
    `Total abundance` = `Total abundance`,
    `log mu` = `log mu`,
    `log skew` = `log skew`
  )
}

metrics <- dat_long %>%
  group_by(sd, status) %>%
  group_modify(~ compute_metrics(.x)) %>% 
  pivot_longer(cols = -c(sd, status), names_to = "metric", values_to = "value") %>% 
  glimpse()

metrics_long_scaled <- metrics %>%
  group_by(metric) %>%
  mutate(value_scaled = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

comparison <- metrics %>%
  pivot_wider(names_from = status, values_from = value) %>%
  mutate(
    ratio = `No-take` / Fished,
    percent_change = (ratio - 1) * 100,
    ratio_label = case_when(
      percent_change > 0 ~ paste0("+", round(percent_change, 0), "%"),
      percent_change < 0 ~ paste0(round(percent_change, 0), "%"),
      TRUE ~ "No change"
    )
  )

metrics_with_ratios <- metrics_long_scaled %>%
  left_join(comparison %>% dplyr::select(sd, metric, ratio_label), by = c("sd", "metric")) %>% 
  glimpse()

# Species rank plot prep
dat_long <- dat %>%
  pivot_longer(
    cols = starts_with("Apogonidae"):last_col(),  # adjust if needed
    names_to = "species",
    values_to = "abundance"
  )
sad_summary <- dat_long %>%
  group_by(sd, status, species) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop")
sad_summary <- sad_summary %>%
  filter(total_abundance > 0)
sad_ranked <- sad_summary %>%
  group_by(sd, status) %>%
  arrange(desc(total_abundance)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

## Plot -----------------------------------------------------------------------

tiles <- ggplot(metrics_with_ratios, aes(x = status, y = metric, fill = value_scaled)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 4)), color = "white", size = 5) +  # Original values
  facet_wrap(~ sd) +
  scale_fill_gradient(low = colour_palette[3], high = colour_palette[2]) +
  custom_theme +
  labs(
    title = "",
    x = "Status", y = "Metric"
  ) +
  theme(strip.text = element_text(size = 12)) +
  geom_label(
    data = subset(metrics_with_ratios, status == "No-take"),
    aes(label = ratio_label),
    vjust = -0.5,
    fill = "white",        # Background of the label
    color = "red",         # Text color
    label.size = 0.5,      # Thickness of label border
    label.r = unit(0.1, "lines"),  # Border radius
    size = 4,
    fontface = "italic"
  )

species_rank <- ggplot(sad_ranked, aes(x = rank, y = total_abundance, color = status, group = interaction(status, sd))) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, NA)) +
  scale_color_manual(values = colour_palette[c(1, 3)]) + 
  labs(
    x = "Species Rank",
    y = "log(total abundance)",
    color = "Status"
  ) +
  facet_wrap(~ sd, labeller = function(x) "") +
custom_theme +
  theme(
    legend.position = "right",
    strip.text = element_blank()
  )

library(patchwork)
(tiles / species_rank) +
  plot_annotation(
    tag_levels = 'A',                     
    tag_prefix = "",
    tag_suffix = ".",
    tag_sep = " ",
    )




### Species accumulation curve ------------------------------------------------

file_analysis_box <- "QGIS layers/polys/com_comp_analysis_limits.shp"
file_all_maxn     <- "data/tidy/2024_geographe_all_tidy_maxn.csv"

# Box around East Geographe, we're only using this area to compare community.
crs <- st_crs(4326)
box <- st_read(file_analysis_box) %>%
  st_make_valid() %>%
  st_transform(crs)
plot(box)


# MaxN data

maxn <- read.csv(file_all_maxn) %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  glimpse()
maxn <- st_as_sf(maxn, coords = c("longitude", "latitude"), crs = 4326); plot(maxn)

# Select only points within community composition box
species <- st_intersection(maxn, box); plot(maxn)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

# 2. Split by sd group
sd_list <- split(species, list(species$sd, species$status), drop = TRUE)

# 3. Function for species accumulation curve per sd group
accum_list <- lapply(names(sd_list), function(sd_name) {
  sd_data <- sd_list[[sd_name]]
  
  # Sample x species abundance matrix
  mat <- sd_data %>%
    st_drop_geometry() %>%
    dplyr::select(opcode, fullspp, MaxN) %>%
    group_by(opcode, fullspp) %>%
    summarise(MaxN = sum(MaxN), .groups = "drop") %>%
    pivot_wider(names_from = fullspp,
                values_from = MaxN,
                values_fill = list(MaxN = 0),
                values_fn = sum)
  
  # Set rownames and drop opcode column
  rownames(mat) <- mat$opcode
  mat <- mat[, !names(mat) %in% "opcode"]
  
  # Flatten any list-columns & force numeric
  mat[] <- lapply(mat, function(x) as.numeric(unlist(x)))
  
  # Confirm structure before passing to vegan
  stopifnot(all(sapply(mat, is.numeric)))
  
  # Run species accumulation
  spec_accum <- specaccum(mat, method = "random", permutations = 100)
  
  # Format for plotting
  data.frame(
    Sites = spec_accum$sites,
    Richness = spec_accum$richness,
    SD = spec_accum$sd,
    sd_group = sd_name
  )
})

# Combine all groups
accum_data <- bind_rows(accum_list) %>%
  tidyr::separate(sd_group, into = c("design", "status"), sep = "\\.")

accum_data$design <- ifelse(accum_data$design == "preferential", "clustered", accum_data$design)

# Plot
sac <- ggplot(accum_data, aes(x = Sites, y = Richness)) +
  geom_line(aes(color = design, linetype = status), size = 1) +
  #facet_grid(~design) +
  geom_ribbon(aes(
    ymin = Richness - SD,
    ymax = Richness + SD,
    fill = design,
    group = interaction(design, status)
  ), alpha = 0.2, color = NA) +
  
  # Custom colors and line types
  scale_color_manual(values = c(
    "clustered" = colour_palette[4],
    "spatially balanced" = colour_palette[6]
  )) +
  scale_fill_manual(values = c(
    "clustered" = colour_palette[4],
    "spatially balanced" = colour_palette[6]
  )) +
  scale_linetype_manual(values = c(
    "Fished" = "dotted",
    "No-take" = "solid"
  )) +
  
  labs(
    title = "",
    x = "Number of samples",
    y = "Cumulative Species Richness",
    color = "Design",
    fill = "Design",
    linetype = "Status"
  ) +
  custom_theme +
  theme(legend.position = "top")

### END ###


glimpse(sad_ranked)

species_rank <- ggplot(sad_ranked, aes(x = rank, y = total_abundance, color = status, group = interaction(status, sd))) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, NA)) +
  scale_color_manual(values = colour_palette[c(1, 3)]) + 
  labs(
    x = "Species Rank",
    y = "log(total abundance)",
    color = "Status"
  ) +
  facet_wrap(~ sd, labeller = function(x) "") +
  custom_theme +
  theme(
    legend.position = "right",
    strip.text = element_blank()
  )


sort(unique(sad_ranked$species))
highlight_species <- subset(sad_ranked, species %in% c(
  "Sparidae.Chrysophrys.auratus",
  "Glaucosomatidae.Glaucosoma.hebraicum",
  "Carangidae..spp",
  "Carangidae.Seriola.hippos",
  "Cheilodactylidae.Nemadactylus.valenciennesi",
  "Sillaginidae.Sillaginodes.punctatus"
)) %>% 
  mutate(label_spp = dplyr::case_when(
    species == "Sparidae.Chrysophrys.auratus" ~ "yijarup",
    species == "Glaucosomatidae.Glaucosoma.hebraicum" ~ "dhubitj",
    species == "Carangidae..spp" ~ "trevally",
    species == "Carangidae.Seriola.hippos" ~ "yellowtail",
    species == "Cheilodactylidae.Nemadactylus.valenciennesi" ~ "blue morwong",
    species == "Sillaginidae.Sillaginodes.punctatus" ~ "king george whiting"
  ))

# Create the plot
ggplot(sad_ranked, aes(x = rank, y = total_abundance, color = status, group = interaction(status, sd))) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, NA)) +
  scale_color_manual(values = colour_palette[c(1, 3)]) + 
  labs(
    x = "Species Rank",
    y = "log(total abundance)",
    color = "Status"
  ) +
  facet_wrap(~ sd, labeller = function(x) "") +
  custom_theme +
  theme(
    legend.position = "right",
    strip.text = element_blank()
  ) +

  geom_text(
    data = highlight_species,
    aes(
      x = rank + 10,
      y = total_abundance + 10,
      label = label_spp,
      color = status,
      hjust = 0,
      angle = 45
    ),
    size = 4,
    inherit.aes = FALSE
  ) +
  
  # Add connecting line
  geom_segment(
    data = highlight_species,
    aes(
      x = rank,
      xend = rank +10,
      y = total_abundance,
      yend = total_abundance + 10,
      color = status
    ),
    linewidth = 0.5,
    inherit.aes = FALSE
  )
