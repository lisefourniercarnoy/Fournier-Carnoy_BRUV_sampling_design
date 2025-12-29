# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (B. Compare communities detected)
# Data:    2024 BRUV data (MEGLAB and DBCA)
# Task:    Calculate metrics of community-level response to NTZ protection
# Author:  Lise Fournier-Carnoy
# Date:    December 2024

# -----------------------------------------------------------------------------

# Important: this is an archive.

rm(list=ls())

library(tidyverse)
library(sf)
library(vegan)
library(e1071) # for skewness
library(purrr)
library(patchwork) # to organise plots

colour_palette <- readLines("chapter_colours.txt"); colour_palette <- eval(parse(text = colour_palette))
source("custom_theme.R")

study_site <- "waatern" # waatu or waatern

## Load data ------------------------------------------------------------------

file_analysis_box <- paste0("QGIS layers/clean/", study_site, "_community_composition_analysis_limits.shp")
file_maxn_wide    <- paste0("data/tidy/C01_2024_", study_site, "_all_tidy_maxn.rds")

dat <- readRDS(file_maxn_wide) %>% 
  st_drop_geometry() %>% 
  glimpse()

dat_long <- dat %>%
  st_drop_geometry() %>%
  select(opcode, fullspp, maxn, sd, status)

# Pivot wider: samples × species
dat_wide <- dat_long %>%
  pivot_wider(names_from = fullspp, values_from = maxn, values_fill = 0)

dat_grouped <- dat_wide %>%
  group_by(sd, status) %>%
  nest()





## SPECIES RANK ---------------------------------------------------------------

sad_summary <- dat %>%
  group_by(sd, status, fullspp) %>%
  summarise(total_abundance = sum(maxn, na.rm = TRUE), .groups = "drop")

sad_summary <- sad_summary %>%
  filter(total_abundance > 0)
sad_ranked <- sad_summary %>%
  group_by(sd, status) %>%
  arrange(desc(total_abundance)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

waatern_species_rank <- ggplot(sad_ranked, aes(x = rank, y = total_abundance, color = status, group = interaction(status, sd))) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, NA)) +
  scale_color_manual(values = colour_palette[c(1, 3)]) + 
  labs(
    x = "Species Rank",
    y = "log(total abundance)",
    color = ""
  ) +
  facet_wrap(~ sd, labeller = function(x) "") +
  custom_theme +
  theme(
    legend.position = "right",
    strip.text = element_blank()
  )

## METRIC TILES ---------------------------------------------------------------

metrics <- dat_long %>%
  group_by(sd, status) %>%
  summarise(
    log_mu = mean(log1p(maxn), na.rm = TRUE),             # log1p avoids log(0) issues
    log_skew = skewness(log1p(maxn), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = -c(sd, status), names_to = "metric", values_to = "value")

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
metrics_with_ratios <- metrics_with_ratios %>%
  mutate(
    label_col = if_else(str_detect(ratio_label, "^\\+") & metric == 'log_mu', "darkgreen", "darkred")
  )
metrics_with_ratios$sd <- ifelse(metrics_with_ratios$sd == "preferential", "Targeted", "Spatially balanced")

waatern_tiles <- ggplot(metrics_with_ratios, aes(x = status, y = metric, fill = value_scaled)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 4)), color = "white", size = 5) +  # Original values
  facet_wrap(~ sd) +
  scale_fill_gradient(low = colour_palette[3], high = colour_palette[2]) +
  custom_theme +
  labs(
    title = "",
    x = "", y = "Metric"
  ) +
  theme(strip.text = element_text(size = 12)) +
  geom_label(
    data = metrics_with_ratios %>% filter(status == "No-take"),
    aes(label = ratio_label, color = I(label_col)),  # Force actual color strings
    vjust = -0.5,
    fill = "white",
    label.size = 0.5,
    label.r = unit(0.1, "lines"),
    size = 4,
    fontface = "italic",
    show.legend = FALSE
  )


## SAC ------------------------------------------------------------------------
crs <- st_crs(4326)
box <- st_read(file_analysis_box) %>%
  st_make_valid() %>%
  st_transform(crs)
plot(box$geometry)

maxn <- readRDS(file_maxn_wide) %>%
  glimpse()
if (study_site == "waatern") {st_crs(maxn) <- crs}

maxn <- st_as_sf(maxn, coords = c("longitude", "latitude"), crs = crs) %>% st_intersection(box); plot(maxn$geometry, add=T)

# split by sampling design group
sd_list <- split(maxn, list(as.factor(maxn$sd), as.factor(maxn$status)), drop = TRUE)

# run the function for species accumulation curve per sd group
accum_list <- lapply(names(sd_list), function(sd_name) {
  sd_data <- sd_list[[sd_name]]
  
  # sample x species abundance matrix
  mat <- sd_data %>%
    st_drop_geometry() %>%
    dplyr::select(opcode, fullspp, maxn) %>%
    group_by(opcode, fullspp) %>%
    summarise(maxn = sum(maxn), .groups = "drop") %>%
    pivot_wider(names_from = fullspp,
                values_from = maxn,
                values_fill = list(maxn = 0),
                values_fn = sum)
  
  # drop opcode column
  mat <- mat[, !names(mat) %in% "opcode"]
  
  # flatten any list-columns & force numeric
  mat[] <- lapply(mat, function(x) as.numeric(unlist(x)))
  
  # confirm structure before passing to vegan
  stopifnot(all(sapply(mat, is.numeric)))
  
  # run species accumulation
  spec_accum <- specaccum(mat, method = "random", permutations = 100)
  
  # format for plotting
  data.frame(
    Sites = spec_accum$sites,
    Richness = spec_accum$richness,
    SD = spec_accum$sd,
    sd_group = sd_name
  )
})

# combine all groups
accum_data <- bind_rows(accum_list) %>%
  tidyr::separate(sd_group, into = c("design", "status"), sep = "\\.")

accum_data$design <- ifelse(accum_data$design == "preferential", "Targeted", "Spatially balanced")

# Plot the SAC
waatern_sac <- ggplot(accum_data, aes(x = Sites, y = Richness)) +
  geom_line(aes(color = design, linetype = status), size = 1) +
  geom_ribbon(aes(
    ymin = Richness - SD,
    ymax = Richness + SD,
    fill = design,
    group = interaction(design, status)
  ), alpha = 0.2, color = NA) +
  
  # Custom colors and line types
  scale_color_manual(values = c(
    "Targeted" = colour_palette[4],
    "Spatially balanced" = colour_palette[6]
  )) +
  scale_fill_manual(values = c(
    "Targeted" = colour_palette[4],
    "Spatially balanced" = colour_palette[6]
  )) +
  scale_linetype_manual(values = c(
    "Fished" = "dotted",
    "No-take" = "solid"
  )) +
  
  labs(
    title = "",
    x = "Number of samples",
    y = "Cumulative\nSpecies Richness",
    color = "",
    fill = "",
    linetype = ""
  ) +
  custom_theme +
  theme(legend.position = "top")


## WAATU ----------------------------------------------------------------------

study_site <- "waatu" # waatu or waatern

## Load data ------------------------------------------------------------------

file_analysis_box <- paste0("QGIS layers/clean/", study_site, "_community_composition_analysis_limits.shp")
file_maxn_wide    <- paste0("data/tidy/C01_2024_", study_site, "_all_tidy_maxn.rds")

dat <- readRDS(file_maxn_wide) %>% 
  st_drop_geometry() %>% 
  glimpse()

dat_long <- dat %>%
  st_drop_geometry() %>%
  select(opcode, fullspp, maxn, sd, status)

# Pivot wider: samples × species
dat_wide <- dat_long %>%
  pivot_wider(names_from = fullspp, values_from = maxn, values_fill = 0)

dat_grouped <- dat_wide %>%
  group_by(sd, status) %>%
  nest()





## SPECIES RANK ---------------------------------------------------------------

sad_summary <- dat %>%
  group_by(sd, status, fullspp) %>%
  summarise(total_abundance = sum(maxn, na.rm = TRUE), .groups = "drop")

sad_summary <- sad_summary %>%
  filter(total_abundance > 0)
sad_ranked <- sad_summary %>%
  group_by(sd, status) %>%
  arrange(desc(total_abundance)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

waatu_species_rank <- ggplot(sad_ranked, aes(x = rank, y = total_abundance, color = status, group = interaction(status, sd))) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, NA)) +
  scale_color_manual(values = colour_palette[c(1, 3)]) + 
  labs(
    x = "Species Rank",
    y = "log(total abundance)",
    color = ""
  ) +
  facet_wrap(~ sd, labeller = function(x) "") +
  custom_theme +
  theme(
    legend.position = "right",
    strip.text = element_blank()
  )

## METRIC TILES ---------------------------------------------------------------

metrics <- dat_long %>%
  group_by(sd, status) %>%
  summarise(
    log_mu = mean(log1p(maxn), na.rm = TRUE),             # log1p avoids log(0) issues
    log_skew = skewness(log1p(maxn), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = -c(sd, status), names_to = "metric", values_to = "value")

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
metrics_with_ratios <- metrics_with_ratios %>%
  mutate(
    label_col = if_else(str_detect(ratio_label, "^\\+") & metric == 'log_mu', "darkgreen", "darkred")
  )
metrics_with_ratios$sd <- ifelse(metrics_with_ratios$sd == "preferential", "Targeted", "Stratified<br>spatially<br>balanced")
waatu_tiles <- ggplot(metrics_with_ratios, aes(x = status, y = metric, fill = value_scaled)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 4)), color = "white", size = 4) +  # Original values
  facet_wrap(~ sd) +
  scale_fill_gradient(low = colour_palette[3], high = colour_palette[2]) +
  custom_theme +
  labs(
    title = "",
    x = "", y = "Metric"
  ) +
  theme(strip.text = element_text(size = 11)) +
  geom_label(
    data = metrics_with_ratios %>% filter(status == "No-take"),
    aes(label = ratio_label, color = I(label_col)),  # Force actual color strings
    vjust = -0.5,
    fill = "white",
    label.size = 0.5,
    label.r = unit(0.1, "lines"),
    size = 4,
    fontface = "italic",
    show.legend = FALSE
  )


## SAC ------------------------------------------------------------------------
crs <- st_crs(4326)
box <- st_read(file_analysis_box) %>%
  st_make_valid() %>%
  st_transform(crs)
plot(box$geometry)

maxn <- readRDS(file_maxn_wide) %>%
  glimpse()
if (study_site == "waatern") {st_crs(maxn) <- crs}

maxn <- st_as_sf(maxn, coords = c("longitude", "latitude"), crs = crs) %>% st_intersection(box); plot(maxn$geometry, add=T)

# split by sampling design group
sd_list <- split(maxn, list(as.factor(maxn$sd), as.factor(maxn$status)), drop = TRUE)

# run the function for species accumulation curve per sd group
accum_list <- lapply(names(sd_list), function(sd_name) {
  sd_data <- sd_list[[sd_name]]
  
  # sample x species abundance matrix
  mat <- sd_data %>%
    st_drop_geometry() %>%
    dplyr::select(opcode, fullspp, maxn) %>%
    group_by(opcode, fullspp) %>%
    summarise(maxn = sum(maxn), .groups = "drop") %>%
    pivot_wider(names_from = fullspp,
                values_from = maxn,
                values_fill = list(maxn = 0),
                values_fn = sum)
  
  # drop opcode column
  mat <- mat[, !names(mat) %in% "opcode"]
  
  # flatten any list-columns & force numeric
  mat[] <- lapply(mat, function(x) as.numeric(unlist(x)))
  
  # confirm structure before passing to vegan
  stopifnot(all(sapply(mat, is.numeric)))
  
  # run species accumulation
  spec_accum <- specaccum(mat, method = "random", permutations = 100)
  
  # format for plotting
  data.frame(
    Sites = spec_accum$sites,
    Richness = spec_accum$richness,
    SD = spec_accum$sd,
    sd_group = sd_name
  )
})

# combine all groups
accum_data <- bind_rows(accum_list) %>%
  tidyr::separate(sd_group, into = c("design", "status"), sep = "\\.")

accum_data$design <- ifelse(accum_data$design == "preferential", "Targeted", "Spatially balanced")

# Plot the SAC
waatu_sac <- ggplot(accum_data, aes(x = Sites, y = Richness)) +
  geom_line(aes(color = design, linetype = status), size = 1) +
  geom_ribbon(aes(
    ymin = Richness - SD,
    ymax = Richness + SD,
    fill = design,
    group = interaction(design, status)
  ), alpha = 0.2, color = NA) +
  
  # Custom colors and line types
  scale_color_manual(values = c(
    "Targeted" = colour_palette[4],
    "Spatially balanced" = colour_palette[6]
  )) +
  scale_fill_manual(values = c(
    "Targeted" = colour_palette[4],
    "Spatially balanced" = colour_palette[6]
  )) +
  scale_linetype_manual(values = c(
    "Fished" = "dotted",
    "No-take" = "solid"
  )) +
  
  labs(
    title = "",
    x = "Number of samples",
    y = "Cumulative\nSpecies Richness",
    color = "",
    fill = "",
    linetype = ""
  ) +
  custom_theme +
  theme(legend.position = "top")



## both study sites -----------------------------------------------------------

# run the above script first
sac <- (
  waatu_sac + theme(legend.position = "right") |
    waatern_sac + theme(legend.position = "right",
                        axis.title.y = element_blank())
) +
  plot_layout(guides = "collect", axis_titles = "collect")
  #plot_annotation(title =)
  theme(legend.position = "right") &
  theme(plot.margin = margin(rep(0, 4)))
sac 


tiles <- (
  waatu_tiles + theme(axis.title.x = element_blank()) |
    waatern_tiles + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank())) +
  plot_layout(axis_titles = "collect")
tiles




rank <- (
  waatu_species_rank + theme(legend.position = "right") |
    waatern_species_rank + theme(legend.position = "right",
                                 axis.title.y = element_blank())
) +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "right") &
  theme(plot.margin = margin(rep(0, 4)))  # Set position for the shared legend

all <- (tiles / rank)

all

save_plot("plots/assemblage_metrics.png", all, base_width = 7, base_height = 6)


### test  -----

# testing package 'rich' from Rossi 2011
#install.packages("C:\\Users\\localuser\\Downloads\\rich_1.0.1.tar.gz", repos = NULL, type = "source")
library(rich)
study_site <- "waatern" 
# the tests require each group to be in a separate matrix, with rows as samples and columns as species
file_analysis_box <- paste0("QGIS layers/clean/", study_site, "_community_composition_analysis_limits.shp")
file_maxn_wide    <- paste0("data/tidy/C01_2024_", study_site, "_all_tidy_maxn.rds")
box <- st_read(file_analysis_box) %>%
  st_make_valid() %>%
  st_transform(crs)
plot(box$geometry)
dat <- readRDS(file_maxn_wide) %>% 
  #st_drop_geometry() %>% 
  glimpse()
dat_long <- dat %>%
  #st_drop_geometry() %>%
  select(opcode, fullspp, maxn, sd, status)
dat_wide <- dat_long %>%
  pivot_wider(names_from = fullspp, values_from = maxn, values_fill = 0)
plot(box$geometry)
st_crs(dat_wide) <- st_crs(box)
dat_wide <- dat_wide %>% st_intersection(box); plot(dat_wide$geometry, add=T)

# separate each group - starting with just sampling design
sb_f <- ungroup(dat_wide) %>% 
  filter(sd == "spatially balanced" & status == "Fished") %>% 
  dplyr::select(-c(opcode, sd, status, id)) %>% 
  st_drop_geometry()
sb_f <- as.matrix(sb_f)

t_f <- ungroup(dat_wide) %>% 
  filter(sd == "preferential" & status == "Fished") %>% 
  dplyr::select(-c(opcode, sd, status, id)) %>% 
  st_drop_geometry()
t_f <- as.matrix(t_f)

sb_nt <- ungroup(dat_wide) %>% 
  filter(sd == "spatially balanced" & status != "Fished") %>% 
  dplyr::select(-c(opcode, sd, status, id)) %>% 
  st_drop_geometry()
sb_nt <- as.matrix(sb_nt)

t_nt <- ungroup(dat_wide) %>% 
  filter(sd == "preferential" & status != "Fished") %>% 
  dplyr::select(-c(opcode, sd, status, id)) %>% 
  st_drop_geometry()
t_nt <- as.matrix(t_nt)

test_dat <- list("sb_f" = sb_f, "sb_nt" = sb_nt, "t_f" = t_f, "t_nt" = t_nt)
summary(test_dat)
names(test_dat)

# targeted observed species richness
test <- rich(matrix = test_dat$t_f, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished targeted samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$t_f), "fished targeted samples :", test$mr)

test <- rich(matrix = test_dat$t_nt, nrandom = 499, verbose = T)
cat("observed cumulative species richness in no-take targeted samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$t_nt), "no-take targeted samples :", test$mr)


# spatially balanced observed species richness
test <- rich(matrix = test_dat$sb_f, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished spatially balanced samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$sb_f), "targeted samples :", test$mr)

test <- rich(matrix = test_dat$sb_nt, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished spatially balanced samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$sb_nt), "targeted samples :", test$mr)


# stat test to see whether the number of species is sig. different between groups
c2cv(com2 = test_dat$t_f,
     com1 = test_dat$t_nt,
     nrandom = 999, 
     verbose = FALSE
     )
c2cv(com1 = test_dat$sb_f, 
     com2 = test_dat$sb_nt,
     nrandom = 999, 
     verbose = FALSE
)

curve_sbf  <- rarc(matrix = test_dat$sb_f, nrandom = 99)$out  %>% mutate(group = "SB - Fished")
curve_sbnt <- rarc(matrix = test_dat$sb_nt, nrandom = 99)$out %>% mutate(group = "SB - No-take")
curve_tf   <- rarc(matrix = test_dat$t_f, nrandom = 99)$out   %>% mutate(group = "Targeted - Fished")
curve_tnt  <- rarc(matrix = test_dat$t_nt, nrandom = 99)$out  %>% mutate(group = "Targeted - No-take")

# Combine all into one dataframe
curve_all <- bind_rows(curve_sbf, curve_sbnt, curve_tf, curve_tnt)

ggplot(curve_all, aes(x = sample, y = mean.richness, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lb.richness, ymax = ub.richness), alpha = 0.2, color = NA) +
  labs(
    x = "Number of Samples",
    y = "Species Richness",
    color = "Sampling Design",
    fill = "Sampling Design"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


## JUST SMAPLIGN DESIGN


# separate each group - starting with just sampling design
sb <- ungroup(dat_wide) %>% 
  filter(sd == "spatially balanced") %>% 
  dplyr::select(-c(opcode, sd, status, id)) %>% 
  st_drop_geometry()
sb <- as.matrix(sb)

t <- ungroup(dat_wide) %>% 
  filter(sd == "preferential") %>% 
  dplyr::select(-c(opcode, sd, status, id)) %>% 
  st_drop_geometry()
t <- as.matrix(t)

test_dat <- list("sb" = sb, "t" = t)
summary(test_dat)
names(test_dat)

# targeted observed species richness
test <- rich(matrix = test_dat$t, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished targeted samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$t), "fished targeted samples :", test$mr)


# spatially balanced observed species richness
test <- rich(matrix = test_dat$sb, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished spatially balanced samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$sb), "targeted samples :", test$mr)




c2cv(com1 = test_dat$t, 
     com2 = test_dat$sb,
     nrandom = 999, 
     verbose = FALSE
)


curve_sb  <- rarc(matrix = test_dat$sb, nrandom = 99)$out  %>% mutate(group = "spatially balanced")
curve_t <- rarc(matrix = test_dat$t, nrandom = 99)$out %>% mutate(group = "targeted")

# Combine all into one dataframe
curve_all <- bind_rows(curve_sb, curve_t)

ggplot(curve_all, aes(x = sample, y = mean.richness, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lb.richness, ymax = ub.richness), alpha = 0.2, color = NA) +
  labs(
    x = "Number of Samples",
    y = "Species Richness",
    color = "Sampling Design",
    fill = "Sampling Design"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

## WAATu TEST

study_site <- "waatu" 
# the tests require each group to be in a separate matrix, with rows as samples and columns as species
file_analysis_box <- paste0("QGIS layers/clean/", study_site, "_community_composition_analysis_limits.shp")
file_maxn_wide    <- paste0("data/tidy/C01_2024_", study_site, "_all_tidy_maxn.rds")
box <- st_read(file_analysis_box) %>%
  st_make_valid() %>%
  st_transform(crs)
plot(box$geometry)
dat <- readRDS(file_maxn_wide) %>% 
  #st_drop_geometry() %>% 
  glimpse()
dat_long <- dat %>%
#  st_drop_geometry() %>%
  select(opcode, fullspp, maxn, sd, status)
dat_wide <- dat_long %>%
  pivot_wider(names_from = fullspp, values_from = maxn, values_fill = 0)
plot(box$geometry)
st_crs(dat_wide) <- 4326
dat_wide <- st_as_sf(dat_wide, coords = c("longitude", "latitude"), crs = 4326) %>% st_intersection(box); plot(dat_wide$geometry, add=T)

# separate each group - starting with just sampling design
sb_f <- ungroup(dat_wide) %>% 
  filter(sd == "spatially balanced" & status == "Fished") %>% 
  dplyr::select(-c(opcode, sd, status)) %>% 
  st_drop_geometry()
sb_f <- as.matrix(sb_f)

t_f <- ungroup(dat_wide) %>% 
  filter(sd == "preferential" & status == "Fished") %>% 
  dplyr::select(-c(opcode, sd, status)) %>% 
  st_drop_geometry()
t_f <- as.matrix(t_f)

sb_nt <- ungroup(dat_wide) %>% 
  filter(sd == "spatially balanced" & status != "Fished") %>% 
  dplyr::select(-c(opcode, sd, status)) %>% 
  st_drop_geometry()
sb_nt <- as.matrix(sb_nt)

t_nt <- ungroup(dat_wide) %>% 
  filter(sd == "preferential" & status != "Fished") %>% 
  dplyr::select(-c(opcode, sd, status)) %>% 
  st_drop_geometry()
t_nt <- as.matrix(t_nt)

test_dat <- list("sb_f" = sb_f, "sb_nt" = sb_nt, "t_f" = t_f, "t_nt" = t_nt)
summary(test_dat)
names(test_dat)

# targeted observed species richness
test <- rich(matrix = test_dat$t_f, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished targeted samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$t_f), "fished targeted samples :", test$mr)

test <- rich(matrix = test_dat$t_nt, nrandom = 499, verbose = T)
cat("observed cumulative species richness in no-take targeted samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$t_nt), "no-take targeted samples :", test$mr)


# spatially balanced observed species richness
test <- rich(matrix = test_dat$sb_f, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished spatially balanced samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$sb_f), "targeted samples :", test$mr)

test <- rich(matrix = test_dat$sb_nt, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished spatially balanced samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$sb_nt), "targeted samples :", test$mr)




c2cv(com2 = test_dat$t_f, # targeted
     com1 = test_dat$t_nt,
     nrandom = 999, 
     verbose = FALSE
)
c2cv(com1 = test_dat$sb_f, # spatially balanced
     com2 = test_dat$sb_nt,
     nrandom = 999, 
     verbose = FALSE
)

c2cv(com2 = test_dat$t_f, # fished
     com1 = test_dat$sb_f,
     nrandom = 999, 
     verbose = FALSE
)
c2cv(com1 = test_dat$t_nt, # no take
     com2 = test_dat$sb_nt,
     nrandom = 999, 
     verbose = FALSE
)

curve_sbf  <- rarc(matrix = test_dat$sb_f, nrandom = 99)$out  %>% mutate(group = "SB - Fished")
curve_sbnt <- rarc(matrix = test_dat$sb_nt, nrandom = 99)$out %>% mutate(group = "SB - No-take")
curve_tf   <- rarc(matrix = test_dat$t_f, nrandom = 99)$out   %>% mutate(group = "Targeted - Fished")
curve_tnt  <- rarc(matrix = test_dat$t_nt, nrandom = 99)$out  %>% mutate(group = "Targeted - No-take")

# Combine all into one dataframe
curve_all <- bind_rows(curve_sbf, curve_sbnt, curve_tf, curve_tnt)

ggplot(curve_all, aes(x = sample, y = mean.richness, color = group, fill = group)) +
  geom_line(size = 1) +
  #geom_ribbon(aes(ymin = lb.richness, ymax = ub.richness), alpha = 0.2, color = NA) +
  labs(
    x = "Number of Samples",
    y = "Species Richness",
    color = "Sampling Design",
    fill = "Sampling Design"
  ) +
  theme_minimal() +
  theme(legend.position = "top")



## JUST SMAPLIGN DESIGN


# separate each group - starting with just sampling design
sb <- ungroup(dat_wide) %>% 
  filter(sd == "spatially balanced") %>% 
  dplyr::select(-c(opcode, sd, status)) %>% 
  st_drop_geometry()
sb <- as.matrix(sb)

t <- ungroup(dat_wide) %>% 
  filter(sd == "preferential") %>% 
  dplyr::select(-c(opcode, sd, status)) %>% 
  st_drop_geometry()
t <- as.matrix(t)

test_dat <- list("sb" = sb, "t" = t)
summary(test_dat)
names(test_dat)

# targeted observed species richness
test <- rich(matrix = test_dat$t, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished targeted samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$t), "fished targeted samples :", test$mr)


# spatially balanced observed species richness
test <- rich(matrix = test_dat$sb, nrandom = 499, verbose = T)
cat("observed cumulative species richness in fished spatially balanced samples is:", test$cr)
cat("observed mean value of species richness over the", nrow(test_dat$sb), "targeted samples :", test$mr)







# rarefaction
# if the number of samples greatly differs, you can rescale the larger dataset to the smaller's size.
nrow(test_dat$sb)
nrow(test_dat$t)

c2cv(com1 = test_dat$t, 
     com2 = test_dat$sb,
     nrandom = 999, 
     verbose = FALSE
)


curve_sb  <- rarc(matrix = test_dat$sb, nrandom = 99)$out  %>% mutate(group = "spatially balanced")
curve_t <- rarc(matrix = test_dat$t, nrandom = 99)$out %>% mutate(group = "targeted")

# Combine all into one dataframe
curve_all <- bind_rows(curve_sb, curve_t)

ggplot(curve_all, aes(x = sample, y = mean.richness, color = group, fill = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lb.richness, ymax = ub.richness), alpha = 0.2, color = NA) +
  labs(
    x = "Number of Samples",
    y = "Species Richness",
    color = "Sampling Design",
    fill = "Sampling Design"
  ) +
  theme_minimal() +
  theme(legend.position = "top")
