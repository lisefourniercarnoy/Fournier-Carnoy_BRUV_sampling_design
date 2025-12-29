# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (B. Compare communities detected)
# Data:    2024 BRUV data (MEGLAB and DBCA)
# Task:    Make CAP plots
# Author:  Lise Fournier-Carnoy
# Date:    December 2024

# -----------------------------------------------------------------------------

rm(list=ls())

## Load libraries -------------------------------------------------------------

library(vegan) # for community composition
library(tidyverse) # for data manipulation
library(mvabund) # for community analysis
library(sf) # for manipulating shapefiles
library(patchwork)  # For combining plots
library(ggtext) # for element_markdown

colour_palette <- readLines("chapter_colours.txt"); colour_palette <- eval(parse(text = colour_palette))
source("custom_theme.R")


# Important: files here were generated with PRIMER.

## WAATERN --------------------------------------------------------------------

# load data
study_site <- "waatern" # waatu or waatern

spp_sd <- read.csv(paste0("data/PRIMER data/", study_site, "/Pearson correlation of CAP (sd) with community (CAP only).csv")) %>% glimpse()
spp_status <- read.csv(paste0("data/PRIMER data/", study_site, "/Pearson correlation of CAP (status) with community (CAP only).csv")) %>% glimpse()

CAP_sd <- read.csv(paste0("data/PRIMER data/", study_site, "/CAP_by_sd_from_Data1.csv")) %>% glimpse()
CAP_status <- read.csv(paste0("data/PRIMER data/", study_site, "/CAP_by_status_from_Data3.csv")) %>% glimpse()

# obtain ubiquitous species
dat <- readRDS("data/tidy/C01_2024_waatern_all_tidy_maxn.rds") %>% ungroup()
n_samples <- length(unique(dat$opcode))

ubi <- dat %>%
  filter(maxn > 0) %>%
  dplyr::select(opcode, fullspp) %>%
  distinct() %>%  # Unique species per sample
  count(fullspp, name = "n_present") %>%
  mutate(prop = n_present / n_samples) %>%
  filter(prop >= 0.10) %>%
  arrange(desc(prop))

ubi <- ubi %>%
  mutate(species = sub("^[^ ]+\\s+", "", fullspp)) # remove family name
ubi$species <- paste0("<i>", ubi$species, "</i>") # format
ubi

## CAP plot -------------------------------------------------------------------

# select ubiquitous spp.
spp_clean <- spp_sd %>% # get only the main species driving the community composition
  filter(!is.na(CAP1)) %>%
  arrange(desc(abs(CAP1))) %>%
  #slice_head(n = 15) %>% 
  mutate(
    species = str_replace(taxa, "^[^.]+\\.", ""),    # Remove family (everything up to first dot)
    species = str_replace_all(species, "\\.", " "),  # Replace dots with spaces
    species = paste0("<i>", species, "</i>")         # Wrap in <i> tags
  )
spp_clean <- spp_clean[spp_clean$species %in% ubi$species & abs(spp_clean$CAP1) > 0.4,]

# split species into top positive and negative associations
spp_pos <- spp_clean %>% filter(CAP1 > 0) %>% arrange(desc(CAP1))
spp_neg <- spp_clean %>% filter(CAP1 < 0) %>% arrange(CAP1)

# 1. Combine all CAP1 values to find shared x-axis limits
x_limits <- range(
  c(spp_pos$CAP1, spp_neg$CAP1, CAP_status$CAP1),
  na.rm = TRUE
)

p_top <- ggplot(spp_neg, aes(y = reorder(species, CAP1), x = 0, xend = CAP1, yend = reorder(species, CAP1))) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), color = colour_palette[1], lwd = 1.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(x_limits) +
  labs(x = "", y = "") +
  custom_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_markdown())

CAP_sd$sd <- ifelse(CAP_sd$sd == "preferential", 'targeted', CAP_sd$sd)
CAP_sd$sd <- paste0("<b>", CAP_sd$sd, "</b>")

p_mid <- ggplot(CAP_sd, aes(x = CAP1, y = sd)) +
  geom_point(aes(color = sd, shape = status), size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(x_limits) +
  scale_color_manual(values = c(
    "<b>spatially balanced</b>" = colour_palette[1],
    "<b>targeted</b>" = colour_palette[4]
  )) +
  labs(x = "", y = "") +
  custom_theme +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_markdown())

p_bot <- ggplot(spp_pos, aes(y = reorder(species, CAP1), x = 0, xend = CAP1, yend = reorder(species, CAP1))) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), color = colour_palette[2], lwd = 1.25) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(x_limits) +
  labs(x = "CAP axis (sampling design)", y = "") +
  theme(axis.text.y = element_markdown()) +
  custom_theme

# Remove plot margins inside each subplot
p_top_clean <- p_top + theme( plot.margin = unit(c(0, 0, 0, 0), "pt"), panel.spacing = unit(0, "pt"))
p_mid_clean <- p_mid + theme(plot.margin = unit(c(0, 0, 0, 0), "pt"), panel.spacing = unit(0, "pt"))
p_bot_clean <- p_bot + theme(plot.margin = unit(c(0, 0, 0, 0), "pt"), panel.spacing = unit(0, "pt"))

# Combine with no spacing between panels
combined_plot <- p_top_clean / p_mid_clean / p_bot_clean +
  plot_layout(heights = c(0.5, 0.2, 0.7)) +
  plot_annotation(
    theme = theme(
      plot.margin = unit(c(10, 10, 10, 10), "pt")  # Outer margin only
    )
  )

combined_plot

ggsave("plots/waatern_CAP_results.png", plot = combined_plot, width = 8, height = 5, dpi = 300)



## WAATU ----------------------------------------------------------------------

study_site <- "waatu" # waatu or waatern

spp_sd <- read.csv(paste0("data/PRIMER data/", study_site, "/Pearson correlation of CAP (sd) with community (CAP only).csv")) %>% glimpse()
spp_status <- read.csv(paste0("data/PRIMER data/", study_site, "/Pearson correlation of CAP (status) with community (CAP only).csv")) %>% glimpse()

CAP_sd <- read.csv(paste0("data/PRIMER data/", study_site, "/CAP_by_sd_from_Data1.csv")) %>% glimpse()
CAP_status <- read.csv(paste0("data/PRIMER data/", study_site, "/CAP_by_status_from_Data3.csv")) %>% glimpse()


spp_clean <- spp_sd %>%
  filter(!is.na(CAP1)) %>%
  arrange(desc(abs(CAP1))) %>%
  slice_head(n = 15)
spp_clean$taxa <- sub("^[^.]+\\.([^.]+)\\.(.+)$", "\\1 \\2", spp_clean$taxa)
spp_clean$taxa <- paste0("<i>", spp_clean$taxa, "</i>")


## obtain species present in >10% of samples
dat <- readRDS("data/tidy/C01_2024_waatu_all_tidy_maxn.rds") %>% ungroup()
n_samples <- length(unique(dat$opcode))
ubi <- dat %>%
  filter(maxn > 0) %>%
  dplyr::select(opcode, fullspp) %>%
  distinct() %>%  # Unique species per sample
  count(fullspp, name = "n_present") %>%
  mutate(prop = n_present / n_samples) %>%
  filter(prop >= 0.10) %>%
  arrange(desc(prop))
ubi$spp <- gsub(" ", ".", ubi$fullspp) # format to fit axes dataset below

# merge the two datasets by 'opcode'
CAP_merged <- CAP_sd %>%
  dplyr::select(opcode, CAP1_sd = CAP1, sd) %>%
  left_join(CAP_status %>% dplyr::select(opcode, status, CAP1_status = CAP1), by = "opcode") %>% 
  mutate(sd = as.factor(sd)) %>% 
  glimpse()
CAP_merged$sd <- ifelse(CAP_sd$sd == "preferential", 'targeted', CAP_sd$sd)
CAP_merged$status <- ifelse(CAP_sd$status == "No-take", 'no-take', "fished")


spp_axes <- left_join(spp_sd %>% rename(CAP1_sd = CAP1),
                      spp_status %>% rename(CAP1_status = CAP1),
                      by = "taxa"
) %>% drop_na(CAP1_sd, CAP1_status) %>%
  mutate(
    CAP1_sd_scaled = scale(CAP1_sd)[, 1],
    CAP1_status_scaled = scale(CAP1_status)[, 1]
  ) %>%
  mutate(
    magnitude = sqrt(CAP1_sd_scaled^2 + CAP1_status_scaled^2)
  ) %>%
  slice_max(order_by = magnitude, n = 15)
spp_axes <- spp_axes[spp_axes$taxa %in% ubi$spp,] # select ubiquitous spp
spp_axes$taxa <- sub("^[^.]+\\.([^.]+)\\.(.+)$", "\\1 \\2", spp_axes$taxa) # format

spp_axes <- spp_axes %>%
  mutate(
    genus = sub("^(\\S+)\\s.*", "\\1", taxa),
    species = sub("^\\S+\\s+(.*)", "\\1", taxa),
    taxa_label = paste0("<i>", genus, "</i><br><i>", species, "</i>")
  )


arrow_scale = 0.1

waatu_CAP <- ggplot(CAP_merged, aes(x = CAP1_sd, y = CAP1_status, color = sd, shape = status)) +
  geom_point(size = 3) +
  geom_segment(
    data = spp_axes,
    aes(
      x = 0, y = 0,
      xend = CAP1_sd_scaled * arrow_scale,
      yend = CAP1_status_scaled * arrow_scale
    ),
    arrow = arrow(length = unit(0.2, "cm")), color = colour_palette[3], inherit.aes = FALSE
  ) +
  geom_richtext(
    data = spp_axes,
    aes(
      x = CAP1_sd_scaled * arrow_scale,
      y = CAP1_status_scaled * arrow_scale,
      label = taxa_label
    ),
    fill = NA, label.color = NA,  # remove box and outline
    size = 3, vjust = 1, hjust = 0.75, inherit.aes = FALSE
  ) +
  labs(
    x = "CAP axis (sampling design)",
    y = "CAP axis (status)",
    shape = "",
    title = paste0("CAP Scores with Species Correlation Vectors - ", study_site),
    col = ""
  ) +
  scale_colour_manual(
    values = c(
      "spatially balanced" = colour_palette[1],
      "targeted" = colour_palette[4]
    )
  ) +
  custom_theme +
  theme(
    plot.title = element_blank(),
    legend.position = "bottom"
  )
waatu_CAP
ggsave("plots/waatu_CAP_results.png", plot = waatu_CAP, width = 7, height = 5.5, dpi = 300)


### END ###