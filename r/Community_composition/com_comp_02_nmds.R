# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (B. Compare communities detected)
# Data:    2024 BRUV data (MEGLAB and DBCA)
# Task:    Make NMDS
# Author:  Lise Fournier-Carnoy
# Date:    December 2024

# -----------------------------------------------------------------------------

# -- The aim here is to investigate the effect of Sampling Design on community
# -- composition. For this, we will visualise the difference with NMDS, and ana-
# -- lyse it using manyglm.

# -----------------------------------------------------------------------------


# Clear memory
rm(list=ls())
set.seed(12345)

## Load libraries -------------------------------------------------------------

library(vegan) # for community composition
library(ggplot2) # for plots
library(tidyverse) # for data manipulation
library(mvabund) # for community analysis
library(sf) # for manipulating shapefiles

colour_palette <- readLines("chapter_colours.txt"); colour_palette <- eval(parse(text = colour_palette))
source("custom_theme.R")


## Load data and transform into OTU format ----------------------------------

# Files used in this script
file_analysis_box <- "QGIS layers/polys/com_comp_analysis_limits.shp"
file_all_maxn     <- "data/tidy/2024_geographe_all_tidy_maxn.csv"
file_habitat      <- "data/fullspp_habitat_table_ref.csv"

# Box around East Geographe, we're only using this area to compare community.
crs <- st_crs(4326)
box <- st_read(file_analysis_box) %>%
  st_make_valid() %>%
  st_transform(crs)
plot(box)


# MaxN data

maxn <- read.csv(file_all_maxn) %>%
  pivot_wider(id_cols = c(opcode, latitude, longitude, depth, status, sd),
              names_from = fullspp, values_from = MaxN, values_fill = 0) %>%
  filter(!is.na(longitude) & !is.na(latitude)) %>%
  glimpse()
maxn <- st_as_sf(maxn, coords = c("longitude", "latitude"), crs = 4326); plot(maxn)

# Select only points within community composition box
com_comp <- st_intersection(maxn, box); plot(com_comp)
saveRDS(com_comp, 'data/rmd/community_composition_analysis_maxn.rds') # for rmarkdown

# Species-habitat lookup table
hab <- read.csv(file_habitat) %>% 
  glimpse()

## PCoA -----------------------------------------------------------------------

# Separate environmental and abundance data
community_data <- as.data.frame(com_comp) %>%
  ungroup() %>%
  dplyr::select(!1:4 & !id & !geometry) %>% # Select only species columns
  mutate(across(everything(), ~replace_na(., 0))) %>%  # Replace NAs with zeroes
  glimpse()

environmental_data <- as.data.frame(com_comp) %>%
  ungroup() %>%
  dplyr::select(sd,
                #depth,
                status) %>%
  mutate(sd = as.factor(sd),
         status = as.factor(status)) %>%
  glimpse()


#### View difference between sampling designs via PCoA ------------------------

# Load necessary packages
library(vegan)
library(ggplot2)
library(ggrepel)

# Perform Bray-Curtis distance
dist_matrix <- vegdist(community_data, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)  # k = 2 for two principal axes
pcoa_scores <- as.data.frame(pcoa_result$points); head(pcoa_scores)

pcoa_scores$sd <- com_comp$sd; pcoa_scores$status <- com_comp$status; pcoa_scores$opcode <- com_comp$opcode

head(pcoa_scores)

# Remove species with zero variance
species_data <- community_data[, apply(community_data, 2, sd) > 0]

# Check alignment of samples between species data and PCoA scores
rownames(species_data) <- rownames(pcoa_scores)  # Ensure rows match

# Calculate correlation between species data and PCoA axes
correlations <- cor(species_data, pcoa_scores[1:2])

# Create a data frame for species vectors
species_vectors <- data.frame(
  species = colnames(species_data),
  cor_PC1 = correlations[, 1],
  cor_PC2 = correlations[, 2]
)

# Set a threshold for correlation to filter major species (e.g., |correlation| > 0.5)
threshold <- 0.35
major_species <- species_vectors[apply(abs(species_vectors[, c("cor_PC1", "cor_PC2")]), 1, max) > threshold, ]
major_species$species <- ifelse(grepl("^(spp|sp1|sp|\\.spp)$", sub("^[^\\.]+\\.", "", major_species$species)),
                                major_species$species,
                                sub("^[^\\.]+\\.", "", major_species$species))
major_species$species <- sub("\\.", "\n", major_species$species); major_species$species <- sub("\\.", "", major_species$species)

head(major_species)

library(ggrepel)

saveRDS(pcoa_scores, file = "data/rmd/pcoa_scores.rds") # for rmarkdown

ggplot(pcoa_scores, aes(x = V1, y = V2)) +
  geom_point(aes( colour = sd, shape = status), size = 3) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_colour_manual(values = colour_palette[c(4, 6)]) +
  theme_minimal() +
  # geom_segment(data = major_species, aes(x = 0, y = 0, xend = cor_PC1, yend = cor_PC2),
  #              arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "red") +  # Add vectors
  # geom_text_repel(data = major_species, aes(x = cor_PC1, y = cor_PC2, label = paste("italic('", species, "')", sep = "")),
  #                 size = 3, box.padding = 1, point.padding = 0,
  #                 max.overlaps = 20, segment.color = "red",
  #                 vjust = -3, hjust = 1.5, parse = TRUE) +  # Italicize species names
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.box.spacing = unit(0.5, "cm"),
    legend.key.height = unit(1, "lines"),
    legend.key.width = unit(1, "lines"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.spacing.x = unit(0.5, "cm"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16)
  ) +
  ggtitle("PCoA with Major Species Vectors") +
  xlab("PC1") +
  ylab("PC2")



## manyglm() ------------------------------------------------------------------

# method from: Wang, Y., Naumann, U., Wright, S. T., & Warton, D. I. (2012). Mvabund-an R package for model-based analysis of multivariate abundance data. Methods Ecol Evol 3: 471b


# Create a list of response and predictor variables
data_list <- list()
data_list$community <- as.matrix(community_data)
data_list$status <- as.factor(com_comp$status)
data_list$sd <- as.factor(com_comp$sd)
#data_list$depth <- as.numeric(com_comp$depth)
summary(data_list)
attach(data_list)

# Convert community data to a mvabund object
data_abund <- mvabund(data_list$community)
plot(data_abund ~ sd, col = as.numeric(status))

# Test predictors to the community data
geo <- manyglm(community ~ sd * status,
               data = data_list, family = "negative.binomial"
               )
?manyglm
# Check assumptions

# -- Mean variance (random cloud of points around zero, no pattern) and log lin-
# -- -earity (residuals not making a u-shape or other irregular shape)
plot.manyglm(geo)


# -- Normality (no big deviation from the qqline)
par(mfrow = c(2,2))
qqnorm(residuals(geo)[,1], main ='sampling design'); qqline(residuals(geo)[,1], col = "red")
qqnorm(residuals(geo)[,2], main ='status'); qqline(residuals(geo)[,2], col = "red")
qqnorm(residuals(geo)[,3], main = 'design:status'); qqline(residuals(geo)[,3], col = "red")


# Check outputs
result <- anova.manyglm(geo, p.uni = 'adjusted') # takes like 20 min
capture.output(result, file = "outputs/Community_comparison/community_composition_analysis_manyglm_result.txt")

### TRY####

species_hab <- hab %>%
  dplyr::select(fullspp, reef, seagrass, sand, pelagic)

# Step 2: Align species names in community matrix to those in habitat lookup
colnames(community_data) <- gsub("\n", ".", colnames(community_data))  # if needed

# Filter to shared species
shared_species <- intersect(colnames(community_data), species_hab$fullspp)

community_filtered <- community_data[, shared_species]
habitat_filtered <- species_hab %>% filter(fullspp %in% shared_species)

# Step 3: Convert abundance to presence/absence
community_pa <- (community_filtered > 0) * 1  # presence = 1, absence = 0

# Convert to matrix for multiplication
comm_mat <- as.matrix(community_pa)

# Build habitat matrix (species × habitat), aligned
hab_mat <- as.matrix(habitat_filtered[match(colnames(comm_mat), habitat_filtered$fullspp),
                                      c("reef", "seagrass", "sand", "pelagic")])

# Matrix multiplication → site × habitat matrix (species richness per habitat)
site_habitat <- as.data.frame(comm_mat %*% hab_mat)

# PCoA (if not already run)
dist_matrix <- vegdist(community_data, method = "bray")
pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
pcoa_scores <- as.data.frame(pcoa_result$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")

# Correlate habitat richness with PCoA axes
habitat_pcoa_correlation <- cor(site_habitat, pcoa_scores, method = "pearson")

# Prepare for plotting
habitat_vectors <- as.data.frame(habitat_pcoa_correlation)
habitat_vectors$habitat <- rownames(habitat_vectors)

# Plot with presence-based habitat vectors
library(ggplot2)
library(ggrepel)

ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(colour = com_comp$sd, shape = com_comp$status), size = 3) +
  geom_segment(data = habitat_vectors,
               aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2, color = habitat),
               arrow = arrow(length = unit(0.3, "cm")), size = 1.2) +
  geom_text_repel(data = habitat_vectors,
                  aes(x = PCoA1, y = PCoA2, label = habitat),
                  size = 4, fontface = "bold") +
  theme_minimal() +
  ggtitle("PCoA with Habitat Association Vectors (Presence-Based)")

