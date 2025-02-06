# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Compare Lengths detected)
# Data:    2014 MEGlab BRUV data. (TO BE CHANGED TO 2024 DATA)
# Task:    Extract the mean and SE of 1000 sampling designs
# Author:  Lise Fournier-Carnoy
# Date:    September 2024

# -----------------------------------------------------------------------------

# Status:  Starting...

# -----------------------------------------------------------------------------

# Load libraries
library(raster)         # Dealing with rasters
library(sp)             # Dealing with shapefiles
library(sf)             # Dealing with shapefiles
library(terra)          # Dealing with shapefiles
library(tidyverse)      # Data manipulation
library(ggstatsplot)    # For plotting with stats integrated
library(FSSgam)         # For GAMs?
library(mgcv)           # For GAMs

# Clear memory
rm(list=ls())
set.seed(12345)

## Files used in this script --------------------------------------------------

file_sim_SZ       <- "data/spatial/shapefiles/simulated_SZ.shp"
file_samp_area    <- "data/spatial/shapefiles/sampling_area_deeper_than_7m.rds"

file_mature_pred  <- "outputs/Length_comparison/predicted_abundance_rasters/mature_predicted_abundance_SZ_increase.rds"

file_spabal_val   <- "outputs/Length_comparison/sampling_designs/spatially_balanced_extracted_values.rds"
file_pref_val     <- "outputs/Length_comparison/sampling_designs/preferential_extracted_values.rds"


## Load files -----------------------------------------------------------------

sim_sz <- st_read(file_sim_SZ) %>% # Simulated SZ
  st_make_valid() %>%
  st_transform(4326); plot(sim_sz)

samp_area <- readRDS(file_samp_area) %>% # Sampling area to crop true abundance raster
  vect(); plot(samp_area)

spabal <- readRDS(file_spabal_val) %>%
  glimpse()
pref <- readRDS(file_pref_val) %>%
  glimpse()

dat <- rbind(spabal, pref) %>% # Full dataframe, with spatially balanced and preferential drops
  na.omit() %>%
  glimpse()


## Normal distribution point --------------------------------------------------

# For each point sampled from the fake abundance map, I want a random abundance within the normal distribution of moean +/- SE

# For this I need a loop that creates the normal distribution using mean and SE
# Then randomly selects an abundance value within it.
# This is done for each point in all sampling designs.

# Create a new column to store the randomly selected values of abundance
dat$random_abundance_before <- NA
dat$random_abundance_after <- NA

for(i in 1:nrow(dat)) {

  ### Before SZ ---------------------------------------------------------------

  # Get the mean and standard deviation for the current row
  mean_value_before <- dat$fit_value_before[i]
  sd_value_before <- dat$se_value_before[i]

  # Generate a random value from the normal distribution
  random_abundance_before <- rnorm(1, mean = mean_value_before, sd = sd_value_before)

  # Replace negative values with 0
  random_abundance_before <- pmax(random_abundance_before, 0)

  # Assign the random value to the new column
  dat$random_abundance_before[i] <- random_abundance_before

  ### Repeat for after SZ -----------------------------------------------------

  # Get the mean and standard deviation for the current row
  mean_value_after <- dat$fit_value_after[i]
  sd_value_after <- dat$se_value_after[i]

  # Generate a random value from the normal distribution
  random_abundance_after <- rnorm(1, mean = mean_value_after, sd = sd_value_after)

  # Replace negative values with 0
  random_abundance_after <- pmax(random_abundance_after, 0)

  # Assign the random value to the new column
  dat$random_abundance_after[i] <- random_abundance_after
}

# View the updated dataframe with the random points
glimpse(dat)
saveRDS(dat, 'data/rmd/SD_points.rds') # for Rmarkdown

## Analyse SZ detection ability -----------------------------------------------

# From there I want to check whether each sampling design is able to detect the
# difference between inside and outside the SZ (which we simulated earlier)

# I will make a ratio of inside/outside abundances for each sampling designs
# (if I have 100 simulations of preferential and 100 simulations of spatially
# balanced, I should get 200 ratios)

# In 05_simulate_SZ, I multiplied the modelled abundance in the SZ by 1.8 (80%
# increase), so I will check how many of the simulations fall within 10% of that
# actual inside/outside ratio.

ratio_after <- dat %>%
  group_by(design_id, SD, in_SZ) %>%
  summarize(mean_random_abundance_after = mean(random_abundance_after, na.rm = TRUE)) %>%
  spread(key = in_SZ, value = mean_random_abundance_after) %>%
  mutate(abundance_ratio = `TRUE` / `FALSE`) %>%
  dplyr::select(SD, design_id, abundance_ratio) %>%
  glimpse()

# How many simulations give you the 1.8x increase (within 10% of 1.8x)
abundance_increase_in_SZ <- 1.8
upper_bound <- abundance_increase_in_SZ+(abundance_increase_in_SZ*0.1) # upper bound of acceptable ratio
lower_bound <- abundance_increase_in_SZ-(abundance_increase_in_SZ*0.1) # lower bound of acceptable ratio

percentage_after <- ratio_after %>%
  filter(abundance_ratio >= lower_bound & abundance_ratio <= upper_bound) %>%  # Filter for the desired range
  group_by(SD) %>%
  summarise(percentage_within_bounds = n() / nrow(ratio_after %>% filter(SD == first(SD))) * 100)  # Calculate the percentage for each SD



### BACI comparison -----------------------------------------------------------

# I can't really know whether the sampling design detects the SZ effect or just
# a difference in habitat that happens to be inside a SZ. For example, my fake-
# SZ has a larger reef area (complex strata) than the non-SZ bits outside it,
# and there may be more fish there naturally.
# If we test the ability of sampling designs to detect difference in abundance
# before and after the SZ (BACI design) we can actually say whether the SZ
# effect is detected correctly.

# So. I need to compare the ratio of ~before~ and the ratio of ~after~. There
# should definitely be a difference.


ratio_before <- dat %>%
  group_by(design_id, SD, in_SZ) %>%
  summarize(mean_random_abundance_before = mean(random_abundance_before, na.rm = TRUE)) %>%
  spread(key = in_SZ, value = mean_random_abundance_before) %>%
  mutate(abundance_ratio = `TRUE` / `FALSE`) %>%
  dplyr::select(SD, design_id, abundance_ratio) %>%
  glimpse()

percentage_before <- ratio_before %>%
  filter(abundance_ratio >= lower_bound & abundance_ratio <= upper_bound) %>%
  group_by(SD) %>%
  summarise(percentage_within_bounds = n() / nrow(ratio_before %>% filter(SD == first(SD))) * 100) %>%  # Calculate the percentage for each SD
  glimpse() # like 4% of the simulations have a 1.8x ratio, which is normal for the ~before~

# View the difference
percentage_before # ~4% of simulations detect a 1.8x difference inside v. outside SZ before the SZ increase (nice and low)
percentage_after # ~12% of simulations detect a 1.8x difference inside v. outside SZ after the SZ increase(should not be that low)

# This result seems to suggest that the detection of SZ effects occurs in only
# 12% of sampling, regardless of the design.

# We've looked at the whole area, inside and outside SZ, but is there a
# difference of 1.8x before and after, only in SZ?

dat_in_SZ <- dat %>% # Re-organising the dataframe so it's easier to work with.
  filter(in_SZ = TRUE) %>%
  pivot_longer(
    cols = c(random_abundance_before, random_abundance_after),
    names_to = "BA",
    values_to = "random_abundance"
  ) %>%
  mutate(BA = recode(BA, "random_abundance_before" = "before", "random_abundance_after" = "after")) %>%
  dplyr::select(c(ID_1, design_id, SD, BA, random_abundance)) %>%
  glimpse()

ratio_SZ <- dat_in_SZ %>%
  group_by(design_id, SD, BA) %>%
  summarize(mean_random_abundance = mean(random_abundance, na.rm = TRUE)) %>%
  spread(key = BA, value = mean_random_abundance) %>%
  mutate(abundance_ratio = `before` / `after`) %>%
  dplyr::select(SD, design_id, abundance_ratio) %>%
  glimpse()

percentage <- ratio_SZ %>%
  filter(abundance_ratio >= 0.7 & abundance_ratio <= 0.9) %>% # Filter for the desired range
  group_by(SD) %>%
  summarise(percentage_within_bounds = n() / nrow(ratio_SZ %>% filter(SD == first(SD))) * 100) # Calculate the percentage for each SD

percentage
saveRDS(percentage, 'data/rmd/SD_performance.rds')

# Okay so. Comparing each sampling design's ability to detect SZ effect,
# preferential design detects 0.7-0.9 difference 30% of the time, and spatially balanced 22.9% of the time.



### END ###
