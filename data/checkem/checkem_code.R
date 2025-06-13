# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. length model)
# Data:    2024 MEGlab Geographe Bay BRUV habitat data.
# Task:    Use CheckEM to check the exported EventMeasure length and MaxN data.
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    February 2025

# -----------------------------------------------------------------------------

# -- The following code is copy-pasted, then commented from https://globalarchivemanual.github.io/CheckEM/articles/r-workflows/check-fish-point.html#save-the-checked-data


# Libraries ----

# install.packages('remotes')
library('remotes')
options(timeout=9999999)
# remotes::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)
library(tidyverse)
library(googlesheets4)
library(sf)
library(terra)
library(here)

name <- "2024-04_Geographe_stereo-BRUVs"
file_loc <- "data/raw/2024_commonwealth/20250225_1230/"



# Metadata ----

metadata <- read_metadata(here::here(file_loc), method = "BRUVs") %>% # Change here to "DOVs"
  dplyr::select(campaignid, sample, status, longitude_dd, latitude_dd, date_time, location, site, depth_m, successful_count, successful_length, successful_habitat_forwards, successful_habitat_backwards) %>%
  glimpse()
saveRDS(metadata, file = here::here(paste0(file_loc,
                                           name, "_Metadata.rds")))

metadata <- readRDS(here::here(paste0(file_loc,
                                      name, "_Metadata.rds")))
metadata_sf <- st_as_sf(metadata, coords = c("longitude_dd", "latitude_dd"), crs = 4326)
metadata_sf <- st_as_sf(metadata, coords = c("longitude_dd", "latitude_dd"), crs = 4326)
regions <- st_as_sf(CheckEM::aus_regions, crs = st_crs(4326))


regions <- st_transform(regions, 4326) %>%
  dplyr::select(REGION)

metadata <- st_join(metadata_sf, regions, join = st_nearest_feature) %>%
  dplyr::rename(marine_region = REGION) %>%
  dplyr::mutate(sample = as.character(sample)) %>%
  as.data.frame() %>%
  dplyr::select(-c(geometry)) %>%
  glimpse()


# Load MaxN points & lengths ----

points <- read_points(here::here(file_loc)) %>%
  glimpse()
# Only run this if there is data in the counts data frame
if(nrow(points) > 1){
  maxn <- points %>%
    dplyr::group_by(campaignid, sample, filename, periodtime, frame, family, genus, species) %>% # If you have MaxN'd by stage (e.g. Adult, Juvenile) add stage here
    dplyr::mutate(number = as.numeric(number)) %>%
    dplyr::summarise(maxn = sum(number)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(campaignid, sample, family, genus, species) %>%
    dplyr::slice(which.max(maxn)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(maxn)) %>%
    dplyr::select(-frame) %>%
    tidyr::replace_na(list(maxn = 0)) %>%
    dplyr::mutate(maxn = as.numeric(maxn)) %>%
    dplyr::filter(maxn > 0) %>%
    dplyr::inner_join(metadata, by = join_by(campaignid, sample)) %>%
    dplyr::filter(successful_count %in% c("Yes")) %>%
    dplyr::filter(maxn > 0) %>%
    dplyr::select(campaignid, sample, family, genus, species, maxn) %>%
    dplyr::glimpse()
}

em_length3dpoints <- read_em_length(here::here(file_loc)) %>%
  dplyr::select(-c(comment))%>% # there is a comment column in metadata, so you will need to remove this column from EM data
  dplyr::inner_join(metadata, by = join_by(sample, campaignid)) %>%
  dplyr::filter(successful_length %in% "Yes") %>%
  # dplyr::rename(length_mm = length) %>%
  glimpse()
gen_length <- read_gen_length(here::here(file_loc)) %>%
  dplyr::full_join(metadata, by = join_by(campaignid, sample)) %>%
  dplyr::filter(successful_length %in% "Yes") %>%
  glimpse()
# If only EventMeasure data then length only includes Length and 3D points data
# If only Generic data then length only includes generic length data
# If both exist, then length includes both Length and 3D points and generic length data
length <- bind_rows(get0("em_length3dpoints"), get0("gen_length")) # this works even if you only have one type of data
count <- maxn %>%
  dplyr::mutate(family = ifelse(family %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(family))) %>%
  dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
  dplyr::mutate(species = ifelse(species %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "spp", as.character(species))) %>%
  dplyr::select(campaignid, sample, family, genus, species, maxn) %>%
  tidyr::complete(nesting(campaignid, sample), nesting(family, genus, species)) %>%
  tidyr::replace_na(list(maxn = 0)) %>%
  group_by(campaignid, sample, family, genus, species) %>%
  dplyr::summarise(count = sum(maxn)) %>%
  ungroup() %>%
  mutate(scientific = paste(family, genus, species, sep = " "))%>%
  dplyr::select(campaignid, sample, scientific, count)%>%
  spread(scientific, count, fill = 0)
count_families <- maxn %>%
  dplyr::mutate(scientific = paste(family, genus, species, sep = " ")) %>%
  filter(!(family %in% "Unknown")) %>%
  dplyr::select(c(family, genus, species, scientific)) %>%
  distinct()

complete_count <- count %>%
  pivot_longer(names_to = "scientific", values_to = "count",
               cols = 3:ncol(.)) %>%
  inner_join(count_families, by = c("scientific")) %>%
  full_join(metadata)%>%
  glimpse()
complete_length <- length %>%
  dplyr::mutate(family = ifelse(family %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(family))) %>%
  dplyr::mutate(genus = ifelse(genus %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "Unknown", as.character(genus))) %>%
  dplyr::mutate(species = ifelse(species %in% c("NA", "NANA", NA, "unknown", "", NULL, " ", NA_character_), "spp", as.character(species))) %>%
  dplyr::filter(!family %in% "Unknown")%>%
  # First make one row for every length measurement
  dplyr::mutate(number = as.numeric(number)) %>%
  tidyr::uncount(number) %>%
  dplyr::mutate(number = 1) %>%
  # Add in missing samples
  dplyr::right_join(metadata) %>%
  dplyr::filter(successful_length %in% "Yes") %>%
  # Complete the data (add in zeros for every species)
  dplyr::select(campaignid, sample, family, genus, species, length_mm, number, any_of(c("range", "rms", "precision"))) %>% # this will keep EM only columns
  tidyr::complete(nesting(campaignid, sample), nesting(family, genus, species)) %>%
  replace_na(list(number = 0)) %>%
  ungroup() %>%
  dplyr::filter(!is.na(number)) %>%
  dplyr::mutate(length_mm = as.numeric(length_mm)) %>%
  left_join(., metadata) %>%
  glimpse()

number_of_samples <- metadata %>%
  dplyr::distinct(campaignid, sample)

message(paste(nrow(number_of_samples), "unique samples in the metadata"))
duplicate_samples <- metadata %>%
  dplyr::group_by(campaignid, sample) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  dplyr::filter(n > 1)
message(paste(nrow(duplicate_samples), "samples duplicated in the metadata"))
metadata_samples <- metadata %>%
  dplyr::select(campaignid, sample, dplyr::any_of(c("opcode", "period")),
                successful_count, successful_length) %>%
  distinct()

samples <- maxn %>%
  distinct(campaignid, sample)

# Check for missing counts ----

missing_count <- anti_join(metadata_samples, samples, by = join_by(campaignid, sample))

message(paste(nrow(missing_count), "samples in the metadata missing count data"))
# 2 unsuccessful counts in the metadata, so this is fine

missing_metadata <- anti_join(samples, metadata_samples, by = join_by(campaignid, sample))
message(paste(nrow(missing_metadata), "samples in count data missing metadata"))
# no missing metadata in count data, all good


# checking missing lengths ----

metadata_samples <- metadata %>%
  dplyr::select(campaignid, sample, dplyr::any_of(c("opcode", "period")),
                successful_count, successful_length) %>%
  distinct()

samples <- length %>%
  distinct(campaignid, sample)

missing_length <- anti_join(metadata_samples, samples, by = join_by(campaignid, sample))
message(paste(nrow(missing_length), "samples in metadata missing length data"))
# 15 videos not lengthable, so this is fine.

missing_metadata <- anti_join(samples, metadata_samples, by = join_by(campaignid, sample))
message(paste(nrow(missing_metadata), "samples in length data missing metadata"))
# no metadata missing in lengths, all good


# Check for missing periods ----
periods <- read_periods(here::here(file_loc)) %>%
  glimpse()
periods_without_end <- periods %>%
  dplyr::filter(has_end == 0)
message(paste(nrow(periods_without_end), "periods without an end"))
# no missing periods, all good

metadata_samples <- metadata %>%
  dplyr::select(campaignid, sample, dplyr::any_of(c("opcode", "period")), successful_count, successful_length) %>%
  dplyr::distinct() %>%
  dplyr::mutate(sample = as.factor(sample))

periods_samples <- periods %>%
  dplyr::select(campaignid, sample, dplyr::any_of(c("opcode", "period"))) %>%
  distinct()

missing_periods <- anti_join(metadata_samples, periods_samples) %>%
  dplyr::select(!sample)
message(paste(nrow(missing_periods), "samples missing period"))
# 2 missing periods = 2 not maxn-able, all good


# Check for points outside periods ----

points_outside_periods <- points %>%
  dplyr::filter(period %in% c("NA", NA, NULL, "")) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, number, frame)
message(paste(nrow(points_outside_periods), "points outside a period"))
glimpse(points_outside_periods)
# no points outside periods, all good


# Check for lengths outside periods ----

lengths_outside_periods <- em_length3dpoints %>%
  dplyr::filter(period %in% c("NA", NA, NULL, "")) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, number)
message(paste(nrow(lengths_outside_periods), "lengths/3D points outside period"))
glimpse(lengths_outside_periods)
# these are sync points (no genus or families)
glimpse(lengths_outside_periods[!is.na(lengths_outside_periods$species),])
# these are actual species 3d-pointed or lengthed outside of periods


# Check for wrong period duration ----
period_length <- 60 # in minutes
periods_wrong <- periods %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), time_start, time_end, has_end) %>%
  dplyr::distinct() %>%
  dplyr::mutate(period_time = round(time_end - time_start)) %>%
  dplyr::filter(!period_time %in% period_length)
message(paste(nrow(periods_wrong), "periods not", period_length, "minutes long"))
glimpse(periods_wrong)
# if the wrong periods are around 60 minutes it's all good

# Count total fish MaxNed and lengthed ----
total_count <- sum(complete_count$count)
message(paste(total_count, "fish counted in the count data"))
total_length <- sum(complete_length$number)
message(paste(total_length, "fish counted in the length data"))


# Check for points without a number ----
points_without_number <- points %>%
  filter(number %in% c("NA", NA, 0, NULL, "", " "))
message(paste(nrow(points_without_number), "points in the _Points.txt file that do not have a number"))
# no points without number, all good


# Check for points without a number ----
lengths_without_number <- em_length3dpoints %>%
  filter(number %in% c("NA", NA, 0, NULL, "", " "))
message(paste(nrow(lengths_without_number), "lengths or 3D points in the EMObs that do not have a number"))
glimpse(lengths_without_number[!is.na(lengths_without_number$species),])
# these are all sync points (no points here are species)


# check for species synonyms ----
synonyms_in_count <- dplyr::left_join(complete_count, CheckEM::aus_synonyms) %>%
  dplyr::filter(!is.na(genus_correct)) %>%
  dplyr::mutate('old name' = paste(family, genus, species, sep = " ")) %>%
  dplyr::mutate('new name' = paste(family_correct, genus_correct, species_correct, sep = " ")) %>%
  dplyr::select('old name', 'new name') %>%
  dplyr::distinct()
message(paste(nrow(synonyms_in_count), "synonyms used in the count data"))
# no synonyms used, all good

synonyms_in_length <- dplyr::left_join(complete_length, CheckEM::aus_synonyms) %>%
  dplyr::filter(!is.na(genus_correct)) %>%
  dplyr::mutate('old name' = paste(family, genus, species, sep = " ")) %>%
  dplyr::mutate('new name' = paste(family_correct, genus_correct, species_correct, sep = " ")) %>%
  dplyr::select('old name', 'new name') %>%
  dplyr::distinct()
message(paste(nrow(synonyms_in_length), "synonyms used in the length data"))
# no synonyms used, all good


# check for region and species ----
count_species_not_observed_region <- complete_count %>%
  dplyr::distinct(campaignid, sample, family, genus, species, marine_region, count) %>%
  dplyr::anti_join(., expand_life_history(CheckEM::australia_life_history), by = c("family", "genus", "species", "marine_region")) %>%
  dplyr::filter(count > 0) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, marine_region) %>%
  dplyr::distinct() %>%
  dplyr::rename('marine region not observed in' = marine_region) %>%
  dplyr::semi_join(., CheckEM::australia_life_history, by = c("family", "genus", "species"))
message(paste(nrow(count_species_not_observed_region), "species not observed in the region before"))
glimpse(count_species_not_observed_region)

length_species_not_observed_region <- complete_length %>%
  dplyr::distinct(campaignid, sample, family, genus, species, marine_region, number) %>%
  dplyr::anti_join(., expand_life_history(CheckEM::australia_life_history), by = c("family", "genus", "species", "marine_region")) %>%
  dplyr::filter(number > 0) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, marine_region) %>%
  dplyr::distinct() %>%
  dplyr::rename('marine region not observed in' = marine_region) %>%
  dplyr::semi_join(., CheckEM::australia_life_history, by = c("family", "genus", "species"))
message(paste(nrow(length_species_not_observed_region), "species not observed in the region before"))
# s. obtusata waiting on confirmation from museum
# R. australiae is fine
# A. aurita is fine
# S. apama is fine
# C. kumu is fine
# S. lewini is fine


# Check for life history list ----

count_species_not_in_list <- complete_count %>%
  dplyr::anti_join(., CheckEM::australia_life_history, by = c("family", "genus", "species")) %>%
  dplyr::filter(count > 0) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species) %>%
  dplyr::distinct()
message(paste(nrow(count_species_not_in_list), "species not in chosen life history list"))
glimpse(count_species_not_in_list)
# sp1s and sp2s, all good

length_species_not_in_list <- complete_length %>%
  dplyr::anti_join(., CheckEM::australia_life_history, by = c("family", "genus", "species")) %>%
  dplyr::filter(number > 0) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species) %>%
  dplyr::distinct()
message(paste(nrow(length_species_not_in_list), "species not in chosen life history list"))
glimpse(length_species_not_in_list)
#sp1s and sp2s all good

# checking lengths ----
incorrect_lengths <- left_join(complete_length, create_min_max(CheckEM::australia_life_history, minimum = 0.15, maximum = 0.85)) %>%
  dplyr::filter(length_mm < min_length_mm | length_mm > max_length_mm) %>%
  mutate(reason = ifelse(length_mm < min_length_mm, "too small", "too big")) %>%
  dplyr::select(campaignid, sample, family, genus, species, length_mm, min_length_mm, max_length_mm, length_max_mm, reason, any_of(c("em_comment", "frame_left")), length_max_mm) %>%
  mutate(difference = ifelse(reason %in% c("too small"), (min_length_mm - length_mm), (length_mm - max_length_mm))) %>%
  dplyr::mutate(percent_of_fb_max = (length_mm/length_max_mm)*100) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, length_mm, min_length_mm, max_length_mm, length_max_mm, reason, any_of(c("em_comment", "frame_left")), difference, percent_of_fb_max)
too_small <- incorrect_lengths %>%
  dplyr::filter(reason %in% "too small")

too_big <- incorrect_lengths %>%
  dplyr::filter(reason %in% "too big")
message(paste(nrow(too_small), "lengths are too small"))
glimpse(too_small) # any that are over 15% need to be checked
message(paste(nrow(too_big), "lengths are too big"))
glimpse(too_big) # check tail end of each spp
# These have all been checked. 

rms_limit <- 20 # in mm
over_rms <- complete_length %>%
  dplyr::filter(as.numeric(rms) > rms_limit)
message(paste(nrow(over_rms), "lengths over RMS limit"))
# all good

precision_limit <- 10 # in %
over_precision <- complete_length %>%
  dplyr::filter(as.numeric(precision) > precision_limit)
message(paste(nrow(over_precision), "lengths over precision limit"))
# idk what that means


# save complete counts & lengths ----

saveRDS(complete_count,
        file = here::here(paste0("data/tidy/",
                                 name, "_complete-count.rds")))
saveRDS(complete_length,
        file = here::here(paste0("data/tidy/",
                                 name, "_complete-length.rds")))
