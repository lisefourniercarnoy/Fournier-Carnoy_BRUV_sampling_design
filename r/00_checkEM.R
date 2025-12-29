# -----------------------------------------------------------------------------

# Project: BRUV sampling design comparison
# Data:    BRUV habitat & abundance data.
# Task:    Check the raw EventMeasure files
# Author:  Lise Fournier-Carnoy / from https://globalarchivemanual.github.io/CheckEM/articles/r-workflows/check-fish-point.html
# Date:    June 2025

# -----------------------------------------------------------------------------

# Status: 

# -----------------------------------------------------------------------------

library('remotes')
options(timeout = 9999999)
# remotes::install_github("GlobalArchiveManual/CheckEM")
library(CheckEM)
library(tidyverse)
library(googlesheets4)
library(sf)
library(terra)
library(here)

# Clear memory
rm(list=ls())

name <- "2024_waatern_BRUVs"
tidy_filename <- "2024_waatern_commonwealth/" # where to save the tidy data THIS SHOULD BE THE SAME IN raw/ AND tidy/

## Format metadata ------------------------------------------------------------

metadata <- read_metadata(here::here(paste0("data/raw/", tidy_filename)), method = "BRUVs") %>% # Change here to "DOVs"
  dplyr::select(campaignid, sample, status, longitude_dd, latitude_dd, date_time, location, site, depth_m, successful_count, successful_length, successful_habitat_forwards, successful_habitat_backwards) %>%
  glimpse()
saveRDS(metadata, file = here::here(paste0("data/tidy/", tidy_filename, name, "_metadata.rds")))


## Format points --------------------------------------------------------------

points <- read_points(here::here(paste0("data/raw/", tidy_filename))) %>%
  glimpse()

if(nrow(points) > 1){
  maxn_points <- points %>% 
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
#maxn <- bind_rows(get0("maxn_points"), get0("maxn_counts")) # this works even if you only have one type of data

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
summary(count$`Pempherididae Parapriacanthus elongatus`)
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


## Format lengths -------------------------------------------------------------


# currently waatu lengths not done, run when they're checked
em_length3dpoints <- read_em_length(here::here(paste0("data/raw/", tidy_filename))) %>%                   
  dplyr::select(-c(comment)) %>% # there is a comment column in metadata, so you will need to remove this column from EM data
  dplyr::inner_join(metadata, by = join_by(sample, campaignid)) %>%
  dplyr::filter(successful_length %in% "Yes") %>%
  # dplyr::rename(length_mm = length) %>%
  glimpse() 

length <- bind_rows(get0("em_length3dpoints"), get0("gen_length")) # this works even if you only have one type of data

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


## Quality control checks -----------------------------------------------------

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
missing_count <- anti_join(metadata_samples, samples, by = join_by(campaignid, sample))
message(paste(nrow(missing_count), "samples in the metadata missing count data"))
glimpse(missing_count)


metadata_samples <- metadata %>%
  dplyr::select(campaignid, sample, dplyr::any_of(c("opcode", "period")), 
                successful_count, successful_length) %>%
  distinct()
samples <- length %>%
  distinct(campaignid, sample)
missing_length <- anti_join(metadata_samples, samples, by = join_by(campaignid, sample))
message(paste(nrow(missing_length), "samples in metadata missing length data"))


missing_metadata <- anti_join(samples, metadata_samples, by = join_by(campaignid, sample))
message(paste(nrow(missing_metadata), "samples in length data missing metadata"))

periods <- read_periods(here::here("r-workflows/data/raw/")) %>%
  glimpse()
periods_without_end <- periods %>%
  dplyr::filter(has_end == 0)
message(paste(nrow(periods_without_end), "periods without an end"))

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


points_outside_periods <- points %>%
  dplyr::filter(period %in% c("NA", NA, NULL, "")) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, number, frame)
message(paste(nrow(points_outside_periods), "points outside a period"))
glimpse(points_outside_periods)


lengths_outside_periods <- em_length3dpoints %>%
  dplyr::filter(period %in% c("NA", NA, NULL, "")) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, number)
message(paste(nrow(lengths_outside_periods), "lengths/3D points outside period"))
glimpse(lengths_outside_periods)


period_length <- 60 # in minutes
periods_wrong <- periods %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), time_start, time_end, has_end) %>%
  dplyr::distinct() %>%
  dplyr::mutate(period_time = round(time_end - time_start)) %>%
  dplyr::filter(!period_time %in% period_length)
message(paste(nrow(periods_wrong), "periods not", period_length, "minutes long"))


total_count <- sum(complete_count$count)
message(paste(total_count, "fish counted in the count data"))

total_length <- sum(complete_length$number)
message(paste(total_length, "fish counted in the length data"))


points_without_number <- points %>%
  filter(number %in% c("NA", NA, 0, NULL, "", " "))
message(paste(nrow(points_without_number), "points in the _Points.txt file that do not have a number"))
glimpse(points_without_number)

lengths_without_number <- em_length3dpoints %>%
  filter(number %in% c("NA", NA, 0, NULL, "", " "))
message(paste(nrow(lengths_without_number), "lengths or 3D points in the EMObs that do not have a number"))

synonyms_in_count <- dplyr::left_join(complete_count, CheckEM::aus_synonyms) %>%
  dplyr::filter(!is.na(genus_correct)) %>%
  dplyr::mutate('old name' = paste(family, genus, species, sep = " ")) %>%
  dplyr::mutate('new name' = paste(family_correct, genus_correct, species_correct, sep = " ")) %>%
  dplyr::select('old name', 'new name') %>%
  dplyr::distinct()
message(paste(nrow(synonyms_in_count), "synonyms used in the count data"))


synonyms_in_length <- dplyr::left_join(complete_length, CheckEM::aus_synonyms) %>%
  dplyr::filter(!is.na(genus_correct)) %>%
  dplyr::mutate('old name' = paste(family, genus, species, sep = " ")) %>%
  dplyr::mutate('new name' = paste(family_correct, genus_correct, species_correct, sep = " ")) %>%
  dplyr::select('old name', 'new name') %>%
  dplyr::distinct()

incorrect_lengths <- left_join(complete_length, create_min_max(CheckEM::australia_life_history, minimum = 0.15, maximum = 0.85)) %>%
  dplyr::filter(length_mm < min_length_mm | length_mm > max_length_mm) %>%
  mutate(reason = ifelse(length_mm < min_length_mm, "too small", "too big")) %>%
  dplyr::select(campaignid, sample, family, genus, species, length_mm, min_length_mm, max_length_mm, length_max_mm, reason, any_of(c("em_comment", "frame_left")), length_max_mm) %>%
  mutate(difference = ifelse(reason %in% c("too small"), (min_length_mm - length_mm), (length_mm - max_length_mm))) %>%
  dplyr::mutate(percent_of_fb_max = (length_mm/length_max_mm)*100) %>%
  dplyr::left_join(metadata) %>%
  dplyr::select(campaignid, dplyr::any_of(c("opcode", "period")), family, genus, species, length_mm, min_length_mm, max_length_mm, length_max_mm, reason, any_of(c("em_comment", "frame_left")), difference, percent_of_fb_max)

## Save the checked, clean data -----------------------------------------------
saveRDS(complete_count,
        file = here::here(paste0("data/tidy/00_",
                                 name, "_tidy_counts.rds")))
saveRDS(complete_length,
        file = here::here(paste0("data/tidy/00_",
                                 name, "_tidy_lengths.rds")))

