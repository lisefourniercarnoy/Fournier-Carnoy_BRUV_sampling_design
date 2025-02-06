#####################################################################
# Project: API query for Australian BRUV synthesis
# Data:    2019 and 2024 BRUV Syntheses
# Task:    Use GlobalArchive API to query data and save as an RDS
# Author:  Brooke Gibbons
# Date:    Feb 2024
#####################################################################

# Load libraries needed -----
library(httr)
library(tidyverse)
library(RJSONIO)
library(devtools)
devtools::install_github("GlobalArchiveManual/CheckEM") # If there has been any updates to the package then CheckEM will install
library(CheckEM)

username <- "public"
password <- "sharedaccess"
synthesis_id <- 14 # For Geographe bay fish and habitat synthesis

# First, API call to the GlobalArchive species list to join to the count and length data ----
species_list <- CheckEM::ga_api_species_list(username, password)

# API call for metadata ----
metadata <- CheckEM::ga_api_metadata(username, password, synthesis_id = synthesis_id)

# API call for count data ----
count <- CheckEM::ga_api_count(username, password, synthesis_id = synthesis_id) %>%
  dplyr::select(sample, family, genus, species, count) %>%
  glimpse()

# API call for length data ----
length <- CheckEM::ga_api_length(username, password, synthesis_id = synthesis_id) %>%
  dplyr::select(sample, family, genus, species, length, number) %>%
  glimpse()

# Save data ----
saveRDS(metadata, "data/raw/geographe_metadata.RDS")
saveRDS(count, "data/raw/geographe_count.RDS")
saveRDS(length, "data/raw/geographe_length.RDS")

# Add in zeros where species are not present in the count data ----
# NOTE: this creates a very large file, it also takes quite a while to run

count_metadata <- metadata %>%
  dplyr::filter(successful_count %in% TRUE)

length(unique(count_metadata$sample)) # 297 successful samples
length(unique(count$sample)) # 297 samples in the count, will need to add in the zeros where no fish were observed

count_wide <- count %>%
  dplyr::full_join(count_metadata) %>%
  dplyr::filter(successful_count %in% TRUE) %>% # now has correct number of samples
  dplyr::select(campaign, sample, family, genus, species, count) %>%
  tidyr::complete(nesting(campaign, sample), nesting(family, genus, species)) %>%
  tidyr::replace_na(list(count = 0)) %>%
  dplyr::group_by(campaign, sample, family, genus, species) %>%
  dplyr::summarise(count = sum(count)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(scientific = paste(family, genus, species, sep = " "))%>%
  dplyr::select(campaign, sample, scientific, count)%>%
  spread(scientific, count, fill = 0)

length(unique(count_wide$sample))

count_families <- count %>%
  dplyr::mutate(scientific = paste(family, genus, species, sep = " ")) %>%
  dplyr::ungroup() %>%
  dplyr::select(c(family, genus, species, scientific)) %>%
  dplyr::distinct()

complete_count <- count_wide %>%
  pivot_longer(names_to = "scientific", values_to = "count", cols = 3:ncol(.)) %>%
  dplyr::inner_join(count_families, by = c("scientific")) %>%
  dplyr::full_join(count_metadata) %>%
  dplyr::filter(successful_count %in% TRUE) %>%
  glimpse()

# Save complete count (zeros added where a species is not present)
saveRDS(complete_count, "data/raw/geographe_complete_count.RDS")

# Add in zeros where species are not present in the count data ----
# NOTE: this creates a very large file, it also takes quite a while to run

length_metadata <- metadata %>%
  dplyr::filter(successful_length %in% TRUE)

length(unique(length_metadata$sample)) # 239 successful samples
length(unique(length$sample)) # 154 samples in the length, will need to add in the zeros where no fish were observed

complete_length <- length %>%
  dplyr::mutate(number = as.numeric(number)) %>%
  dplyr::filter(!is.na(number)) %>%
  tidyr::uncount(number) %>%
  dplyr::mutate(number = 1) %>%
  dplyr::full_join(length_metadata) %>%
  dplyr::filter(successful_length %in% TRUE) %>%
  dplyr::select(campaign, sample, family, genus, species, length, number) %>%
  tidyr::complete(nesting(campaign, sample), nesting(family, genus, species)) %>%
  replace_na(list(number = 0)) %>%
  ungroup() %>%
  dplyr::mutate(length_mm = as.numeric(length)) %>%
  full_join(length_metadata) %>%
  dplyr::filter(!is.na(number)) %>%
  dplyr::filter(successful_length %in% TRUE) %>%
  glimpse()

length(unique(complete_length$sample))

# Save complete count (zeros added where a species is not present)
saveRDS(complete_length, "data/raw/geographe_complete_length.RDS")
