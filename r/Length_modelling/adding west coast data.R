#####################################################################
# Project: Pilbara Marine conservation Partnership
# Data:    Campaigns from the PMCP, from the 2024 Aus BRUV synthesis
# Task:    Use GlobalArchive API to query data and filter, add habitat annotations
# Author:  Claude spencer
# Date:    May 2025
#####################################################################

# Clear environment
rm(list = ls())

# Install CheckEM package ----
options(timeout = 9999999) # the package is large, so need to extend the timeout to enable the download.
remotes::install_github("GlobalArchiveManual/CheckEM") # If there has been any updates to the package then CheckEM will install, if not then this line won't do anything

# Load libraries needed -----
library(CheckEM)
library(httr)
library(tidyverse)
library(RJSONIO)
library(devtools)
library(leaflet)
library(arrow)

# Set your API token to access GlobalArchive data shared with you ----
# It is extremely important that you keep your API token out of your scripts, and github repository!
# This function will ask you to put your API token in the console
# It will then create a folder in your project folder called "secrets" and saves your API token to use in the functions later
# The function adds the token into the .gitignore file so it will never be put in version control with Git
CheckEM::ga_api_set_token()

# Load the saved token
token <- readRDS("secrets/api_token.RDS")

# Load the metadata, count and length ----
# This way does not include the zeros where a species isn't present - it returns a much smaller dataframe
CheckEM::ga_api_all_data(synthesis_id = "57",
                         token = token,
                         dir = "data/raw/2024_swc_state/",
                         include_zeros = FALSE)

benthos <- CheckEM::ga_api_habitat(token = token,
                                   synthesis_id = "57")






## checking which campaigns there are -----------------------------------------

metadata <- readRDS("data/raw/2024_swc_state/metadata.rds") %>% glimpse()
benthos <- readRDS("data/raw/2024_swc_state/benthos_summarised.rds") %>% glimpse()
count <- readRDS("data/raw/2024_swc_state/count.rds") %>% glimpse()
length <- readRDS("data/raw/2024_swc_state/length.rds") %>% glimpse()
write.csv(metadata, "QGIS layers/state_BRUVs_to_check.csv")
#count has 2024_02
count_dat <- left_join(count, metadata, by = "sample_url")

hab_dat <- left_join(benthos, metadata, by = "sample_url")
length_dat <- left_join(length, metadata, by = "sample_url")
