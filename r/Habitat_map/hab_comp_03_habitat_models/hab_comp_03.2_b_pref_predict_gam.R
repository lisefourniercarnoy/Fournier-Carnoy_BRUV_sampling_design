# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (C. Habitat map comp.)
# Data:    2024 DBCA BRUV habitat data.
# Task:    Predict habitats using GAMs.
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  Got predictions for spa-bal, models needs tweaking + needs preferential data.

# -----------------------------------------------------------------------------


# Load libraries
library(reshape2)
library(mgcv)
library(ggplot2)
library(viridis)
library(terra)
library(predicts)
library(tidyverse)
library(sf)
library(nlraa)
library(patchwork)

rm(list = ls())

name <- "2024_geographe"
sampling_design <- "_pref"

# SPATIALLY BALANCED ----

## Load & format the data ----

habi <- readRDS(paste0("data/tidy/", name, sampling_design, "_tidy_habitat.rds")) %>%
  glimpse()
  #filter(OpCode != "12RC9") # Also removing one OpCode with point count issues.


preds <- readRDS(paste0("data/spatial/rasters/", name, sampling_design, "_bathymetry-derivatives.rds"))
plot(preds)

preddf <- as.data.frame(preds, xy = TRUE, na.rm = TRUE)


## Make and select top models for each habitat ----
#Not sure my model selection is correct ?

# Rock
m_rock <- gam(cbind(Consolidated, total_pts - Consolidated) ~ # Literally none are sig. ???
                s(aspect, k = 5, bs = "cr")  +
                s(detrended, k = 5, bs = "cr")  +
                s(gadepth, k = 5, bs = "cr") +
                s(roughness, k = 5, bs = "cr"),
              data = habi, method = "REML", family = binomial("logit"))
summary(m_rock)
plot(m_rock, pages = 1, residuals = T, cex = 5)

# Macroalgae
m_macro <- gam(cbind(Macroalgae, total_pts - Macroalgae) ~
                 s(aspect, k = 5, bs = "cr")  +
                 s(detrended, k = 5, bs = "cr")  +
                 s(gadepth, k = 5, bs = "cr") +
                 s(roughness, k = 5, bs = "cr"),
               data = habi, method = "REML", family = binomial("logit"))
summary(m_macro)
plot(m_macro, pages = 1, residuals = T, cex = 5)

# Seagrass
m_seagrass <- gam(cbind(Seagrasses, total_pts - Seagrasses) ~
                    s(aspect, k = 5, bs = "cr")  +
                    #s(detrended, k = 5, bs = "cr")  # not significant
                    #s(gadepth, k = 5, bs = "cr") +  # not significant
                    s(roughness, k = 5, bs = "cr"),
                  data = habi, method = "REML", family = binomial("logit"))
summary(m_seagrass)
plot(m_seagrass, pages = 1, residuals = T, cex = 5)

# Sand
m_sand <- gam(cbind(Unconsolidated, total_pts - Unconsolidated) ~
                s(aspect, k = 5, bs = "cr")  +
                s(detrended, k = 5, bs = "cr")  +
                s(gadepth, k = 5, bs = "cr") +
                s(roughness, k = 5, bs = "cr"),
              data = habi, method = "REML", family = binomial("logit"))
summary(m_sand)
plot(m_sand, pages = 1, residuals = T, cex = 5)


## predict, rasterise and plot ----
preddf <- cbind(preddf,
                "pmacro" = predict(m_macro, preddf, type = "response", se.fit = T),
                "prock" = predict(m_rock, preddf, type = "response", se.fit = T),
                "psand" = predict(m_sand, preddf, type = "response", se.fit = T),
                "pseagrass" = predict(m_seagrass, preddf, type = "response", se.fit = T)) %>%
  glimpse()

prasts <- rast(preddf %>% dplyr::select(x, y, pmacro.fit, prock.fit, psand.fit, pseagrass.fit),
               crs = crs(preds))
par(mfrow = c(2, 3)); plot(prasts)
summary(prasts)

# Need to fix this bit at a later stage

# Calculate MESS and mask predictions
# xy <- habi %>%
#   dplyr::select(longitude , latitude) %>%
#   glimpse()
#
# resp.vars <- c("pmacro", "prock", "psand", "pseagrass", "pinverts", "preef")
#
# for (i in 1:length(resp.vars)) {
#   print(resp.vars[i])
#   mod <- get(str_replace_all(resp.vars[i], "p", "m_"))
#
#   temppred <- preddf %>%
#     dplyr::select(x, y, resp.vars[i]) %>%
#     rast(crs = "epsg:4326")
#
#   dat <- terra::extract(subset(preds, names(mod$model)[2:length(names(mod$model))]), xy) %>%
#     dplyr::select(-ID)
#   messrast <- predicts::mess(subset(preds, names(mod$model)[2:length(names(mod$model))]), dat) %>%
#     terra::clamp(lower = -0.01, values = F)
#   messrast <- terra::crop(messrast, temppred)
#   temppred_m <- terra::mask(temppred, messrast)
#
#
#   if (i == 1) {
#     preddf_m <- as.data.frame(temppred_m, xy = T)
#   }
#   else {
#     preddf_m <- as.data.frame(temppred_m, xy = T) %>%
#       full_join(preddf_m)
#   }
# }
#
# glimpse(preddf_m)

# Categorise by dominant tag
preddf$dom_tag <- apply(preddf %>% dplyr::select(pmacro.fit, prock.fit, psand.fit, pseagrass.fit), 1,
                        FUN = function(x){names(which.max(x))})
unique(preddf$dom_tag)

saveRDS(preddf, paste0("outputs/Habitat_comparison/", name, sampling_design, "_predicted-habitat.rds"))
