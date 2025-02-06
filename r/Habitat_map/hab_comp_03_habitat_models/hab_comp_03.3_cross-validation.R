# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (C. Habitat map comp.)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Check model validity
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  Pearson's coef = NA... wtf

# -----------------------------------------------------------------------------


# Load libraries
library(tidyverse)
library(raster)
library(terra)
library(sp)
library(sf)
sf_use_s2(T)
library(ggplot2)
library(RStoolbox)
library(blockCV)
library(ggnewscale)
library(scales)
library(mgcv)
library(pROC)
library(caret)

# Set your study name
name <- "2024_geographe"
sampling_design <- "_spabal"

# Load habitat data
dat <- readRDS(paste0("data/tidy/", name, sampling_design, "_tidy_habitat.rds")) %>%
  pivot_longer(cols = c("Macroalgae", "Unconsolidated", "Consolidated", "Seagrasses", "Stony.corals"),
               values_to = "number", names_to = "taxa") %>%
  dplyr::mutate(depth = as.double(depth) * -1) %>%
  dplyr::select(taxa, total_pts, number, gadepth, aspect, roughness, detrended) %>%
  filter(!(number > total_pts)) %>% # Removing any invalid rows that have more taxa counts than total points
  glimpse()



# Read in and crop the rasters to your study extent
# the rasters should be the same as the as the values in the shapefile

stack <- readRDS(paste(paste0('data/spatial/rasters/', name, sampling_design),
                       'bathymetry-derivatives.rds', sep = "_")) %>%
  project("EPSG:9473")  # Meter projection, otherwise the blocks can't compute

n_sample <- 1000

rpc = rasterPCA(stack, nComp = 1 , spca = TRUE, nSamples = n_sample) # first PCA component of the raster stack to estimate spatial AC

poly_sf <- st_sf(geometry = st_as_sfc(st_bbox(stack))) # Convert the polygon to sf for later
aoi_r <- as(object = poly_sf, Class = "Spatial") # an sp polygon of your study areas
boss_sf <- dat %>%   # any csv with lat lon that has your samples and rasters extracted at those points
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = "epsg:4326") %>%   ###### can use 'coords = c("x", "y")' instead of wkt if there is no geometry colomn
  st_transform(9473)
?cv_spatial_autocor
# AC range estimation
range1 <- cv_spatial_autocor(
  r = rpc$map,                                                            # For some reason even if you input spatraster the output is rasterbrick
  num_sample = n_sample,
  progress = F,
  plot = T
) # the blocks keep getting smaller ??

sb1 = cv_spatial(x = boss_sf,
                 r = rpc$map,
                 size = range1$range,
                 k = 5,
                 selection = "random",
                 iteration = 100,
                 seed = 1,
                 biomod2 = TRUE,
                 hexagon = FALSE)

# creating a dataframe for the RF model

mod_df = as_Spatial(boss_sf) %>%
  data.frame() %>%
  add_column(block = sb1$folds_ids) %>%
  drop_na() ## drop any nas in your covariates

# # Loop through cross validation for each taxa (habitat) and each fold
habi <- mod_df %>%
  dplyr::select(-OpCode) %>%
  dplyr::rename(x = coords.x1, y = coords.x2) %>%
  glimpse()

resp.vars <- unique(habi$taxa)
blocks <- unique(habi$block)

preds <- readRDS(paste(paste0('data/tidy/', name, sampling_design),
                     'metadata-bathymetry-derivatives.rds', sep = "_"))
preddf <- as.data.frame(preds, xy = TRUE, na.rm = TRUE)
wgscrs <- "+proj=longlat +datum=WGS84 +south"

# Set the models manually here
resp.vars

# Loop through each taxa
# Train GAM off 4 folds, test against last fold
# Absolute value of the distance from the observed to predicted data
for (i in 1:length(resp.vars)) {
  use.dat <- habi[habi$taxa == resp.vars[i],]
  use.dat   <- as.data.frame(use.dat)
  use.dat$observed <- use.dat$number/use.dat$total_pts
  print(resp.vars[i])

  for (b in 1:length(blocks)) {
    print(blocks[b])
    train.dat <- use.dat %>% dplyr::filter(!block == blocks[b])
    test.dat  <- use.dat %>% dplyr::filter(block == blocks[b])

    # Macroalgae
    mod1 <- gam(cbind(number, total_pts - number) ~
                  s(gadepth,    k = 5, bs = "cr")  +
                  s(aspect,     k = 5, bs = "cr")  +
                  s(roughness,  k = 5, bs = "cr")  +
                  s(detrended,  k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    # Unconsolidated
    mod2 <- gam(cbind(number, total_pts - number) ~
                  s(gadepth,    k = 5, bs = "cr")  +
                  s(aspect,     k = 5, bs = "cr")  +
                  s(roughness,  k = 5, bs = "cr")  +
                  s(detrended,  k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    # Consolidated
    mod3 <- gam(cbind(number, total_pts - number) ~
                  s(gadepth,    k = 5, bs = "cr")  +
                  s(aspect,     k = 5, bs = "cr")  +
                  s(roughness,  k = 5, bs = "cr")  +
                  s(detrended,  k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    # Seagrasses
    mod4 <- gam(cbind(number, total_pts - number) ~
                  s(gadepth,    k = 5, bs = "cr")  +
                  s(aspect,     k = 5, bs = "cr")  +
                  s(roughness,  k = 5, bs = "cr")  +
                  s(detrended,  k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    # Stony.corals
    mod5 <- gam(cbind(number, total_pts - number) ~
                  s(gadepth,    k = 5, bs = "cr")  +
                  s(aspect,     k = 5, bs = "cr")  +
                  s(roughness,  k = 5, bs = "cr")  +
                  s(detrended,  k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    mod <- list(mod1, mod2, mod3, mod4, mod5)

    modpred <- cbind(preddf,
                     "predicted" = predict(mod[[i]], preddf, type = "response"))
    modpredr <- rast(modpred %>% dplyr::select(longitude, latitude, predicted))
    habi_sp <- vect(test.dat, geom = c("x", "y"), crs = wgscrs)
    habi_df   <- cbind(test.dat, terra::extract(modpredr, habi_sp)) %>%
      dplyr::filter(!is.na(predicted),
                    !is.na(observed))
    if (b == 1) {
      pearsons <- data.frame(taxa = resp.vars[i],
                             fold = blocks[b],
                             pcc = cor(habi_df$observed, habi_df$predicted, method = "pearson"))

      plot.dat.b <- habi_df
    }
    else {
      pearsons <- data.frame(taxa = resp.vars[i],
                             fold = blocks[b],
                             pcc = cor(habi_df$observed, habi_df$predicted, method = "pearson")) %>%
        bind_rows(pearsons)

      plot.dat.b <- habi_df %>%
        bind_rows(plot.dat.b)
    }
  }
  if (i == 1) {
    pearsons_table <- pearsons

    plot.dat <- plot.dat.b
  }
  else {
    pearsons_table <- pearsons %>%
      bind_rows(pearsons_table)

    plot.dat <- plot.dat.b %>%
      bind_rows(plot.dat)
  }
}

pearsons_final <- pearsons_table %>%
  dplyr::mutate(fold = as.character(fold)) %>%
  group_by(taxa) %>%
  group_modify(~ add_row(.x, pcc = mean(.$pcc), fold = "mean")) %>%
  ungroup() %>%
  pivot_wider(names_from = taxa, values_from = pcc) %>%
  glimpse()

# Save output
write.csv(pearsons_final, paste0("output/habitat/", name, "_pearsons-coefficient.csv"),
          row.names = F)

# Loop through each taxa
# Train GAM off 4 folds, test against last fold
# Accuracy - number of True positives / total observations

for (i in 1:length(resp.vars)) {
  use.dat <- habi[habi$taxa == resp.vars[i],]
  use.dat   <- as.data.frame(use.dat)
  print(resp.vars[i])

  for (b in 1:length(blocks)) {
    train.dat <- use.dat %>% dplyr::filter(!block == blocks[b])
    test.dat  <- use.dat %>% dplyr::filter(block == blocks[b])
    print(blocks[b])

    # Macroalgae
    mod1 <- gam(cbind(value, total.pts - value) ~
                  s(detrended,     k = 5, bs = "cr")  +
                  s(UCUR, k = 5, bs = "cr") +
                  s(Z, k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    # Unconsolidated
    mod2 <- gam(cbind(value, total.pts - value) ~
                  s(detrended,     k = 5, bs = "cr")  +
                  s(roughness, k = 5, bs = "cr") +
                  s(Z, k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    # Consolidated
    mod3 <- gam(cbind(value, total.pts - value) ~
                  s(detrended,     k = 5, bs = "cr")  +
                  s(SLA, k = 5, bs = "cr") +
                  s(Z, k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    # Seagrasses
    mod4 <- gam(cbind(value, total.pts - value) ~
                  s(detrended,     k = 5, bs = "cr")  +
                  s(UCUR, k = 5, bs = "cr") +
                  s(Z, k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    # Stony.corals
    mod5 <- gam(cbind(value, total.pts - value) ~
                  s(SLA, k = 5, bs = "cr") +
                  s(SST, k = 5, bs = "cr") +
                  s(Z, k = 5, bs = "cr"),
                data = train.dat, method = "REML", family = binomial("logit"))

    mod <- list(mod1, mod2, mod3, mod4, mod5)

    if (b == 1 & i == 1) {
      modpred <- cbind(preddf, "predicted" = predict(mod[[i]], preddf, type = "response")) %>%
        dplyr::mutate(block = blocks[b]) %>%
        dplyr::select(x, y, predicted, block) %>%
        dplyr::mutate(taxa = resp.vars[i])

    }
    else {
      modpred <-  cbind(preddf, "predicted" = predict(mod[[i]], preddf, type = "response")) %>%
        dplyr::mutate(block = blocks[b],
                      taxa = resp.vars[i]) %>%
        dplyr::select(x, y, predicted, block, taxa) %>%
        rbind(modpred)
    }
  }
}

# Save out so you don't need to keep running loop
saveRDS(modpred, "data/tidy/sw-network_modpreds.rds")                           # Save out to stop re running mega loop -- ignored

# Categorise raw data into habitat classes
max.dat <- habi %>%
  group_by(block) %>%
  slice(which.max(number/total_pts)) %>%
  ungroup() %>%
  dplyr::select(x, y, number, taxa, block) %>%
  dplyr::rename(observed = taxa) %>%
  glimpse()

# Categorise model predictions into habitat classes for each fold
modpred <- readRDS("data/tidy/sw-network_modpreds.rds")

# Make a raster of predicted values for each taxa per fold
for (b in 1:5) {
  print(b)
  rastdat <- modpred %>%
    dplyr::filter(block == b)
  for (i in 1:length(resp.vars)) {
    print(resp.vars[i])
    tempdat <- rastdat %>% dplyr::filter(taxa %in% resp.vars[i]) %>% dplyr::select(x, y, predicted)
    colnames(tempdat)[3] = paste(b, resp.vars[i], sep = "_")
    taxarast <- rast(tempdat)
    if (b == 1 & i == 1) {
      testrast <- taxarast
    }
    else {
      testrast <- rast(list(testrast,
                            taxarast))
    }}}

names(testrast)
plot(testrast)

# Convert raw data into spatial vector file
max.datv <- vect(max.dat, geom = c("x", "y"))
plot(testrast[[1]])
plot(max.datv, add = T) # Aligns correctly

# Load in predicted habitat from 04_predict-gam
predhabdf <- readRDS(paste0("output/habitat/", name, "_predicted-habitat.rds")) %>%   # Ignored
  dplyr::select(x, y, dom_tag_1) %>%
  dplyr::rename(dom_tag = dom_tag_1)

roc.dat <- cbind(max.dat, terra::extract(testrast, max.datv)) %>%
  dplyr::select(-c(ID, block)) %>%
  pivot_longer(cols = c(starts_with(c("1_", "2_", "3_", "4_", "5_"))),
               names_to = "block", values_to = "predicted") %>%
  separate(block, into = c("block", "taxa"), sep = "_") %>%
  pivot_wider(names_from = "taxa", values_from = "predicted") %>%
  dplyr::filter(!is.na(sand)) %>%
  glimpse()

roc.dat$predicted <- apply(roc.dat %>% dplyr::select(sand, inverts, macroalgae, rock, seagrass), 1,
                           FUN = function(x){names(which.max(x))})
unique(roc.dat$predicted)


# Loop through each spatial block
# Generate kappa and auc
for (i in 1:5) {
  test.dat <- roc.dat %>%
    dplyr::filter(block == blocks[i]) %>%
    dplyr::mutate(observed = as.factor(observed),
                  predicted = as.factor(predicted))
  levels(test.dat$predicted) <- levels(test.dat$observed)

  test.matrix <- as.matrix(test.dat %>%
                             dplyr::select(sand, rock, macroalgae, inverts, seagrass))
  conf = confusionMatrix(test.dat$predicted, test.dat$observed)
  kapp = conf$overall[1]

  roc.hab <- multiclass.roc(response = test.dat$observed, predictor = test.matrix) # Second should be numeric - predicted value

  if (i == 1) {
    kappa.dat <- data.frame(kappa = kapp,
                            auc = roc.hab$auc) %>%
      dplyr::mutate(block = blocks[i])
  }
  else {
    kappa.dat <- data.frame(kappa = kapp,
                            auc = roc.hab$auc) %>%
      dplyr::mutate(block = blocks[i]) %>%
      rbind(kappa.dat)
  }

}

kappa.dat <- kappa.dat %>%
  dplyr::mutate(block = as.character(block),
                auc = as.numeric(auc)) %>%
  add_row(block = "mean",
          kappa = mean(.$kappa),
          auc = mean(.$auc))

write.csv(kappa.dat, paste0("output/habitat/", name, "_kappa-auc.csv"),
          row.names = F)
