# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (B. Compare communities detected)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Mess around and see
# Author:  Lise Fournier-Carnoy
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  Trying different predictors of presence/absence

# -----------------------------------------------------------------------------

library(dismo)

# Clear memory
rm(list=ls())


# Get Presence/Absence data ----

data <- readRDS("data/raw/geographe_complete_count.RDS") %>%
  filter(scientific == "Sparidae Chrysophrys auratus") %>%
  dplyr::select(longitude, latitude, count) %>%
  glimpse()

length <- readRDS("data/raw/geographe_complete_length.RDS") %>%
  filter(genus == "Chrysophrys" & !is.na(length_mm)) %>%
  dplyr::select(longitude, latitude, length_mm) %>%
  glimpse()

presdata <- data %>%
  filter(count > 0) %>%
  dplyr::select(longitude, latitude)


# Get environmental data ----

# Cropping to extent of choice
extent <- as(extent(114.9, 115.7, -33.7, -33.2), 'SpatialPolygons')
crs(extent) <- "+proj=longlat +datum=WGS84 +no_defs"

## Habitat ----

hab <- readRDS("data/spatial/rasters/wadandi_predicted-habitat.RDS")
hab <- raster::crop(hab, extent)
reef <- raster::crop(raster(hab[[1]]), extent); writeRaster(hab[[1]], "data/spatial/rasters/geographe_reef.tif", overwrite = T)
invert <- raster::crop(raster(hab[[2]]), extent); writeRaster(hab[[2]], "data/spatial/rasters/geographe_invert.tif", overwrite = T)
seagrass <- raster::crop(raster(hab[[3]]), extent); writeRaster(hab[[3]], "data/spatial/rasters/geographe_seagrass.tif", overwrite = T)
macro <- raster::crop(raster(hab[[4]]), extent); writeRaster(hab[[4]], "data/spatial/rasters/geographe_macro.tif", overwrite = T)
rock <- raster::crop(raster(hab[[5]]), extent); writeRaster(hab[[5]], "data/spatial/rasters/geographe_rock.tif", overwrite = T)
sand <- raster::crop(raster(hab[[6]]), extent); writeRaster(hab[[6]], "data/spatial/rasters/geographe_sand.tif", overwrite = T)


## Just bathymetry ----

bathy <- raster("data/spatial/rasters/GB-SW_250mBathy.tif")

extent <- as(extent(114.9, 115.7, -33.7, -33.2), 'SpatialPolygons')
crs(extent) <- "+proj=longlat +datum=WGS84 +no_defs"
bathy <- raster::crop(bathy, extent)
plot(bathy)


## Bathymetry derivatives ----

bathy_d <- readRDS("data/spatial/rasters/2014_geographe_bathymetry-derivatives.RDS")
bathy_d <- raster::crop(bathy_d, extent)
plot(bathy_d)
gadepth <- raster::crop(raster(bathy_d[[1]]), extent); writeRaster(bathy_d[[1]], "data/spatial/rasters/geographe_gadepth.tif", overwrite = T)
aspect <- raster::crop(raster(bathy_d[[2]]), extent); writeRaster(bathy_d[[2]], "data/spatial/rasters/geographe_aspect.tif", overwrite = T)
roughness <- raster::crop(raster(bathy_d[[3]]), extent); writeRaster(bathy_d[[3]], "data/spatial/rasters/geographe_roughness.tif", overwrite = T)
detrended <- raster::crop(raster(bathy_d[[4]]), extent); writeRaster(bathy_d[[4]], "data/spatial/rasters/geographe_detrended_bathy.tif", overwrite = T)


## Making a file stack for all environmental data ----

#Re
rasters <- list(reef = reef, invert = invert, seagrass = seagrass, macro = macro, rock = rock, sand = sand, aspect = aspect, gadepth = gadepth, roughness = roughness, detrended = detrended)
resampled_rasters <- lapply(rasters, function(r) resample(r, bathy))
list2env(resampled_rasters, envir = .GlobalEnv)

predictors <- stack(bathy, gadepth, aspect, roughness, detrended,
                    reef, invert, seagrass, macro, rock, sand) # Doing this otherwise the stack isn't cropped

names(predictors) <- c("bathy", "depth", "aspect", "roughness", "detrended_bathy",
                       "reef", "invert", "seagrass", "macro", "rock", "sand") # Add names for whatever variable you're adding
predictors

plot(predictors[[1]])
points(data, col = "blue")


# Extract env. data to occurrence points

presval <- extract(predictors, presdata) # Environmental values for where fish was found

# Question below: the guide I used simulated absence data (randomPoints), but I have absence data... Which do I use?
backgr <- randomPoints(predictors, 100)
backgr <- data %>%
  filter(count == 0) %>%
  dplyr::select(longitude, latitude)

absval <- extract(predictors, backgr) # Environmental values for where fish was not found.

pb <- c(rep(1, nrow(presval)), rep(0, nrow(absval))) # Makes a vector of presence (1, from our data) or absence (0, from the random points) ...
sdmdata <- data.frame(cbind(pb, rbind(presval, absval))) # ... which we can then use to put together bathymetry data of presence points and absence points
head(sdmdata)


# Model fitting ----

m1 <- glm(pb ~ bathy + reef + invert + seagrass + macro + rock + sand, data = sdmdata)
summary(m1) # glm can take presence/absence data so pb is used

# Model prediction ----

p <- predict(predictors, m1)
plot(p)

# This plot is the predicted distribution of Pink Snapper according to the predictors we have.


