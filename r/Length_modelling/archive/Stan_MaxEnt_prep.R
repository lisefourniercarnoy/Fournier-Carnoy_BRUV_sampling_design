# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (A. Simulated sizes)
# Data:    Spatial data
# Task:    Adapt ENVT5566 labs to Geographe data
# Author:  Lise Fournier-Carnoy / Adapted from ENVT5566 Lab b
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  Just started

# -----------------------------------------------------------------------------


# Clear memory
rm(list=ls())

require(raster)

## Species data (spatially balanced) ----

# Load presence data
data <- readRDS("data/raw/geographe_complete_count.RDS") %>%
  filter(scientific == "Sparidae Chrysophrys auratus") %>%
  dplyr::select(longitude, latitude, count)
data

presdata <- data %>%
  filter(count > 0) %>%
  dplyr::select(longitude, latitude)

# We want to subset the data to train the model.
samp <- sort(sample(nrow(presdata), nrow(presdata)*.7)) # We'll chose a random 70% of the data.
train <- presdata[samp,] # Create a training dataset from the 70% of the data
test <- presdata[-samp,] # And a "test" dataset from the remaining 30% of the data

# MaxEnt needs the dataset to be a certain format.
pink <- train; pink$Species <- "Chrysophrys_auratus" # Add a species name
pink <- pink[, c(3,1:2)] # Reorder the columns so that species is the first, then longitude, then latitude.
plot(pink)

# Save train and test data
write.csv(pink, file="C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/Pink_Maxent.csv", row.names = FALSE)
write.csv(test, file="C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/Pink_Maxent_test_data.csv", row.names = FALSE)


# ------------------------------------------------------------------------------

## Environmental data ----

# Rugosity Variables from Depth ----

require(raster)

# Select depth data
extent <- as(extent(114.9, 115.7, -33.7, -33.2), 'SpatialPolygons')
dat <- raster("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/GB-SW_250mBathy.tif")
dat <- raster::crop(dat, extent)
plot(dat)

# Convert into an .asc file, needed for MaxEnt
writeRaster(dat, filename="C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/depth.asc",
            datatype='FLT4S', format="ascii", overwrite=TRUE)

## some kernel sizes what ever you want here as odd number
## but large Kernels run slowly

msize <- c(3,5,9)

## some kernel functions
krange <- function(x){
  abs(min(x) - max(x))}

kfun <- c(krange,sd)
## function names as the functions can be passed as a string
kfunname <- c("krange","sd")
## counter to loop through function names
counter <- 1

for (j in kfun){
  for(i in msize){
    # run the kernel
    dattemp <- focal(dat,w=matrix(1,i,i),fun=j)
    #file names
    filename1 <- paste0("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/kernel_", kfunname[counter], "_pix", i, ".asc")
    filename2 <- paste0("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/kernel_", kfunname[counter], "_pix", i, "_residuals_from_depth.asc")
    # write out the files
    writeRaster(dattemp,filename=filename1,datatype='FLT4S',format="ascii",overwrite=TRUE)
    writeRaster(dat-dattemp,filename=filename2,datatype='FLT4S',format="ascii",overwrite=TRUE)
  }
  counter <- counter + 1
}

neighbours <- c(4,8)

for(i in neighbours) {

  # filenames
  fnslope <- paste0("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/kernel_slope_pix_", i, ".asc")
  fnaspect <- paste0("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/kernel_aspect_pix_", i, ".asc")
  # topographic postion index
  fntpi <- paste0("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/kernel_tpi_pix_", i, ".asc")
  # caculate measures
  try(tslp <- terrain(dat, opt='slope', unit='radians', neighbors=i))
  try(tasp <- terrain(dat, opt='aspect', unit='radians', neighbors=i))
  ttpi <- terrain(dat, opt='TPI', unit='radians', neighbors=i)
  #write files
  try(writeRaster(tslp,filename=fnslope,datatype='FLT4S',format="ascii",overwrite=TRUE))
  try(writeRaster(tasp,filename=fnaspect,datatype='FLT4S',format="ascii",overwrite=TRUE))
  writeRaster(ttpi,filename=fntpi,datatype='FLT4S',format="ascii",overwrite=TRUE)
}

# Have a look at the slope and aspect, calculated from depth
par(mfrow = c(2, 2))
plot(tslp); plot(tasp); plot(ttpi)


# Habitat ----

macro <- raster("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/geographe_macro.tif")
invert <- raster("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/geographe_invert.tif")
reef <- raster("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/geographe_reef.tif")
rock <- raster("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/geographe_rock.tif")
sand <- raster("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/geographe_sand.tif")
seagrass <- raster("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/geographe_seagrass.tif")

habitat <- stack(reef, invert, seagrass, macro, rock, sand) # Doing this otherwise the stack isn't cropped


# Convert into an .asc file, needed for MaxEnt
writeRaster(macro, filename="C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/macro.asc",
            datatype='FLT4S', format="ascii", overwrite=TRUE)
writeRaster(invert, filename="C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/invert.asc",
            datatype='FLT4S', format="ascii", overwrite=TRUE)
writeRaster(reef, filename="C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/reef.asc",
            datatype='FLT4S', format="ascii", overwrite=TRUE)
writeRaster(rock, filename="C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/rock.asc",
            datatype='FLT4S', format="ascii", overwrite=TRUE)
writeRaster(sand, filename="C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/sand.asc",
            datatype='FLT4S', format="ascii", overwrite=TRUE)
writeRaster(seagrass, filename="C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/seagrass.asc",
            datatype='FLT4S', format="ascii", overwrite=TRUE)

## some kernel sizes what ever you want here as odd number
## but large Kernels run slowly

msize <- c(3,5,9)

## some kernel functions
krange <- function(x){
  abs(min(x) - max(x))}

kfun <- c(krange,sd)
## function names as the functions can be passed as a string
kfunname <- c("krange","sd")
## counter to loop through function names
counter <- 1

for (layer_index in 1:nlayers(habitat)) {
  # Extract the current raster layer
  current_raster <- habitat[[layer_index]]

  # Iterate over kernel functions and sizes
  for (j in kfun) {
    for (i in msize) {
      # Run the kernel function on the current raster layer
      dattemp <- focal(current_raster, w = matrix(1, i, i), fun = j)

      # Generate file names
      filename1 <- paste0("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/kernel_", kfunname[counter], "_pix", i, "_layer", names(habitat[[layer_index]]), ".asc")
      filename2 <- paste0("C:/Users/localuser/Documents/Fournier-Carnoy_Geographe/data/spatial/rasters/MaxEnt_files/kernel_", kfunname[counter], "_pix", i, "_layer", names(habitat[[layer_index]]), "_residuals_from_hab.asc")

      # Write out the files
      writeRaster(dattemp, filename = filename1, datatype = 'FLT4S', format = "ascii", overwrite = TRUE)
      writeRaster(current_raster - dattemp, filename = filename2, datatype = 'FLT4S', format = "ascii", overwrite = TRUE)
    }
    counter <- counter + 1
  }
  counter <- 1  # Reset counter for the next layer
}

# Have a look at the habitat rasters
par(mfrow = c(2, 2))
plot(tslp); plot(tasp); plot(ttpi)

# Before running MaxEnt, open the .asc files from a habitat file, anc copy-paste the NCOLS, NROWS, XLLCORNER, YLLCORNER, and CELLSIZE lines into every output from the depth derivatives.
# I don't know why, but if not, MaxEnt thinks that the depth files and the habitat files are different dimentions... Despite cropping them earlier


