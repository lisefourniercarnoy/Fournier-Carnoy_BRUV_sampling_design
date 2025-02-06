# Following https://rspatial.org/raster/sdm/1_sdm_introduction.html as an initial guide

remove(list=ls()) # Clean workspace env.

library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(sf)
library(raster)
library(dismo)


# Get Presence/Absence data ----

data <- readRDS("data/raw/geographe_complete_count.RDS") %>%
  filter(scientific == "Sparidae Chrysophrys auratus") %>%
  dplyr::select(longitude, latitude, count)
data

presdata <- data %>%
  filter(count > 0) %>%
  dplyr::select(longitude, latitude)

# Get environmental data ----

# Cropping Bathymetry to extent of choice

bathy <- raster("data/spatial/rasters/GB-SW_250mBathy.tif")

extent <- as(extent(114.9, 115.7, -33.7, -33.2), 'SpatialPolygons')
crs(extent) <- "+proj=longlat +datum=WGS84 +no_defs"
bathy <- raster::crop(bathy, extent)

plot(bathy)

# Making a file stack for all environmental data

path <- file.path(paste0("data/spatial/rasters"))
files <- list.files(path, pattern='tif$', full.names=TRUE) # This makes a list of all file of .tif format in the folder specified above ...
predictors <- stack(files) # ... And this makes a stack of them

predictors <- stack(bathy) # Doing this otherwise the stack isn't cropped

names(predictors) <- c("bathy") # Add names for whatever variable you're adding
predictors

plot(predictors, add = T)
points(data, col = "blue")


# Extract env. data to occurrence points

presval <- extract(predictors, presdata) # Bathymetry values for where fish was found

backgr <- randomPoints(predictors, 500)
absval <- extract(predictors, backgr) # Bathymetry value for where fish was not found (using random points)

pb <- c(rep(1, nrow(presval)), rep(0, nrow(absval))) # Makes a vector of presence (1, from our data) or absence (0, from the random points) ...
sdmdata <- data.frame(cbind(pb, rbind(presval, absval))) # ... which we can then use to put together bathymetry data of presence points and absence points
head(sdmdata)


# Model fitting ----

m1 <- glm(pb ~ bathy, data = sdmdata)
summary(m1) # glm can take presence/absence data so pb is used

bc <- bioclim(predictors, presdata) # bioclim uses presence-only data
response(bc) # Distribution of predicted values across my env. variable according to bioclim model


# Model prediction ----

p <- predict(predictors, m1)
plot(p)

# Okay so based on just a bathymetry variable, predict() gives me a map of suitability for my species
# With just bathymetry it looks like a big chunk of the region is suitable, but hopefully habitat makes it more interesting.


# Model evaluation ----

# We need to check the quality of the predictions.
# The Area Under the Receiver Operator Curve (AUROC or AUC for short) is a measure of rank-correlation.
# In unbiased data, high AUC = predicted presence and observed presence match up very well.
# AUC of 0.5 = predicted presence is a random guess.


par(mfrow=c(1, 2))
plot(sort(absval), col='red', pch=21)
points(sort(presval), col='blue', pch=24)
legend(200, 0.75 * min(absval,presval), c('absence', 'presence'),
       pch=c(21,24), col=c('red', 'blue'))
comb <- c(absval, presval)
group <- c(rep('presence', length(presval)), rep('absence', length(absval)))
boxplot(comb~group, col=c('red', 'blue'))


e <- evaluate(p = as.vector(presval), a = as.vector(absval))
e # Here we have the AUC that tells us about the difference between presence and absence data groups.
# AUC is ~0.66 which is not very good but that's expected with only bathymetry as a variable.

par(mfrow=c(1, 2))
density(e) ; boxplot(e, col=c('blue', 'red'))

samp <- sample(nrow(sdmdata), round(0.75 * nrow(sdmdata))) # Sub-sampling a random portion of our presence/absence data
traindata <- sdmdata[samp,] # Selecting the rows that correspond to the random selection
traindata <- traindata[traindata[,1] == 1, 1:2] # Selecting presence point within presence/absence and their corresponding env. variable values. Change 1:2 to howver many column you have ?
testdata <- sdmdata[-samp,] # Getting the inverse of the samp selection ?
bc <- bioclim(traindata)
e <- evaluate(testdata[testdata==1,], testdata[testdata==0,], bc)
e

par(mfrow=c(1, 1))
plot(e, 'ROC') # Okay this is saying AUC is very high, basically our model is really good.. dodgy

# Just above we separated the data into 2 groups, randomly.
# It's better to partition your data using kfold():
# What it does is assign a value 1:k to the rows of your presence matrix

pres <- sdmdata[sdmdata[,1] == 1, 1:2]
back <- sdmdata[sdmdata[,1] == 0, 1:2]

k <- 5 # Number of partitions
group <- kfold(presdata, k) # Partitioning

e <- list() # What this does is make k models (5 models) using the subsampled presence/absence data
for (i in 1:k) {
  train <- presdata[group != i,]
  test <- presdata[group == i,]
  bc <- bioclim(train)
  e[[i]] <- evaluate(p=test, a=train, bc)
}

auc <- sapply(e, function(x){x@auc})
auc # And we get 5 AUC values ...
mean(auc)

sapply( e, function(x){ threshold(x)['spec_sens'] } ) # 'maximum of the sum of the sensitivity (true positive rate) and specificity (true negative rate)' idk what that is exactly

# One problem with AUC is that it is relative: depending on the spatial extent used to select background points.
# Generally, bigger extents = bigger AUC. SoAUCs are biased and can't be compared.
# We can remove this 'spatial sorting bias' (SSB, diff. between testing-presence to training-presence and distance between testing-absence and training-absence) through 'point-wise distance sampling'

nr <- nrow(data)
s <- sample(nr, 0.25 * nr)
pres_train <- data[-s, ]
pres_test <- data[s, ]
nr <- nrow(backgr)
set.seed(9)
s <- sample(nr, 0.25 * nr)
back_train <- backgr[-s, ]
back_test <- backgr[s, ]

sb <- ssb(pres_test, back_test, pres_train)
sb[,1] / sb[,2] # This is a measure of SSB. The closer to 1, the less SSB there is. The closer to 0, the stronger the SSB.
# Using the full data, SSB is strong with a value of 0.205.
# Below we'll reduce that SSB.

i <- pwdSample(pres_test, back_test, pres_train, n=1, tr=0.1)
pres_test_pwd <- pres_test[!is.na(i[,1]), ]
back_test_pwd <- back_test[na.omit(as.vector(i)), ]
sb2 <- ssb(pres_test_pwd, back_test_pwd, pres_train)
sb2[1]/ sb2[2] # SSB much reduced, value very close to 1.


bc <- bioclim(predictors, pres_train[, 1:2])
evaluate(bc, p=pres_test[, 1:2], a=back_test, x=predictors) # AUC with no reduction of SSB is ~0.65
evaluate(bc, p=pres_test_pwd[, 1:2], a=back_test_pwd, x=predictors) # AUC with reduction of SSB is ~0.51

# This SSB reduction avoids false positive interpretation of variables' association with species presence.
#woohoo
