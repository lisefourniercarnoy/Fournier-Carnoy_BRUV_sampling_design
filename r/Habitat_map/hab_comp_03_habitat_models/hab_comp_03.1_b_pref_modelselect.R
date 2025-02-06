# -----------------------------------------------------------------------------

# Project: Geographe Bay BRUV sampling design comparison (C. Habitat map comp.)
# Data:    2014 MEGlab BRUV habitat data. (TO BE CHANGED TO 2024 DATA)
# Task:    Select models to predict habitats.
# Author:  Lise Fournier-Carnoy / adapted from Claude Spencer
# Date:    August 2024

# -----------------------------------------------------------------------------

# Status:  Spatially balanced models all good, need preferential data.

# -----------------------------------------------------------------------------


# Load libraries
library(tidyverse)
library(mgcv)
library(MuMIn)
library(FSSgam)
library(CheckEM)
library(patchwork)
library(ggplot2)

rm(list=ls())

name <- "2024_geographe"
sampling_design <- "_pref"


# PREFERENTIAL ----

## Load & format the data ----

habi <- readRDS(paste0("data/tidy/", name, sampling_design, "_tidy_habitat.rds")) %>%
  pivot_longer(cols = c("Macroalgae", "Unconsolidated", "Consolidated", "Seagrasses"),
               values_to = "number", names_to = "taxa") %>%
  dplyr::mutate(depth = as.double(depth) * -1) %>%
  filter(longitude > 115.3) %>%
  glimpse()

ggplot(data = habi, aes(x = depth, y = gadepth)) +
  geom_point() +
  geom_smooth(method = "lm")
plot(habi$longitude, habi$latitude)
habi <- habi[!(habi$number > habi$total_pts | habi$number < 0 | habi$total_pts < 0), ] # Removing invalid rows that were making the GAMs messed up.


## Set predictor variables ----

names(habi)
pred.vars <- c("gadepth", "aspect", "roughness", "detrended")

# Check for correlation of predictor variables- remove anything highly correlated (>0.95)
round(cor(habi[ , pred.vars]), 2) # All good, no high correlations
plot_transformations(pred.vars = pred.vars, dat = habi)

# Check to make sure Response vector has not more than 80% zeros
unique.vars = unique(as.character(habi$taxa))

unique.vars.use = character()
for(i in 1:length(unique.vars)){
  temp.dat = habi[which(habi$taxa == unique.vars[i]),]
  if(length(which(temp.dat$taxa == 0))/nrow(temp.dat)<0.9){
    unique.vars.use = c(unique.vars.use,unique.vars[i])}
}

unique.vars.use # All remain


## Run the full subset model selection ----
outdir    <- ("outputs/Habitat_comparison/preferential/")
resp.vars <- unique.vars.use
out.all   <- list()
var.imp   <- list()


## Loop through the FSS function for each Abiotic taxa ----
?generate.model.set
for(i in 1:length(resp.vars)){
  print(resp.vars[i])
  use.dat <- habi[habi$taxa == resp.vars[i],]
  use.dat   <- as.data.frame(use.dat)
  Model1  <- gam(cbind(number, (total_pts - number)) ~ # Success and failure counts for each habitat
                   s(gadepth, bs = 'cr', k = 3),
                 family = binomial("logit"),  data = use.dat)

  model.set <- generate.model.set(use.dat = use.dat,
                                  test.fit = Model1,
                                  pred.vars.cont = pred.vars,
                                  cyclic.vars = c("aspect"),
                                  k = 5,
                                  cov.cutoff = 0.7,
                                  max.predictors = 3
  )
  out.list <- fit.model.set(model.set,
                            max.models = 600,
                            parallel = T)
  names(out.list)

  out.list$failed.models # examine the list of failed models
  mod.table <- out.list$mod.data.out  # look at the model selection table
  mod.table <- mod.table[order(mod.table$AICc), ]
  mod.table$cumsum.wi <- cumsum(mod.table$wi.AICc)
  out.i     <- mod.table[which(mod.table$delta.AICc <= 2), ]
  out.all   <- c(out.all, list(out.i))
  var.imp   <- c(var.imp, list(out.list$variable.importance$aic$variable.weights.raw))



  # plot the best models
  for(m in 1:nrow(out.i)){
    best.model.name <- as.character(out.i$modname[m])

    png(file = paste(outdir, m, resp.vars[i], "mod_fits.png", sep = ""))
    if(best.model.name != "null"){
      par(mfrow = c(3, 1), mar = c(9, 4, 3, 1))
      best.model = out.list$success.models[[best.model.name]]
      plot(best.model, all.terms = T, pages = 1, residuals = T, pch = 16)
      mtext(side = 2, text = resp.vars[i], outer = F)}
    dev.off()
  }
}


## Model fits and importance ----

names(out.all) <- resp.vars
names(var.imp) <- resp.vars
all.mod.fits <- list_rbind(out.all, names_to = "response")
all.var.imp  <- do.call("rbind", var.imp)
write.csv(all.mod.fits[ , -2], file = paste0(outdir, name, sampling_design, "_abiotic_all.mod.fits.csv"))
write.csv(all.var.imp,         file = paste0(outdir, name, sampling_design, "_abiotic_all.var.imp.csv"))

