# Initial messing around with BRUV Geographe data

remove(list=ls()) # Clean workspace env.

library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(sf)
library(readr)
library(patchwork)


# Setting up files ----

## State Metadata ----

state22md <- read.table("data/tidy/2022-03_NgariCapes.MP.Monitoring_stereoBRUVs_tidy_Metadata", sep = ",", header = T)


## State Lengths ----

state22l <- read.table("data/tidy/2022-03_NgariCapes.MP.Monitoring_stereoBRUVs_tidy_Lengths.csv", sep = ",", header = T)

state22l$fullspp <- paste(state22l$Genus, state22l$Species)

names(state22l)[1] <- paste("opcode") # Renaming the column so it's the same as the metadata
state22l <- merge(state22l, state22md[,c("opcode", "status", "longitude_dd", "latitude_dd")], by = "opcode") # Adding the status (Fished/SanctuaryZone) to the lengths data


## State MaxN ----

state22c <- read.table("data/tidy/2022-03_NgariCapes.MP.Monitoring_stereoBRUVs_tidy_Points.csv", sep = ",", header = T)

state22c$fullspp <- paste(state22c$Genus, state22c$Species)

state22maxn <- state22c %>%
  group_by(Frame, fullspp, opcode) %>%
  summarise(Count = n()) %>%
  ungroup() # Counting each frame's number of each species == MaxN


## Federal Metadata ----

fed14md <- readRDS("data/raw/geographe_metadata.RDS")


## Federal Lengths ----

fed14l <- readRDS("data/raw/geographe_complete_length.RDS")

fed14l$fullspp <- paste(fed14l$genus, fed14l$species)


## Federal MaxN ----

fed14c <- readRDS("data/raw/geographe_complete_count.RDS")

fed14c$fullspp <- paste(fed14c$genus, fed14c$species)

fed14maxn <- fed14c # Federal data is already in MaxN format


# Initial Exploration ----

## Fish Length ----

# I want to compare the distributions and the sensitivity of BS and P sampling on fish length
# I spent hours on sommething i won't need in the end. cuz im not comparing real data for length

### Visualisation ----

# I want a histogram per species, with length bins. Overlay one SD with the other to see whether there's a diff.
# To do that I need to clean up the bloody dataframes

SB_l <- fed14l %>%
  dplyr::select(sample, length_mm, longitude, latitude, status, fullspp) %>%
  filter(!is.na(length_mm)) %>%
  mutate(sample = gsub("[^0-9]", "", sample), # Simplify the opcode column
         sd = "SB") %>% # Identify which SD this is
  rename(opcode = sample, length = length_mm) # Renaming column names to rbind with P_l


P_l <- state22l %>%
  dplyr::select(opcode, Length, fullspp, status, longitude_dd, latitude_dd) %>%
  rename(longitude = longitude_dd, latitude = latitude_dd, length = Length) %>% # Renaming column names to rbind with SB_l
  mutate(sd = "P") # Identify which SD this is

full_l <- rbind(SB_l, P_l)


# Selecting species that have enough lengths

species_counts <- full_l %>%
  group_by(fullspp) %>%
  summarise(count = n()) %>%
  ungroup()

sp <- species_counts %>%
  group_by(fullspp) %>%
  filter(all(count >= 40)) %>%
  pull(fullspp) %>%
  unique()

plot_list <- vector("list", length(sp))
for (i in seq_along(sp)) { # For all the species in the selected species vector ...

  p <- ggplot(full_l[full_l$fullspp == sp[i],], aes(x = length, fill = sd)) + # ... Make a plot ...
    geom_histogram(data=subset(full_l[full_l$fullspp == sp[i],], sd == 'SB'), fill = "red", alpha = 0.2) +
    geom_histogram(data=subset(full_l[full_l$fullspp == sp[i],], sd == 'P'), fill = "blue", alpha = 0.2) +
    facet_wrap(~ factor(fullspp), scales = "free_x") +
    theme_minimal()


  plot_list[[i]] <- p # ... And store it in an object.

}

combined_plot <- wrap_plots(plot_list, ncol = 5) # Combine all plots into a single layout
print(combined_plot)


### Stats tests ?? ----


## Community Composition ----

### pref NMDS ----

# P design first
# Okay I want to see whether the community of fish detected is different in SZ and fished areas
# First I need to format my data so that rows are the drops and columns are species with counts

SZ_opcodes <- state22md %>% # Finding the opcodes that are in SZs
  filter(status == "No-take") %>%
  pull(opcode)

df_matrix <- state22maxn %>% # Making an OTU dataframe.
  pivot_wider(names_from = fullspp, values_from = Count, values_fill = 0) %>%
  group_by(opcode) %>%
  summarise_all(.funs = sum, na.rm = TRUE) %>%
  ungroup()

df_matrix <- as.data.frame(df_matrix) %>% # Adding fished/unfished status from opcodes.
  mutate(status = ifelse(opcode %in% SZ_opcodes, "SZ", "Fished"))
rownames(df_matrix) <- df_matrix$opcode

df_matrix <- df_matrix[, 3:122] # Removes the first columns ('Frame' and 'opcode") which we don't need anymore.


# Now doing the NMDS

dist_matrix <- vegdist(df_matrix[, 1:119], method = "bray")
nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100) # Stress is below 0.2 - all good
nmds_coords <- data.frame(nmds_result$points)
status <- df_matrix$status

nmds_coords$site <- rownames(nmds_coords) # Adding a status column to see difference between community in fished/SZ sites
nmds_coords$status <- ifelse(nmds_coords$site %in% SZ_opcodes, "SZ", "fished")


ggplot(nmds_coords, aes(x = MDS1, y = MDS2, colour = status)) +
  geom_point(size = 3) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_minimal()


#### pref PERMANOVA ----

permanova_result <- adonis2(dist_matrix ~ status, data = df_matrix)
print(permanova_result)

### Combined dataset ----

# I want to PERMANOVA from a combined dataset


SB_c <- fed14maxn %>%
  dplyr::select(sample, fullspp, count, status) %>%
  filter(!is.na(count)) %>%
  mutate(sample = gsub("[^0-9]", "", sample), # Simplify the opcode column
         sd = "SB") %>% # Identify which SD this is
  rename(opcode = sample) # Renaming column names to rbind with P_c.

P_c <- state22maxn %>%
  dplyr::select(opcode, fullspp, Count) %>%
  rename(count = Count) %>% # Renaming column names to rbind with SB_c.
  mutate(sd = "P",
         status = ifelse(opcode %in% SZ_opcodes, "SZ", "Fished"), # Identify which SD this is.
         count = as.numeric(count))

full_c <- rbind(SB_c, P_c)

full_c$count <- as.integer(full_c$count)

df_matrix <- full_c %>% # Making an OTU dataframe.
  pivot_wider(names_from = fullspp, values_from = count, values_fill = list(count = 0)) %>%
  group_by(opcode) %>%
  summarise_all(.funs = sum, na.rm = TRUE) %>%
  ungroup()

df_matrix <- as.data.frame(df_matrix) %>% # Adding fished/unfished status from opcodes.
  mutate(status = ifelse(opcode %in% SZ_opcodes, "SZ", "Fished"))
rownames(df_matrix) <- df_matrix$opcode

df_matrix <- df_matrix[, 3:122] # Removes the first columns ('Frame' and 'opcode") which we don't need anymore.


# Now doing the NMDS

dist_matrix <- vegdist(df_matrix[, 1:119], method = "bray")
nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100) # Stress is below 0.2 - all good
nmds_coords <- data.frame(nmds_result$points)
status <- df_matrix$status

nmds_coords$site <- rownames(nmds_coords) # Adding a status column to see difference between community in fished/SZ sites
nmds_coords$status <- ifelse(nmds_coords$site %in% SZ_opcodes, "SZ", "fished")


ggplot(nmds_coords, aes(x = MDS1, y = MDS2, colour = status)) +
  geom_point(size = 3) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_minimal()


#### PERMANOVA ----

permanova_result <- adonis2(dist_matrix ~ status, data = df_matrix)
print(permanova_result)



### Spatially Balanced Sampling Design ----

species_counts <- fed14l %>%
  group_by(fullspp) %>%
  summarise(count = n()) %>%
  ungroup()

filtered_species <- species_counts %>%
  group_by(fullspp) %>%
  filter(all(count >= 10)) %>%
  pull(fullspp) %>%
  unique() # Filtering species that have n > 10 in both Sanctuary Zone and Fished areas

filtered_df <- fed14l %>%
  filter(fullspp %in% filtered_species) # Filtering out

SZ_longitudes <- c(115.4117, 115.422, 115.3901, 115.4144, 115.104, 115.2039)
filtered_df$status <- ifelse(filtered_df$longitude %in% SZ_longitudes, "SZ", "fished")
# There are 6 drops that are in the (future) Sanctuary Zones. They are marked as SZ here

ggplot(filtered_df) +
  geom_boxplot(aes(y = length, fill = as.factor(status)), varwidth = T) +
  facet_wrap(~fullspp, scales = "free")

#### NMDS ----

# Copied and slightly modified from above
# Okay I want to see whether the community of fish detected is different in SZ and fished areas
# First I need to format my data so that rows are the drops and columns are species with counts

fed14maxn <- fed14maxn %>%
  mutate(sample = as.factor(sample))

df_matrix <- fed14maxn %>%
  pivot_wider(names_from = fullspp, values_from = count, values_fill = 0) %>% # Makes each species a new column with a count for that species in that opcode
  group_by(sample) %>%
  summarise_all(.funs = sum, na.rm = TRUE) %>% # Sums up each opcode's count for each species
  ungroup()

df_matrix <- df_matrix %>%
  mutate(status = ifelse(longitude %in% SZ_longitudes, "SZ", "Fished"))

df_matrix <- as.data.frame(df_matrix)
rownames(df_matrix) <- df_matrix$sample # Changes the name of the rows to the opcode
df_matrix <- df_matrix[, 31:150] # Removes the first column ('Frame') which we don't need anymore
names(df_matrix)

# Now doing the NMDS

dist_matrix <- vegdist(df_matrix[1:120], method = "bray")
nmds_result <- metaMDS(dist_matrix, k = 2, trymax = 100) # Stress is below 0.2 - all good
nmds_coords <- data.frame(nmds_result$points)

nmds_coords$site <- rownames(nmds_coords) # Adding a status column to see difference between community in fished/SZ sites
nmds_coords$status <- ifelse(nmds_coords$site %in% SZ_opcodes, "SZ", "fished")


ggplot(nmds_coords, aes(x = MDS1, y = MDS2, colour = status)) +
  geom_point(size = 3) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_minimal()
