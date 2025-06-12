
### BACI comparison -----------------------------------------------------------

# I can't really know whether the sampling design detects the SZ effect or just
# a difference in habitat that happens to be inside a SZ. For example, my fake-
# SZ has a larger reef area (complex strata) than the non-SZ bits outside it,
# and there may be more fish there naturally.
# If we test the ability of sampling designs to detect difference in abundance
# before and after the SZ (BACI design) we can actually say whether the SZ
# effect is detected correctly.

# So. I need to compare the ratio of ~before~ and the ratio of ~after~. There
# should definitely be a difference.

names(dat)
ratio_before <- dat %>%
  group_by(design_id, SD, in_SZ) %>%
  summarize(mean_random_abundance_before = mean(random_abundance_before, na.rm = TRUE)) %>%
  spread(key = in_SZ, value = mean_random_abundance_before) %>%
  mutate(abundance_ratio = `TRUE` / `FALSE`) %>%
  dplyr::select(SD, design_id, abundance_ratio) %>%
  glimpse()

ggplot(ratio_before, aes(x = abundance_ratio, fill = SD)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 20) +
  scale_fill_manual(values = c("spabal" = "#fff2cc", "pref" = "#f1c232", "clump" = "#7f6000")) +
  labs(title = "Abundance Ratio Distributions BEFORE", x = "abundance ratio BEFORE", y = "Count") +
  theme_minimal()



percentage_before <- ratio_before %>%
  filter(abundance_ratio >= lower_bound & abundance_ratio <= upper_bound) %>%
  group_by(SD) %>%
  summarise(
    percentage_within_bounds = (n() / nrow(ratio_before %>% filter(SD == first(SD)))) * 100
  ) %>%
  glimpse()


# View the difference
percentage_before # no pref or clump design detects inside/outside abundance ratios over the 1.8 limit, but 3% of spabal design do.
percentage_after # ~15% of spabal and pref designs detect 1.8x abundance increase, and 21% of clump designs do.

# This result seems to suggest that the detection of SZ effects occurs in only
# 15-20% of sampling, regardless of the design.

# We've looked at the whole area, inside and outside SZ, but is there a
# difference of 1.8x before and after, only in SZ?

dat_in_SZ <- dat %>% # Re-organising the dataframe so it's easier to work with.
  filter(in_SZ = TRUE) %>%
  pivot_longer(
    cols = c(random_abundance_before, random_abundance_after),
    names_to = "BA",
    values_to = "random_abundance"
  ) %>%
  mutate(BA = recode(BA, "random_abundance_before" = "before", "random_abundance_after" = "after")) %>%
  dplyr::select(c(ID_1, design_id, SD, BA, random_abundance)) %>%
  glimpse()

ratio_SZ <- dat_in_SZ %>%
  group_by(design_id, SD, BA) %>%
  summarize(mean_random_abundance = mean(random_abundance, na.rm = TRUE)) %>%
  spread(key = BA, value = mean_random_abundance) %>%
  mutate(abundance_ratio = `before` / `after`) %>%
  dplyr::select(SD, design_id, abundance_ratio) %>%
  glimpse()

percentage <- ratio_SZ %>%
  filter(abundance_ratio >= 0.7 & abundance_ratio <= 0.9) %>% # Filter for the desired range
  group_by(SD) %>%
  summarise(percentage_within_bounds = n() / nrow(ratio_SZ %>% filter(SD == first(SD))) * 100) # Calculate the percentage for each SD

ggplot(ratio_SZ, aes(x = abundance_ratio, fill = SD)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 20) +
  scale_fill_manual(values = c("spabal" = "#fff2cc", "pref" = "#f1c232", "clump" = "#7f6000")) +
  labs(title = "Abundance Ratio Distributions AFTER", x = "abundance ratio AFTER", y = "Count") +
  theme_minimal()


percentage
saveRDS(percentage, paste0('data/rmd/SD_performance', spp, '.rds'))

# Okay so. Comparing each sampling design's ability to detect SZ effect,

# for pink snapper, clump design is best at 33.2%, then spabal at 20.6%, then pref at 18.2%
