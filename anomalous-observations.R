# Anomalous observations script
# Originally developed by A. Rosemartin and T. Crimmins
# Adapted from https://github.com/alyssarosemartin/spring-media/anomalous_obs
# by E. Zylstra
# 30 Jan 2025

library(rnpn)
library(dplyr)
library(ggplot2)

# This code identifies anomalous observations of plants in the current year 
# compared to the long term record (2009-XXXX). First observation date in 
# current year is compared to the distribution of dates from prior years. Early 
# (lower than Xth percentile) and outliers (earlier than 1.5x the IQR) are 
# flagged and plotted

# Need individual phenometrics for current year
# Need individual phenometrics for prior years (for current year species)

# Set parameters --------------------------------------------------------------#

# Current year
year <- 2025

# Prior years
prior <- 2009:2024

# Define early observations (below X percentile)
lowerq <- 0.05

# Colors...

# Helper calls ----------------------------------------------------------------#

# Get species IDs
species <- npn_species()

# Get phenophase classes
phenoclasses <- npn_pheno_classes() %>% data.frame()
head(phenoclasses, 13)
  # 1: initial shoot/leaf growth
  # 3: leaves/needles
  # 6: flowers/cones
  # 7: open flowers/cones 

# Aquire, format phenometric data in current year -----------------------------#

current <- npn_download_individual_phenometrics(
  request_source = 'erinz', 
  years = year,
  pheno_class_ids = c(1, 3, 6, 7), 
  additional_fields = c("observed_status_conflict_flag",
                        "species_functional_type")
)

# Convert -9999's to NA
current <- current %>%
  mutate(across(c(elevation_in_meters, numdays_since_prior_no, 
                  numdays_until_next_no), 
         ~ ifelse(. == -9999, NA, .))) %>%
  mutate(across(c(state, observed_status_conflict_flag), 
                ~ ifelse(. == "-9999", NA, .)))
  


