# Anomalous observations script
# Originally developed by A. Rosemartin and T. Crimmins
# Adapted from https://github.com/alyssarosemartin/spring-media/anomalous_obs
# by E. Zylstra
# 31 Jan 2025

library(rnpn)
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(elevatr)
library(tigris)
options(tigris_class = "sf")

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
prior_years <- 2009:(year - 1)

# Define early observations (below X percentile)
lowerq <- 0.05

# Colors...

# Get states layer from tigris package ----------------------------------------#

states <- states()
states <- vect(states)
states <- terra::project(states, "epsg:4326")
states <- states[, "STUSPS"]

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

# Note: for this script, not going to differentiate among phenophases in 
# each phenophase class. Mostly phenophases differ among species or functional
# type, but there are also some old phenophases that are no longer used.

# Aquire, format phenometric data in current year -----------------------------#

current <- npn_download_individual_phenometrics(
  request_source = 'erinz', 
  years = year,
  pheno_class_ids = c(1, 3, 6, 7), 
  additional_fields = c("observed_status_conflict_flag",
                        "species_functional_type")
)

# Simplify data
current <- current %>%
  select(site_id, latitude, longitude, elevation_in_meters, state, 
         common_name, species_id, species_functional_type, individual_id, 
         phenophase_description, pheno_class_id, first_yes_year, first_yes_day, 
         numdays_since_prior_no, observed_status_conflict_flag) %>%
  rename(lat = latitude,
         lon = longitude, 
         elev = elevation_in_meters, 
         func_type = species_functional_type,
         id = individual_id, 
         phenophase = phenophase_description,
         year = first_yes_year,
         first_yes = first_yes_day,
         prior_no = numdays_since_prior_no,
         flag = observed_status_conflict_flag) %>%
  mutate(across(c(elev, prior_no), 
         ~ ifelse(. == -9999, NA, .))) %>%
  mutate(across(c(state, flag), 
                ~ ifelse(. == "-9999", NA, .))) %>%
  mutate(state = ifelse(state == "", NA, state)) %>%
  data.frame()
  
# Remove observations of 'ohi'a lehua 
current <- current %>%
  filter(common_name != "'ohi'a lehua")

# Remove observations with conflicts or those without a prior no
current <- current %>%
  filter(is.na(flag)) %>%
  filter(!is.na(prior_no))

# Keep only first observation of individual for a phenophase in a year
current <- current %>%
  arrange(id, pheno_class_id, year, first_yes) %>%
  distinct(id, pheno_class_id, year, .keep_all = TRUE)

# Some observations missing state ID. Will fill in using the tigris package
states48 <- state.abb[! state.abb %in% c("AK", "HI")]
state_fill <- filter(current, is.na(state)) %>%
  select(site_id, lon, lat) %>%
  distinct() 
state_fillv <- vect(state_fill, crs = "epsg:4326")
state_new <- terra::extract(states, state_fillv)
state_fill <- cbind(state_fill, state_new = state_new$STUSPS)
current <- current %>%
  left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
  mutate(state = ifelse(!is.na(state), state, state_new)) %>%
  select(-state_new) %>%
  # Restrict observations to the continental US
  filter(state %in% states48) 

# Some observations missing elevation. Will fill in using the elevatr package
elev_fill <- filter(current, is.na(elev)) %>%
  distinct(site_id, lat, lon) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
elev_fill <- get_elev_point(locations = elev_fill, src = "epqs")
elev_fill <- data.frame(elev_fill) %>%
  mutate(elev_new = round(elevation)) %>%
  select(site_id, elev_new)
current <- current %>%
  left_join(elev_fill, by = "site_id") %>%
  mutate(elev = ifelse(!is.na(elev), elev, elev_new)) %>%
  select(-elev_new)

# Check which phenophases, functional groups have sufficient data
count(current, pheno_class_id)
count(current, pheno_class_id, func_type)
  # Cactus, deciduous broadleaf, deciduous conifer, drought deciduous broadleaf,
  # evergreen broadleaf, evergreen conifer, forb, graminoid, 
  # semi-evergreen broadleaf, semi-evergreen forb

# For now, pick a phenophase class to focus on #################################
pheno_class <- 7 

# Extract first observations of the year and species list
current_ph <- current %>%
  filter(pheno_class_id == pheno_class) %>%
  arrange(id, pheno_class_id, first_yes) %>%
  distinct(id, pheno_class_id, .keep_all = TRUE)
current_spp <- unique(current_ph$species_id)

# Aquire, format phenometric data in previous years ---------------------------#

prior <- npn_download_individual_phenometrics(
  request_source = 'erinz', 
  years = prior_years,
  species_ids = current_spp,
  pheno_class_ids = pheno_class, 
  additional_fields = c("observed_status_conflict_flag",
                        "species_functional_type")
)

# Simplify data
prior <- prior %>%
  select(site_id, latitude, longitude, elevation_in_meters, state, 
         common_name, species_id, species_functional_type, individual_id, 
         phenophase_description, pheno_class_id, first_yes_year, first_yes_doy, 
         numdays_since_prior_no, observed_status_conflict_flag) %>%
  rename(lat = latitude,
         lon = longitude, 
         elev = elevation_in_meters, 
         func_type = species_functional_type,
         id = individual_id, 
         phenophase = phenophase_description,
         year = first_yes_year,
         first_yes = first_yes_doy,
         prior_no = numdays_since_prior_no,
         flag = observed_status_conflict_flag) %>%
  mutate(across(c(elev, prior_no), 
                ~ ifelse(. == -9999, NA, .))) %>%
  mutate(across(c(state, flag), 
                ~ ifelse(. == "-9999", NA, .))) %>%
  mutate(state = ifelse(state == "", NA, state)) %>%
  data.frame()

# Remove observations with conflicts or those without a prior no
prior <- prior %>%
  filter(is.na(flag)) %>%
  filter(!is.na(prior_no))

# Keep only first observation of individual for a phenophase in a year
prior <- prior %>%
  arrange(id, pheno_class_id, year, first_yes) %>%
  distinct(id, pheno_class_id, year, .keep_all = TRUE)

# Fill in state or elevation where needed
state_fill <- filter(prior, is.na(state)) %>%
  select(site_id, lon, lat) %>%
  distinct() 
state_fillv <- vect(state_fill, crs = "epsg:4326")
state_new <- terra::extract(states, state_fillv)
state_fill <- cbind(state_fill, state_new = state_new$STUSPS)
prior <- prior %>%
  left_join(select(state_fill, site_id, state_new), by = "site_id") %>%
  mutate(state = ifelse(!is.na(state), state, state_new)) %>%
  select(-state_new) %>%
  # Restrict observations to the continental US
  filter(state %in% states48) 

# Some observations missing elevation. Will fill in using the elevatr package
elev_fill <- filter(prior, is.na(elev)) %>%
  distinct(site_id, lat, lon) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
elev_fill <- get_elev_point(locations = elev_fill, src = "epqs")
elev_fill <- data.frame(elev_fill) %>%
  mutate(elev_new = round(elevation)) %>%
  select(site_id, elev_new)
prior <- prior %>%
  left_join(elev_fill, by = "site_id") %>%
  mutate(elev = ifelse(!is.na(elev), elev, elev_new)) %>%
  select(-elev_new)

# Compare current year observations with distribution of first observation dates
# in prior years --------------------------------------------------------------#

# Want a minimum number of years/observations in prior years to create
# a reasonable distribution.

# What is our spatial unit for comparison?
  # The site where the 2025 plant is located?
  # XX-mi radius around 2025 site (also within xx m elevation?)
  # State the 2025 site is located in?

# Start with an analysis by state, since that's how it was done previously ----#
current_ph <- current_ph %>%
  mutate(state_spp = paste0(state, "_", common_name))
prior <- prior %>%
  mutate(state_spp = paste0(state, "_", common_name))
prior_st <- prior %>%
  filter(state_spp %in% unique(current_ph$state_spp))
state_spp_n <- prior_st %>%
  group_by(state_spp) %>%
  summarize(n_obs = n(),
            n_indiv = n_distinct(id),
            n_years = n_distinct(year),
            .groups = "keep") %>%
  data.frame()
# Only evaluate 2025 anomalies for state-species combinations that have 8+
# years of data and 20+ observations
min_yrs <- 8
min_obs <- 20
state_spp_n <- state_spp_n %>%
  filter(n_years >= 8 & n_obs >= 20)

prior_st <- prior_st %>%
  filter(state_spp %in% state_spp_n$state_spp)
quantiles_st <- prior_st %>%
  group_by(state_spp) %>%
  summarize(q0.05 = quantile(first_yes, 0.05),
            q0.25 = quantile(first_yes, 0.25),
            q0.75 = quantile(first_yes, 0.75),
            q0.95 = quantile(first_yes, 0.95),
            IQR = IQR(first_yes)) %>%
  data.frame()

# PICK UP HERE






# Evaluate how much data we have at each site:
  site_spp <- prior %>%
    group_by(site_id, lat, lon, elev, common_name) %>%
    summarize(n_obs = n(),
              n_indiv = n_distinct(id),
              n_years = n_distinct(year),
              .groups = "keep") %>%
    mutate(spp_site = paste0(common_name, "_", site_id)) %>%
    data.frame()
  
  yrs10 <- site_spp$spp_site[site_spp$n_years >= 10]
  length(yrs10)
  length(yrs10) / nrow(site_spp)
  # Few site/species have >= 10 years
  
  obs10yrs5 <- site_spp$spp_site[site_spp$n_years >= 5 & site_spp$n_obs >= 10]
  length(obs10yrs5)
  length(obs10yrs5) / nrow(site_spp)
  # Few site/species have >= 10 obs over 5+ years
  
  sum(paste0(current_ph$common_name, "_", current_ph$site) %in% yrs10)
  # 6 (of 28) site-spp that have 10+ years of prior data
  sum(paste0(current_ph$common_name, "_", current_ph$site) %in% obs10yrs5)
  # 12 (of 28) site-spp that have 10+ obs over 5+ years of prior data

  
# For each 2025 location, see how much data we have for a 50 km radius 
# and within 

# First one: 
  current1 <- current_ph[1,]
  
locs <- prior1 %>%
  filter(common_name == current1$common_name) %>%
  select(site_id, lat, lon) %>%
  distinct()

dist <- terra::distance(x = as.matrix(current1[,c("lon", "lat")]),
                        y = as.matrix(locs[,c("lon", "lat")]), 
                        lonlat = TRUE,
                        unit = "km")
locs50 <- locs[which(dist <= 50), ] %>%
  mutate(within50 = 1)

# Want to also filter by elevation, but will need to fill in missing values
# at some sites

prior50 <- prior1 %>%
  left_join(select(locs50, site_id, within50), by = "site_id") %>%
  filter(within50 == 1) %>%
  filter(common_name == current1$common_name)

boxplot(prior50$first_yes)








  
  
  
  
  





