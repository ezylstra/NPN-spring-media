# Anomalous observations script
# Originally developed by A. Rosemartin and T. Crimmins
# Adapted from https://github.com/alyssarosemartin/spring-media/anomalous_obs
# 3 Feb 2025

library(rnpn)
library(dplyr)
library(stringr)
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
##TODO: Save layer in project?

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
pheno_class_ids <- c(1, 3, 6, 7)

# Note: for this script, not going to differentiate among phenophases in 
# each phenophase class. Mostly phenophases differ among species or functional
# type, but there are also some old phenophases that are no longer used.

# Aquire, format phenometric data in current year -----------------------------#

current_dl <- npn_download_individual_phenometrics(
  request_source = 'erinz', 
  years = year,
  pheno_class_ids = pheno_class_ids, 
  additional_fields = c("observed_status_conflict_flag",
                        "species_functional_type")
)

# Simplify data (remove unnecessary columns; replace -9999 with NA)
current <- current_dl %>%
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

# Some observations missing state ID. Will fill in using states layer from the 
# tigris package
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
##TODO: move later in script since only needed if using elevation to filter
## observations in prior years?
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

# Extract first observations of the year and species list
current <- current %>%
  arrange(id, pheno_class_id, first_yes) %>%
  distinct(id, pheno_class_id, .keep_all = TRUE)

# Aquire, format phenometric data in previous years ---------------------------#

for (pheno_class in pheno_class_ids) {
  spp <- unique(current$species_id[current$pheno_class_id == pheno_class])
  dl_temp <- npn_download_individual_phenometrics(
               request_source = 'erinz', 
               years = prior_years,
               species_ids = spp,
               pheno_class_ids = pheno_class, 
               additional_fields = c("observed_status_conflict_flag",
                                     "species_functional_type")
              )
  if (pheno_class == pheno_class_ids[1]) {
    prior_dl <- dl_temp
  } else {
    prior_dl <- rbind(prior_dl, dl_temp)
  }
  rm(dl_temp)
}

# Simplify data (remove unnecessary columns; replace -9999 with NA)
prior <- prior_dl %>%
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

# Remove observations with conflicts or those without a prior no in last 30 days
prior <- prior %>%
  filter(is.na(flag)) %>%
  filter(!is.na(prior_no)) %>%
  filter(prior_no <= 30) %>%
  select(-flag)

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
# for plants in the same state ------------------------------------------------#

current <- current %>%
  mutate(state_spp = paste0(state, "_", common_name),
         state_spp_ph = paste0(state_spp, "_", pheno_class_id))
prior <- prior %>%
  mutate(state_spp = paste0(state, "_", common_name),
         state_spp_ph = paste0(state_spp, "_", pheno_class_id))

# Identify early or outlier observations based on previous observations in 
# same state
prior_st <- prior %>%
  filter(state_spp_ph %in% unique(current$state_spp_ph))
state_spp_n <- prior_st %>%
  group_by(state_spp, pheno_class_id, state_spp_ph) %>%
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
  filter(n_years >= min_yrs & n_obs >= min_obs)

prior_st <- prior_st %>%
  filter(state_spp_ph %in% state_spp_n$state_spp_ph)
quantiles_st <- prior_st %>%
  group_by(state_spp_ph) %>%
  summarize(min = min(first_yes),
            q0.05 = quantile(first_yes, 0.05),
            q0.25 = quantile(first_yes, 0.25),
            q0.75 = quantile(first_yes, 0.75),
            q0.95 = quantile(first_yes, 0.95),
            IQR = IQR(first_yes)) %>%
  data.frame() %>%
  mutate(whisker = q0.25 - 1.5 * IQR,
         outlier_threshold = ifelse(whisker < min, min, whisker))
current_st <- current %>%
  left_join(quantiles_st, by = "state_spp_ph") %>%
  mutate(early = ifelse(first_yes < q0.05, 1, 0),
         outlier = ifelse(first_yes < outlier_threshold, 1, 0)) %>%
  # Remove state_spp_ph combinations that don't have enough data in prior years
  filter(!is.na(early))

# Look at what's left
current_st %>%
  select(pheno_class_id, state_spp, func_type, first_yes, min, 
         outlier_threshold, q0.25, q0.75, early, outlier) %>%
  filter(early == 1 | outlier == 1) %>%
  arrange(pheno_class_id, func_type, state_spp, first_yes)

prior %>%
  filter(state_spp == "CA_coyotebrush" & pheno_class_id == 7) %>%
  ggplot(aes(first_yes)) +
    geom_histogram(bins = 60)
# Coyotebrush histogram highlights why boxplots might be deceiving. It's a 
# multimodal distribution, with the highest frequency in fall, but a distinct
# peak at the beginning of the year. 

# Extract just early/outlier observations
current_steo <- current_st %>%
  filter(early == 1 | outlier == 1) %>%
  select(common_name, state, state_spp, pheno_class_id, state_spp_ph, func_type, 
         first_yes, min, q0.05, q0.25, q0.75, early, outlier) %>%
  mutate(panel = paste0(common_name, " (", state, ")"),
         eo = ifelse(outlier == 1, "Outlier", "Early"),
         eo = factor(eo, levels = c("Early", "Outlier")))

# Loop through phenophase classes
for (phc in pheno_class_ids) {
  current_steo_p <- filter(current_steo,  pheno_class_id == phc)
  current_eo_stspp <- sort(unique(current_steo_p$state_spp))
  current_eo_spp <- str_split_i(current_eo_stspp, "_", 2)
  current_eo_st <- str_split_i(current_eo_stspp, "_", 1)
  
  # Put prior data into better format for plotting
  for (i in 1:length(current_eo_stspp)) {
    prior_spp <- prior %>%
      filter(common_name == current_eo_spp[i] & pheno_class_id == phc) %>%
      select(common_name, state_spp, state, func_type, year, first_yes) %>%
      mutate(group = "All") %>%
      mutate(panel = paste0(current_eo_spp[i], " (", current_eo_st[i], ")"))
    prior_stspp <- prior %>%
      filter(state_spp == current_eo_stspp[i] & pheno_class_id == phc) %>%
      select(common_name, state_spp, state, func_type, year, first_yes) %>%
      mutate(group = "State") %>%
      mutate(panel = paste0(current_eo_spp[i], " (", current_eo_st[i], ")"))
    prior_plot_new <- rbind(prior_spp, prior_stspp)
    if (i == 1) {
      prior_plot <- prior_plot_new
    } else {
      prior_plot <- rbind(prior_plot, prior_plot_new)
    }
  }
  
  phc_name <- paste0(phenoclasses$name[phc], ", first observation of the year")
  
  plot_temp <- prior_plot %>%
    ggplot(aes(x = first_yes, fill = group)) +
    geom_histogram(bins = 60, position = "identity") +
    scale_fill_manual(values = c("gray", "steelblue3")) +
    facet_wrap(~panel, scales = "free_y", ncol = 2) +
    geom_vline(data = current_steo_p,
               aes(xintercept = first_yes, linetype = eo), 
               color = "red") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Day of year", y = "Count", title = phc_name,
         fill = "Prior observations", linetype = paste0(year, " observations")) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1, "cm"))
  assign(paste0("plot_", phc), plot_temp)
}

for (phc in pheno_class_ids) {
  print(get(paste0("plot_", phc)))
}

# Compare current year observations with distribution of first observation dates
# for plants within XX-km radius (and within xx m elevation?) -----------------#

radius <- 50

current_orig <- current

rm(quantiles_r)
rm(quantiles_temp)
rm(prior_temp)
rm(prior_plot_r)
current <- current_orig
current$n_obs_r <- NA
current$n_indiv_r <- NA
current$n_yrs_r <- NA

for (i in 1:nrow(current)) {
  current1 <- current[i, ]
  # Get locations for observations of that species, phenophase
  locs <- prior %>%
    filter(common_name == current1$common_name, 
           pheno_class_id == current1$pheno_class_id) %>%
    select(site_id, lat, lon) %>%
    distinct()
  # Calculate distance between 2025 plant and all others
  dists <- terra::distance(x = as.matrix(current1[,c("lon", "lat")]),
                           y = as.matrix(locs[,c("lon", "lat")]), 
                           lonlat = TRUE,
                           unit = "km")
  locs_r <- locs[which(dists <= radius), ] %>%
    mutate(withinr = 1)
  
  prior_r <- prior %>%
    left_join(select(locs_r, site_id, withinr), by = "site_id") %>%
    filter(withinr == 1) %>%
    filter(common_name == current1$common_name) %>%
    filter(pheno_class_id == current1$pheno_class_id)
  
  current$n_obs_r[i] <- nrow(prior_r)
  current$n_indiv_r[i] <- n_distinct(prior_r$id)
  current$n_yrs_r[i] <- n_distinct(prior_r$year)
  
  if (current$n_obs_r[i] <= min_obs | current$n_yrs_r[i] < min_yrs) {next}
  
  quantiles_temp <- data.frame(indiv = paste0(current1$common_name, "_", current1$id), 
                               common_name = current1$common_name,
                               site_id = current1$site_id,
                               lon = current1$lon,
                               lat = current1$lat, 
                               state = current1$state,
                               id = current1$id, 
                               pheno_class_id = current1$pheno_class_id,
                               first_yes = current1$first_yes,
                               min = min(prior_r$first_yes),
                               q0.05 = quantile(prior_r$first_yes, 0.05),
                               q0.25 = quantile(prior_r$first_yes, 0.25),
                               q0.75 = quantile(prior_r$first_yes, 0.75),
                               q0.95 = quantile(prior_r$first_yes, 0.95),
                               IQR = IQR(prior_r$first_yes)) %>%
    mutate(whisker = q0.25 - 1.5 * IQR,
           outlier_threshold = ifelse(whisker < min, min, whisker),
           early = ifelse(first_yes < q0.05, 1, 0),
           outlier = ifelse(first_yes < outlier_threshold, 1, 0))
  
  if (exists("quantiles_r")) {
    quantiles_r <- rbind(quantiles_r, quantiles_temp)
  } else {
    quantiles_r <- quantiles_temp
  }
  
  if (quantiles_temp$early == 1 | quantiles_temp$outlier == 1) {
    # Create prior dataframe for ggplot
    prior_temp <- prior_r %>%
      select(common_name, func_type, pheno_class_id, site_id, lat, lon, state,
             year, first_yes) %>%
      mutate(indiv = quantiles_temp$indiv) %>%
      mutate(panel = paste0(indiv, " (", quantiles_temp$state, ")"))
      
    if (exists("prior_plot_r")) {
      prior_plot_r <- rbind(prior_plot_r, prior_temp)
    } else {
      prior_plot_r <- prior_temp
    }
  }
}  

# Extract just early/outlier observations
quantiles_reo <- quantiles_r %>%
  filter(early == 1 | outlier == 1) %>%
  select(indiv, common_name, state, site_id, pheno_class_id, 
         first_yes, min, q0.05, q0.25, q0.75, early, outlier) %>%
  mutate(panel =  paste0(indiv, " (", state, ")"),
         eo = ifelse(outlier == 1, "Outlier", "Early"),
         eo = factor(eo, levels = c("Early", "Outlier")))

# Loop through phenophase classes
for (phc in pheno_class_ids) {
  
  quantiles_reo_p <- filter(quantiles_reo,  pheno_class_id == phc)
  prior_plot_r_p <- filter(prior_plot_r, pheno_class_id == phc)

  phc_name <- paste0(phenoclasses$name[phc], ", first observation of the year")
  
  plot_temp <- prior_plot_r_p %>%
    ggplot(aes(x = first_yes)) +
    geom_histogram(bins = 60, fill = "steelblue3") +
    facet_wrap(~panel, scales = "free_y", ncol = 2) +
    geom_vline(data = quantiles_reo_p,
               aes(xintercept = first_yes, linetype = eo), 
               color = "red") +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(x = "Day of year", y = "Count", title = phc_name,
         linetype = paste0(year, " observations")) +
    theme_bw() +
    theme(legend.position = "bottom")
  assign(paste0("plot_r_", phc), plot_temp)
}

for (phc in pheno_class_ids) {
  print(get(paste0("plot_r_", phc)))
}



# Not sure if we'll use stuff below -------------------------------------------#

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
  