# Spring rankings, by county
# Originally developed by A. Rosemartin and T. Crimmins
# Adapted from https://github.com/alyssarosemartin/spring-media/tree/main/earliest-spring-ranking
# 10 Feb 2025

library(rnpn)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(exactextractr)

# library(tigris) 
  # Only needed if us_counties shapefile not already present in resources folder

# This code extracts the average day of year of leaf out and bloom by county 
# (weighted average of pixels that fall within the county), enabling us
# to rank the top 10 earliest springs in this county and assess how the 
# current year compares.

# Helper calls ----------------------------------------------------------------#

# layers <- npn_get_layer_details()
# filter(layers, name == "si-x:average_leaf_prism")
# filter(layers, name == "si-x:average_bloom_prism")

# Using annual rasters with DOY that requirements for the first leaf or bloom
# Spring Index were met, averaged for 3 species, based on PRISM data

# Note: will obtain rasters for current year from: 
# https://data.usanpn.org/geoserver-request-builder/
# Selecting WCS layer: Spring Indices, Current Year - First Leaf- Spring Index
# Selecting WCS layer: Spring Indices, Current Year - first Bloom - Spring Index
#TODO: Check that this is correct!
# Will save current-year raster as: ncep_yyyymmdd_si-x_average_leaf/bloom.tif

# Set parameters --------------------------------------------------------------#

# Current year
year <- 2025

# First year to include in rankings
year1 <- 1981

# Prior years
prior_years <- year1:(year - 1)

# Location to save rasters
rast_folder <- "resources/rasters/"

# Logical to indicate whether script should be re-run to calculate current year
# index values as of today's date (if not done already)
today <- FALSE

# Minimum proportion of county that has reached spring index threshold 
min_prop <- 0.90
# If less than (min_prop * 100)% of county has reached spring index threshold, 
# we will not calculate a county mean.

# Get US county layer ---------------------------------------------------------#

# File location
counties_shp <- "resources/us_counties.shp"

# Download first and save to file, if not done already
if (!file.exists(counties_shp)) {
  library(tigris)
  us_counties <- tigris::counties(cb = TRUE)
  us_counties <- vect(us_counties)
  writeVector(us_counties, filename = counties_shp)
}

# Load file and limit to lower 48 states
counties <- vect(counties_shp)
states48 <- state.abb[! state.abb %in% c("AK", "HI")]
states48 <- c(states48, "DC")
counties <- terra::subset(counties, counties$STUSPS %in% c(states48))

# Get raster layers with annual leaf/bloom index values -----------------------#

index <- c("leaf", "bloom")

for (ind in index) {
  
  ind_folder <- paste0(rast_folder, ind)
  ind_files <- paste0(ind_folder, "/", prior_years, "_si-x_", ind, "_prism.tif")
  
  # If one or more index files for prior years doesn't exist, download them all
  if (!all(file.exists(ind_files))) {
    for (i in 1:length(prior_years)) {
      npn_download_geospatial("si-x:average_leaf_prism",
                              date = paste0(prior_years[i], "-01-01"),
                              output_path = ind_files[i])
    }
  }
  
  # If no current year files exist, will produce a warning. 
  # Need to download file from geoserver request builder first
  current_files <- list.files(ind_folder, pattern = "ncep_", full.names = TRUE)
  if (length(current_files) == 0) {
    warning(paste0("Current year ", ind, " index file has not been downloaded.",
                   " Download from geoserver request builder first."))
  }
  
  # If we need today's raster and that hasn't been downloaded yet, do that first
  today_flag <- 0
  today_file <- paste0(ind_folder, "/ncep_", 
                       str_remove_all(Sys.Date(), "-"), 
                       "_si-x_average_", ind, ".tif")
  if (today & !file.exists(today_file)) {
    today_flag <- 1
    warning(paste0("Current year ", ind, " index file has not been downloaded today.",
                   " Download from geoserver request builder first."))
  }
  
  # Resample current year raster to have the same geometry as rasters from prior
  # years if that hasn't been done already.
  current_file <- paste0(ind_folder, "/", year, "_si-x_", ind, "_prism.tif")
  if (!file.exists(current_file) | today_flag == 1) {
    current_files_orig <- list.files(ind_folder, 
                                     pattern = "ncep_", 
                                     full.names = TRUE)
    mostcurrent_file <- last(sort(current_files_orig))
    current_rast <- rast(mostcurrent_file)
    
    prior_rast <- rast(last(ind_files))
    
    # Resample current year raster to match geometry of prior years
    current_rast <- resample(current_rast, prior_rast, method = "bilinear")
    names(current_rast) <- paste0(year, "_si-x_", ind, "_prism")
    
    # Save new resampled raster
    writeRaster(current_rast, filename = current_file, overwrite = TRUE)
  }
}

# Calculate annual index values by county -------------------------------------#

# For each year, we will calculate the mean DOY for each county by calculating 
# the mean value of raster cells in each county, weighted by the proportion of 
# each raster cell that is within the county boundary. 

# Did some checks with terra and exactextractr packages
# terra::extract() with na.rm = TRUE and exact = TRUE produces same results as
# exactextractr::exact_extract() with fun = "mean".
# Functions in exactextractr require sf vector objects, which is a pain, but 
# run much faster than terra functions, so we'll stick with exactextractr

# After an initial run, realized that some county values were extreme (eg, more
# than 30 days early) and high resolution spring maps on NPN website didn't
# have values anywhere near that. This seems to be occurring when we have large
# counties where only a few cells have reached the threshold and so most raster
# cells are NA. Given this, we will only calculate county means when most or all
# of the raster cells have reached the threshold. 

for (ind in index) {
  # Create multi-layer raster
  ann_files <- list.files(paste0(rast_folder, ind), full.names = TRUE)
  ann_files <- str_subset(ann_files, pattern = "ncep_", negate = TRUE)
  mrast <- rast(ann_files)
  
  names(mrast) <- paste0("y", str_sub(names(mrast), 1, 4))
  ind_doy <- exactextractr::exact_extract(x = mrast,
                                          y = sf::st_as_sf(counties),
                                          fun = c("mean", "count"),
                                          append_cols = c("STATEFP", "COUNTYFP",
                                                          "STUSPS", "NAME"))
  # The count function provides the sum of all cell coverage fractions, ignoring
  # NA values. Counts should be the same for all prior years. We will use the 
  # ratio of current/prior coverages to determine whether enough of the county
  # has reached the spring threshold value that it's worth calculating the mean.
  
  cmean <- paste0("mean.y", year)
  ccount <- paste0("count.y", year)
  ind_doy$frac <- ind_doy[ ,ccount] / ind_doy[ ,"count.y2024"]
  # Set county mean to NA if < min_prop of county has reached threshold
  ind_doy[, cmean] <- ifelse(ind_doy$frac < min_prop, NA, ind_doy[, cmean])
  ind_doy <- ind_doy %>%
    mutate(across(everything(), ~ifelse(is.nan(.), NA, .))) %>%
    select(STATEFP, COUNTYFP, STUSPS, NAME, starts_with("mean.y")) %>%
    rename(state = STUSPS, 
           county = NAME) %>%
    rename_with(~ str_replace(.x, pattern = "mean.y", replacement = "doy_"),
                matches("mean.y"))

  # Save to file
  doy_file <- paste0("output/county-", ind, "-index-", year1, "-", year, ".csv")
  write.csv(ind_doy,
            file = doy_file, 
            row.names = FALSE)
  rm(ind_doy)
}

# Summarize/visualize spring index values -------------------------------------#

# Now we have dataframes with the mean DOY leaf/bloom spring indices were
# reached in each county and year (value = NA if 90% of county has not reached
# the first leaf/bloom).

# For each county: we can list the top 10 earliest/latest years, and note which
# counties have current year in top 10 lists.

# For counties where spring index has been reached in current year: we can 
# calculate anomalies (ie, the number of days earlier/later the current year is 
# than historical average). This is similar to the spring index maps on the NPN 
# website, but aggregated by county.

# Visualize leaf spring index values ------------------------------------------#

# Read csv with annual leaf index values
leaf_doy_file <- paste0("output/county-leaf-index-", year1, "-", year, ".csv") 
leaf <- read.csv(leaf_doy_file, colClasses = c(STATEFP = "character",
                                               COUNTYFP = "character"))

# Calculate anomalies
prior_yr_cols <- paste0("doy_", prior_years)
leaf$mean <- apply(leaf[, prior_yr_cols], 1, mean)
leaf$median <- apply(leaf[, prior_yr_cols], 1, median)
leaf <- leaf %>%
  mutate(across(.cols = starts_with("doy_"),
                .fns = ~.x - mean,
                .names = "{sub('doy', 'anom', col)}"))

# Add current year anomalies to county polygons
current_anoms <- paste0("anom_", year)
counties2 <- merge(counties, 
                   select(leaf, STATEFP, COUNTYFP, contains(current_anoms)),
                   by.x = c("STATEFP", "COUNTYFP"))

# Limits for nice colors in map
limit <- max(abs(data.frame(counties2[,current_anoms])), na.rm = TRUE)
limit_round <- ifelse(limit < 10, 5, floor(limit / 10) * 10)
limits <- limit * c(-1, 1)

# Create map
leaf_a_map <- ggplot(counties2) +
  geom_spatvector(aes(fill = get(current_anoms)), color = NA) +
  scale_fill_distiller(palette = "Spectral", limit = limits, direction = 1,
                       breaks = c(-limit_round, 0, limit_round),
                       labels = c(paste0(limit_round, " days early"), "Average",
                                  paste0(limit_round, " days late"))) +
  labs(title = paste0("Leaf spring index, ", year, " anomalies")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.key.width = unit(3, "cm"), 
        legend.ticks.length = unit(c(-.15, 0), "cm"), 
        legend.ticks = element_line(color = "black"))
leaf_a_map

# Create dataframe with spring rankings for prior years (not including current)
# I think there might be ties we'll need to deal with...

prior_cols <- paste0("doy_", prior_years)
leafr <- leaf %>%
  select(-contains("anom")) %>%
  select(-contains(paste0("doy_", year))) %>%
  pivot_longer(cols = all_of(prior_cols),
               names_to = "year",
               values_to = "doy") %>%
  mutate(year = as.numeric(str_remove(year, "doy_"))) %>%
  group_by(STATEFP, COUNTYFP, state, county, mean, median) %>%
  mutate(rank = rank(doy)) %>%
  data.frame()
leaf_rank <- leafr %>%
  group_by(STATEFP, COUNTYFP, state, county, mean, median) %>%
  summarize(early1 = year[rank == 1],
            early2 = year[rank == 2], 
            .groups = "keep") %>%
  data.frame()
  
