# Spring indices, by county
# Originally developed by A. Rosemartin and T. Crimmins
# Adapted from https://github.com/alyssarosemartin/spring-media/tree/main/earliest-spring-ranking
# 13 Feb 2025

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

# This code calculates the mean day of year that the leaf/bloom spring index is
# reached in each US county

# Set parameters --------------------------------------------------------------#

# Current year
year <- 2025

# First year to download layers for
year1 <- 1981

# Prior years (for data download)
prior_years <- year1:(year - 1)

# Location to save rasters
rast_folder <- "resources/rasters/"

# Logical to indicate whether script should be re-run to calculate current year
# index values as of today's date (if not done already)
rerun_today <- TRUE

# Minimum proportion of county that has reached spring index threshold 
min_prop <- 0.90
# If less than (min_prop * 100)% of county has reached spring index threshold, 
# we will not calculate a county mean.

# Identify whether to restrict output to current year + previous 40 years
# (this would be similar to approach NPN has used previously, particularly for
# return intervals)
prior40 <- TRUE

# Identify whether to create a formatted table and save a csv file with leaf 
# index and/or bloom index values 
leaf_table <- TRUE
bloom_table <- TRUE

# For csv file output: Select whether to have one row per county 
# (columns = years; county_row = TRUE) or one column per county (rows = years; 
# county_row = FALSE)
county_row <- TRUE

# Helper calls ----------------------------------------------------------------#

layers <- npn_get_layer_details()

# filter(layers, name == "si-x:average_leaf_prism")
# filter(layers, name == "si-x:average_bloom_prism")
  # Annual rasters with DOY that requirements for the first leaf/bloom index 
  # were met, averaged for 3 species (PRISM data)

# filter(layers, name == "si-x:average_leaf_ncep")
# filter(layers, name == "si-x:average_bloom_ncep")
  # Raster with DOY leaf/bloom index reached in current year (NCEP data). Can
  # download any day up to current or up to 6 days forecasted in advance.

# filter(layers, name == "si-x:leaf_anomaly")
# filter(layers, name == "si-x:bloom_anomaly")
  # Raster with anomalies for current year (DOY - 30-year mean DOY; NCEP data). 
  # Could use this as a gut-check for anomalies we calculate at the county level

# filter(layers, name == "si-x:leaf_return_interval")
# filter(layers, name == "si-x:bloom_return_interval")
  # Raster with return intervals for current year based on data over previous 
  # 40 years (NCEP data). Could use this as a gut-check for return intervals we 
  # calculate at the county level

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
  
  # Download rasters for prior years 
  # If one or more index files for prior years doesn't exist, download them all
  if (!all(file.exists(ind_files))) {
    for (i in 1:length(prior_years)) {
      npn_download_geospatial("si-x:average_leaf_prism",
                              date = paste0(prior_years[i], "-01-01"),
                              output_path = ind_files[i])
    }
  }
  
  # Download rasters for current year. Only do so if any current-year file 
  # doesn't exist or if rerun_today == TRUE and layer for yesterday (most
  # current non-forecasted layer) hasn't already been downloaded
  ncep_doy_files <- list.files(ind_folder, pattern = "ncep_", full.names = TRUE)
  ncep_doy_yest <- paste0(ind_folder, "/ncep_", str_remove_all(today() - 1, "-"), 
                          "_si-x_average_", ind, ".tif")
  need_yest <- !file.exists(ncep_doy_yest)
  if (length(ncep_doy_files) == 0 | (rerun_today & need_yest)) {
    npn_download_geospatial("si-x:average_leaf_ncep",
                            date = today() - 1, 
                            output_path = ncep_doy_yest)
    
    # Resample current year raster to have the same geometry as rasters from 
    # prior years (based on PRISM data).
    ncep_doy_rast <- rast(ncep_doy_yest) 
    prioryr_rast <- rast(last(ind_files))
    current_rast <- resample(ncep_doy_rast, prioryr_rast, method = "bilinear")
    names(current_rast) <- paste0(year, "_si-x_", ind, ".tif")
    
    # Save new resampled raster
    writeRaster(current_rast, 
                filename = paste0(ind_folder, "/", year, "_si-x_", ind, ".tif"),
                overwrite = TRUE)
    
    rm(ncep_doy_rast, prioryr_rast, current_rast)
  }
}
  
# Calculate annual index values by county -------------------------------------#

# For each year, we will calculate the mean DOY for each county by calculating 
# the mean value of raster cells in each county, weighted by the proportion of 
# each raster cell that is within the county boundary. 

# Did some checks with terra and exactextractr packages
# terra::extract() with na.rm = TRUE and exact = TRUE produces same results as
# exactextractr::exact_extract() with fun = "mean".
# Functions in exactextractr require sf vector objects, which is a pain but 
# run much faster than terra functions, so we'll stick with exactextractr

# After an initial run, realized that some county values were extreme (eg, more
# than 30 days early) and high resolution spring maps on NPN website didn't
# have values anywhere near that. This seems to be occurring when we have large
# counties where only a few cells have reached the threshold, Given this, we 
# will only calculate county means when most or all of the raster cells have 
# reached the threshold. 

for (ind in index) {
  
  # Create multi-layer raster
  ann_files <- list.files(paste0(rast_folder, ind), full.names = TRUE)
  ann_files <- str_subset(ann_files, pattern = "ncep_", negate = TRUE)
  mrast <- rast(ann_files)
  names(mrast) <- paste0("y", str_sub(names(mrast), 1, 4))
  
  # Create dataframe with annual means for each county
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

# For each county we can provide tables with annual values, including those 
# from the current year.

# For counties where spring index has been reached in current year: we can 
# calculate anomalies (ie, the number of days earlier/later the current year is 
# than historical average). We can also calculate the return interval,
# or the frequency over the last 40 years with which spring has been as early or 
# late as that in the current year. We can "gut check" our calculations of 
# anomalies and return intervals at the county level with high resolution maps
# provided on NPN website.

# Format and save tables with raw data ----------------------------------------#

# Identify whether to restrict output to current year + previous 40 years
# (this would be similar to approach NPN has used previously, particularly for
# return intervals)
prior40 <- TRUE

# Select whether to have one row per county (columns = years; county_row = TRUE)
# or one column per county (rows = years; county_row = FALSE)
county_row <- TRUE

# Identify whether to create table for leaf index and/or bloom index
leaf_table <- TRUE
bloom_table <- TRUE

if (leaf_table) {
  
  # Read csv with annual leaf index values
  leaf_doy_file <- paste0("output/county-leaf-index-", year1, "-", year, ".csv") 
  leaf <- read.csv(leaf_doy_file, colClasses = c(STATEFP = "character",
                                                 COUNTYFP = "character"))
  
  # Create columns to uniquely identify each county
  leaf <- leaf %>%
    mutate(FIPS = paste0(STATEFP, COUNTYFP),
           name_FIPS = paste0(state, "_", county, "_", FIPS)) %>%
    relocate(c(FIPS, name_FIPS)) %>%
    arrange(FIPS)

  if (prior40) {
    earliest_yr <- year - 40
    cols_to_remove <- paste0("doy_", year1:(earliest_yr - 1))
    leaf <- leaf %>%
      select(-all_of(cols_to_remove))
  }
  
  if (!county_row) {
    doys <- t(as.matrix(select(leaf, contains("doy")))) %>%
      data.frame()
    colnames(doys) <- leaf$name_FIPS
    doys <- doys %>%
      mutate(year = as.numeric(str_remove(rownames(doys), "doy_"))) %>%
      relocate(year)
    rownames(doys) <- NULL
    leaf <- doys
  }

  # Write to csv file
  # Filename will begin with "OUTPUT" and will contain today's date
  earliest_yr <- ifelse(earliest_yr > year1, earliest_yr, year1)
  leaf_out <- paste0("output/OUTPUT-county-leaf-index-", earliest_yr, "-", 
                     str_remove_all(today(), "-"), ".csv")
  write.csv(leaf, leaf_out, row.names = FALSE)

}

if (bloom_table) {
  
  # Read csv with annual bloom index values
  bloom_doy_file <- paste0("output/county-bloom-index-", year1, "-", year, ".csv") 
  bloom <- read.csv(bloom_doy_file, colClasses = c(STATEFP = "character",
                                                 COUNTYFP = "character"))
  
  # Create columns to uniquely identify each county
  bloom <- bloom %>%
    mutate(FIPS = paste0(STATEFP, COUNTYFP),
           name_FIPS = paste0(state, "_", county, "_", FIPS)) %>%
    relocate(c(FIPS, name_FIPS)) %>%
    arrange(FIPS)
  
  if (prior40) {
    earliest_yr <- year - 40
    cols_to_remove <- paste0("doy_", year1:(earliest_yr - 1))
    bloom <- bloom %>%
      select(-all_of(cols_to_remove))
  }
  
  if (!county_row) {
    doys <- t(as.matrix(select(bloom, contains("doy")))) %>%
      data.frame()
    colnames(doys) <- bloom$name_FIPS
    doys <- doys %>%
      mutate(year = as.numeric(str_remove(rownames(doys), "doy_"))) %>%
      relocate(year)
    rownames(doys) <- NULL
    bloom <- doys
  }
  
  # Write to csv file
  # Note: Appending today's date to file name 
  # Filename will begin with "OUTPUT" and will contain today's date
  earliest_yr <- ifelse(earliest_yr > year1, earliest_yr, year1)
  bloom_out <- paste0("output/OUTPUT-county-bloom-index-", earliest_yr, "-", 
                      str_remove_all(today(), "-"), ".csv")
  write.csv(bloom, bloom_out, row.names = FALSE)
  
}



# Visualize spring leaf index anomalies in each county ------------------------#
# And see how they compare to NPN gridded product

# Read csv with annual leaf index values
leaf_doy_file <- paste0("output/county-leaf-index-", year1, "-", year, ".csv") 
leaf <- read.csv(leaf_doy_file, colClasses = c(STATEFP = "character",
                                               COUNTYFP = "character"))

# Restrict values to last 40 years + current
earliest_yr <- year - 40
cols_to_remove <- paste0("doy_", year1:(earliest_yr - 1))
leaf <- leaf %>%
  select(-all_of(cols_to_remove))

# Calculate current-year anomalies (based on 30-year means)
yr30_cols <- paste0("doy_", 1991:2020)
leaf$mean <- apply(leaf[, yr30_cols], 1, mean)
leaf <- leaf %>%
  mutate(anom = doy_2025 - mean)

# Add current year anomalies to county polygons
counties2 <- merge(counties, 
                   select(leaf, STATEFP, COUNTYFP, anom),
                   by.x = c("STATEFP", "COUNTYFP"))

# Limits for nice colors in map
limit <- max(abs(data.frame(counties2[,"anom"])), na.rm = TRUE)
limit_round <- ifelse(limit < 10, 5, floor(limit / 10) * 10)
limits <- limit * c(-1, 1)

# Create map
leaf_a_map <- ggplot(counties2) +
  geom_spatvector(aes(fill = anom), color = NA) +
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

# Download gridded layer to compare
anom_grid <- npn_download_geospatial("si-x:leaf_anomaly",
                                     date = today() - 1)
anom_grid <- rast(anom_grid)
names(anom_grid) <- "anom"
freq(anom_grid)

leaf_ag_map <- ggplot() +
  geom_spatraster(data = anom_grid, aes(fill = anom)) +
  scale_fill_distiller(palette = "Spectral", direction = 1,
                       limit = c(-20, 10)) 
leaf_ag_map
leaf_a_map

# Visualize spring leaf index return intervals in each county -----------------#
# And see how they compare to NPN gridded product

leafl <- leaf %>%
  select(STATEFP, COUNTYFP, state, county, contains("doy_")) %>%
  rename(currentyr = all_of(paste0("doy_", year))) %>%
  pivot_longer(cols = contains("doy_"),
               names_to = "year",
               values_to = "doy") %>%
  mutate(year = as.numeric(str_remove(year, "doy_"))) %>%
  group_by(STATEFP, COUNTYFP, state, county, currentyr) %>%
  summarize(mean40 = mean(doy),
            median40 = median(doy),
            n_earlier = length(doy[doy < currentyr]),
            n_later = length(doy[doy > currentyr]),
            .groups = "keep") %>%
  data.frame()

# If it's the latest observation over the 41 years, then the return interval
# should be 41 years.
# If there's only one year with a later observation, then the return interval
# should be 20.5 (41/2) years since something this extreme happened 2x in 41 yrs
# If it falls right in the middle, then the return interval shoudl be 2 years 
# (41/20), right?


  # How did they calculate it?
  returns <- npn_download_geospatial("si-x:leaf_return_interval",
                                     date = "2025-02-12")
  returns <- rast(returns)
  ann_files <- list.files(paste0(rast_folder, "leaf"), full.names = TRUE)
  ann_files <- str_subset(ann_files, pattern = "ncep_", negate = TRUE)
  ann_files <- str_subset(ann_files, pattern = "1981|1982|1983|1984", 
                          negate = TRUE)
  mrast <- rast(ann_files)
  names(mrast) <- paste0("y", str_sub(names(mrast), 1, 4))
  point_ext <- terra::extract(x = mrast, y = data.frame(x = -111, y = 32), 
                              raw = TRUE)
  priors <- point_ext[2:41]
  summary(priors)
  point_ext[42]
  
  return_ext <- terra::extract(x = returns, y = data.frame(x = -111, y = 32), 
                               raw = TRUE)
  # If late, I think it should be:
  41 / sum(priors > point_ext[42])
  # If early:
  41 / sum(priors < point_ext[42])


  