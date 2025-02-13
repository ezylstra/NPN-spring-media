# Spring rankings, by county
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

# This code extracts day of year of leaf out and bloom by county 
# (weighted average of pixels that fall within the county) and ranks how the 
# current spring compares to those in previous years.

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

# Set parameters --------------------------------------------------------------#

# Current year
year <- 2025

# First year to download layers for
year1 <- 1981

# Prior years
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
  # Wasn't sure what to do with ties (see WV: Wood for an example), so for now,
  # put more recent year first
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
  arrange(STATEFP, COUNTYFP, rank, desc(year)) %>%
  data.frame()

leafrr <- leafr %>%
  group_by(STATEFP, COUNTYFP, state, county, mean, median) %>%
  summarize(early01 = year[1],
            early02 = year[2],
            early03 = year[3],
            early04 = year[4],
            early05 = year[5],
            early06 = year[6],
            early07 = year[7],
            early08 = year[8],
            early09 = year[9],
            early10 = year[10],
            late01 = year[length(prior_years)],
            late02 = year[length(prior_years) - 1],
            late03 = year[length(prior_years) - 2],
            late04 = year[length(prior_years) - 3],
            late05 = year[length(prior_years) - 4],
            late06 = year[length(prior_years) - 5],
            late07 = year[length(prior_years) - 6],
            late08 = year[length(prior_years) - 7],
            late09 = year[length(prior_years) - 8],
            late10 = year[length(prior_years) - 9],
            .groups = "keep") %>%
  data.frame()

# Create indicators to see whether recent years in top 10 earliest
leafrr$e2024 <- apply(leafrr[, grepl("early", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2024), 1, 0))
leafrr$e2023 <- apply(leafrr[, grepl("early", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2023), 1, 0))
leafrr$e2022 <- apply(leafrr[, grepl("early", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2022), 1, 0))
leafrr$e2021 <- apply(leafrr[, grepl("early", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2021), 1, 0))
leafrr$e2020 <- apply(leafrr[, grepl("early", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2020), 1, 0))

# Create indicators to see whether recent years in top 10 latest
leafrr$l2024 <- apply(leafrr[, grepl("late", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2024), 1, 0))
leafrr$l2023 <- apply(leafrr[, grepl("late", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2023), 1, 0))
leafrr$l2022 <- apply(leafrr[, grepl("late", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2022), 1, 0))
leafrr$l2021 <- apply(leafrr[, grepl("late", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2021), 1, 0))
leafrr$l2020 <- apply(leafrr[, grepl("late", colnames(leafrr))], 1, 
                      function(x) ifelse(any(x == 2020), 1, 0))

leafrr$e20_sum <- apply(leafrr[, paste0("e202", 0:4)], 1, sum)
leafrr$l20_sum <- apply(leafrr[, paste0("l202", 0:4)], 1, sum)
# Average number of 2020 years in the top 10 earliest or latest
summary(leafrr[, c("e20_sum", "l20_sum")])
# Number of 2020 years in top 10 earliest or latest
count(leafrr, e20_sum)
count(leafrr, l20_sum)
# Recent years more likely to be in early than late.

# Could also present ranked years in table where each column is a county
leafrr_t <- leafrr %>%
  # Using FIPS codes because there are a few county names that are duplicates
  # (Baltimore county and city; 24005, 24510)
  mutate(fips = paste0("fp_", STATEFP, COUNTYFP)) %>%
  select(fips, 
         paste0("early0", 1:9), "early10", 
         paste0("late0", 1:9), "late10") %>%
  pivot_longer(cols = -fips, 
               names_to = "rank") %>%
  # group_by(state_co) %>%
  # mutate(row = row_number()) %>%
  pivot_wider(names_from = fips,
              values_from = value) %>%
  # select(-row) %>%
  data.frame()

# Create dataframe with spring rankings that include the current year
# If ties, put more recent year first
year_cols <- paste0("doy_", c(prior_years, year))
leafr_current <- leaf %>%
  select(-contains("anom")) %>%
  filter(!is.na(get(paste0("doy_", year)))) %>%
  pivot_longer(cols = all_of(year_cols),
               names_to = "year",
               values_to = "doy") %>%
  mutate(year = as.numeric(str_remove(year, "doy_"))) %>%
  group_by(STATEFP, COUNTYFP, state, county, mean, median) %>%
  mutate(rank = rank(doy)) %>%
  arrange(STATEFP, COUNTYFP, rank, desc(year)) %>%
  data.frame()
leafrr_current <- leafr_current %>%
  group_by(STATEFP, COUNTYFP, state, county, mean, median) %>%
  summarize(earliest = ifelse(year[1] == 2025, 1, 0),
            top5 = ifelse(any(year[1:5] == 2025), 1, 0),
            top10 = ifelse(any(year[1:10] == 2025), 1, 0),
            .groups = "keep") %>%
  data.frame() %>%
  mutate(rank = case_when(
    earliest == 1 ~ "Earliest", 
    top5 == 1 ~ "Top 5", 
    top10 == 1 ~ "Top 10", 
    .default = "Not early"
  )) %>%
  mutate(rank = factor(rank, levels = c("Not early", "Top 10", 
                                        "Top 5", "Earliest")))
filter(leafrr_current, top10 == 1)

# Create map
counties3 <- left_join(counties, 
                       select(leafrr_current, STATEFP, COUNTYFP, rank),
                       by = c("STATEFP", "COUNTYFP"))

rank_map <- ggplot(counties3) +
  geom_spatvector(aes(fill = rank), color = NA) +
  scale_fill_manual(values = c("Not early" = "white", 
                               "Top 10" = "yellow",
                               "Top 5" = "orange",
                               "Earliest" = "red"), 
                    breaks = c("Not early", "Top 10", "Top 5"),
                    drop = TRUE) +
  labs(title = paste0("Leaf spring index, ", year, " ranking"),
       fill = "Rank")
rank_map

  