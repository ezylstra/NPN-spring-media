# Spring rankings, by county
# Originally developed by A. Rosemartin and T. Crimmins
# Adapted from https://github.com/alyssarosemartin/spring-media/tree/main/earliest-spring-ranking
# 7 Feb 2025

library(rnpn)
library(dplyr)
library(lubridate)
library(stringr)
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

layers <- npn_get_layer_details()
filter(layers, name == "si-x:average_leaf_prism")

# Using annual rasters with DOY that requirements for the first leaf Spring 
# Index were met, averaged for 3 species, based on PRISM data

# Get raster layers with annual spring index values for first leaf ------------#

# Location to save rasters
rast_folder <- "resources/rasters/"

# Get layers for all years except current
for (yr in 1981:2024) {
  output_path <- paste0(rast_folder, "leaf/", yr, "_si-x_leaf_prism.tif")
  npn_download_geospatial("si-x:average_leaf_prism",
                          date = paste0(yr, "-01-01"),
                          output_path = output_path)
}

# Obtaining raster for current year from https://data.usanpn.org/geoserver-request-builder/
  # Selecting WCS layer: Spring Indices, Current Year - first Leaf- Spring Index
# TODO: Check that this is correct!
current_leaf_files <- list.files(paste0(rast_folder, "leaf"),
                                 pattern = "ncep_",
                                 full.names = TRUE)
mostcurrent_leaf_file <- last(sort(current_leaf_files))
l2025 <- rast(mostcurrent_leaf_file)

# Resample current year raster to match resolution of prior years
l2025 <- resample(l2025, l2024, method = "bilinear")

# Save new resampled raster
current_file <- paste0(rast_folder, "leaf/2025_si-x_leaf_prism.tif")
writeRaster(l2025, filename = current_file, overwrite = TRUE)

# Create multi-layer raster
leaf_files <- list.files(paste0(rast_folder, "leaf"), full.names = TRUE)
leaf_files <- str_subset(leaf_files, pattern = "ncep_", negate = TRUE)
leaf <- rast(leaf_files)

# County boundaries -----------------------------------------------------------#

# File location
counties_shp <- "resources/us_counties.shp"

# Download first and save to file, if not done already
if (!file.exists(counties_shp)) {
  library(tigris)
    us_counties <- tigris::counties(cb = TRUE)
    us_counties <- vect(us_counties)
    writeVector(us_counties, filename = counties_shp)
}

counties <- vect(counties_shp)
states48 <- state.abb[! state.abb %in% c("AK", "HI")]
states48 <- c(states48, "DC")
counties <- terra::subset(counties, counties$STUSPS %in% c(states48))

# Combine ---------------------------------------------------------------------#

# For each year, will calculate the mean DOY for each county by calculating the
# mean value of raster cells in each county, weighted by the proportion of each
# raster cell that is within the county boundary. 

  # # Check to make sure terra::extract() is working the way I think it should  
  # and produces the same results as exactextractr::exact_extract(). 
  # testr <- leaf[[44]]
  # teste <- terra::extract(x = testr,
  #                         y = counties,
  #                         fun = "mean",
  #                         na.rm = TRUE,
  #                         exact = TRUE)
  # wgts <- terra::extract(x = testr,
  #                        y = counties,
  #                        exact = TRUE)
  # # plot(counties[2131,])
  # # plot(testr, add = TRUE)
  # # plot(counties[2131,], add = TRUE)
  # wgts_sub <- wgts[wgts$ID == 2131,]
  # wgts_sub$values <- wgts_sub[,2] * wgts_sub[,3]
  # 
  # testee <- exactextractr::exact_extract(x = testr,
  #                                        y = sf::st_as_sf(counties),
  #                                        fun = "mean")
  # sum(wgts_sub$values) / sum(wgts_sub[,3])
  # teste[2131, ]
  # testee[2131]
  # Everything works as expected, but exact_extract() is much faster so we'll
  # use that moving forward

names(leaf) <- paste0("y", str_sub(names(leaf), 1, 4))
leaf_doy <- exactextractr::exact_extract(x = leaf,
                                         y = sf::st_as_sf(counties),
                                         fun = "mean")
names(leaf_doy) <- str_replace_all(names(leaf_doy), "mean.y", "doy_")
leaf_doy <- leaf_doy %>%
  mutate(across(everything(), ~ifelse(is.nan(.), NA, .))) %>%
  mutate(STATEFP = counties$STATEFP,
         COUNTYFP = counties$COUNTYFP,
         state = counties$STUSPS,
         county = counties$NAME) %>%
  select(STATEFP, COUNTYFP, state, county, paste0("doy_", 1981:2025))

# Write annual leaf index values to file
leaf_doy_file <- "output/county-leaf-index-1981-2025.csv"
write.csv(leaf_doy,
          file = leaf_doy_file, 
          row.names = FALSE)

# Read csv with annual leaf index values
leaf_doy <- read.csv(leaf_doy_file, 
                     colClasses = c(rep("character", 4), 
                                    rep("numeric", length(1981:2025))))

# What to do with this?
  # For each county, list of top 10 earliest/latest years
  # For those counties where spring index has been reached in 2025:
    # Number of days earlier/later than historical average (like the spring
    # index maps on NPN website, but aggregated by county)
    # Identify those in top 10/5/1 earliest/latest

# Calculate anomalies and add 2025 values to SpatVector
prior_yr_cols <- paste0("doy_", 1981:2024)
leaf_doy$mean <- apply(leaf_doy[, prior_yr_cols], 1, mean)
leaf_doy$median <- apply(leaf_doy[, prior_yr_cols], 1, median)
leaf_doy <- leaf_doy %>%
  mutate(across(.cols = starts_with("doy_"),
                .fns = ~.x - mean,
                .names = "{.col}_anom"))

counties2 <- merge(counties, select(leaf_doy, STATEFP, COUNTYFP, doy_2025_anom),
                   by.x = c("STATEFP", "COUNTYFP"))

limit <- max(abs(counties2$doy_2025_anom)) * c(-1, 1)
ggplot(counties2) +
  geom_spatvector(aes(fill = doy_2025_anom), color = NA) +
  scale_fill_distiller(palette = "Spectral", limit = limit, direction = 1)

# Comparing some of the output with what's on NPN website for spring indices
# (The difference here is that we're aggregating things by county)
leaf_doy25 <- leaf_doy %>%
  select(STATEFP, COUNTYFP, state, county, doy_2025, mean, doy_2025_anom)

leaf_doy25 %>%
  filter(state == "TX", 
         county %in% c("Hidalgo", "Cameron", "San Patricio", "Winkler"))
leaf_doy25 %>%
  filter(state == "OK", county %in% c("Pittsburg"))
# Some county values are extreme (on the early side), but looking closely that's
# occurring because we have big counties in the West and we're getting a value
# for any county that has one more more cells that aren't NA. It's really
# skewing things when only a few cells in a large county have reached the 
# threshold. 

# Probably should go back and only calculate county-level values when all or 
# most of the area has reached the index.
leaf <- rast(leaf_files)
names(leaf) <- paste0("y", str_sub(names(leaf), 1, 4))
leaf_doy <- exactextractr::exact_extract(x = leaf,
                                         y = sf::st_as_sf(counties),
                                         c("mean", "count"))
# The count function provides the sum of all cell coverage fractions, ignoring
# NA values. Counts should be the same for all prior years. Could use the ratio
# for 2025/prior values to see whether the 2025 mean is based on too few cells

# Random test
all.equal(leaf_doy$count.y2016, leaf_doy$count.y2024)

leaf_doy <- leaf_doy %>%
  mutate(frac2025 = count.y2025 / count.y2024)
hist(leaf_doy$frac2025[leaf_doy$frac2025 > 0], breaks = 30)
leaf_doy <- leaf_doy %>%
  mutate(mean.y2025.new = ifelse(frac2025 < 0.90, NA, mean.y2025)) %>%
  mutate(across(everything(), ~ifelse(is.nan(.), NA, .))) %>%
  mutate(STATEFP = counties$STATEFP,
         COUNTYFP = counties$COUNTYFP,
         state = counties$STUSPS,
         county = counties$NAME)

# Calculate anomalies and add 2025 values to SpatVector
prior_yr_cols <- paste0("mean.y", 1981:2024)
leaf_doy$mean <- apply(leaf_doy[, prior_yr_cols], 1, mean)
leaf_doy$median <- apply(leaf_doy[, prior_yr_cols], 1, median)
leaf_doy <- leaf_doy %>%
  mutate(across(.cols = starts_with("mean.y"),
                .fns = ~.x - mean,
                .names = "{.col}_anom"))

counties3 <- merge(counties, select(leaf_doy, STATEFP, COUNTYFP, mean.y2025.new_anom),
                   by.x = c("STATEFP", "COUNTYFP"))

limit <- max(abs(counties3$mean.y2025.new_anom)) * c(-1, 1)
ggplot(counties3) +
  geom_spatvector(aes(fill = mean.y2025.new_anom), color = NA) +
  scale_fill_distiller(palette = "Spectral", limit = limit, direction = 1)
# OK
  




