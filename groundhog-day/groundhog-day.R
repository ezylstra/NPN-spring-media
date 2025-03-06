# Groundhog day script
# Adapted from https://github.com/alyssarosemartin/spring-media/blob/main/groundhog/groundhog.R#L1
# 31 Jan 2025

library(here)
library(rnpn)
library(dplyr)
library(lubridate)
library(ggplot2)
library(terra)
library(tidyterra)
library(ggimage)

# Set parameters --------------------------------------------------------------#

# Year
year <- 2025

# Define breaks that will be used to delineate areas cooler/warmer than 30-year
# normals. Breaks = AGDD anomolies on Groundhog day (deg F).
breaks <- c(-500, -20, 20, 500)

# Colors for map (blue = cooler, offwhite = normal, orange = warmer)
cols <- c("#6fa8d6", "#f7f0da", "#fa824d")
# plot(1:3, rep(1, 3), pch = 19, cex = 5, col = cols)

# Helper calls ----------------------------------------------------------------#

layers <- npn_get_layer_details()
# layers$abstract

phenoclasses <- npn_pheno_classes() %>% data.frame()
# phenoclasses

# Aquire, format phenometric data ---------------------------------------------#

# Download individual phenometrics (first yeses) so far this year: 
# breaking leaf bud / initial growth (class 1) and open flowers (class 7)
df <- npn_download_individual_phenometrics(
  request_source = 'erinz', 
  years = year,
  pheno_class_ids = c(1, 7)
)

# Remove observations outside the continental US (for now, assuming
# that if state is missing, location is in lower 48)
states48 <- state.abb[! state.abb %in% c("AK", "HI")]
df48 <- df %>%
  filter(state %in% c(states48, "-9999"))

# Simplify first observation data
df48 <- df48 %>%
  select(site_id, latitude, longitude, elevation_in_meters, state, 
         common_name, individual_id, phenophase_description, pheno_class_id,
         first_yes_doy, numdays_since_prior_no, last_yes_doy, 
         numdays_until_next_no) %>%
  rename(lat = latitude,
         lon = longitude, 
         elev = elevation_in_meters, 
         id = individual_id, 
         phenophase = phenophase_description,
         first_yes = first_yes_doy,
         prior_no = numdays_since_prior_no,
         last_yes = last_yes_doy,
         next_no = numdays_until_next_no) %>%
  mutate(prior_no = ifelse(prior_no == -9999, NA, prior_no),
         next_no = ifelse(next_no == -9999, NA, next_no),
         elev = ifelse(elev == -9999, NA, elev),
         state = ifelse(state == "-9999", NA, state)) %>%
  data.frame()

# Look at species included
count(df48, common_name)
# Delete any 'ohi'a lehua records
df48 <- df48 %>%
  filter(common_name != "'ohi'a lehua")

# Delete any observations after groundhog day (in case we're running this 
# after the fact)
groundhog_doy <- yday(ymd(paste0(year, "-02-02")))
df48 <- df48 %>%
  filter(first_yes <= groundhog_doy)

# Distribution of first observation dates
ggplot(df48) +
  geom_histogram(aes(x = first_yes)) +
  facet_wrap(~pheno_class_id, ncol = 1)

# How many observations and individuals?
df48 %>%
  group_by(pheno_class_id) %>%
  summarize(n_obs = n(),
            n_indiv = n_distinct(id))

# Keep just one "first yes" for each plant and phenophase class 
df48 <- df48 %>%
  arrange(pheno_class_id, common_name, id, first_yes) %>%
  distinct(pheno_class_id, id, .keep_all = TRUE)

# How many observations, and what proportion have prior no?
df48 %>%
  group_by(pheno_class_id) %>%
  summarize(n_obs = n(),
            with_prior_no = sum(is.na(prior_no))) %>%
  data.frame() %>%
  mutate(prop_prior_no = with_prior_no / n_obs)

leaf <- subset(df48, pheno_class_id == 1)
flower<- subset(df48, pheno_class_id == 7)

# Aquire forecasted AGDD anomalies --------------------------------------------#

# Acquire raster forecast anomaly for Feb 2 (32 base)
layers %>%
  filter(name == "gdd:agdd_anomaly") %>%
  select(name, title, abstract, dimension.name)

gh_day <- paste0(year, "-02-02")
groundhog <- npn_download_geospatial(coverage_id = "gdd:agdd_anomaly", 
                                     date = gh_day)  
# plot(groundhog)
# hist(groundhog)

# Create map ------------------------------------------------------------------#

# Note: flower and leaf icons already downloaded into resources folder
# (leaf.svg, flower.ico)

# Convert to spatraster and bin data
ghd <- rast(groundhog)
ghd_class <- terra::classify(x = ghd, 
                             rcl = breaks, 
                             include.lowest = TRUE)
freq(ghd_class)
plot(ghd_class)

ghd_plot <- ggplot() +
  geom_spatraster(data = ghd_class, maxcell = Inf) +
  scale_fill_manual(values = cols, na.value = "transparent") +
  geom_image(data = leaf, 
             aes(x = lon, y = lat, image = here("icons", "leaf.png"))) +
  geom_image(data = flower, 
             aes(x = lon, y = lat, image = here("icons", "flower.ico"))) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) 
ghd_plot

ghd_filename <- here("groundhog-day/output",
                     paste0("groundhog-day", year, ".png"))

# Save to file (commented out for now, so we don't accidentally overwrite anything)
# ggsave(filename = ghd_filename,
#        plot = ghd_plot,
#        width = 6,
#        height = 5,
#        units = "in",
#        dpi = 600,
#        bg = "transparent")
