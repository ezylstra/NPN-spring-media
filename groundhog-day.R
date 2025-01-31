# Groundhog day script
# Adapted from https://github.com/alyssarosemartin/spring-media/blob/main/groundhog/groundhog.R#L1
# 31 Jan 2025

library(rnpn)
library(dplyr)
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
         first_yes_day, numdays_since_prior_no, last_yes_day, 
         numdays_until_next_no) %>%
  rename(lat = latitude,
         lon = longitude, 
         elev = elevation_in_meters, 
         id = individual_id, 
         phenophase = phenophase_description,
         first_yes = first_yes_day,
         prior_no = numdays_since_prior_no,
         last_yes = last_yes_day,
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

# Note: flower and leaf icons already downloaded into icons folder
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
  geom_image(data = leaf, aes(x = lon, y = lat, image = "icons/leaf.png")) +
  geom_image(data = flower, aes(x = lon, y = lat, image = "icons/flower.ico")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) 
ghd_plot

ghd_filename <- paste0("output/groundhog-day-", year, ".png")
ggsave(filename = ghd_filename,
       plot = ghd_plot,
       width = 6,
       height = 5,
       units = "in",
       dpi = 600,
       bg = "transparent")



# Map using leaflet (OLD) -----------------------------------------------------#

# library(leaflet)
# library(leaflet.extras)
# library(mapview)
# library(webshot2)
# library(htmlwidgets)
# 
# # Prep point data viz, by creating and sizing leaf and flower icons
# my_icons <- iconList(
#   ileaf <- makeIcon(iconUrl = "https://www.freeiconspng.com/uploads/leaf-icon-16.png",
#                     iconWidth = 20, iconHeight = 20),
#   iflower <- makeIcon(iconUrl = "https://www.freeiconspng.com/uploads/purple-flower-png-17.png",
#                       iconWidth = 20, iconHeight = 20)
# )
# 
# # Prep raster viz with  bins and palette (original, with 2 warm bins)
# # bins <- c(-500, -20, 20, 200, 500) # Why create a break at +200?
# # pal <- colorBin(c("#6fa8d6", "#f7f0da", "#fa824d", "#f44f07"), 
# #                 domain = groundhog, bins = bins,  na.color = "transparent")
# 
# # Prep raster viz with  bins and palette (new, with 1 warm bin)
# pal <- colorBin(cols,
#                 domain = groundhog, 
#                 bins = breaks,  
#                 na.color = "transparent")
# 
# # Create map
# m <- leaflet(leaf) %>% 
#   addRasterImage(groundhog, colors = pal, opacity = 1) %>%
#   addMarkers(lng = ~leaf$lon, lat = ~leaf$lat, icon = ileaf) %>%
#   addMarkers(lng = ~flower$lon, lat = ~flower$lat, icon = iflower) %>%
#   setMapWidgetStyle(list(background= "transparent"))
# m
# 
# saveWidget(m, "temp.html", selfcontained = FALSE)
# webshot("temp.html", file = "leaflet_map.png", cliprect = "viewport")
