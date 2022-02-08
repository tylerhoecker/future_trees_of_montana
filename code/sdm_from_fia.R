# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
# This script will work best if you download/clone the entire GitHub repo where it lives: 
# https://github.com/tylerhoecker/future_trees_of_montana/

# Otherwise, the script assumes it is stored in a directory called 'code' and that
# there is another directory at the same level called 'data'. 
# The root directory (`future_trees_of_montana`) could be named anything 
# Example:
# - future_trees_of_montana/code/sdm_from_fia.R
# - future_trees_of_montana/data/

# Install missing packages and load them
using <- function(...) {
  libs <- unlist(list(...))
  req <- unlist(lapply(libs,require,character.only=TRUE))
  need <- libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

using('tidyverse','sf','tmap','rFIA','USAboundaries','raster','terra','GGally')

# ------------------------------------------------------------------------------
# Download FIA data - only run once
# ------------------------------------------------------------------------------
# FIA guide will be helpful: https://www.fia.fs.fed.us/library/database-documentation/current/ver80/FIADB%20User%20Guide%20P2_8-0.pdf
# FIA data are organized by state, so let's specify a few likely states to pull from
fia_states <- c('ID','MT') 

# Use `getFIA` to download data... did not make object, just saving to data directory
# This package give us the option to save downloaded datasets in a specified directory
# This saves time of downloading the data every time
# We are downloading all of the TABLES, later we can read in specific TABLES
getFIA(states = fia_states, dir = '../data') # Note that '../' means go up one directory

# Explore the structure of FIA data that is automatically printed after we download
# Which tables should we focus on?

# ------------------------------------------------------------------------------
# Read in and Explore FIA data
# ------------------------------------------------------------------------------
# Now read the data table into memory
fia_df <- readFIA(dir = '../data', 
                  states = fia_states, 
                  tables = c('SEEDLING','PLOT')) 

# Access the tables we read in to memory.
# When FIA data is loaded into R with readFIA or getFIA, those data are stored 
# in an FIA.Database object. An FIA.Database object is essentially just a list, 
# and users can access individual tables with $, [''], and [['']] operators
fia_df[['SEEDLING']] %>% 
  glimpse()

# Save this table as a dataframe
fia_df_seed <- fia_df[['SEEDLING']]

# Download the species reference table from here: https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html
# Note that you don't include 'REF_' at the beginning of the table name, and use 'ref' for states
# We just saved this to our data directory, but didn't read it in
getFIA(states = 'ref', tables = c('SPECIES'), dir = '../data')

# Let's read this in as a CSV
spp_ref <- read_csv('../data/REF_SPECIES.csv')
# Make a simplified df with code and common
spp_names <- dplyr::select(spp_ref, SPCD, COMMON_NAME)

# Which species are present as seedlings in these states?
# Let's see which are more or less abundant
left_join(fia_df_seed, spp_names) %>% 
  group_by(SPCD, COMMON_NAME) %>% 
  tally(sort = TRUE) %>% 
  View()


# ------------------------------------------------------------------------------
#  Create presence-absence data for a species of interest
# ------------------------------------------------------------------------------
# Let's focus on Utah juniper (Juniperus osteosperma) SPCD = 65

# Create a simplified dataframe of JUOS presence or absence at each plot
fia_juos <- fia_df_seed %>% 
  # Create a special indicator column - this will be response in binomial model
  mutate(JUOS_PRES = ifelse(SPCD == 65, 1, 0)) %>% 
  # Group by PLT_CN (bit confusing... but PLT_CN will link us to CN in the PLOT table, which is the unique visit ID)
  group_by(PLOT, PLT_CN) %>% 
  # Max of presence indicator will summarize whether JUOS present in each plot
  summarise(JUOS_PRES = as.factor(max(JUOS_PRES)))


# ------------------------------------------------------------------------------
#  Map plots with seedlings of species of interest
# ------------------------------------------------------------------------------
# Get the location information associated with these plots
juos_plot_coords <- fia_df[['PLOT']] %>% 
  dplyr::select(CN, PLOT, INVYR, LAT, LON) %>% 
  left_join(fia_juos, ., by = c('PLT_CN' = 'CN', 'PLOT'='PLOT'))

# Make this a spatial object
juos_sf  <- st_as_sf(juos_plot_coords, coords = c('LON','LAT'), crs = 4326) # coords must be in x, y order!

# Get US state boundaries to make maps more interpretable
state_bounds <- us_states(map_date = NULL, resolution = c("low"), states = NULL)

# Map... looks reasonable?
# Learn more about tmap: https://github.com/r-tmap/tmap
tm_shape(state_bounds, bbox = juos_sf) +
  tm_borders() +
  tm_fill(col = 'grey90') +
  tm_shape(juos_sf) +
  tm_dots(col = 'JUOS_PRES',
          size = 0.08,
          palette = c("grey10", "red"), 
          alpha = 0.5) 


# ------------------------------------------------------------------------------
#  Extract and explore climate data for this species
# ------------------------------------------------------------------------------
# Use the `getdata` function in raster to download data from WorldClim
# WorldClim provides 19 'biologically meaningful' variables: https://www.worldclim.org/data/bioclim.html
bio_hist <- raster::getData('worldclim', var = 'bio', res = 2.5, path = '../data/')

# Once downloaded, read in like this (I changed the name of the directory manually)
bio_hist <- list.files('../data/wc2-5/', 
                       pattern = '*.bil$', 
                       full.names = T) %>% 
  stack()

juos_bio <- terra::extract(bio_hist, juos_sf)

# Select some climate variables... this is important and there are lots of ways to do it
# This is just a simplified example, picking a priori and checking for autocorrelation
juos_bio %>% 
  ggcorr(label = TRUE)


juos_hist_sf <- juos_bio %>% 
  as.data.frame() %>% 
  dplyr::select(temp_min = bio6, temp_seas = bio4, precip_dry = bio17) %>% 
  cbind(juos_sf, .) %>% 
  mutate(temp_min = temp_min,
         temp_seas = temp_seas)

# Map... looks reasonable?
climate_map <- 
  tm_shape(state_bounds, bbox = juos_hist_sf) +
  tm_fill(col = 'grey60') +
  tm_shape(juos_hist_sf) +
  tm_dots(col = 'precip_dry', 
          size = 0.09, 
          palette = 'BrBG',
          midpoint = NA, 
          style = 'cont') +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(juos_sf, JUOS_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6) 
climate_map

# ------------------------------------------------------------------------------
#  Fit a model
# ------------------------------------------------------------------------------
# This section is extremely bare-bones. Here is the opportunity for students to 
# use what they know from stats courses etc. to fit and evaluate a model following
# the lastest and greatest best practices.

# Make a non-spatial dataframe 
juos_hist_df <- juos_hist_sf %>% 
  st_drop_geometry() %>% 
  mutate(presence = as.numeric(as.character(JUOS_PRES)))


# Fit a model
juos_glm <- glm(presence ~ temp_min + temp_seas + precip_dry, 
               data = juos_hist_df,
               family = binomial(link = "logit"))

summary(juos_glm)

# Needs work...
ggplot(juos_hist_df, aes(x = temp_seas, y = presence)) +
  geom_point() +
  geom_smooth(method = 'glm', 
              se = F, 
              method.args = list(family = binomial(link = "logit")))

# Predict historical across space
hist_preds <- bio_hist[[c('bio6','bio4','bio17')]]
names(hist_preds) <- c('temp_min','temp_seas','precip_dry')
hist_preds <- terra::crop(hist_preds, juos_sf)
hist_space <- terra::predict(hist_preds, model = juos_glm, type = 'response')

hist_map <- 
  tm_shape(hist_space, bbox = juos_hist_sf, raster.downsample = FALSE) +
  tm_raster(palette = 'YlOrRd', style = 'cont') +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(juos_sf, JUOS_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6) 
hist_map

# ------------------------------------------------------------------------------
#  Project the model into the future
# ------------------------------------------------------------------------------

# Download future projected climate data
bio_cmip5 <- raster::getData('CMIP5', 
                            var = 'bio', 
                            res = 2.5, 
                            model = 'HD',
                            rcp = 60,
                            year = 50,
                            path = '../data/')

# Once downloaded, read in like this (I changed the name of the directory manually)
bio_cmip5 <- list.files('../data/cmip5/2_5m/', 
                       pattern = '*.tif$', 
                       full.names = T) %>% 
  stack() %>% 
  # Rename these so they are the same as the model
  `names<-` (names(bio_hist))


juos_cmip5 <- terra::extract(bio_cmip5, juos_sf)

cmip5_preds <- bio_cmip5[[c('bio6','bio4','bio17')]]
names(cmip5_preds) <- c('temp_min','temp_seas','precip_dry')
cmip5_preds <- terra::crop(cmip5_preds, juos_sf)
cmip5_space <- terra::predict(cmip5_preds, model = juos_glm, type = 'response')

cmip5_map <- 
  tm_shape(cmip5_space, bbox = juos_cmip5_sf) +
  tm_raster(palette = 'YlOrRd', style = 'cont') +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(juos_sf, JUOS_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6) 
cmip5_map


# Map difference between historical and future
sdm_difference <- cmip5_space - hist_space

diff_map <- 
  tm_shape(cmip5_space, bbox = juos_cmip5_sf) +
  tm_raster(palette = 'cividis', style = 'cont') +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(juos_sf, JUOS_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6) 
diff_map











