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
# Pick out your focal species - do some a priori learnin'
# ------------------------------------------------------------------------------
foc_sp <- 'Picea pungens'

# Example:
# Picea pungens (blue spruce): https://pubs.usgs.gov/pp/p1650-a/pages/pipushrt.pdf 

# Look at these maps from the USFS to determine which states data you will likely need
# Group 1
# Juniperus osteosperma (Utah juniper): https://pubs.usgs.gov/pp/p1650-a/pages/juostrim.pdf
# Group 2
# Quercus gambellii (Gambel oak): https://pubs.usgs.gov/pp/p1650-a/pages/gambtrim.pdf 
# Group 3
# Pinus edulis (pinyon pine): https://pubs.usgs.gov/pp/p1650-a/pages/piedtrim.pdf 


# ------------------------------------------------------------------------------
# Download FIA data - only run once
# ------------------------------------------------------------------------------
# FIA guide will be helpful: https://www.fia.fs.fed.us/library/database-documentation/current/ver80/FIADB%20User%20Guide%20P2_8-0.pdf
# FIA data are organized by state, so put in the ones where your species occurs:
fia_states <- c('WY','UT','CO','NM') 

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

# Get information about species names...
# Download the species reference table from here: https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html
# Note that you don't include 'REF_' at the beginning of the table name, and use 'ref' for states
# We just saved this to our data directory, but didn't read it in
getFIA(states = 'ref', tables = c('SPECIES'), dir = '../data')

# Let's read this in as a CSV
spp_ref <- read_csv('../data/REF_SPECIES.csv')
# Make a simplified df with code and common name, species, genus
spp_names <- dplyr::select(spp_ref, SPCD, COMMON_NAME, GENUS, SPECIES) %>% 
  mutate(latin_name = paste(GENUS, SPECIES))
# Match our focal species latin name to it's code 
foc_spcd <- filter(spp_names, latin_name == foc_sp)[['SPCD']]

# Check how many records of this species are present in the FIA databases from these states 
left_join(fia_df_seed, spp_names) %>% 
  group_by(SPCD, COMMON_NAME, latin_name) %>% 
  tally(sort = TRUE) %>% 
  View()


# ------------------------------------------------------------------------------
#  Create presence-absence data for a species of interest
# ------------------------------------------------------------------------------
# Recall that we just saved the FIA species code (SPCD) of our focal species
foc_sp
foc_spcd

# Create a simplified dataframe of foc_sp presence or absence at FIA plots
fia_foc_sp <- fia_df_seed %>% 
  # Each row indicates a subplot where the species is present
  # Create a presence/absense indicator column - this will be response in binomial model
  mutate(foc_sp_PRES = ifelse(SPCD == foc_spcd, 1, 0)) %>% 
  # Group by PLT_CN (bit confusing... but PLT_CN will link us to CN in the PLOT table, which is the unique visit ID)
  group_by(PLOT, PLT_CN) %>% 
  # Max of presence indicator will summarize whether foc_sp present in each plot
  summarise(foc_sp_PRES = as.factor(max(foc_sp_PRES)))


# ------------------------------------------------------------------------------
#  Map plots with seedlings of species of interest
# ------------------------------------------------------------------------------

# Get the location information associated with these plots
foc_sp_plot_coords <- fia_df[['PLOT']] %>% 
  dplyr::select(CN, PLOT, INVYR, LAT, LON) %>% 
  # Join our presence/absense data to a list of plots and their coordinates
  # Why left_join?
  left_join(fia_foc_sp, ., by = c('PLT_CN' = 'CN', 'PLOT'='PLOT')) 

# Make this a spatial object
foc_sp_sf  <- st_as_sf(foc_sp_plot_coords, coords = c('LON','LAT'), crs = 4326) # coords must be in x, y order!

# Get US state boundaries to make maps more interpretable
state_bounds <- us_states(map_date = NULL, resolution = c("low"), states = NULL)

# Map... looks reasonable?
# Check against the USFS maps...link above
tm_shape(state_bounds, bbox = foc_sp_sf) +
  tm_borders() +
  tm_fill(col = 'grey90') +
  tm_shape(foc_sp_sf) +
  tm_dots(col = 'grey30') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(col = 'red', border.col = 'black', size = 0.3, shape = 21, alpha = 0.7)

# Learn more about tmap: https://github.com/r-tmap/tmap
tm_shape(state_bounds, bbox = foc_sp_sf) +
  tm_borders() +
  tm_fill(col = 'grey90') +
  tm_shape(foc_sp_sf) +
  tm_dots(col = 'foc_sp_PRES', 
          palette = c('grey30','red')) 
  
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

# This can be slow, you may decide to save this as a Rdata object for quick access later
foc_sp_bio <- terra::extract(bio_hist, foc_sp_sf)

# Select some climate variables... this is important and there are lots of ways to do it
# This is just a simplified example, picking a priori and checking for autocorrelation
# See function ggpairs for pairwise scatterplots... might be useful... could cbind presence/absence too
foc_sp_bio %>% 
  ggcorr(label = TRUE)


foc_sp_hist_sf <- foc_sp_bio %>% 
  as.data.frame() %>% 
  dplyr::select(temp_min = bio6, temp_seas = bio4, precip_dry = bio17) %>% 
  cbind(foc_sp_sf, .) %>% 
  mutate(temp_min = temp_min,
         temp_seas = temp_seas)

# Map... looks reasonable?
climate_map <- 
  tm_shape(state_bounds, bbox = foc_sp_hist_sf) +
  tm_fill(col = 'grey60') +
  tm_shape(foc_sp_hist_sf) +
  tm_dots(col = 'precip_dry', 
          size = 0.09, 
          palette = 'BrBG',
          midpoint = NA, 
          style = 'cont') +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6) 
climate_map

# ------------------------------------------------------------------------------
#  Fit a model
# ------------------------------------------------------------------------------
# This section is extremely bare-bones. Here is the opportunity for students to 
# use what they know from stats courses etc. to fit and evaluate a model following
# the lastest and greatest best practices.

# Make a non-spatial dataframe 
foc_sp_hist_df <- foc_sp_hist_sf %>% 
  st_drop_geometry() %>% 
  mutate(presence = as.numeric(as.character(foc_sp_PRES)))


# Fit a model
foc_sp_glm <- glm(presence ~ temp_min + temp_seas + precip_dry, 
                data = foc_sp_hist_df,
                family = binomial(link = "logit"))

summary(foc_sp_glm)

# Needs work... but shouled evaluate response curves
ggplot(foc_sp_hist_df, aes(x = temp_seas, y = presence)) +
  geom_point() +
  geom_smooth(method = 'glm', 
              se = F, 
              method.args = list(family = binomial(link = "logit")))

# Predict historical across space
hist_preds <- bio_hist[[c('bio6','bio4','bio17')]]
names(hist_preds) <- c('temp_min','temp_seas','precip_dry')

# Add Montana state boundary to crop the future prediction
pred_bounds <- us_states(map_date = NULL, resolution = c("low"), states = c('ID','MT',fia_states))

hist_preds <- terra::crop(hist_preds, pred_bounds)
hist_space <- terra::predict(hist_preds, model = foc_sp_glm, type = 'response')

hist_map <- 
  tm_shape(hist_space, bbox = pred_bounds, raster.downsample = FALSE) +
  tm_raster(title = 'Historical probability',
            palette = 'PuRd', style = 'cont') +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
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


foc_sp_cmip5 <- terra::extract(bio_cmip5, foc_sp_sf)

cmip5_preds <- bio_cmip5[[c('bio6','bio4','bio17')]]
names(cmip5_preds) <- c('temp_min','temp_seas','precip_dry')
cmip5_preds <- terra::crop(cmip5_preds, pred_bounds)
cmip5_space <- terra::predict(cmip5_preds, model = foc_sp_glm, type = 'response')

cmip5_map <- 
  tm_shape(cmip5_space, bbox = pred_bounds) +
  tm_raster(title = 'Future probability',
            palette = 'PuRd', style = 'cont') +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6) 
cmip5_map

# Arrow between historical and future to get a eyeball sense of differences

# Map difference between historical and future
sdm_difference <- cmip5_space - hist_space

diff_map <- 
  tm_shape(sdm_difference, bbox = pred_bounds) +
  tm_raster(title = 'Difference in probability',
            palette = 'RdGy', 
            style = 'cont', 
            breaks = seq(-.05, .05, by=.01)) +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6) 
diff_map











