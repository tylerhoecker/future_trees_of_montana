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

using('tidyverse','sf','tmap','rFIA','USAboundaries','raster','terra','GGally', 'rgdal', 'gridExtra','e1071', 'ROCR', 'gbm','dismo', 'PresenceAbsence')




setwd("TYLERS WORKING DIRECTORY!!!!")
# ------------------------------------------------------------------------------
# Pick out your focal species - do some a priori learnin'
# ------------------------------------------------------------------------------
#foc_sp <- 'Juniperus osteosperma'
#foc_sp <- 'Quercus gambelii'
foc_sp <- 'Pinus edulis'


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
fia_states <- c('ID','WY','UT','CO','NM', 'AZ', 'NV','WA','OR','MT','CA') # all western states


# Use `getFIA` to download data... did not make object, just saving to data directory
# This package give us the option to save downloaded datasets in a specified directory
# This saves time of downloading the data every time
# We are downloading all of the TABLES, later we can read in specific TABLES

getFIA(states = fia_states, dir = '../data', tables = c('SEEDLING','PLOT')) # Note that '../' means go up one directory

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
getFIA(states = 'ref', tables = c('SPECIES'), dir = '/data')

# Let's read this in as a CSV
spp_ref <- read_csv('/data/REF_SPECIES.csv')
# Make a simplified df with code and common name, species, genus
spp_names <- dplyr::select(spp_ref, SPCD, COMMON_NAME, GENUS, SPECIES) %>% 
  mutate(latin_name = paste(GENUS, SPECIES))
# Match our focal species latin name to it's code 
foc_spcd <- filter(spp_names, latin_name == foc_sp)[['SPCD']]

# Check how many records of this species are present in the FIA databases from these states 
left_join(fia_df_seed, spp_names) %>% 
  group_by(SPCD, COMMON_NAME, latin_name) %>% 
  tally(sort = TRUE)


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
  summarise(foc_sp_PRES = max(foc_sp_PRES))

# ------------------------------------------------------------------------------
#  Map plots with seedlings of species of interest
# ------------------------------------------------------------------------------

# Get the location information associated with these plots
foc_sp_plot_coords <- fia_df[['PLOT']] %>% 
  dplyr::select(CN, PLOT, INVYR, LAT, LON) %>% 
  drop_na(LON) %>% #get rid of any NA in lat or long columns
  drop_na(LAT) %>%
  #Join our presence/absence data to a list of plots and their coordinates
  #Why left_join?
  full_join(fia_foc_sp, ., by = c('PLT_CN' = 'CN', 'PLOT'='PLOT')) %>% 
  
  mutate(foc_sp_PRES = ifelse(is.na(foc_sp_PRES),0, foc_sp_PRES)) %>%
  ungroup() %>% #ungroup by 'plot' and 'plot_cn' 
  group_by(foc_sp_PRES) %>% # group by presence and absences
  slice_sample(n = sum(fia_foc_sp$foc_sp_PRES)) %>% #determined tot. presence and randomly sampled an equivalent subset of absences 
  group_by(PLOT, PLT_CN) #regroup by 'plot' and 'plot_cn'... not sure if necessary


# Make this a spatial object
foc_sp_sf  <- st_as_sf(foc_sp_plot_coords, coords = c('LON','LAT'), crs = 4326) # coords must be in x, y order!
#saveRDS(foc_sp_sf, file = "foc_sp_sf.rds")


# Get US state boundaries to make maps more interpretable
state_bounds <- us_states(map_date = NULL, resolution = c("low"), states = c(fia_states))

missoula_coords<- data.frame(latitude = c(46.8721), longitude = c(-113.9940))
loc_missoula <- st_as_sf(missoula_coords, coords = c('longitude','latitude'), crs = 4326)
# Map... looks reasonable?

# Check against the USFS maps...link above
tm_shape(state_bounds, bbox = state_bounds) +
  tm_borders() +
  tm_fill(col = 'grey90') +
  tm_shape(foc_sp_sf) +
  tm_dots(col = 'grey30') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(col = 'red', border.col = 'black', size = 0.3, shape = 21, alpha = 0.7)+
  tm_shape(loc_missoula) + 
  tm_dots(size = 0.2, shape = 4, border.lwd = 3)

# Learn more about tmap: https://github.com/r-tmap/tmap
tm_shape(state_bounds, bbox = state_bounds) +
  tm_borders() +
  tm_fill(col = 'grey90') +
  tm_shape(foc_sp_sf) +
  tm_dots(col = 'foc_sp_PRES', 
          palette = c('grey30','red'))+
  #  tm_shape(ch)+
  #  tm_borders()+
  tm_shape(loc_missoula) + 
  tm_dots(size = 0.2, shape = 4, border.lwd = 3)


# ------------------------------------------------------------------------------
#  Extract and explore climate data for this species
# ------------------------------------------------------------------------------
# Use the `getdata` function in raster to download data from WorldClim
# WorldClim provides 19 'biologically meaningful' variables: https://www.worldclim.org/data/bioclim.html
bio_hist <- raster::getData('worldclim', var = 'bio', res = 2.5, path = '../data')

# Once downloaded, read in like this (I changed the name of the directory manually)
bio_hist <- list.files('../data/wc2-5/', 
                       pattern = '*.bil$', 
                       full.names = T) %>% 
  stack()

# This can be slow, you may decide to save this as a Rdata object for quick access later

foc_sp_bio <- terra::extract(bio_hist, foc_sp_sf)
#saveRDS(foc_sp_bio, file = "foc_sp_bio.rds")

# Select some climate variables... this is important and there are lots of ways to do it
# This is just a simplified example, picking a priori and checking for autocorrelation
# See function ggpairs for pairwise scatterplots... might be useful... could cbind presence/absence too

bio_order <- c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11',
               'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')
foc_sp_bio <- foc_sp_bio[, bio_order]
foc_sp_bio %>% 
  ggcorr(label = TRUE)


#create list of climate variables
clim.var <- list('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11',
                 'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')

foc_sp_hist_sf <- foc_sp_bio %>% 
  as.data.frame() %>% 
  dplyr::select(as.character(clim.var)) %>% 
  cbind(foc_sp_sf, .)  


# Map... looks reasonable?
climate_map <- 
  tm_shape(state_bounds, bbox = state_bounds) +
  tm_fill(col = 'grey60') +
  tm_shape(foc_sp_hist_sf) +
  tm_dots(col = 'black', 
          size = 0.09, 
          palette = 'BrBG',
          midpoint = NA, 
          style = 'cont') +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6)+ 
  tm_shape(loc_missoula) + 
  tm_dots(size = 0.2, shape = 4, border.lwd = 3)
climate_map

# ------------------------------------------------------------------------------
#  Fit a model
# ------------------------------------------------------------------------------

# Make a non-spatial dataframe 
foc_sp_hist_df <- foc_sp_hist_sf %>% 
  st_drop_geometry() %>% 
  mutate(presence = as.numeric(as.character(foc_sp_PRES)))

#saveRDS(foc_sp_hist_df, file = "foc_sp_hist_df.rds")

# split data into 80%(train) + 20%(valid)
dt <- sort(sample(nrow(foc_sp_hist_df), nrow(foc_sp_hist_df)*.8))
train<-foc_sp_hist_df[dt,]
valid<-foc_sp_hist_df[-dt,]

#full model formula

response <- 'presence'

all.var <- c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11',
             'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')

# =======================
# Boosted Regression Tree - FULL
# =======================


# Model requires tuning of the learning rate (lr) for each species to reach target of ~1000 trees. The lr's 
# for each species in the below ifelse statement were determined through trial and error.
lr <- ifelse(foc_sp == 'Juniperus osteosperma', 0.03,
             ifelse(foc_sp == 'Quercus gambelii', 0.06,
                    ifelse(foc_sp == 'Pinus edulis', 0.06, 'na'))) 



# FULL
gbm1.step <-gbm.step(data=train, gbm.x = all.var, gbm.y = response, 
                     family = "bernoulli", tree.complexity = 6,
                     learning.rate = lr, bag.fraction = 0.5) 

gbm1.step
summary(gbm1.step)

## look at partial dependency plots and relative importance

gbm.plot(gbm1.step, plot.layout= c(5,4), rug=TRUE, smooth = TRUE,common.scale=TRUE, write.title = FALSE)


##accuracy metrics 
pred <- predict(gbm1.step,valid,type="response")

mean(foc_sp_hist_df$presence, na.rm=T)
pred.nom <- ifelse(pred>0.05,1,0)

con.tab <- table(valid$presence,pred.nom)
con.tab

acc.metrics <- classAgreement(con.tab)
acc.metrics

pred.rocr <- prediction(pred,valid$presence)
perf <- performance(pred.rocr,"acc")
plot(perf)

perf <- performance(pred.rocr,"tpr","fpr")
plot(perf,colorize=T)

perf.auc <- performance(pred.rocr,"auc")
AUC<- slot(perf.auc,"y.values")
AUC #predictive skill

dt <- data.frame(cbind(foc_sp_hist_df$PLT_CN,foc_sp_hist_df$foc_sp_PRES,pred))
names(dt)<-c("plot","con_pa","pred")
dt$con_pa<-as.numeric(dt$con_pa)
dt$pred<-as.numeric(dt$pred)
acc.metrics2 <- optimal.thresholds(dt,opt.methods=c(3,4,6,7))


# ------------------------------------------------------------------------------
#  Project the model to current distro
# ------------------------------------------------------------------------------

# Predict historical across space
hist_preds <- bio_hist[[c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11',
                          'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')]]

# Add Montana state boundary to crop the future prediction
pred_bounds <- us_states(map_date = NULL, resolution = c("low"), states = c(fia_states))

hist_preds <- terra::crop(hist_preds, pred_bounds)

hist_space <- terra::predict(hist_preds, model = gbm1.step, type = 'response')

maxKappa<-acc.metrics2$pred[2] #probability threshold to make binary based on 'maxKappa'

hist_binary <- reclassify(hist_space, as.matrix(c(0,maxKappa,0,maxKappa,1,1)))

hist_con_map <- 
  tm_shape(hist_space, bbox = state_bounds, raster.downsample = FALSE) +
  tm_raster(title = 'Historical Climate Niche',
            palette = 'BuGn', 
            style = 'cont',
            breaks = seq(0, 1, by=.1)) +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.1, shape = 1, border.lwd = 1.6)+
  tm_shape(loc_missoula) + 
  tm_dots(size = 0.2, shape = 4, border.lwd = 3)


hist_bin_map <- 
  tm_shape(hist_binary, bbox = state_bounds, raster.downsample = FALSE) +
  tm_raster(title = 'Historical Climate Niche',
            palette = c('white','green'), 
            style = 'cat',
            breaks = seq(0, 1, by=.1)) +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.1, shape = 1, border.lwd = 1.6)+
  tm_shape(loc_missoula) + 
  tm_dots(size = 0.2, shape = 4, border.lwd = 3) 

hist_con_map
hist_bin_map

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
                             path = '/data')



# Once downloaded, read in like this (I changed the name of the directory manually)
bio_cmip5 <- list.files('/data/cmip5/2_5m/', 
                        pattern = '*.tif$', 
                        full.names = T) %>% 
  stack() %>% 
  # Rename these so they are the same as the model
  `names<-` (names(bio_hist))


foc_sp_cmip5 <- terra::extract(bio_cmip5, foc_sp_sf)

cmip5_preds <- bio_cmip5[[c('bio1','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','bio10','bio11',
                            'bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19')]]

#names(cmip5_preds) <- c('temp_min','temp_seas','precip_dry')
cmip5_preds <- terra::crop(cmip5_preds, pred_bounds)
cmip5_space <- terra::predict(cmip5_preds, model = gbm1.step, type = 'response')

cmip5_binary<- reclassify(cmip5_space, as.matrix(c(0,maxKappa,0,maxKappa,1,1)))

cmip5_cont_map <- 
  tm_shape(cmip5_space, bbox = state_bounds) +
  tm_raster(title = 'Future Climate Niche',
            #palette = c('white','green'), 
            #style = 'cat',
            palette = 'BuGn', 
            style = 'cont',
            breaks = seq(0, 1, by=.1)) +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6)+
  tm_shape(loc_missoula) + 
  tm_dots(size = 0.2, shape = 4, border.lwd = 3)

cmip5_bin_map <- 
  tm_shape(cmip5_binary, bbox = state_bounds) +
  tm_raster(title = 'Future Climate Niche',
            palette = c('white','green'), 
            style = 'cat',
            #palette = 'BuGn', 
            #style = 'cont',
            breaks = seq(0, 1, by=.1)) +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6)+
  tm_shape(loc_missoula) + 
  tm_dots(size = 0.2, shape = 4, border.lwd = 3)

cmip5_cont_map
cmip5_bin_map


# Arrow between historical and future to get a eyeball sense of differences

# Map difference between historical and future
sdm_difference_con <- cmip5_space - hist_space

diff_con_map <- 
  tm_shape(sdm_difference_con, bbox = pred_bounds) +
  tm_raster(title = 'Difference in probability',
            palette = "BrBG", 
            style = 'cont',
            breaks = seq(-1, 1, by=.2)) +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6)+
  tm_shape(loc_missoula) + 
  tm_dots(size = 0.2, shape = 4, border.lwd = 3)

sdm_difference_bin <- cmip5_binary - hist_binary

diff_bin_map <- 
  tm_shape(sdm_difference_bin, bbox = pred_bounds) +
  tm_raster(title = 'Difference in geog. suitability',
            palette = c('red','white','green'), 
            style = 'cat', 
            breaks = seq(-1, 1, by=.2)) +
  tm_shape(state_bounds) +
  tm_borders(col = 'black') +
  tm_shape(filter(foc_sp_sf, foc_sp_PRES == 1)) +
  tm_dots(size = 0.12, shape = 1, border.lwd = 1.6)+
  tm_shape(loc_missoula) + 
  tm_dots(size = 0.2, shape = 4, border.lwd = 3)

diff_con_map
diff_bin_map

setwd("C:/Users/sv181974/Box/Arboretum_project/figures")

writeRaster(hist_space, paste0(foc_sp, "_hist_space.tif"))
writeRaster(hist_binary, paste0(foc_sp, "_hist_binary.tif"))
writeRaster(cmip5_space, paste0(foc_sp, "_cmip5_space.tif"))
writeRaster(cmip5_binary, paste0(foc_sp, "_cmip5_binary.tif"))
writeRaster(sdm_difference_con, paste0(foc_sp, "_sdm_difference_con.tif"))
writeRaster(sdm_difference_bin, paste0(foc_sp, "_sdm_difference_bin.tif"))




