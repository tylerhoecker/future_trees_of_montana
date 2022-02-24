# ------------------------------------------------------------------------------
# This script extracts TerraClim data for a set of points. User species the 
# variables and the product. It's not particularly efficient... could be much
# faster if parallelized...
# ------------------------------------------------------------------------------
library(ncdf4)

# Insert this code after you've already created a spatial point dataframe 
# (must be sf for this to work as-is)
foc_sp_sf

# From the sf object, extract coordinates and transform into a list
point_sample <- foc_sp_sf %>% 
  st_transform(., crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  dplyr::select(Longitude = X, Latitude = Y) %>% 
  # Optionally, randomly sample a subset of these data - I don't recommend doing 
  # this completely randomly b/c you'll lose a bunch of presences. Could first
  # use `group_by(foc_sp_PRES)` and then `slice_sample` so you pick a subset of 
  # size n from every group (in this case, two groups: 0,1).
  # slice_sample(n = 10) %>% 
  # Split the dataframe into a list, to facilitate running function over list
  split(.,1:nrow(.))

# This function based on code from TerraClimate website: 
# http://www.climatologylab.org/uploads/2/2/1/3/22133936/read_terraclimate_point.r
terraclim_dl <- function(coords, variable, product){
  
  
  # enter in longitude, latitude here
  coords <- c(coords[1,"Longitude"], coords[1,"Latitude"])
  
  # Some progress info
  print(paste0('Running: ',product,'; ',variable,'; ',coords[1],',',coords[2]))
  
  # enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html
  # Updated this part to use climatologies with URL from here: http://thredds.northwestknowledge.net:8080/thredds/catalog/TERRACLIMATE_ALL/summaries/catalog.html
  # '#fillmismatch' to fix netCDF update that broke code
  base_url <- paste0('http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/summaries/')
  prod_url <- paste0(base_url, paste0(product,'_',variable,'.nc#fillmismatch'))
  nc <- nc_open(prod_url)
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  flat = match(abs(lat - coords[2]) < 1/48, 1)
  latindex = which(flat %in% 1)
  flon = match(abs(lon - coords[1]) < 1/48, 1)
  lonindex = which(flon %in% 1)
  start <- c(lonindex, latindex, 1)
  count <- c(1, 1, -1)
  
  # read in the full period of record using aggregated files
  data <- as.numeric(ncvar_get(nc, varid = variable, start = start, count))
  nc_close(nc)
  
  # ----------------------------------------------------------------------------
  # EDIT HERE (and in `result`) TO CALCULATE DERIVATIVES OF VARIABLE
  # ----------------------------------------------------------------------------
  # Calculate mean of `variable` for summer months only 
  summer_mean <- mean(data[5:7])
  
  # Put all info indo dataframe
  result <- data.frame('lon' = coords[1],
                       'lat' = coords[2],
                       'period' = product,
                       'variable' = variable,
                       'stat' = 'summer_mean',
                       'value' = summer_mean)
  return(result)
}

# ------------------------------------------------------------------------------
# EDIT HERE THE VARIABLES AND TIME PERIODS/PRODUCTS OF INTEREST
# ------------------------------------------------------------------------------
variables <- list('aet','def') 
products <- list('TerraClimate19812010','TerraClimate2C') 

# Run this set of nested functions over the 3 lists (variables, products, coordinates)
terra_clim_df <- 
  map_df(products, function(product){
    map_df(variables, function(variable){
      map_df(point_sample, function(coords){
        terraclim_dl(coords,variable,product)
      })
    })
  })

# Suggest saving these data in some format so you don't have to re-run!
saveRDS(terra_clim_df, 'terraclim.R')
