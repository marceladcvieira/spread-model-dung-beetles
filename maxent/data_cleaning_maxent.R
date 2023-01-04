#### Cleaning raw GBIF data ####

# install.packages("dplyr")
library(readr)
library("dplyr") 
library(tidyverse)
library(CoordinateCleaner)
library(tibble)

# set directory

setwd("C:/")  ## add file location

occurrence_alexis_gbif <- read_csv("occurrence_alexis_gbif.csv")
occurrence_binodis_gbif <- read_csv("occurrence_binodis_gbif.csv")
occurrence_taurus_gbif <- read_csv("occurrence_taurus_gbif.csv")
occurrence_intermedius_gbif <- read_csv("occurrence_intermedius_gbif.csv")

## select columns of interest and remove duplicates
!duplicated(occurrence_intermedius_gbif) #TRUE = unique, FALSE = duplicate

intermedius_gbif<- occurrence_intermedius_gbif %>%
  select(institutionCode,year,decimalLatitude,decimalLongitude,
         species, level0Name) %>% distinct(decimalLatitude,decimalLongitude, .keep_all = TRUE)%>% 
  arrange(desc(year))

pallipes_gbif<- occurrence_pallipes_gbif %>%
  select(institutionCode,year,decimalLatitude,decimalLongitude,
         species, level0Name)%>% distinct(decimalLatitude,decimalLongitude, .keep_all = TRUE)%>% 
  arrange(desc(year))

alexis_gbif<- occurrence_alexis_gbif %>%
  select(institutionCode,year,decimalLatitude,decimalLongitude,
         species, level0Name) %>% distinct(decimalLatitude,decimalLongitude, .keep_all = TRUE)%>% 
  arrange(desc(year))

binodis_gbif<- occurrence_binodis_gbif %>%
  select(institutionCode,year,decimalLatitude,decimalLongitude,
         species, level0Name)%>% distinct(decimalLatitude,decimalLongitude, .keep_all = TRUE)%>% 
  arrange(desc(year))

taurus_gbif<- occurrence_taurus_gbif %>%
  select(institutionCode,year,decimalLatitude,decimalLongitude,
         species, level0Name) %>% distinct(decimalLatitude,decimalLongitude, .keep_all = TRUE)%>% 
  arrange(desc(year))

# For more info see Penny Edwards 2007
# intermedius origin Africa
# imported from South Africa 
intermedius_gbif_origin<- intermedius_gbif %>%
  filter(level0Name=="Angola"|level0Name=="Botswana"|
           level0Name=="Democratic Republic of the Congo"|
           level0Name=="Ghana"|level0Name=="Kenya"|
           level0Name=="Malawi"|level0Name=="Mozambique"|
           level0Name=="Namibia"|level0Name=="Nigeria"|
           level0Name=="Rwanda"|level0Name=="South Africa"|
           level0Name=="Swaziland"|level0Name=="Tanzania"|
           level0Name=="Zambia"|level0Name=="Zimbabwe") %>% 
  rename(latitude = decimalLatitude,
         longitude = decimalLongitude) %>%
relocate(species,longitude,latitude,year,country=level0Name)

write.csv(intermedius_gbif_origin,'C:/.csv')  ## add path


# alexis origin Africa and soutern Europe
alexis_gbif_origin <- alexis_gbif %>%
  filter(level0Name=="Algeria"|level0Name=="Angola"|
           level0Name=="Botswana"|
           level0Name=="Central African Republic"|
           level0Name=="Democratic Republic of the Congo"|
           level0Name=="Egypt"|level0Name=="Ethiopia"|
           level0Name=="France"|level0Name=="Ghana"|
           level0Name=="Kenya"|level0Name=="Lesotho"|
           level0Name=="Malawi"|level0Name=="Mozambique"|
           level0Name=="Namibia"|level0Name=="Nigeria"|
           level0Name=="Rwanda"|level0Name=="Senegal"|
           level0Name=="South Africa"|level0Name=="Swaziland"|
           level0Name=="Tanzania"|level0Name=="Zambia"|
           level0Name=="Zimbabwe") %>% 
  rename(latitude = decimalLatitude,
         longitude = decimalLongitude)%>%
  relocate(species,longitude,latitude,year,country=level0Name)

write.csv(alexis_gbif_origin,'C:/.csv')  ## add path

# binodis origin South Africa
binodis_gbif_origin <- binodis_gbif %>%
  filter(level0Name=="Democratic Repblic of the Congo"|
           level0Name=="Lesotho"|level0Name=="Zimbabwe"|
           level0Name=="South Africa") %>% 
  rename(latitude = decimalLatitude,
         longitude = decimalLongitude)%>%
  relocate(species,longitude,latitude,year,country=level0Name)

write.csv(binodis_gbif_origin,'C:/.csv')  ## add path


# taurus origin Europe, North Africa, Middle East, Spain, Portugal, France, Italy, Greece, Morocco and Turkey
taurus_gbif_origin <-taurus_gbif %>%
  filter(level0Name=="Algeria"|level0Name=="Austria"|
           level0Name=="Croatia"|level0Name=="Czech Republic"|
           level0Name=="France"|level0Name=="Germany"|
           level0Name=="Greece"|level0Name=="Hungary"|
           level0Name=="Israel"|level0Name=="Italy"|
           level0Name=="Morocco"|level0Name=="Netherlands"|
           level0Name=="Poland"|level0Name=="Portugal"|
           level0Name=="Romania"|level0Name=="Spain"|
           level0Name=="Tunisia"|level0Name=="Turkey"|
           level0Name=="Turkmenistan") %>% 
  rename(latitude = decimalLatitude,
         longitude = decimalLongitude)%>%
  relocate(species,longitude,latitude,year,country=level0Name)

write.csv(taurus_gbif_origin,'C:/.csv')  ## add path


# Use package CoordinateCleaner to automatically remove spatial outliers so that the model fitting and selection process can be done on records that seem to be more likely
# DBs are getting moved around to help with dung and then there are times they may be introduced, a presence recorded, but then not persist? 
# rm(list=ls())

# save botht's possible to explore the contribution of those outliers to a final model too,

remove_outliers_save<-function(occ.path, species.name){

  occspecies<-read.csv(occ.path) [,c("species",	"longitude","latitude")]

  outliers<- cc_outl(
    occspecies,
    lon = "longitude",
    lat = "latitude",
    species = "species",
    method = "quantile",
    mltpl = 5,
    tdi = 1000,
    value = "clean",
    sampling_thresh = 0,
    verbose = TRUE,
    min_occs = 7,
    thinning = FALSE,
    thinning_res = 0.5
  )
  
  write.csv(outliers,paste0("C://", species.name,"_clean_pl_up_to_2020_outliers.csv"))   ## add path
  print(paste('Finished remove outliers for', species.name))
}

# binodis
remove_outliers_save (occ.path="C:/.csv", species.name = "binodis")  ## add path

# intermedius
remove_outliers_save (occ.path= "C:/.csv", species.name ="intermedius" )  ## add path

# alexis
remove_outliers_save (occ.path="C:/.csv" , species.name ="alexis" )  ## add path

# taurus
remove_outliers_save (occ.path="C:/.csv" , species.name ="taurus" )  ## add path


# Next step before adding data into Maxent is to subsample (delete points in the same pixel) and create buffer around occ points (cleaned, removed outliers and subsampled)

# creates data for maxent
# resample the presence locations 
# and crop the environmental data so it has a buffer around presence locations

library(raster)
library(sp)
library(rgdal)
library(dismo)
library(dplyr)
library(rgeos)
library(sf)


sp.file.name <- 'pallipes' # CHANGE species name 

save.buffer.points <- TRUE  # TRUE calculates buffer & and FALSE doesn't

# get predictor variables
load.climate.data <- function(dir.name) {
  pred_files <- list.files(dir.name, '\\.tif$', full.names=TRUE)
  asc_pred_files <- list.files(dir.name, '\\.asc$', full.names=TRUE)
  pred_files <- c(pred_files, asc_pred_files)
  predictors <- raster::stack(pred_files)
  
  return(predictors)
}
predictors <- load.climate.data("C:/") ## add file location

plot(predictors[[1]])

# load presence data
d <- read.csv(file = paste0("C://", sp.file.name, "_clean_pl_up_to_2020_outliers", ".csv"), stringsAsFactors = FALSE)  ## add path
head(d)

# subsample presence data
subsample.data <- function(d, reso, pred.raster.data, do.plot=FALSE) {
  if (!is.data.frame(d) || ncol(d) != 2) {
    stop('d should be a dataframe with two columns (Longitude first then Latitude).') 
  }
  if (length(names(pred.raster.data)) > 1) {
    stop('pred.raster.data should be a single raster file, not a raster stack.')
  }
  if (class(pred.raster.data) !='RasterLayer') {
    stop('pred.raster.data should be a RasterLayer')
  }
  r <- pred.raster.data
  print(paste('resolution of predictor raster data is',
              res(r)[1],
              'while the resolution of subsampled data will be',
              reso))
  if (round(reso, 4) < round(res(r)[1], 4)) {
    stop(paste('resolution of subsampled data is less than the resolution of predictor raster data.
               Note that Maxent will automatically subsample this further to the resolution of predictor raster data.'))
  }
  res(r) <- reso # degree
  out <- gridSample(d[c('longitude', 'latitude')], r, n=1)
  if (do.plot) {
    plot(r, xlim=c(-10,50), ylim=c(-40,40))
    points(d)
    points(out, cex=1, col='red', pch='x')
  }
  return(out)
}

species.name <- d[1,1]
d <- subsample.data(
  d[,c("longitude","latitude")],                   # 30 seconds -> reso =  0.0083333334 in decimal degrees
  reso=0.00833333333333333, # 2.5 arc minutes -> reso = 0.0416666667 in decimal degrees
  predictors[[1]],
  do.plot=TRUE
)

if (save.buffer.points) {
  # because the below functions require huge amounts of temporary storage, set the temp folder to somewhere with lots of storage
  # space
  rasterOptions(tmpdir="C:/temp") 
  
  # get the buffer mask
  buffer.size <- 2 # in units of raster data 2 = 2 degrees
  buffer.pnts <- d
  colnames(buffer.pnts) <- c('X', 'Y')
  coordinates(buffer.pnts) <- ~ X + Y
  crs(buffer.pnts) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  buff.mask <- gBuffer(buffer.pnts, width=buffer.size)
  
  # crop climate data using the buffer mask
  out <- raster::crop(predictors, raster::extent(buff.mask))
  out <- raster::mask(out, buff.mask)
  
  
  # plot to make sure it looks correct
  plot(out[[1]])
  
  # save stuff
  for (i in 1:length(names(out))) {
    layer.name <- names(out)[i]
    writeRaster(out[[i]], paste0("C:/", layer.name, "_buffer2d_", sp.file.name, ".asc")) ## add file location
  }
}

# add species name back to d
d <- cbind(species.name, d)
write.csv(d, paste0('C:/', sp.file.name, '.csv'))  ## add file location
