# combine carrying capacity data with maxent output to incorporate
# climate effects on dung beetle distribution

library(raster)
library(sf)
library(gdalUtils)

format_and_save_data <- function(maxent.path, carrying.cap.path, species.name) {
  # read in data
  maxent.output <- raster(maxent.path,
                          crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
  carrying.cap <- raster(carrying.cap.path)
  
  # convert maxent output to the correct crs
  maxent.output <- projectRaster(maxent.output, crs = crs(carrying.cap))
  
  # check crs of both
  crs(maxent.output)
  crs(carrying.cap)
  
  # convert maxent output to the correct spatial resolution
  raster.to.be.resampled <- maxent.output
  size.to.resample.to <- carrying.cap
  
  r_resam <- resample(raster.to.be.resampled, size.to.resample.to, method='bilinear')
  
  new.maxent <- r_resam
  # values(new.maxent)[values(new.maxent) < 0.1] <- 0
  
  # set any na values to 0
  values(carrying.cap)[is.na(values(carrying.cap))] <- 0
  values(new.maxent)[is.na(values(new.maxent))] <- 0
  
  writeRaster(carrying.cap, paste0(species.name,  '_carrying_capacity_wa.tif'), overwrite=TRUE)
  writeRaster(new.maxent,  paste0(species.name,  '_probability_of_survival_wa.tif'), overwrite=TRUE)
  
  par(mfrow=c(1, 2))
  plot(carrying.cap)
  plot(new.maxent)
  par(mfrow=c(1,1))
  
}

## Create maxent and carrying capacity outputs for each species

# without competition
# binodis
format_and_save_data(
  # UPDATE VALUES HERE \/
  maxent.path="Inputs120721/Onthophagus_binodis_Australia all bio clips_general.asc",
  carrying.cap.path="Inputs120721/carrying_capacity_o_binodis3_grid_5km.tif",
  species.name="o_binodis"
)

# intermedius
format_and_save_data(
  maxent.path="Inputs120721/Euoniticellus_intermedius_Australia all bio clips_general.asc",
  carrying.cap.path="Inputs120721/carrying_capacity_e_intermedius3_grid_5km.tif",
  species.name="e_intermedius"
)

# pallipes
format_and_save_data(
  maxent.path="Inputs120721/Euoniticellus_pallipes_Australia all bio clips_general.asc",
  carrying.cap.path="Inputs120721/carrying_capacity_e_intermedius3_grid_5km.tif",
  species.name="e_pallipes"
)

# taurus
format_and_save_data(
  maxent.path="Inputs120721/Onthophagus_taurus_Australia all bio clips_general.asc",
  carrying.cap.path="Inputs120721/carrying_capacity_o_binodis3_grid_5km.tif",
  species.name="o_taurus"
)

# alexis
format_and_save_data(
  maxent.path="Inputs120721/Onitis_alexis_Australia all bio clips_general_03.11.21.asc",
  carrying.cap.path="Inputs120721/carrying_capacity_o_alexis3_grid_5km.tif",
  species.name="o_alexis"
)


# lets plot these outputs to make sure they look right
par(mfrow=c(1, 3))
plot(raster('o_binodis_carrying_capacity_wa.tif'))
plot(raster('o_binodis_probability_of_survival_wa.tif'))
par(mfrow=c(1, 1))

par(mfrow=c(1, 3))
plot(raster('o_taurus_carrying_capacity_wa.tif'))
plot(raster('o_taurus_probability_of_survival_wa.tif'))
par(mfrow=c(1, 1))
