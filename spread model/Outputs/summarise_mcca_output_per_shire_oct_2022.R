library(sf)
library(terra)
library(tidyr)
library(progress)
library(stars)
library(exactextractr)

#' Summarises continuous rasters by polygons
#' 
#' @param raster Raster* object with data to be extracted
#' @param polygon sf object with polygons
#' @param summary_functions vector of Characters specifying function, or supply a custom function
#' ("min"|"max"|"count"|"sum"|
#'  "mean"|"median"|"quantile"|"mode"|"majority"|
#'  "minority"|"variety"|"variance|stdev|
#'  coefficient_of_variation|weighted_mean|weighted_sum")
#'  See https://cran.r-project.org/web/packages/exactextractr/exactextractr.pdf
#'  for more options
#' @export
summarise_continuous_raster_in_polygon <- function(
  r,
  polygon,
  summary_functions
) {
  # remove empty polygon geometries
  polygon <- polygon[!st_is_empty(polygon), ]

  # convert polygon to the crs of the raster file if necessary
  poly_crs <- st_crs(polygon)
  r_crs <- st_crs(r)
  if (poly_crs != r_crs) {
    # print(
    #   paste(
    #     'Converting polygon crs (',
    #     poly_crs$epsg,
    #     ') to crs of raster file (',
    #     r_crs$epsg,
    #     ') ...'
    #   )
    # )
    conv_polygon <- st_transform(polygon, r_crs)
  } else {
    conv_polygon <- polygon
  }
  
  # extract and sum raster values per polygon
  # see https://isciences.gitlab.io/exactextractr/articles/vig1_population.html
  # print('Extracting values from raster file... \n')
  extracted <-
    exact_extract(x = r,
                  y = conv_polygon,
                  fun = summary_functions,
                  progress=FALSE)

  return(extracted)
}

# shapefile with polygons of shires
polygon <- st_read('Reprojected_census1981_WA.gpkg')
polygon <- polygon[!st_is_empty(polygon),]

# build the table we want to create
output.dir <- c(
  'D:/Marcela/ca-spread-model-dung-beetle/Outputs/Output_parameters_alexis_calibration_flight_minus10_and_updated_maxnalates_Oct2022_time_2022-10-17 19-59-25/',
  'D:/Marcela/ca-spread-model-dung-beetle/Outputs/Output_parameters_taurus_calibration_flight_minus10_updated_maxnalates_time_2022-08-03 16-30-59/',
  'D:/Marcela/ca-spread-model-dung-beetle/Outputs/Output_parameters_intermedius_calibration_flight_minus10_time_2022-07-13 08-22-31/',
  'D:/Marcela/ca-spread-model-dung-beetle/Outputs/Output_parameters_binodis_calibration_flight_minus10_time_2022-07-04 07-15-40/'
)
time.steps <- 1:20
mc.runs <- 1:1000   ## before 1:50 ## before 1: 1000
LGA.codes <- polygon$LGA_CODE_1981
output <- crossing(output.dir, time.steps, mc.runs, LGA.codes)
output$abundance <- 0

pb <- progress_bar$new(
  format = " Summarising mcca output [:bar] :percent eta: :eta",
  total=nrow(output), clear=FALSE
)

for (od in output.dir) {
  for (t in time.steps) {
    for (run in mc.runs) {
      r.file.path <- paste0(od, 'abundance_Run', run, '_', t, '.img')
      if (file.exists(r.file.path)) {
        # load the raster file with information we want to extract
        r <- rast(r.file.path)
        
        
        #####
        # new way
        sum <- summarise_continuous_raster_in_polygon(
          r,
          polygon,
          summary_functions=c('sum')
        )

        if (!all(output[index,]$LGA.codes == polygon$LGA_CODE_1981)) {
          stop('LGA codes are mixed up. The below assignment will be incorrect')
        }        
          index <- (
            output$output.dir==od
            & output$time.steps==t
            & output$mc.runs==run
          )
          output$abundance[index] <- sum
        
          pb$tick(length(sum))
          
        #####
        
        # old way
        # loop through each shire and extract information from the raster file: r
        # for (LGA.code in polygon$LGA_CODE_1981) {
        #   shire <- polygon[polygon$LGA_CODE_1981==LGA.code,]
        # 
        #   values <- terra::extract(r, vect(shire))
        #   total <- sum(values$Layer_1)
        # 
        #   index <- (
        #     output$output.dir==od
        #     & output$time.steps==t
        #     & output$mc.runs==run
        #     & output$LGA.codes==LGA.code
        #   )
        #   output$abundance[index] <- total
        # 
        #   # plot(r)
        #   # plot(shire$geom, add=T)
        #   # Sys.sleep(0.1)
        # 
        #   pb$tick()
        # }
          
        #####
      }
    }
  }
}

write.csv(output, 'summary_mcca_output_May_2022_sensitivity_intermedius.csv')
