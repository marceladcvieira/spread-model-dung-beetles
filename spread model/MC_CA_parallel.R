#-----------------------------------------------------------------------------------------
# Name:         MC_CA_parallel.r
# Purpose:      MonteCarlo Cellular Automata Simulation test file. Stand-alone script.
# Author:       Adapted and extended by Jake Manger (based on original model developed by Francesco Tonini), 
# Email:        jakesmanger@gmail.com
# Created:      10/07/2020
# Copyright:    (c) 2020 by Jake Manger. See https://github.com/f-tonini/Termite-Dispersal-Simulation for original model (c) 2013 by Francesco Tonini.
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-3.6.1 64-bit version(http://www.r-project.org/)
#-----------------------------------------------------------------------------------------
# - added parallel processing
# - added human-assisted spread
# - added varying number of alates produced depending on colony age
# - for code simplicity, modified to work with only vector files .shp and .csv
# - fixed issues with simulation area
# - improved plotting (includes habitat block layer, human assisted spread and alate flights)
#   and plots termite biology growth/spread functions
# - added simulation parameter files
# - added precalculation of values to speed up processing time
# - made to work for R-3.6.1
#
# works correctly with: 
# - vector .shp habitat suitability layer in projected coordinates in meters
# - .csv X Y starting infestations layer in projected coordinates in meters
# - .csv X Y building locations in projected coordinates in meters
# - .csv NowCol_table.csv
#-----------------------------------------------------------------------------------------

# clear everything
rm(list=ls())

## SETUP
#setwd("D:/Marcela/ca-spread-model-dung-beetle/")

library(parallel)
library(raster)
library(rgdal)

# parameter_file <- readline(prompt="Enter parameter filename (e.g. parameters_test.R):\n")
parameter_file <- 'parameters_intermedius.R'

source(paste0('parameter_files/', parameter_file))
if ((doPlot || doPause) && parallel.proc) {stop('cannot plot and pause with parallel processing. Set parallel.proc to FALSE if you want to plot and pause.')}
source('myfunctionsMC_CA_parallel.R')


## set Path to folders in which you want to save all your vector & raster files
simulation_identifier = paste0(gsub(".R","", parameter_file),'_time_',gsub(":", '-', Sys.time()))
fOutput <- paste0('Outputs/Output_', simulation_identifier,'/')
# Create a physical copy of the subdirectory folder(s) where you save your output
# If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
dir.create(fOutput, showWarnings = FALSE)


## CREATE MODEL ENV
right <- left + grain * KernelSize_east
top <- bottom + grain * KernelSize_north
print('Extent of simulation area:')
print(paste('Left', left))
print(paste('Right', right))
print(paste('bottom', bottom))
print(paste('top', top))
#Set study area extent
ext <- list(Left = left, Right = right, Bottom = bottom, Top = top)				

#Read look-up table for number of new colonies
nCol_table <- read.table('./NewCol_table.csv', header = T, stringsAsFactors = F, sep=',')
nCol_table <- nCol_table[nCol_table$ratio == SexRatio, ]

## Read in carrying capacity
# read as a shape file (only for 0 and 1 data)
# habitat_suitability <- read.NSHabitat(ext, suitable_habitat, grain, proj_crs)
# read as a raster file (for data with the actual carrying capacity)
habitat_suitability <- raster(CCapacity)
habitat_suitability <- crop.or.extend.habitat.block(ext, habitat_suitability)
Carrying.capacity.lst <- sapply(rasterToPoints(habitat_suitability)[,3], FUN=list)

##Read in probability of survival layer
probability_of_adult_survival <- raster(ProbAdultSurvival)
probability_of_adult_survival <- crop.or.extend.habitat.block(ext, probability_of_adult_survival)
probability_of_adult_survival.lst <- sapply(rasterToPoints(probability_of_adult_survival)[,3], FUN=list)


if (doPlot) {
  plot(habitat_suitability, main = 'Habitat suitability')
  if (doPause) {
    readline(prompt="Here is the habitat suitability. Press [enter] to continue\n")
  }
  plot(probability_of_adult_survival, main = 'Probability of adult survival')
  if (doPause) {
    readline(prompt="Here is the probability of adult survival. Press [enter] to continue\n")
  }
}

## Load data to weight human assisted spread gravity function 
# prop <- generate.property.area(ext, filename = 'properties-coords-utm.csv', grain, min(KernelSize_east, KernelSize_north), proj_crs)
# 
# if (doPause) {
#   plot(prop$n.properties, main = 'Number of properties')
#   readline(prompt="Here are the Number of properties (to weight the human-assisted spread gravity function). Press [enter] to continue\n")
# }


## CREATE SPREAD FUNCTIONS
# natural spread
# make kernal of percentage of alates at distances from each point
alateProbKernel <- generate.alate.kernel(grain, AlateFlight_KernelSize, FUN = "exp", alpha = FlightDistance)   #use "exp" for exponential kernel and "gauss" for gaussian kernel 								 
# human assisted spread
# options(scipen=5)
# humanAssistedProbs <- generate.human.assist.probs(grain, min(KernelSize_east, KernelSize_north), prop$n.properties, psi, omega)


## PLOT GROWTH AND SPREAD FUNCTIONS
if (doPause) {
  plot(alateProbKernel, main = 'Probability of alate flight distance')
  readline(prompt="Here is the probability kernal for alate flights. Press [enter] to continue\n")
  
  temp <- as.data.frame(alateProbKernel, xy = T)
  plot(temp$y[temp$x==0 & temp$y>=0], temp$layer[temp$x==0 & temp$y>=0], main = 'Probability of offspring flight distance',
       xlab='Distance from mating pair', ylab='Proportion of offspring in single direction')
  readline(prompt="Here is the probability kernal in one direction. Press [enter] to continue\n")
  
  sequ <- seq(0, MaxAge)
  plot(sequ, alates.gen.precalculated[sequ+1], type='l', xlab='Colony age (time intervals)\n', ylab='Number of offspring that survive to be reproductive')
  readline(prompt="Here is the number of offspring with age. Press [enter] to continue\n")  
  
  sequ <- seq(0, (MaxAge+1))
  plot(sequ, colony.survival.precalculated[sequ+1], type='l', ylim=c(0, 1), xlab='Colony age', ylab='Survivability')
  readline(prompt="Here is the chance of survival of a reproductive female with age. Press [enter] to continue\n")
  
  # a <- vector()
  # for (i in 1:10000) {
  #   a <- c(a, max.age_fun(0:31))
  # }
  # plot(sort(table(a),decreasing=TRUE))
  # readline(prompt="Here is the most common survivability of colonies of different ages. Press [enter] to continue\n")
}


# for new introductions
time.points.of.introductions <- read.csv(starting_colony_coords)
time.points.of.introductions <- unique(time.points.of.introductions$release_time_int)

# define the main loop for the model
main <- function(run) {
  cat(paste('\n___________________\n','Run',run,'\n___________________\n'))
  ### CREATE DATABASES over time
  temp.dataset <- as.data.frame(matrix(0, ncol=length(seq(start_time,end_time)), nrow=1))
  colnames(temp.dataset) <- seq(start_time,end_time)
  ## General
  area.dataset <- temp.dataset
  n.colonies.dataset <- temp.dataset
  
  save.time_interval.stats <- function(out) {
    n.colonies.dataset[current_time-start_time+1] <<- out$num.colonies
    area.dataset[current_time-start_time+1] <<- out$Tot_area_km2
  }
  
  for (current_time in seq(start_time,end_time)){
    cat(paste('current_time',current_time,'\n___________________\n'))
    if (current_time == start_time) {
      ## let's initiate the colonies in the study area
      # we do this for each MC simulation run, so that ages are random between MC runs.
      # if there is no starting colony locations, can use n_inf_pixel in generate.Area
      # init <- generate.Area(ext, n_inf_pixel = 20, randomise.density=TRUE, grain=grain, KernelSize=AlateFlight_KernelSize, proj_crs=proj_crs)
      init <- generate.Area(ext, filename=starting_colony_coords, introduction.time.point=current_time, randomise.density=FALSE, grain=grain, KernelSize=AlateFlight_KernelSize, proj_crs=proj_crs)
      
      # create init (starting variables)
      simArea <- init$rr 	
      Age.lst <- init$age 
      Colonies.lst <- init$cl
      Alates <- simArea
      Alates[] <- 0
      
    } else {
      # Check for a new introduction event
      if (any(current_time==time.points.of.introductions)) {
        # if new introduction happened, generate a new area and starting variables and add those to the current variables
        starting_colonies <- read.csv(starting_colony_coords)
        starting_colonies <- as.data.frame(starting_colonies)
        starting_colonies <- starting_colonies[starting_colonies$release_time_int==current_time,]
        ##Rasterize points
        r <- simArea
        z <- rasterize(starting_colonies[,1:2], r, background=0, field=starting_colonies[,4], fun='sum')
        rm(r)
        
        # add new colonies to colonies list and age list
        out <- newCol.introduction(l1 = Age.lst, l2 = Colonies.lst, l3 = sapply(rasterToPoints(z)[,3], FUN=list), carrying.cap = Carrying.capacity.lst)		
        rm(z)

        Age.lst <- out$age
        Colonies.lst <- out$cl
        rm(out)
        
        simArea[] <- unlist(Colonies.lst)
        
        gc()
      }
      
      ### COLONY GROWTH
      # Colonies.lst = list of the number of colonies in each cell
      # Age.lst = list with the count of colony ages in each cell
      # if count.ages, e.g. 14 3 0 would be 14 colonies age 0, 3 age 1 and 0 age 2 in an individual cell. If not count.ages, e.g. 14 3 0 would be the ages of three colonies in an individual cell
      
      ##Increase age by 1 in those cells with at least one colony inside
      if (count.ages) {
        Age.lst <- lapply(Age.lst, FUN=function(x){ if (!is.null(x)) x <- c(0, x) })
        ##Remove colonies over MaxAge or that died randomly
        Age.lst <- lapply(seq_along(Age.lst), FUN=function(i) {max.age_fun(Age.lst[[i]], probability_of_adult_survival.lst[[i]])})
        # whenever a single element within age list has no more colonies inside
        # update the corresponding list element in colony list and set it to 0
        # also remove any colonies from the Colonies.lst that were over the max age (if they are no longer in Age.lst)
        Colonies.lst <- mapply(FUN=function(x,y){ if(is.null(x)) {y <- 0} else {y <- sum(x)}; return(list(y)) }, Age.lst, Colonies.lst) 
      } else {
        Age.lst <- lapply(Age.lst, FUN=function(x){ if (!is.null(x)) x <- x + 1 })
        ##Remove colonies over MaxAge or that died randomly
        Age.lst <- lapply(seq_along(Age.lst), FUN=function(i) {max.age_fun(Age.lst[[i]], probability_of_adult_survival.lst[[i]])})
        # whenever a single element within age list has no more colonies inside
        # update the corresponding list element in colony list and set it to 0
        # also remove any colonies from the Colonies.lst that were over the max age (if they are no longer in Age.lst)
        Colonies.lst <- mapply(FUN=function(x,y){ if(is.null(x)) {y <- 0} else {y <- length(x)}; return(list(y)) }, Age.lst, Colonies.lst) 
      }
      
      # if all colonies died
      if ( !any(unlist(Colonies.lst) > 0) ) {
        print(paste('No colonies survived at time', current_time))
        simArea[] <- unlist(Colonies.lst)
        
        out <- calculate.time_interval.stats(Colonies.lst, simArea, colAgeSwarmers, fOutput, run, current_time, grain)	
        save.time_interval.stats(out)
        
        ##EXIT the inner loop and go to the next MC simulation run
        break
      }
      
      ## Generate alates
      # Generate total number of offspring per cell by using the age of each existing colony
      Alates.lst <- lapply(Age.lst, FUN=alates.gen_fun)
      Alates.count <- lapply(Alates.lst, FUN=function(x){if(is.null(x))x <- 0 else sum(x)})
      # spread the alates from their colonies location according to the alateProbKernel
      Alates[] <- unlist(Alates.count)
      Alates <- round(focal(Alates, w=as.matrix(alateProbKernel), fun=sum, pad=TRUE, padValue=0))
      # remove alates from unsuitable habitat
      Alates[habitat_suitability[] == 0] <- 0
      # use look-up table now to decide how many new colonies based on alates number in each cell
      NewCol <- Alates
      NewCol[] <- sapply(Alates[], FUN=function(x) newCol.gen(x, tab=nCol_table))
      
      print(paste(sum(NewCol[]), 'new colonies generated via alate flight (ignoring overpopulation)'))
      if (doPlot) {
        plot(Alates, main = 'Number of Alates')
        if (doPause) {
          readline(prompt="Here are the number of alates in each cell. Press [enter] to continue\n")
          
        }
        plot(NewCol, main = 'New colonies from alate flight (ignoring overpopulation)')
        if (doPause) {
          readline(prompt="Here are the new colonies spread via alate flight (ignoring overpopulation). Press [enter] to continue\n")
        }
      }
      
      # # Generate human assisted neotenics
      # Neotenics <- generate.neotenics(Age.lst, prop$n.properties)
      # 
      # # add human-assisted dispersed neotenics with new colonies from alates to NewCol raster
      # NewCol[] <- NewCol[] + Neotenics
      
      NewCol.lst <- sapply(rasterToPoints(NewCol)[,3], FUN=list)	
      
      # add new colonies to colonies list and age list
      out <- newCol.addAge(l1 = Age.lst, l2 = Colonies.lst, l3 = NewCol.lst, carrying.cap = Carrying.capacity.lst)		
      Age.lst <- out$age
      Colonies.lst <- out$cl
      
      rm(NewCol.lst)
      
      simArea[] <- unlist(Colonies.lst)
      gc()
    }
    
    if (doPlot) {
      plot(habitat_suitability, legend = F)
      
      plot(simArea > 0, main = paste("Run ", run, "\n", current_time, sep=''), add = T, alpha=0.7)
      
      transp <- gray.colors(12, alpha = 0)
      
      if (sum(Alates[])>0) {
        micolor <- rev(gray.colors(12, alpha = 0.55))
        micolor[1] <- transp[1]
        plot(Alates, add = T, legend = F, col=micolor)
      }
      
 
      if (doPause) {
        readline(prompt="Plotted colonies, alates and habitat block. Press [enter] to continue\n")
      } 
    }
    
    # calculate area infested, economic impact
    # save raster files of this data and total sum
    # for each time and run
    out <- calculate.time_interval.stats(Colonies.lst, simArea, colAgeSwarmers, fOutput, run, current_time, grain)	
    save.time_interval.stats(out)
    gc()
    cat('\n___________________\n')
  }
  
  return(list(
    n.colonies.dataset=n.colonies.dataset,
    area.dataset=area.dataset
  ))
}

## NOW RUN THE MODEL!
if (parallel.proc) {
  # with parallel processing
  system.time({
    cl <- makeCluster(n.cpu.cores)
  
    ex <- ls(.GlobalEnv)
    clusterExport(cl, ex)
    clusterEvalQ(cl,
                 {
                   library(raster)
                   library(rgdal)
                   library(plyr)
                 })
  
    out_1 <- parLapply(cl, as.list(1:NRuns), main)
  })
  stopCluster(cl)
} else {
  # with a single cpu core (good for debugging)
  system.time({
    out_1 <- lapply(as.list(1:NRuns), main)
  })
}


## save all each times statistics
write.csv(do.call(rbind, lapply(out_1, '[[', 'n.colonies.dataset')),
          paste0(fOutput,'/n_colonies_dataset.csv'))
write.csv(do.call(rbind, lapply(out_1, '[[', 'area.dataset')),
          paste0(fOutput,'/area_dataset.csv'))


## create a summary of occupancy statistics across all MC runs 
summary.stats(do.call(rbind, lapply(out_1, '[[', 'area.dataset')))

## create and save averaged/sum rasters
envelope.raster('occupancy', calculate.percent.categories=TRUE)
envelope.raster('abundance', calculate.percent.categories=TRUE)

