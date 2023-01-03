# maximum density of reproductive pairs in cell (carrying capacity per cell) is defined by the carrying capacity raster (carrying_capacity)
carrying_capacity_location <- 'o_alexis_carrying_capacity_wa.tif'  
probability_of_adult_survival_location <- 'o_alexis_probability_of_survival_wa.tif'

# should have the location and number of starting reproductive females (not the number of introduced individuals)
starting_colony_coords <- 'Inputs/o_alexis_by_six_months.csv'

# doPlot for plotting while running simulation
# doPause for pausing and examining at each step 
# (also makes additional plots for debugging)
doPlot <- FALSE
doPause <- FALSE

# a memory optimisation - if the number of possible individuals in a cell is greater than
# the number of possible ages, then set count.ages to TRUE
count.ages = TRUE

# use parallel processing?
parallel.proc <- TRUE
# the number of cpu cores to use for running the model
# (default is set to detectCores, but it may be faster to use less if using few NRuns)
n.cpu.cores <- 5 # detectCores()

# note should be in meter projected coordinates, so change in pixels is the same length in meters
proj_crs <- '+init=EPSG:3577'

# set extent of study area (in the proj_crs coordinate reference system)
grain <- 5000 #cell resolution in METERS...each cell corresponds to a single mating pair
left = -1890000


bottom = -3920000


KernelSize_east <- 320 #n cells
KernelSize_north <- 492 #n cells


# size of cells used to calculate dispersal distances of flights
AlateFlight_KernelSize <- 41

# AlateFlight_KernelSize <- 201

# number of Monte Carlo runs
NRuns <- 1000
# NRuns <- 100
# first time interval of Simulation 
start_time <- 1
# last time interval of Simulation 
end_time <- 20
# we are using 6 months time intervals

# more conservative parameters
# the age at which colonies start producing first offspring (become mature)
# in the time intervals of the simulation
AgeOfMaturity <-  ## add parameter value

# the maximum life expectancy for a colony
# in the time intervals of the simulation (1 = dies at 2)
MaxAge <-     ## add parameter value

# probability that an offspring survives to find a mate
# (how many brood balls survive)
alateSurvival <-    ## add parameter value
                         

# male-female ratio of colony offspring 
SexRatio <-  ## add parameter value

# mean flight distance (meters)
meanFlightDist <-   ## add parameter value

# max number of offspring produced by each mating pair in a time interval
MaxNAlates <-  ## add parameter value

# max survival probability of a adult reproductive female
maxSurvival <-  ## add parameter value

