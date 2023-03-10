# maximum density of reproductive pairs in cell (carrying capacity per cell) is defined by the carrying capacity raster (carrying_capacity)
CCapacity <- 'carrying capacity/o_binodis_carrying_capacity_wa.tif'   # same as carrying_capacity_e.intermedius_per_gridcell_for_ca_model2.tif
ProbAdultSurvival <- 'prob survival/o_binodis_probability_of_survival_wa.tif'

# should have the location and number of starting reproductive females (not the number of introduced individuals)
starting_colony_coords <- 'release data/o_binodis_by_six_months.csv'

# doPlot for plotting while running simulation
# doPause for pausing and examining at each step 
# (also makes additional plots for debugging)
doPlot <- FALSE
doPause <- FALSE

# a memory optimisation - if the number of possible individuals in a cell is greater than
# the number of possible ages, then set count.ages to TRUE
count.ages = TRUE

# use parallel processing?
parallel.proc <- FALSE
# the number of cpu cores to use for running the model
# (default is set to detectCores, but it may be faster to use less if using few NRuns)
n.cpu.cores <- 5 # detectCores()

# note should be in meter projected coordinates, so change in pixels is the same length in meters
proj_crs <- '+init=EPSG:3577'

# set extent of study area (in the proj_crs coordinate reference system)
grain <- 5000 #cell resolution in METERS...each cell corresponds to a single mating pair
left = -1890000
bottom = -3920000
# left = -1492000
# bottom = -3734000


KernelSize_east <- 320 #n cells
KernelSize_north <- 492 #n cells

# KernelSize_east <- 500 # n cells
# KernelSize_north <- 500 # n cells

# KernelSize_east <- 250 # n cells
# KernelSize_north <- 250 # n cells

# size of cells used to calculate dispersal distances of flights
AlateFlight_KernelSize <- 41

# number of Monte Carlo runs
NRuns <- 1000
# first time interval of Simulation 
start_time <- 1
# last time interval of Simulation 
end_time <- 20
# we are using 1 month time intervals

# more conservative parameters
# the age at which colonies start producing first offspring (become mature)
# in the time intervals of the simulation
AgeOfMaturity <- 1

# the maximum life expectancy for a colony
# in the time intervals of the simulation (1 = dies at 2)
MaxAge <- 1

# probability that an offspring survives to find a mate
# (how many brood balls survive)
EggToAdultSurvival <- 0.5      ### James  Ridsdill-Smith Reproductive Behaviour of Insects page 272 O. binodis
                          

# male-female ratio of colony offspring 
SexRatio <- 0.5

# mean flight distance (meters)
FlightDistance <- 2000 * 0.9 #Fincher et al. 1983 12km/3years for D. gazella
                        # 2km in 6 months    
                        # using 1.8km/time step from calibration process

# max number of offspring produced by each mating pair in a time interval
MaxNOffspring <- 104  #  Ridsdill-Smith, Hall & Craig 1982

# max survival probability of a adult reproductive female
MaxSurvival <- 0.9 


