# ca-spread-model-dung-beetle
This is a cellular automata based spread model, implemented for predicting the spread of introduced dung beetles throughout Australia.

-The model has been tested and works correctly on R version 3.6.1. and R 4.0.2

# Main loop
The model has a main loop which consists of the following:
Add any new introduced individuals -> increase age of individuals by 1 -> remove any individuals greater than max age or that died randomly -> generate offspring and spread -> remove offspring in unlivable habitat (probability of survival of 0) -> select number of offspring to turn into reproductive females -> repeat


# Prerequisites
To get started, you must have R and RStudio installed.

# Setup
1. Clone this repository to a local folder on your computer
2. Open the rstudio project folder 'ca-spread-model-dung-beetle.Rproj'
3. Edit your parameter file in 'parameter_files/parameters_test.R', or whatever parameter file you want to use.
4. Open the main file that runs the model


There are three r scripts for the model:
-	MC_CA_parallel.R � The main file. This has the main loop of the model and logic for running the model
-	myfunctionsMC_CA_parallel.R � The myfunctions file. This has a lot of functions used by the main script above
-	parameters_test.R � The parameter file. This has all the parameters to change the model

There are four input files/folders for the model (all should be in UTM coordinates, have the same names and should be in the same format as the data in this example):
-	starting-colony-coords-utm.csv - has the locations of the starting colonies in the first year of the simulation
-	properties-coords-utm.csv - has the locations of buildings, which is used as the weight for the gravity function in the human-assisted dispersal module
-	NewCol_table.csv � has precalculated values for the number of new colonies to be created when there are X numbers of alates that have landed in a grid cell. You should change this to values that suit your species (e.g. by randomly sampling numbers from an appropriate distribution)
-	Habitat_block_layer � a folder containing a shapefile to designate �uninhabitable areas�. This should be a multipolygon in UTM coordinates

To run the model, simply install the dependencies with 
install.packages(c(�parallel�, �raster�, �rgdal�)) 

and then Source the MC_CA_parallel.R file.

I would recommend that you first run the model with plotting and pausing, so you can see how the model works.
To do this:
change doPlot and doPause variables to TRUE in the parameter file
and change parallel.proc to FALSE in the parameter file.
Change doPlot and doPause back to FALSE when not testing a model and change parallel.proc to TRUE (parallel processing) if you want the model to run a lot faster. Note that parallel processing does not work when doPlot or doPause are TRUE.

To change the model, alter the parameters in the parameter file, or change the inputs to the model. Making additional changes will have to involve adding code to the main script of the model or editing the model�s functions in the myfunctions file. For example, if you needed to change the dispersal distance of the species, you could change the meanFlightDist variable in the parameter file. If, however, you do not want an exponential spread function, you would have to change the code in the generate.alate.kernel() function in the myfunctions file.

The model outputs raster files and csv files with the results of the simulation. There should be summary raster files (mean, sum and % categories, i.e. present in >0% of simulations, present in >25% of simulations, present in >50% of simulations, present in >75% of simulations and present in 100% of simulations) with locations of colonies for each year of the simulation across the number of montecarlo runs, as well as the locations of colonies in each run. You can drag these into QGIS or similar software and use these for analysis. But just a warning, these files can get quite big, so keep track of your storage space.
