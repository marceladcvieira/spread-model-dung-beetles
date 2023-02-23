#-----------------------------------------------------------------------------------------
# Name:         myfunctions_CA_parallel.r
# Purpose:      Modules (functions) called by the main script
# Author:       Adapted and extended by Jake Manger (based on original model developed by Francesco Tonini), 
# Email:        jakesmanger@gmail.com
# Created:      10/07/2020
# Copyright:    (c) 2020 by Jake Manger. See https://github.com/f-tonini/Termite-Dispersal-Simulation for original model (c) 2013 by Francesco Tonini.
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-3.6.1 64-bit version(http://www.r-project.org/)
#-----------------------------------------------------------------------------------------

##HABITAT/BACKGROUND LAYER MODULE
##This module is used to read a background landscape/habitat layer
read.NSHabitat <- function(ext, layer, grain, proj_crs)
{
  cat('Reading unsuitable habitat layer...\n')
  ##Define all accepted file formats (OGR & GDAL)
  OGR_ext <- c('shp')
  
  ##Strip the extension from the file name
  file.ext <- unlist(strsplit(layer, '\\.'))[2]
  
  ##Strip the name from the file name
  file.name <- unlist(strsplit(basename(layer), '\\.'))[1]
  
  ##If the file extension is .shp read the vector file using
  ##the readOGR() function from the rgdal package
  ##Otherwise the routine assumes it is a raster file and uses
  ##the readGDAL() function
  if (any(file.ext %in% OGR_ext)) { 
    
    #For 'rgdal' library all supported OGR data formats are listed under ogrDrivers()
    habitat_block <- readOGR(dsn=dirname(layer), layer = file.name)
    
    ##Check if the background layer is georeferenced
    if (!is.na(proj4string(habitat_block))){
      
      ##If the coord. system of the layer is Geographic ("LAT-LON")
      ##Ask the user to define parameters of a Projected coord. system
      if(substring(proj4string(habitat_block),7,13) == 'longlat'){
        CRS <- proj_crs
        
        #To Change Projection and/or datum use the spTransform() function within package 'rgdal'
        habitat_block <- spTransform(habitat_block,CRS=CRS)
      }
      
    }else{
      CRS <- proj_crs
      
      proj4string(habitat_block) <- CRS
    }
    
    #get bounding box
    left <- habitat_block@bbox[1,1]
    right <- habitat_block@bbox[1,2]
    bottom <- habitat_block@bbox[2,1]
    top <- habitat_block@bbox[2,2]
    
    if (left > ext$Right || right < ext$Left) stop('Habitat_block layer extents are outside simulation')
    if (bottom > ext$Top || top < ext$Bottom) stop('Habitat_block layer extents are outside simulation')
    
    # round habitat layer extents to nearest simulation extents and grains
    # rounds up or down to make sure that the whole extent of the habitat block layer is included
    max.boarder <- 100 * grain
    left <- round.to.arbitrary.num(left, seq(ext$Left-max.boarder, ext$Right+max.boarder, grain), 'up')
    bottom <- round.to.arbitrary.num(bottom, seq(ext$Bottom-max.boarder, ext$Top+max.boarder, grain), 'up')
    right <- round.to.arbitrary.num(right, seq(ext$Left-max.boarder, ext$Right+max.boarder, grain), 'down')
    top <- round.to.arbitrary.num(top, seq(ext$Bottom-max.boarder, ext$Top+max.boarder, grain), 'down')
    
    #create an empty raster container
    r <- raster(xmn=left,xmx=right,ymn=bottom,ymx=top)
    
    h <- (r@extent@xmax - r@extent@xmin) / grain
    if (!(h%%1==0)) {
      new.xmx <- r@extent@xmin + round(h) * grain
      r@extent@xmax <- new.xmx
    }
    v <- (r@extent@ymax - r@extent@ymin) / grain
    if (!(v%%1==0)) {
      new.ymx <- r@extent@ymin + round(v) * grain
      r@extent@ymax <- new.ymx
    }
    res(r) <- grain
    habitat_block <- rasterize(habitat_block, r, field='DN')

    nas <- is.na(habitat_block[])
    if (any(nas)) {
      warning('NAs detected in habitat suitability file.
              Converting these to a suitability of 0.')
      habitat_block[nas] <- 0
    }
    
    plot(habitat_block)
    # update CRS as rasterize removes it
    CRS <- proj_crs
    proj4string(habitat_block) <- CRS
    
  }else{
    stop('not implemented!')
  }	
  
  # crop or extend habitat_block to simulation extent - is required as
  # if extent of habitat block is different to simulation,
  # then the wrong areas can be deemed habitable and uninhabitable
  habitat_block <- crop.or.extend.habitat.block(ext, habitat_block)
  
  cat('Done!\n')
  
  return(habitat_block)
}

round.to.arbitrary.num <- function(number, arbitrary.numbers, round.direction='up') {
  # e.g. usage round.to.arbitrary.num(5000.00000, seq(0, 8000, 155), 'up')
  if (round.direction=='up') {
    add <- 1
  } else if (round.direction=='down') {
    add <- -1
  } else {
    stop('round.direction should be up or down')
  }
  return(arbitrary.numbers[which.min(abs(arbitrary.numbers-number))+add])
}

crop.or.extend.habitat.block <- function(new_ext, old_habitat_block) {
  # crop or extend habitat_block to simulation extent - is required as
  # if extent of habitat block is different to simulation, then the wrong areas can be deemed habitable and uninhabitable
  new_habitat_block <- crop(old_habitat_block, extent(new_ext$Left, new_ext$Right, new_ext$Bottom, new_ext$Top))
  # if extent of simulation is larger than habitat block, then extend the extent
  new_habitat_block <- extend(new_habitat_block, extent(new_ext$Left, new_ext$Right, new_ext$Bottom, new_ext$Top), value=0)
  return(new_habitat_block)
}

generate.alate.kernel <- function(grain, KernelSize, xy = c(0,0), FUN, alpha){
  # Let's generate a Euclidean Distance Kernel (in an automatic manner)
  # and transform it into a Probability Kernel (by choosing a reference prob. distr.)
  # then choose an average distance (e.g. 200meters)
  
  if (KernelSize %% 2 == 0) {
    warning('You specified a moving window with an even number of cells! Size would be increased to nearest odd number of cells!')
    KernelSize <- KernelSize + 1
  }
  
  # create empty raster layer
  r <- raster(xmn= -(grain * KernelSize/2), xmx= grain * KernelSize/2, ymn= -(grain * KernelSize/2), ymx= grain * KernelSize/2)
  # set resolution of raster layer
  res(r) <- grain 
  # make a kernel with values as the distance from the center
  distKernel <- distanceFromPoints(r, xy)
  
  if (FUN == 'exp'){
    c <- 1
    # exponential function with distance
    distKernel[] <- ( c / (2 * alpha * gamma(1/c)) ) * exp( -abs((distKernel[])/alpha)^c  ) 
    # make each pixel a percentage of all pixels
    # each pixel is then the percentage of all alates during a flight
    distKernel[] <- distKernel[] / sum(distKernel[])
  }else if (FUN == 'gauss'){
    c <- 2
    distKernel[] <- ( c / (2 * alpha * gamma(1/c)) ) * exp( -abs(distKernel[]/alpha)^c  ) 
    distKernel[] <- distKernel[] / sum(distKernel[])
  }
  
  return(distKernel)
}


generate.Area <- function(ext, n_inf_pixel = NULL, filename = NULL, introduction.time.point, randomise.density=TRUE, grain, KernelSize, proj_crs){
  
  
  out <- list()
  if (ext$Right - ext$Left < grain * KernelSize | ext$Top - ext$Bottom < grain * KernelSize) stop("The extent of the chosen study area is less than the chosen distKernel size!\nPlease choose larger dimensions!")
  
  if (!is.null(filename)) {
    
    cat('\nReading Input File...\n')
    
    ##Define all accepted file formats
    Extensions <- c('txt','csv')
    
    ##Strip the extension from the file name
    file.ext <- unlist(strsplit(filename, '\\.'))[2]
    
    
    if (file.ext == 'txt') {
      starting_colonies <- as.matrix(read.delim(file = filename, header = TRUE, stringsAsFactors = FALSE)) #For TAB-delimited files
    }else if (file.ext == 'csv'){
      starting_colonies <- as.matrix(read.table(file = filename, header = TRUE, stringsAsFactors = FALSE, sep=',')) #For Comma-delimited files
    }

    # Grab the header row and change it to upper-case, as a default
    header <- colnames(starting_colonies)
    correct_labels <- c('X','Y')
    
    # Check for label inconsistencies and/or errors:
    
    if (header[1]!='X' | header[2]!='Y' | header[3]!='release_time_int' | header[4]!='qty') {
      stop('Column name mismatch in the starting colonies file.\n
           The starting colonies file MUST have 4 columns with X, Y, year and qty in this order.')
    }
    
    if (length(header) != 4) {
      stop('The starting colonies file should have 4 columns with X, Y, year and qty in this order.')
    }
    
    # # Remove all duplicate points from the uploaded file
    # if(any(duplicated(starting_colonies))) {
    #   
    #   exclude <- duplicated(starting_colonies)
    #   starting_colonies <- starting_colonies[-which(exclude),]
    #   cat(paste('\nRemoved', sum(exclude),'Duplicate Records!'))
    #   cat('\n')
    # }
    
    cat('Creating colonies dataset...\n')
    
    ##Convert the input matrix to a dataframe
    starting_colonies <- as.data.frame(starting_colonies)
    starting_colonies <- starting_colonies[starting_colonies$release_time_int==introduction.time.point,]
    ##Rasterize points
    r <- raster()
    extent(r) <- c(ext$Left, ext$Right, ext$Bottom, ext$Top)
    res(r) <- grain
    z <- rasterize(starting_colonies[,1:2], r, background=0, field=starting_colonies[,4], fun='sum')
    # z[z > maxdensity] <- maxdensity
    if (randomise.density) {
      stop('randomise.density is not yet implemented')
      # n_colony_areas <- length(z[z>0])
      # z[z>0] <- sample(1:maxdensity, n_colony_areas, replace=TRUE)
    }
    
    cat('Dataset created!\n')
    
  }else{
    
    z <- raster()
    extent(z) <- c(ext$Left, ext$Right, ext$Bottom, ext$Top)
    res(z) <- grain
    z[] <- 0
    
    a <- rasterToPoints(z)
    cond <- apply(a >= ((ext$Right - ext$Left) / 2) - 1000 & a <= ((ext$Right - ext$Left) / 2) + 1000, 1, sum)
    if (n_inf_pixel > length(cond)) stop('number of initial infested pixels cannot exceed available number of cells!') 
    ss <- sample(which(cond == 2), n_inf_pixel)
    if (!randomise.density) {
      stop('Randomise.density was set to FALSE, this needs to be TRUE if you are creating random initial starting colony data')
    }
    z[ss] <- sample(1:maxdensity, n_inf_pixel, replace=T)
  }
  
  crs(z) <- proj_crs
  
  out$cl <- sapply(rasterToPoints(z)[,3], FUN=list)
  # creates random starting age of colonies from 3
  # (an age that a colony is likely to be first detected) to max age
  # out$age <- lapply(out$cl, FUN=function(x){if(x > 0) sample(3:MaxAge, x, replace=T)})
  # creates random starting age of colonies from 0 to max age
  if (count.ages) {
    out$age <- lapply(out$cl, FUN=function(x){if(x > 0) rand.vect.with.total(0, MaxAge, x)})
  } else {
    out$age <- lapply(out$cl, FUN=function(x){if(x > 0) sample(0:MaxAge, x, replace=T)})
  }
  
  out$rr <- z
  
  return(out)
}

rand.vect.with.total <- function(min, max, total) {
  # generate random numbers
  x <- sample(min:max, total, replace=TRUE)
  # count numbers
  sum.x <- table(x)
  # convert count to index position
  out = vector()
  for (i in 1:length(min:max)) {
    out[i] <- sum.x[as.character(i)]
  }
  out[is.na(out)] <- 0
  return(out)
}

habitat.survival <- function(a, z){
  
  out <- list()
  w <- habitat_block[] == 1
  z[w] <- 0
  
  out$rr <- z	
  out$cl <- sapply(rasterToPoints(z)[,3], FUN=list)
  out$age <- mapply(FUN=function(x,y){if(x == 0) y <- NULL; return(list(y))}, out$cl, a)
  
  return(out)
  
}

# max.age_fun <- function(x){
# 	if (!is.null(x) & any(x > MaxAge)) {
# 	  # get any colonies over max age,
# 	  # then set to null
# 		x <- x[-which(x > MaxAge)]
# 		if(length(x) == 0) x <- NULL
# 	}else{x}
# 
# 	return(x)
# }

max.age_fun <- function(x, probability.of.survival.at.pixel.location){
  if (!is.null(x)) {
    if (count.ages) {
      # remove colonies over max age (the MaxAge+2 index)
      x <- x[1:(MaxAge+1)]
      # see if any colonies died from random causes.
      # if so, remove from x. if all colonies are dead, set
      # list item to NULL.
      # go through each age and calculate which female reproductives die
      for (i in 1:length(x)) {
        rand.death <- runif(x[i])
        survivability <- colony.survival.precalculated[i] * probability.of.survival.at.pixel.location
        colonies.died <- (rand.death >= survivability)
        x[i] <- x[i] - sum(colonies.died)
      }
      # if all ages have 0 reproductive females, set to null
      if (sum(x)==0) x <- NULL
    } else {
      # check for errors
      if (any(x < 0 | x > (MaxAge+1))) {
        stop(paste('x (colony age) is an impossible value:', x,'\n'))
      }
      
      # see if any colonies died from random causes.
      # if so, remove from x. if all colonies are dead, set
      # list item to NULL.
      # also note that colony survival at index max age + 1 is 0
      rand.death <- runif(length(x))
      survivability <- colony.survival.precalculated[x+1] * probability.of.survival.at.pixel.location
      colonies.died <- (rand.death >= survivability)
      x <- x[!colonies.died] 
      
      if(length(x) == 0) x <- NULL
    }
  }
  
  return(x)
}


colony.survival <- function(current_ages) {
  z <- rep(NA, length(current_ages))
  z_ind <- 1
  
  for (current_age in current_ages) {
    # check for errors
    if (current_age > (MaxAge + 1)) {
      stop('current_age can\'t be greater than MaxAge + 1')
    }
    if (current_age < 0) {
      stop('current_age can\'t be less than 0')
    }
    ###
    # fixed regardless of age
    if (current_age <= MaxAge) {
      z[z_ind] <- MaxSurvival
    } else {
      z[z_ind] <- 0
    }
    ###
    # depending on age
    # age_of_max_maturity <- round(AgeOfMaturity + (0.75 * (MaxAge-AgeOfMaturity)))
    # 
    # if (current_age < AgeOfMaturity) {
    #   # if not mature
    #   x <- c(-1, AgeOfMaturity)
    #   y <- c(0, MaxSurvival)
    # 
    #   z[z_ind] <- approx(x, y, xout=current_age)$y
    # } else if (current_age >= age_of_max_maturity) {
    #   # if at full maturity
    #   x <- c(age_of_max_maturity, (MaxAge+1))
    #   y <- c(MaxSurvival, 0)
    # 
    #   z[z_ind] <- approx(x, y, xout=current_age)$y
    # } else {
    #   # if mature, but not less than age_of_max_maturity
    #   z[z_ind] <- MaxSurvival
    # }
    z_ind <- z_ind + 1
  }
  return(z)
}
# precalculate for speed
colony.survival.precalculated <- colony.survival(0:(MaxAge+1))

# old francesco tonini function
# alates.gen.tonini <- function(x, scenario='optimistic'){
#   
#   if (scenario == 'optimistic'){
#     
#     z <- ifelse(x < AgeOfMaturity, 0, 100000 * EggToAdultSurvival)
#     z <- ifelse(x >= AgeOfMaturity & x < 10, 1000 * EggToAdultSurvival, z)
#     z <- ifelse(x >= 10 & x < 15, 10000 * EggToAdultSurvival, z)
#     
#   }else if (scenario == 'pessimistic'){
#     
#     z <- ifelse(x < AgeOfMaturity, 0, 100000 * EggToAdultSurvival)
#     z <- ifelse(x >= AgeOfMaturity & x < AgeOfMaturity * 2, 10000 * EggToAdultSurvival, z)
#     z <- ifelse(x >= AgeOfMaturity * 2 & x < AgeOfMaturity * 3, 50000 * EggToAdultSurvival, z)
#     
#   }
#   
#   return(z)
# }

alates.gen_fun <- function(x) {
  if (!is.null(x)) {
    if (count.ages) {
      return(alates.gen.precalculated[1:length(x)]*x)
    } else {
      # check for errors
      if (any(x < 0 | x > MaxAge)) {
        stop(paste('x (colony age) is an impossible value:', x,'\n'))
      }
      return(alates.gen.precalculated[x+1])
    }
  }
}

alates.gen <- function(current_ages) {
  # the number of offspring produced with mating pair age
  z <- rep(NA, length(current_ages))
  z_ind <- 1
  
  for (current_age in current_ages) {
    # check for errors
    if (current_age > MaxAge) {
      stop('current_age can\'t be greater than MaxAge')
    }
    
    if (current_age < AgeOfMaturity) {
      # if too young to produce alates
      z[z_ind] <- 0
    } else {
      # if can produce offspring
      z[z_ind] <- MaxNOffspring
    }
    
    z[z_ind] <- z[z_ind] * EggToAdultSurvival
    
    z_ind <- z_ind + 1
  }
  return(z)
}
# precalculate for speed
alates.gen.precalculated <- alates.gen(0:MaxAge)

newCol.gen <- function(x, tab){
  if (!is.na(x) & x > 1){
    idx <- sample(1:100, 1)
    int <- findInterval(x, tab[,1])
    sub <- tab[(int - 100 + 1):int,4] # is essentially doing this: nCol_table$NewCol[nCol_table$N==x]
    x <- sub[idx]
  }else{x <- 0}
  return(x)
}

newCol.introduction <- function(l1,l2,l3, carrying.cap){
  # l1 = Age.lst, l2 = Colonies.lst, l3 = NewCol.lst, carrying.cap = carrying capacity of each cell
  
  out <- list()
  for (i in 1:length(l1)){
    maxdensity <- carrying.cap[[i]]
    
    if(l3[[i]]==0) next
    if(l2[[i]] < maxdensity) {
      if (count.ages && is.null(l1[[i]])) {
        l1[[i]] <- rep(0, MaxAge+1)
      }
      if (l2[[i]] + l3[[i]] > maxdensity) {
        if (count.ages) {
          l1[[i]] <- l1[[i]] + rand.vect.with.total(0, MaxAge, maxdensity - l2[[i]])
          l2[[i]] <- sum(l1[[i]])
        } else {
          l1[[i]] <- append(l1[[i]], sample(0:MaxAge, maxdensity - l2[[i]], replace=T)) 
          l2[[i]] <- length(l1[[i]]) 
        }
      }else{
        if (count.ages) {
          l1[[i]] <- l1[[i]] + rand.vect.with.total(0, MaxAge, l3[[i]])
          l2[[i]] <- sum(l1[[i]])
        } else {
          l1[[i]] <- append(l1[[i]], sample(0:MaxAge, l3[[i]], replace=T))
          l2[[i]] <- length(l1[[i]]) 
        }
      }
    }
    
  }
  
  out$age <- l1
  out$cl <- l2
  
  return(out)
}

newCol.addAge <- function(l1,l2,l3, carrying.cap){
  # l1 = Age.lst, l2 = Colonies.lst, l3 = NewCol.lst, carrying.cap = carrying capacity of each cell
  
  out <- list()
  for (i in 1:length(l1)){
    maxdensity <- carrying.cap[[i]]
    if(l3[[i]] == 0) next
    if(l2[[i]] < maxdensity) {
      if (count.ages && is.null(l1[[i]])) {
        l1[[i]] <- rep(0, MaxAge+1)
      }
      if (l2[[i]] + l3[[i]] > maxdensity) {
        if (count.ages) {
          l1[[i]][1] <- l1[[i]][1] + (maxdensity - l2[[i]])
          l2[[i]] <- sum(l1[[i]]) 
        } else {
          l1[[i]] <- append(l1[[i]], rep(0, maxdensity - l2[[i]])) 
          l2[[i]] <- length(l1[[i]]) 
        }
      }else{
        if (count.ages) {
          l1[[i]][1] <- l1[[i]][1] + l3[[i]]
          l2[[i]] <- sum(l1[[i]]) 
        } else {
          l1[[i]] <- append(l1[[i]], rep(0, l3[[i]]))
          l2[[i]] <- length(l1[[i]]) 
        }
      }
    }
  }
  
  out$age <- l1
  out$cl <- l2
  
  return(out)
}

## HUMAN ASSISTED SPREAD MODULE
generate.human.assist.probs <- function(grain, KernelSize, n.properties, psi, omega) {
  if (KernelSize %% 2 == 0) {
    warning('You specified a moving window with an even number of cells! Size would be increased to nearest odd number of cells!')
    KernelSize <- KernelSize + 1
  }
  
  coords <- cbind(coordinates(n.properties), values(n.properties))
  
  gravs <- vector(mode = "list", length = nrow(coords))
  # if the coordinate has a property, generate the distance kernel
  # if the coordinate has no property, leave as NULL
  coords.with.properties <- which(coords[,3]>0)
  num.calculated <- 0
  for (i in coords.with.properties) {
    # calculate distance kernel between properties (distij)
    # calculate euclidean distances between j and all properties
    distances <- pointDistance(coords[i,1:2], coords[coords.with.properties,1:2], lonlat=FALSE)
    rj <- coords[i,3]
    ri <- coords[coords.with.properties,3]
    
    # omega <- 1.2
    # rerun doubling and halving omega 
    # to see what this does to the distance of spread.
    
    lam <- psi*(distances^(-omega))*rj*ri
    
    # set this (i) lam to 0
    lam[coords.with.properties==i] <- 0
    lam <- lam/sum(lam)
    # browser()
    
    # plot(distances, lam)
    
    # plot
    # all.lams <- rep(0, nrow(coords))
    # all.lams[coords.with.properties] <- lam
    # r <- n.properties
    # values(r) <- 0
    # values(r) <- all.lams
    # plot(r)
    
    
    gravs[[i]] <- lam
    
    num.calculated <- num.calculated + 1
    if (num.calculated %% 50 == 0) {
      print(paste('human-assisted spread calculated for', num.calculated, 'of',
                  length(coords.with.properties), 'cells with properties'))  
    }
  }
  
  return(gravs)
}

count.mature.cols <- function(x) {
  if (!is.null(x)) {
    if (count.ages) {
      return(sum(x[(AgeOfMaturity+1):length(x)]))
    } else {
      return(length(x[x>=AgeOfMaturity]))
    }
  } else {
    return(0)
  }
}

generate.neotenics <- function(ages, n.properties) {
  prop.values <- values(n.properties)
  # get the number of mature colonies in each cell
  n.mature.cols.vec <- rep(0, length(ages))
  n.mature.cols.vec <- sapply(ages, FUN=count.mature.cols)
  # remove mature colonies with no buildings
  n.mature.cols.vec[prop.values==0] <- 0
  # get index of mature colonies with buildings
  mat.indx <- which(n.mature.cols.vec > 0)
  # repeat index by the number of mature colonies
  mat.indx <- rep(mat.indx, times=n.mature.cols.vec[n.mature.cols.vec > 0])
  # get probabilities
  probs <- humanAssistedProbs[c(mat.indx)]
  # sum probabilites to get the total probability
  total.prob.with.houses <- Reduce(`+`, humanAssistedProbs[c(mat.indx)])
  total.prob <- rep(0, length(prop.values))
  total.prob[which(prop.values>0)] <- total.prob.with.houses
  # scale the total probability by the chances of a spread event occuring per colony per year
  total.prob <- total.prob * probHumanAssSpread
  # generate random spread based on these probabilities
  rand.spread <- runif(length(total.prob))
  new.colonies <- rep(0, length(total.prob))
  new.colonies[rand.spread <= total.prob & total.prob > 0] <- 1
  
  print(paste(length(new.colonies[new.colonies>0]), 'new colonies generated via human assistance (ignoring overpopulation)'))
  if (doPlot) {
    r <- n.properties
    values(r) <- total.prob
    plot(r, main = 'Probability of human-assisted spread')
    if (doPause) {
      readline(prompt="Here is the probability of human-assisted spread. Press [enter] to continue\n")
      
    }
    values(r) <- new.colonies
    plot(r, main = 'New colonies from human-assisted spread (ignoring overpopulation)')
    if (doPause) {
      readline(prompt="Here are the new colonies spread via human-assistance (ignoring overpopulation). Press [enter] to continue\n")
    }
  }
  return(new.colonies)
}

generate.property.area <- function(ext, filename = NULL, grain, KernelSize, proj_crs){
  # this is used to inform the # humans for the human assisted spread gravity function
  
  out <- list()
  if (ext$Right - ext$Left < grain * KernelSize | ext$Top - ext$Bottom < grain * KernelSize) stop("The extent of the chosen study area is less than the chosen distKernel size!\nPlease choose larger dimensions!")
  
  if (!is.null(filename)) {
    
    cat('\nReading Input File...\n')
    
    ##Define all accepted file formats
    Extensions <- c('txt','csv')
    
    ##Strip the extension from the file name
    file.ext <- unlist(strsplit(filename, '\\.'))[2]
    
    if (file.ext == 'txt') {
      properties <- as.matrix(read.delim(file = filename, header = TRUE, stringsAsFactors = FALSE)) #For TAB-delimited files
    }else if (file.ext == 'csv'){
      properties <- as.matrix(read.table(file = filename, header = TRUE, stringsAsFactors = FALSE, sep=',')) #For Comma-delimited files
    }
    
    ##Remove all duplicate points from the uploaded file
    if(any(duplicated(properties))) {
      
      exclude <- duplicated(properties)
      properties <- properties[-which(exclude),]
      cat(paste('\nRemoved', sum(exclude),'Duplicate Records!'))
      cat('\n')
    }
    
    cat('Creating properties dataset...\n')
    
    ##Convert the input matrix to a dataframe
    properties <- as.data.frame(properties)
    
    ##Rasterize points
    r <- raster()
    extent(r) <- c(ext$Left, ext$Right, ext$Bottom, ext$Top)
    res(r) <- grain
    # counts the number of properties in each grid cell
    n.properties <- rasterize(properties[,1:2], r, background=0, fun='count')
    # sums the area of properties in each grid cell
    coordinates(properties) <- ~X+Y
    # update CRS as rasterize removes it
    CRS <- proj_crs
    proj4string(n.properties) <- CRS
    
    cat('Dataset created!\n')
    
  } else {
    stop('not implemented!')
  }
  
  out$pl <- sapply(rasterToPoints(n.properties)[,3], FUN=list)
  
  out$n.properties <- n.properties
  
  return(out)
}


## time_interval INVASION STATISTICS MODULE
# this module calculates the statistics for the invasion and saves raster data
calculate.time_interval.stats <- function(Colonies.lst, simArea, colAgeSwarmers, fOutput, run, current_time, grain) {
  ### COLONY STATS
  num.colonies <- sum(unlist(Colonies.lst), na.rm=TRUE)
  print(paste('Number of colonies:', num.colonies))
  writeRaster(simArea, filename=paste('./',fOutput,'/abundance_Run',run,'_',current_time,'.img',sep=''), format='HFA', datatype='INT4S', overwrite=TRUE)

  # estimate the approx. area covered by current colonies in Km^2
  # all colonies (detectable via active and passive detection)
  Tot_area <- cellStats(simArea > 0, stat='sum') * grain * grain  #in meters^2
  Tot_area_km2 <- Tot_area/1000000  
  print(paste('Area covered by colonies (Km^2):', Tot_area_km2))
  writeRaster(simArea > 0, filename=paste('./',fOutput,'/occupancy_Run',run,'_',current_time,'.img',sep=''), format='HFA', datatype='INT1U', overwrite=TRUE)
  
  return(list(
    num.colonies=num.colonies,
    Tot_area_km2=Tot_area_km2
  ))
}


##SUMMARY STATISTICS MODULE
##This module is used to calculate required statistics on simulation results
summary.stats <- function(area.dataset)
{
  
  #Now let's have a look at the final stats
  simTimeIntervals <- seq(start_time,end_time)
  
  #Store the information about the area covered after each time step 
  Area_dataset <- data.frame(Time_interval = simTimeIntervals)
  
  SD <- sapply(area.dataset, sd)
  MEAN <- sapply(area.dataset, mean)
  Area_dataset$MEAN <- round(MEAN,2)
  Area_dataset$SD <- round(SD,2)
  Area_dataset$CV <- round(SD / MEAN, 2)
  
  if (nrow(Area_dataset) > 1) {
    
    Area_dataset$areaSpeed <- 0
    Area_dataset$percGrowth <- 0
    
    for (i in 2:nrow(Area_dataset)){
      
      Area_dataset$areaSpeed[i] <- round(MEAN[i] - MEAN[i-1], 2) #Km^2 / time interval				
      Growth <- (MEAN[i] - MEAN[i-1]) / MEAN[i-1]
      Area_dataset$percGrowth[i] <- round(Growth * 100)
      
      if (exists('n_inf_pixel') && n_inf_pixel == 1) {
        
        Area_dataset$EuclDistSpeed[i] <- round(MEAN2[i] - MEAN2[i-1], 2)  #Km / time interval	
        GrowthEucl <- (MEAN2[i] - MEAN2[i-1]) / MEAN2[i-1]
        if (!is.finite(GrowthEucl)) GrowthEucl <- 0
        Area_dataset$EuclDistGrowth[i] <- round(GrowthEucl * 100)
      }
      
    }
    
    colnames(Area_dataset) <- colnames(Area_dataset) <- c('Time intervals','Avg. Area (Km^2)', 'St.Dev. Area (Km^2)',
                                                          'CV Area', 'Avg. Km^2/time interval', 'Avg. Areal Growth (%)')
    var1 <- c('Time intervals', Area_dataset[,1])
    var2 <- c('Avg_Area_Km^2', Area_dataset[,2])
    var3 <- c('St.Dev_Area_Km^2', Area_dataset[,3])
    var4 <- c('CV_Area', Area_dataset[,4])
    var5 <- c('Avg.Km^2/time interval', Area_dataset[,5])
    var6 <- c('Avg.Areal_Growth(%)', Area_dataset[,6])
    
    myArray <- c(var1,var2,var3,var4,var5,var6)	
    
    
  }else{
    
    colnames(Area_dataset) <- c('Time intervals','Avg. Area (Km^2)', 'St.Dev. Area (Km^2)','CV Area')
    var1 <- c('Time intervals', Area_dataset[,1])
    var2 <- c('Avg_Area_Km^2', Area_dataset[,2])
    var3 <- c('St.Dev_Area_Km^2', Area_dataset[,3])
    var4 <- c('CV_Area', Area_dataset[,4])
    
    myArray <- c(var1,var2,var3,var4)
    
  }
  
  #This part is run only in the case of an invasion starting from ONE source
  #So that we can use the average euclidean distance from the center (invasion point)
  #when comparing the expasion rate of simulation Vs. theoretical uniform distribution
  ###if (n_inf_pixel == 1)###
  #Radial Increase general formula: [sqrt(A1)/pi - sqrt(A0)/pi ] / t1- t0	
  
  #Final message of simulation is over
  print('Final statistics saved to main folder!')
  
  write.table(Area_dataset, paste0(fOutput, '_Final_Stats.csv'), row.names=F, sep=',')
}

##OCCUPANCY ENVELOPE MODULE
##This module is used to compute occupancy envelopes across MonteCarlo simulation runs
##and save raster files accordingly
envelope.raster <- function(file.identifier) {
  # set file.identifier for different metrics
  # 'occupancy' is occupancy envelope
  # etc.
  
  for (yrs in start_time:end_time){
    RasDir <- paste0(getwd(), '/', fOutput)
    
    yrs <- 15
    
    #Store full path and name of all raster files corresponding to the current time intervals of the LOOP 
    etiq <- paste0(RasDir,file.identifier,'_Run',seq(NRuns),'_',yrs,'.img')

    # now find the raster files that are missing/empty (these should mean that the invader died out with 0s everywhere)
    n.empty <- table(file.exists(etiq))['FALSE']
    if (is.na(n.empty)) {n.empty <- 0}
    # remove these missing files from etiq
    etiq <- etiq[file.exists(etiq)]
    
    #Use the stack() function of the 'raster' package to stack all layers stored in 'etiq'
    #Note: our rasters only have 1 band...in case you have multispectral images you can add that
    Ras_stack <- stack(etiq)
    # create an empty stack if needed and add these 0s to the stack
    # to make sure that averages and sums are correct
    if (n.empty > 0) {
      empty <- raster(etiq[1])
      values(empty) <- 0
      names(empty) <- 'empty'
      for (i in 1:n.empty) {
        Ras_stack <- stack(Ras_stack, empty)  
      }
    }
    
    #Overlay all rasters of the stack and sum values for each pixel/cell
    Overlap <- overlay(Ras_stack, fun=function(x){sum(x)})
    writeRaster(Overlap,filename=paste0(RasDir,file.identifier,'_',yrs,'_sum','.img'),
                format='HFA', overwrite=TRUE)
    average <- overlay(Ras_stack, fun=function(x){mean(x)})
    writeRaster(average,filename=paste0(RasDir,file.identifier,'_',yrs,'_mean','.img'),
                format='HFA', overwrite=TRUE)

  ##Delete all single raster files (comment this line in case you want to keep all of them)
  #do.call(file.remove,list(list.files(RasDir, pattern=paste0(file.identifier,'_Run'),full.names=TRUE)))
  
  ##Now you have all raster files ready to be visualized in any GIS software, showing areas/pixels occupied
  ##In more than X% of the MonteCarlo simulation runs
  print(paste(file.identifier, 'raster has been saved to the Raster folder'))
  print('Done')
}
}