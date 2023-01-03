# used to work out for every number of offspring and sex ratio,
# what are the number of resulting new mating pairs
library(progress)

xmin <- 0
xmax <- 1000
ymin <- 0
ymax <- 1000

attraction_distance <- 1000

first.lot <- seq(2, 10, 2)
second.lot <- seq(20, 100, 10)
third.lot <- seq(200, 1000, 100)
fourth.lot <- seq(2000, 10000, 1000)
fifth.lot <- seq(20000, 100000, 10000)
sixth.lot <- seq(200000, 1000000, 100000)
seventh.lot <- seq(2000000, 10000000, 1000000)
offspring_numbers <- c(first.lot, second.lot, third.lot, fourth.lot, fifth.lot, sixth.lot, seventh.lot)

reps <- 1:100
percent.female <- 0.5

max.females.males.can.mate.with <- 30

output <- data.frame(
  N=rep(offspring_numbers, each=length(reps)),
  rep_num=rep(reps, times=length(unique(offspring_numbers))),
  NewCol=0
)

pb <- progress_bar$new(total=length(offspring_numbers)*length(reps),format = " [:bar] :percent eta: :eta")

for (n in offspring_numbers) {
  for (rep in reps) {
    # calculate random x and y positions
    x <- runif(n, xmin, xmax)
    y <- runif(n, ymin, ymax)
    sex <- runif(n, 0, 1)
    
    # seperate acccording to male female distribution
    x_male <- x[sex > percent.female]
    y_male <- y[sex > percent.female]
    
    x_female <- x[sex <= percent.female]
    y_female <- y[sex <= percent.female]
    
    if (attraction_distance==1000) {
      n.female.males.mate.with <- length(x_male) * max.females.males.can.mate.with
      n.females <- length(x_female)
      output$NewCol[output$rep_num==rep & output$N==n] <- min(n.female.males.mate.with, n.females)
    } else {
      stop('not yet implemented without an attraction_distance of 1000')
    }

    pb$tick()
  }
}

write.csv(output, 'NewCol_table.csv')