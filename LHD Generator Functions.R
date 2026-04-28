# Code to randomly generate Latin Hypercube Designs on a 2d input space ( lhd_generator() ), 
# and find a better LHD design in terms of run spread across the space ( best_lhd() )

################################################################
# Latin Hypercube Design for Runs
################################################################

# generates 2d LHDs
lhd_generator <- function(n,              # number of points to place on 2d grid
                          x1_lim=c(0,1),  # limits of the x1 coordinate input space
                          x2_lim=c(0,1)   # limits of the x2 coordinate input space 
                          ){
  x_lhd <- cbind("x1"=sample(0:(n-1)), "x2"=sample(0:(n-1))) /n +0.5/n
  
  # Rescale x_lhd
  x_lhd[,"x1"] <- x1_lim[1] + x_lhd[,"x1"]*(x1_lim[2] - x1_lim[1])
  x_lhd[,"x2"] <- x2_lim[1] + x_lhd[,"x2"]*(x2_lim[2] - x2_lim[1])
  for(i in 1:1000){
    mat <- as.matrix(dist(x_lhd)) + diag(10, n)
    closest_runs <- which(mat==min(mat), arr.ind=TRUE)
    ind <- closest_runs[sample(nrow(closest_runs), 1), 1]
    swap_ind <- sample(setdiff(1:n, ind), 1) 
    x_lhd2 <- x_lhd 
    x_lhd2[ind[1], 1] <- x_lhd[swap_ind, 1]
    if(min(dist(x_lhd2)) >= min(dist(x_lhd))-0.00001){
      x_lhd <- x_lhd2
    }
    return(x_lhd)
  }
}

# Generates M 2d LHDs and finds that with the maximum minimum distance 
# between points
best_lhd <- function(n,                 # number of points to place on 2d grid
                     M=50,              # number of LHDs to test over for better maximum minimum point distance
                     x1_lim=c(0,1),     # limits of the x1 coordinate input space
                     x2_lim=c(0,1),     # limits of the x2 coordinate input space
                     print_switch=TRUE  # print the number of times a better LHD is found
                     ){
  x_lhd <- lhd_generator(n, x1_lim=x1_lim, x2_lim=x2_lim)
  switch_count <- 0
  for(m in 1:M){
    x_trial <- lhd_generator(n, x1_lim=x1_lim, x2_lim=x2_lim)
    if(min(dist(x_trial)) > min(dist(x_lhd))-0.00001){     # if min distance between points is same or better, replace with new LHD
      x_lhd <- x_trial
      switch_count <- switch_count + 1
    }
  }
  if(print_switch){
    cat("\nSwitch Count:", switch_count,
        "\nMinimum distance between points", min(dist(x_lhd)))
  }
  return(x_lhd)
}
