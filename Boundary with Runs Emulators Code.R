# (Separate) code to implement the Boundary with Runs Emulator on a 2d input 
# space with a single boundary, 2 perpendicular boundaries or 2 parallel 
# boundaries and

# (Complete) function that will identify relevant emulator (single- or 
# perpendicular-boundary) to use, given the boundary definitions

# Note: calls on global function of interest f()

################################################################
# (Single-) Boundary with Runs Emulator on a 2d input space
################################################################

library(pdist)

single_boundary_with_runs_BL_emulator <- function(xP,                               # emulator prediction points
                                                  xD,                               # the run input locations
                                                  D,                                # the outputs, D = (f(x^1), ..., f(x^n))
                                                  boundary=list("x1"=0, "x2"=NULL), # boundary definition: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL)
                                                  theta=1,                          # the correlation length vector (or scalar)
                                                  sigma=1,                          # the prior standard deviation
                                                  E_f=0                             # the prior expectation of f(x)
                                                  ){
  if(length(theta)==1){
    theta <- rep(theta, 2)
  }
  
  axis <- ifelse(is.null(boundary$x1), 2, 1)
  constant <- ifelse(is.null(boundary$x1), boundary$x2, boundary$x1)
  
  # Incorporating the boundary
  # projections of xP onto K
  E_fx <- rep(E_f, nrow(xP))
  Var_fx <- rep(sigma^2, nrow(xP))
  
  xK <- xP
  xK[,axis] <- rep(constant, nrow(xP))
  aK <- xP[,axis] - rep(constant, nrow(xP))
  # a = (difference between corresponding xP and xK points) NOT rescaled by theta
  
  EK_fx <- E_fx + exp(-(aK/theta[axis])^2) * (f(xK) - E_fx)
  VarK_fx <- Var_fx - sigma^2 * exp(-2*(aK/theta[axis])^2)
  
  # Updated covariance function
  Cov_fx_fxdash <- function(dist_matrix, d=axis) {
    # Cov(f(x), f(x'))
    sigma^2 * exp(-(dist_matrix)^2)
  }
  
  CovK_fx_fxdash <- function(points1, points2=NULL){
    pointsK1 <- points1
    pointsK1[,axis] <- rep(constant, nrow(points1))
    a1 <- points1[,axis] - rep(constant, nrow(points1))
    
    if(!is.null(points2)){
      pointsK2 <- points2
      pointsK2[,axis] <- rep(constant, nrow(points2))
      a2 <- points2[,axis] - rep(constant, nrow(points2))
    }else{
      a2 <- a1
    }
    
    if(is.null(points2)){
      points1_resc <- t(t(points1)/theta)
      pointsK1_resc <- t(t(pointsK1)/theta)
      dist_points1_points2_resc <- as.matrix( dist(points1_resc) )
      dist_pointsK1_pointsK2_resc <- as.matrix( dist(pointsK1_resc) )
    }else{
      points1_resc <- t(t(points1)/theta)
      points2_resc <- t(t(points2)/theta)
      pointsK1_resc <- t(t(pointsK1)/theta)
      pointsK2_resc <- t(t(pointsK2)/theta)
      dist_points1_points2_resc <- as.matrix( pdist(points1_resc, points2_resc) )
      dist_pointsK1_pointsK2_resc <- as.matrix( pdist(pointsK1_resc, pointsK2_resc) )
    }
    
    exp_a1_exp_a2 <- outer(exp(-(a1/theta[axis])^2), exp(-(a2/theta[axis])^2), FUN="*")
    
    return( Cov_fx_fxdash(dist_points1_points2_resc) - exp_a1_exp_a2 * Cov_fx_fxdash(dist_pointsK1_pointsK2_resc) )
  }
  
  # Incorporating the runs xD
  E_D <- rep(E_f, nrow(xD))
  
  xDK <- xD
  xDK[,axis] <- rep(constant, nrow(xDK))
  aD <- xD[,axis] - rep(constant, nrow(xD))
  EK_D <- E_D + exp(-(aD/theta[axis])^2) * (f(xDK) - E_D)
  
  VarK_D <- CovK_fx_fxdash(xD)
  diag(VarK_D) <- diag(VarK_D) + delta
  
  VarK_D_inv <- chol2inv(chol(VarK_D))
  
  CovK_fx_D <- CovK_fx_fxdash(points1=xP, points2=xD)
  covK_fx_D_VarK_D_inv <- CovK_fx_D %*% VarK_D_inv
  
  EKD_fx <- EK_fx + covK_fx_D_VarK_D_inv %*% (D - EK_D)
  VarKD_fx <- VarK_fx - apply(covK_fx_D_VarK_D_inv * CovK_fx_D, 1, sum)
  
  return(cbind("ExpD_f.x."=c(EKD_fx), "VarD_f.x."=VarKD_fx))
}





################################################################
# (Perpendicular-) Boundary with Runs Emulator on a 2d input space
################################################################

library(pdist)

perpendicular_boundary_with_runs_BL_emulator <- function(xP,                              # emulator prediction points
                                                         xD,                              # the run input locations
                                                         D,                               # the outputs, D = (f(x^1), ..., f(x^n))
                                                         boundary=list("x1"=0.5, "x2"=0), # boundary definitions: x1 = boundary$x1 AND x2 = boundary$x2
                                                         theta=1,                         # the correlation length vector
                                                         sigma=1,                         # the prior standard deviation
                                                         E_f=0                            # the prior expectation of f(x)
                                                         ){
  if(length(theta)==1){
    theta <- rep(theta, 2)
  }
  
  constant <- c(boundary$x1, boundary$x2)
  
  # Incorporating the boundary
  # projections of xP onto K
  E_fx <- rep(E_f, nrow(xP))
  Var_fx <- rep(sigma^2, nrow(xP))
  
  xK <- xP
  xK[,1] <- rep(constant[1], nrow(xP))
  aK <- xP[,1] - rep(constant[1], nrow(xP))
  # a = (difference between corresponding xP and xK points) NOT rescaled by theta
  
  xL <- xP
  xL[,2] <- rep(constant[2], nrow(xP))
  aL <- xP[,2] - rep(constant[2], nrow(xP))
  
  xKL <- matrix(rep(constant, times=nrow(xP)), ncol=2, byrow=TRUE)
  
  EKL_fx <- E_fx + exp(-(aK/theta[1])^2) * (f(xK) - E_fx) +
    exp(-(aL/theta[2])^2) * (f(xL) - E_fx) -
    exp(-(aK/theta[1])^2) * exp(-(aL/theta[2])^2) * (f(xKL) - E_fx)
  VarKL_fx <- sigma^2 * (1 - exp(-2*(aK/theta[1])^2)) * (1 - exp(-2*(aL/theta[2])^2))
  
  # Updated covariance function
  Cov_fx_fxdash <- function(dist_matrix, d) {
    # Cov(f(x), f(x'))
    sigma^2 * exp(-(dist_matrix)^2)
  }
  
  CovKL_fx_fxdash <- function(points1, points2=NULL){
    # pointsK1 <- points1
    # pointsK1[,axis] <- rep(constant, nrow(points1))
    aK1 <- points1[,1] - rep(constant[1], nrow(points1))
    aL1 <- points1[,2] - rep(constant[2], nrow(points1))
    
    pointsKL1 <- matrix(rep(constant, times=nrow(points1)), ncol=2, byrow=TRUE)
    
    if(!is.null(points2)){
      # pointsK2 <- points2
      # pointsK2[,axis] <- rep(constant, nrow(points2))
      aK2 <- points2[,1] - rep(constant[1], nrow(points2))
      aL2 <- points2[,2] - rep(constant[2], nrow(points2))
      
      pointsKL2 <- matrix(rep(constant, times=nrow(points2)), ncol=2, byrow=TRUE)
    }else{
      aK2 <- aK1
      aL2 <- aL1
    }
    
    if(is.null(points2)){
      pointsKL1_resc <- t(t(pointsKL1)/theta)
      # dist_points1_points2_resc <- as.matrix( dist(points1_resc) )
      dist_pointsKL1_pointsKL2_resc <- as.matrix( dist(pointsKL1_resc) )
    }else{
      # points1_resc <- t(t(points1)/theta)
      # points2_resc <- t(t(points2)/theta)
      pointsKL1_resc <- t(t(pointsKL1)/theta)
      pointsKL2_resc <- t(t(pointsKL2)/theta)
      # dist_points1_points2_resc <- as.matrix( pdist(points1_resc, points2_resc) )
      dist_pointsKL1_pointsKL2_resc <- as.matrix( pdist(pointsKL1_resc, pointsKL2_resc) )
    }
    
    aK1_aK2 <- outer(aK1, aK2, FUN="-")
    aL1_aL2 <- outer(aL1, aL2, FUN="-")
    
    exp_aK1_exp_aK2 <- outer(exp(-(aK1/theta[1])^2), exp(-(aK2/theta[1])^2), FUN="*")
    exp_aL1_exp_aL2 <- outer(exp(-(aL1/theta[2])^2), exp(-(aL2/theta[2])^2), FUN="*")
    
    return( (exp(-(aK1_aK2/theta[1])^2) - exp_aK1_exp_aK2) * (exp(-(aL1_aL2/theta[2])^2) - exp_aL1_exp_aL2) * Cov_fx_fxdash(dist_pointsKL1_pointsKL2_resc) )
  }
  
  # Incorporating the runs xD
  E_D <- rep(E_f, nrow(xD))
  
  xDK <- xD
  xDK[,1] <- rep(constant[1], nrow(xD))
  aDK <- xD[,1] - rep(constant[1], nrow(xD))
  xDL <- xD
  xDL[,2] <- rep(constant[2], nrow(xD))
  aDL <- xD[,2] - rep(constant[2], nrow(xD))
  xDKL <- matrix(rep(constant, times=nrow(xD)), ncol=2, byrow=TRUE)
  
  EKL_D <- E_D + exp(-(aDK/theta[1])^2) * (f(xDK) - E_D) +
    exp(-(aDL/theta[2])^2) * (f(xDL) - E_D) -
    exp(-(aDK/theta[1])^2) * exp(-(aDL/theta[2])^2) * (f(xDKL) - E_D)
  
  VarKL_D <- CovKL_fx_fxdash(xD)
  diag(VarKL_D) <- diag(VarKL_D) + delta
  
  VarKL_D_inv <- chol2inv(chol(VarKL_D))
  
  CovKL_fx_D <- CovKL_fx_fxdash(points1=xP, points2=xD)
  covKL_fx_D_VarKL_D_inv <- CovKL_fx_D %*% VarKL_D_inv
  
  EKLD_fx <- EKL_fx + covKL_fx_D_VarKL_D_inv %*% (D - EKL_D)
  VarKLD_fx <- VarKL_fx - apply(covKL_fx_D_VarKL_D_inv * CovKL_fx_D, 1, sum)
  
  return(cbind("ExpD_f.x."=c(EKLD_fx), "VarD_f.x."=VarKLD_fx))
}





################################################################
# (Complete) Boundary with Runs Emulator on a 2d input space
################################################################

boundary_with_runs_BL_emulator_fast <- function(xP,                            # emulator prediction points
                                                xD,                            # the run input locations
                                                D,                             # the outputs, D = (f(x^1), ..., f(x^n))
                                                boundary=list("x1"=0, "x2"=0), # boundary definition: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL)
                                                theta=1,                       # the correlation length vector (or scalar)
                                                sigma=1,                       # the prior standard deviation
                                                E_f=0                          # the prior expectation of f(x)
                                                ){
  
  if(!is.null(boundary$x1) & !is.null(boundary$x2)){ # x1=c(), x2=c()
    return(perpendicular_boundary_with_runs_BL_emulator(xP, xD=xD, D=D, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f))
  }else{
    return(single_boundary_with_runs_BL_emulator(xP, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f))
  }
}

