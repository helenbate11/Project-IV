# Code to implement the Boundary with Derivative and Runs Emulator on a 2d input 
# space with a single boundary or 2 perpendicular boundaries

# Note: calls on global functions of interest f() and gradf() (at the relevant boundary)

################################################################
# (Single-) Boundary with Derivative and Runs Emulator on a 2d input space
################################################################

library(pdist)

single_BwD_with_runs_BL_emulator <- function(xP,                               # emulator prediction points
                                             xD,                               # the run input locations
                                             D,                                # the outputs, D = (f(x^1), ..., f(x^n))
                                             boundary=list("x1"=0, "x2"=NULL), # boundary definition: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL)
                                             theta=1,                          # the correlation length vector (or scalar)
                                             sigma=1,                          # the prior standard deviation
                                             E_f=0,                            # the prior expectation of f(x)
                                             E_df=0                            # the prior expectation of df/dx(x) in the direction perpendicular to the boundary
                                             ){
  if(length(theta)==1){
    theta <- rep(theta, 2)
  }
  
  axis <- ifelse(is.null(boundary$x1), 2, 1)
  constant <- ifelse(is.null(boundary$x1), boundary$x2, boundary$x1)
  
  # Incorporating the boundary
  # projections of xP onto K
  E_fx <- rep(E_f, nrow(xP))
  E_dfx <- rep(E_df, nrow(xP))
  Var_fx <- rep(sigma^2, nrow(xP))
  
  xK <- xP
  xK[,axis] <- rep(constant, nrow(xP))
  aK <- xP[,axis] - rep(constant, nrow(xP))
  
  EK_fx <- E_fx + exp(-(aK/theta[axis])^2) * (f(xK) - E_fx) + 
    aK*exp(-(aK/theta[axis])^2) * (gradf(xK)[,axis] - E_dfx)
  VarK_fx <- Var_fx * (1 - exp(-2*(aK/theta[axis])^2) - 2*(aK/theta[axis])^2*exp(-2*(aK/theta[axis])^2))
  
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
    
    a1a2 <- outer(a1, a2, FUN="*")
    exp_a1_exp_a2 <- outer(exp(-(a1/theta[axis])^2), exp(-(a2/theta[axis])^2), FUN="*")
    
    return( Cov_fx_fxdash(dist_points1_points2_resc) - exp_a1_exp_a2 * Cov_fx_fxdash(dist_pointsK1_pointsK2_resc) 
            - (2/theta[axis]^2)*a1a2*exp_a1_exp_a2 * Cov_fx_fxdash(dist_pointsK1_pointsK2_resc) )
  }
  
  # Incorporating the runs xD
  E_f_D <- rep(E_f, nrow(xD))
  E_df_D <- rep(E_df, nrow(xD))
  
  xDK <- xD
  xDK[,axis] <- rep(constant, nrow(xDK))
  aD <- xD[,axis] - rep(constant, nrow(xD))
  EK_D <- E_f_D + exp(-(aD/theta[axis])^2) * (f(xDK) - E_f_D) +
    aD*exp(-(aD/theta[axis])^2) * (gradf(xDK)[,axis] - E_df_D)
  
  VarK_D <- CovK_fx_fxdash(xD)
  diag(VarK_D) <- diag(VarK_D) + delta
  
  VarK_D_inv <- chol2inv(chol(VarK_D))
  
  CovK_fx_D <- CovK_fx_fxdash(points1=xP, points2=xD)
  covK_fx_D_VarK_D_inv <- CovK_fx_D %*% VarK_D_inv
  
  EKD_fx <- EK_fx + covK_fx_D_VarK_D_inv %*% (D - EK_D)
  VarKD_fx <- VarK_fx - apply(covK_fx_D_VarK_D_inv * CovK_fx_D, 1, sum)
  
  return(cbind("ExpKD_f.x."=c(EKD_fx), 
               "VarKD_f.x."=c(VarKD_fx),
               "ExpK_f.x."=c(EK_fx),
               "VarK_f.x."=c(VarK_fx)))
}






################################################################
# (Perpendicular-) Boundary with Derivative and Runs Emulator on a 2d input space
################################################################

perp_BwDR_BL_emulator <- function(xP,                              # emulator prediction points
                                  xD,                              # the run input locations
                                  D,                               # the outputs, D = (f(x^1), ..., f(x^n))
                                  boundary=list("x1"=0.5, "x2"=0), # boundary definitions: x1 = boundary$x1 AND x2 = boundary$x2
                                  theta=1,                         # the correlation length vector (or scalar)
                                  sigma=1,                         # the prior standard deviation
                                  E_f=0,                           # the prior expectation of f(x)
                                  E_df=c(0,0,0)                    # the prior expectation of df/dx(x) in the direction perpendicular to each boundary and with respect to both (df/dx_1, df/dx_2, d^2f/dx_1dx_2)
                                  ){  
  if(length(theta)==1){
    theta <- rep(theta, 2)
  }
  
  constant <- c(boundary$x1, boundary$x2)
  
  # Incorporating the boundary
  # projections of xP onto K and L
  xK <- xP
  xL <- xP
  xKL <- xP
  xK[,1] <- rep(boundary$x1, nrow(xP))  # (x1=constant[1])
  xL[,2] <- rep(boundary$x2, nrow(xP))  # (x1=constant[2])
  xKL <- matrix(c(rep(boundary$x1, nrow(xP)), rep(boundary$x2, nrow(xP))), ncol=2) # projection onto x1, x2 corner
  a <- xP[,1] - rep(boundary$x1, nrow(xP))
  b <- xP[,2] - rep(boundary$x2, nrow(xP))
  
  E_fx <- rep(E_f, nrow(xP))
  E_dfx <- matrix(rep(E_df, nrow(xP)), ncol=3, byrow=TRUE)
  Var_fx <- rep(sigma^2, nrow(xP))
  
  EKL_fx <- E_fx + exp(-(a/theta[1])^2) * (f(xK) - E_fx) + a*exp(-(a/theta[1])^2) * (gradf(xK)[,1] - E_dfx[1]) +
    exp(-(b/theta[2])^2) * (f(xL) - E_fx) + b*exp(-(b/theta[2])^2) * (gradf(xL)[,2] - E_dfx[2]) -
    exp(-(a/theta[1])^2) * exp(-(b/theta[2])^2) * (f(xKL) - E_fx) -
    a * exp(-(a/theta[1])^2) * exp(-(b/theta[2])^2) * (gradf(xKL)[,1] - E_dfx[1]) -
    b * exp(-(a/theta[1])^2) * exp(-(b/theta[2])^2) * (gradf(xKL)[,2] - E_dfx[2]) -
    a * b * exp(-(a/theta[1])^2) * exp(-(b/theta[2])^2) * (gradf(xKL)[,3] - E_dfx[3])
  
  VarKL_fx <- Var_fx * ( 1 - exp(-2*(a/theta[1])^2) - 2*(a/theta[1])^2*exp(-2*(a/theta[1])^2) ) *
    ( 1 - exp(-2*(b/theta[2])^2) - 2*(b/theta[2])^2*exp(-2*(b/theta[2])^2) )
  
  # Updated covariance function
  Cov_fx_fxdash <- function(dist_matrix, d) {
    # Cov(f(x), f(x'))
    sigma^2 * exp(-(dist_matrix)^2)
  }
  
  CovKL_fx_fxdash <- function(points1, points2=NULL){
    
    a1 <- points1[,1] - rep(boundary$x1, nrow(points1))
    b1 <- points1[,2] - rep(boundary$x2, nrow(points1))
    
    pointsKL1 <- matrix(c(rep(boundary$x1, nrow(points1)), rep(boundary$x2, nrow(points1))), ncol=2)
    
    if(!is.null(points2)){
      
      a2 <- points2[,1] - rep(boundary$x1, nrow(points2))
      b2 <- points2[,2] - rep(boundary$x2, nrow(points2))
      
      pointsKL2 <- matrix(c(rep(boundary$x1, nrow(points2)), rep(boundary$x2, nrow(points2))), ncol=2)
    }else{
      a2 <- a1
      b2 <- b1
    }
    
    if(is.null(points2)){
      pointsKL1_resc <- t(t(pointsKL1)/theta)
      dist_pointsKL1_pointsKL2_resc <- as.matrix( dist(pointsKL1_resc) )
    }else{
      pointsKL1_resc <- t(t(pointsKL1)/theta)
      pointsKL2_resc <- t(t(pointsKL2)/theta)
      dist_pointsKL1_pointsKL2_resc <- as.matrix( pdist(pointsKL1_resc, pointsKL2_resc) )
    }
    
    a1_a2 <- outer(a1, a2, FUN="-")
    b1_b2 <- outer(b1, b2, FUN="-")
    
    a1a2 <- outer(a1, a2, FUN="*")
    b1b2 <- outer(b1, b2, FUN="*")
    
    exp_a1_exp_a2 <- outer(exp(-(a1/theta[1])^2), exp(-(a2/theta[1])^2), FUN="*")
    exp_b1_exp_b2 <- outer(exp(-(b1/theta[2])^2), exp(-(b2/theta[2])^2), FUN="*")
    
    return( (exp(-(a1_a2/theta[1])^2) - (1 + (2*a1a2)/(theta[1]^2))*exp_a1_exp_a2) * (exp(-(b1_b2/theta[2])^2) - (1 + (2*b1b2)/(theta[2]^2))*exp_b1_exp_b2) * Cov_fx_fxdash(dist_pointsKL1_pointsKL2_resc) )
  }
  
  # Incorporating the runs xD
  xDK <- xD
  xDL <- xD
  xDK[,1] <- rep(boundary$x1, nrow(xD))
  xDL[,2] <- rep(boundary$x2, nrow(xD))
  xDKL <- matrix(c(rep(boundary$x1, nrow(xD)), rep(boundary$x2, nrow(xD))), ncol=2)
  
  aD <- xD[,1] - rep(boundary$x1, nrow(xD))
  bD <- xD[,2] - rep(boundary$x2, nrow(xD))
  
  E_f_D <- rep(E_f, nrow(xD))
  E_df_D <- matrix(rep(E_df, nrow(xD)), ncol=3, byrow=TRUE)
  Var_f_D <- rep(sigma^2, nrow(xD))
  
  EKL_D <- E_f_D + exp(-(aD/theta[1])^2) * (f(xDK) - E_f_D) + aD*exp(-(aD/theta[1])^2) * (gradf(xDK)[,1] - E_df_D[1]) +
    exp(-(bD/theta[2])^2) * (f(xDL) - E_f_D) + bD*exp(-(bD/theta[2])^2) * (gradf(xDL)[,2] - E_df_D[2]) -
    exp(-(aD/theta[1])^2) * exp(-(bD/theta[2])^2) * (f(xDKL) - E_f_D) -
    aD * exp(-(aD/theta[1])^2) * exp(-(bD/theta[2])^2) * (gradf(xDKL)[,1] - E_df_D[1]) -
    bD * exp(-(aD/theta[1])^2) * exp(-(bD/theta[2])^2) * (gradf(xDKL)[,2] - E_df_D[2]) -
    aD * bD * exp(-(aD/theta[1])^2) * exp(-(bD/theta[2])^2) * (gradf(xDKL)[,3] - E_df_D[3])
  
  # VarKL_D <- Var_f_D * ( 1 - exp(-2*(aD/theta[1])^2) - 2*(aD/theta[1])^2*exp(-2*(aD/theta[1])^2) ) *
  #   ( 1 - exp(-2*(bD/theta[2])^2) - 2*(bD/theta[2])^2*exp(-2*(bD/theta[2])^2) )
  
  VarKL_D <- CovKL_fx_fxdash(xD)
  diag(VarKL_D) <- diag(VarKL_D) + delta
  
  VarKL_D_inv <- chol2inv(chol(VarKL_D))
  
  CovKL_fx_D <- CovKL_fx_fxdash(points1=xP, points2=xD)
  covKL_fx_D_VarKL_D_inv <- CovKL_fx_D %*% VarKL_D_inv
  
  EKLD_fx <- EKL_fx + covKL_fx_D_VarKL_D_inv %*% (D - EKL_D)
  VarKLD_fx <- VarKL_fx - apply(covKL_fx_D_VarKL_D_inv * CovKL_fx_D, 1, sum)
  
  return(cbind("ExpD_f.x."=c(EKLD_fx), "VarD_f.x."=VarKLD_fx))
}
