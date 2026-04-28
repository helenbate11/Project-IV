# Code to implement the Simple Emulator on 1d and 2d input spaces

################################################################
# Simple Emulator on a 1d input space
################################################################

# Simple Emulator for a single input x on a 1d input space
simple_BL_emulator_1d_single <- function(x,        # single emulator prediction point
                                         xD,       # the run input locations
                                         D,        # the outputs, D = (f(x^1), ..., f(x^n))
                                         theta=1,  # the correlation length
                                         sigma=1,  # the prior standard deviation, sigma = sqrt(Var[f(x)])
                                         E_f=0     # the prior expectation of f
                                         ){
  n <- length(xD)
  
  Cov_fx_fxdash <- function(x, xdash) {
    sigma^2 * exp(-(x-xdash)^2 / (theta^2))
  }
  
  E_f <- E_f
  E_D <- rep(E_f, n)
  E_D
  
  Var_D <- matrix(0, nrow=n, ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      Var_D[i,j] <- Cov_fx_fxdash(xD[i], xD[j])
    }
  }
  diag(Var_D) <- diag(Var_D) + delta
  
  E_fx <- E_f
  Var_fx <- sigma^2
  Cov_fx_D <- matrix(0, nrow=1, ncol=n)
  for(j in 1:n){
    Cov_fx_D[1,j] <- Cov_fx_fxdash(x, xD[j])
  }
  
  # Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x)
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% t(Cov_fx_D)
  
  # Return emulator expectation and variance
  return(c("ExpD_f(x)"=ED_fx, "VarD_f(x)"=VarD_fx))
}

# Simple Emulator for on a 1d input space
simple_BL_emulator_1d <- function(xP,       # emulator prediction points
                                  xD,       # the run input locations
                                  D,        # the outputs, D = (f(x^1), ..., f(x^n))
                                  theta=1,  # the correlation length
                                  sigma=1,  # the prior standard deviation, sigma = sqrt(Var[f(x)])
                                  E_f=0     # the prior expectation of f
                                  ){
  # nP <- length(xP)
  # em_out <- matrix(0, nrow=nP, ncol=2, dimnames=list(NULL, c("ExpD_f(x)", "VarD_f(x)")))
  em_out <- t(sapply(xP, simple_BL_emulator_1d_single, xD=xD, D=D, theta=theta, sigma=sigma, 
                     E_f=E_f))
  colnames(em_out) <- c("ExpD_f(x)", "VarD_f(x)")
  return(em_out)
}






################################################################
# Simple Emulator on a 2d input space
################################################################
library(pdist)

simple_BL_emulator_2d_fast <- function(xP,       # emulator prediction points
                                       xD,       # the run input locations
                                       D,        # the run outputs, D
                                       theta=1,  # the correlation length vector
                                       sigma=1,  # the prior standard deviation, sigma
                                       E_f=0     # the prior expectation of f(x)
                                       ){
  nD <- length(D)
  nP <- nrow(xP)
  
  # Rescale each input by theta (works both for different and similar theta)
  xP <- t(t(xP)/theta)
  xD <- t(t(xD)/theta)
  
  # Define covariance structure of f(x): Cov[f(x), f(xdash)] acting on matrix of distances
  Cov_fx_fxdash <- function(dist_matrix) {
    sigma^2 * exp(-(dist_matrix)^2)
  }
  
  E_D <- rep(E_f, nD)
  Var_D <- Cov_fx_fxdash( as.matrix(dist(xD)) )
  E_fx <- rep(E_f, nP)
  Var_fx <- rep(sigma^2, nP)
  diag(Var_D) <- diag(Var_D) + delta
  
  Cov_fx_D <- Cov_fx_fxdash( as.matrix(pdist(xP, xD)) )
  Var_D_inv <- chol2inv(chol(Var_D))
  
  cov_fx_D_Var_D_inv <- Cov_fx_D %*% Var_D_inv
  ED_fx <- E_fx + cov_fx_D_Var_D_inv %*% (D - E_D)
  VarD_fx <- Var_fx - apply(cov_fx_D_Var_D_inv * Cov_fx_D, 1, sum)
  
  return(cbind("ExpD_f.x."=c(ED_fx), "VarD_f.x."=VarD_fx))
}


