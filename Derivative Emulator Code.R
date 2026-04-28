# Code to implement the Simple Emulator on 1d and 2d input spaces

################################################################
# Derivative Emulator on a 1d input space
################################################################

# Derivative Emulator for a single input x on a 1d input space
deriv_BL_emulator_1d_single <- function(x,        # emulator prediction point
                                        xD,       # run input locations xD
                                        D,        # outputs D=(f(x^1), ..., f(x^n))
                                        theta=1,  # correlation lengths
                                        sigma=1,  # prior SD sigma=sqrt(Var[f(x)])
                                        E_f=0,
                                        E_dfdx=0) {
  n <- length(xD)
  
  Cov_fx_fxdash <- function(x, xdash) {
    sigma^2 * exp(-(x-xdash)^2 / (theta^2))
  }
  
  Cov_fx_df_dxdash <- function(x, xdash) {
    2*sigma^2 / (theta^2) * (x-xdash) * exp(-(x-xdash)^2 / (theta^2))
  }
  
  Cov_dfdx_fxdash <- function(x, xdash) {
    - Cov_fx_df_dxdash(x, xdash)
  }
  
  Cov_dfdx_df_dxdash <- function(x, xdash) {
    2*sigma^2 / (theta^2) * (1 - 2/(theta^2) * (x-xdash)^2) * exp(-(x-xdash)^2 / (theta^2))
  }
  
  E_f <- E_f
  E_dfdx <- E_dfdx
  E_D <- c(rep(E_f, n), rep(E_dfdx, n))
  E_D
  
  Var_D <- matrix(0, nrow=2*n, ncol=2*n)
  for(i in 1:n){
    for(j in 1:n){
      Var_D[i,j] <- Cov_fx_fxdash(xD[i], xD[j])
      Var_D[i,n+j] <- Cov_fx_df_dxdash(xD[i], xD[j])
      Var_D[n+i,j] <- Cov_dfdx_fxdash(xD[i], xD[j])
      Var_D[n+i,n+j] <- Cov_dfdx_df_dxdash(xD[i], xD[j])
    }
  }
  diag(Var_D) <- diag(Var_D) + delta
  
  E_fx <- E_f
  Var_fx <- sigma^2
  Cov_fx_D <- matrix(0, nrow=1, ncol=2*n)
  for(j in 1:n){
    Cov_fx_D[1,j] <- Cov_fx_fxdash(x, xD[j])
    Cov_fx_D[1,n+j] <- Cov_fx_df_dxdash(x, xD[j])
  }
  
  Cov_D_fx <- matrix(0, nrow=2*n, ncol=1)
  for(i in 1:n){
    Cov_D_fx[i,1] <- Cov_fx_fxdash(xD[i], x)
    Cov_D_fx[n+i,1] <- Cov_dfdx_fxdash(xD[i], x)
  }
  
  # Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x)
  ED_fx <- E_fx + Cov_fx_D %*% solve(Var_D) %*% (D - E_D)
  VarD_fx <- Var_fx - Cov_fx_D %*% solve(Var_D) %*% Cov_D_fx
  
  # Return emulator expectation and variance
  return(c("ExpD_f(x)"=ED_fx, "VarD_f(x)"=VarD_fx))
}


# Derivative Emulator for on a 1d input space
deriv_BL_emulator_1d <- function(xP,       # emulator prediction points
                                 xD,       # the run input locations
                                 D,        # the outputs, D = (f(x^1), ..., f(x^n))
                                 theta=1,  # the correlation length
                                 sigma=1,  # the prior standard deviation, sigma = sqrt(Var[f(x)])
                                 E_f=0,    # the prior expectation of f
                                 E_dfdx=0  # the prior expectation of df/dx
                                 ){
  
  em_out <- t(sapply(xP, deriv_BL_emulator_1d_single, xD=xD, D=D, theta=theta, sigma=sigma, 
                     E_f=E_f, E_dfdx=E_dfdx))
  colnames(em_out) <- c("ExpD_f(x)", "VarD_f(x)")
  return(em_out)
}






################################################################
# Derivative Emulator on a 2d input space
################################################################
library(pdist)

deriv_BL_emulator_2d_fast <- function(xP,              # emulator prediction points
                                      xD,              # the run input locations
                                      D,               # the outputs, D = (f(x^1), ..., f(x^n))
                                      theta=1,         # the correlation length
                                      sigma=1,         # the prior standard deviation, sigma = sqrt(Var[f(x)])
                                      E_f=0,           # the prior expectation of f
                                      E_dfdx=c(0,0)    # the prior expectation of df/dx, (df/dfx_1, df/dx_2)
                                      ){
  nD <- nrow(xD)
  nP <- nrow(xP)
  
  xP_resc <- t(t(xP)/theta)
  xD_resc <- t(t(xD)/theta)
  
  dist_matrix_xP_xD <- as.matrix(pdist(xP_resc, xD_resc))
  dist_matrix_resc1_xP_xD <- outer(xP_resc[,1], xD_resc[,1], "-")
  dist_matrix_resc2_xP_xD <- outer(xP_resc[,2], xD_resc[,2], "-")
  dist_matrix_xD_xD <- as.matrix(dist(xD_resc))
  dist_matrix_resc1_xD_xD <- outer(xD_resc[,1], xD_resc[,1], "-")
  dist_matrix_resc2_xD_xD <- outer(xD_resc[,2], xD_resc[,2], "-")
  
  Cov_fx_fxdash <- function(dist_matrix) {
    # Cov(f(x), f(x'))
    sigma^2 * exp(-(dist_matrix)^2)
  }
  
  Cov_fx_df_dxdash1 <- function(dist_matrix, dist_matrix_resc1) { 
    # Cov(f(x), df/dx'_1(x'))
    2*sigma^2 / (theta[1]) * dist_matrix_resc1 * exp(-(dist_matrix)^2)
  }
  
  Cov_fx_df_dxdash2 <- function(dist_matrix, dist_matrix_resc2) { 
    # Cov(f(x), df/dx'_2(x'))
    2*sigma^2 / (theta[2]) * dist_matrix_resc2 * exp(-(dist_matrix)^2)
  }
  
  Cov_dfdx1_fxdash <- function(dist_matrix, dist_matrix_resc1) {
    # Cov(df/dx_1(x), f(x'))
    - Cov_fx_df_dxdash1(dist_matrix, dist_matrix_resc1)
  }
  
  Cov_dfdx2_fxdash <- function(dist_matrix, dist_matrix_resc2) {
    # Cov(df/dx_2(x), f(x'))
    - Cov_fx_df_dxdash2(dist_matrix, dist_matrix_resc2)
  }
  
  Cov_dfdx1_df_dxdash1 <- function(dist_matrix, dist_matrix_resc1) {
    # Cov(df/dx_1(x), df/dx'_1(x'))
    2*sigma^2 / (theta[1]^2) * (1 - 2*dist_matrix_resc1^2) * exp(-(dist_matrix)^2)
  }
  
  Cov_dfdx2_df_dxdash2 <- function(dist_matrix, dist_matrix_resc2) {
    # Cov(df/dx_2(x), df/dx'_2(x'))
    2*sigma^2 / (theta[2]^2) * (1 - 2*dist_matrix_resc2^2) * exp(-(dist_matrix)^2)
  }
  
  Cov_dfdx_df_dxdash <- function(dist_matrix, dist_matrix_resc1, dist_matrix_resc2) {
    # Cov(df/dx_1(x), df/dx'_2(x')),    Cov(df/dx_2(x), df/dx'_1(x'))
    (-4)*sigma^2 / (theta[1] * theta[2]) * dist_matrix_resc1 * dist_matrix_resc2 * exp(-(dist_matrix)^2)
  }
  
  E_f <- E_f
  E_dfdx1 <- E_dfdx[1]
  E_dfdx2 <- E_dfdx[2]
  E_D <- c(rep(E_f, nD), rep(E_dfdx1, nD), rep(E_dfdx2, nD))
  
  Var_D <- rbind(
    cbind(Cov_fx_fxdash(dist_matrix_xD_xD), 
          Cov_fx_df_dxdash1(dist_matrix_xD_xD, dist_matrix_resc1_xD_xD), 
          Cov_fx_df_dxdash2(dist_matrix_xD_xD, dist_matrix_resc2_xD_xD)), 
    cbind(Cov_dfdx1_fxdash(dist_matrix_xD_xD, dist_matrix_resc1_xD_xD), 
          Cov_dfdx1_df_dxdash1(dist_matrix_xD_xD, dist_matrix_resc1_xD_xD), 
          Cov_dfdx_df_dxdash(dist_matrix_xD_xD, dist_matrix_resc1_xD_xD, dist_matrix_resc2_xD_xD)),
    cbind(Cov_dfdx2_fxdash(dist_matrix_xD_xD, dist_matrix_resc2_xD_xD), 
          Cov_dfdx_df_dxdash(dist_matrix_xD_xD, dist_matrix_resc1_xD_xD, dist_matrix_resc2_xD_xD),
          Cov_dfdx2_df_dxdash2(dist_matrix_xD_xD, dist_matrix_resc2_xD_xD)))
  diag(Var_D) <- diag(Var_D) + delta
  
  Var_D_inv <- chol2inv(chol(Var_D))
  
  # Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of f(x)
  E_fx <- rep(E_f, nP)
  Var_fx <- rep(sigma^2, nP)
  Cov_fx_D <- cbind(Cov_fx_fxdash(dist_matrix_xP_xD), 
                    Cov_fx_df_dxdash1(dist_matrix_xP_xD, dist_matrix_resc1_xP_xD), 
                    Cov_fx_df_dxdash2(dist_matrix_xP_xD, dist_matrix_resc2_xP_xD))
  
  cov_fx_D_Var_D_inv <- Cov_fx_D %*% Var_D_inv
  
  ED_fx <- E_fx + cov_fx_D_Var_D_inv %*% (D - E_D)
  VarD_fx <- Var_fx - apply(cov_fx_D_Var_D_inv * Cov_fx_D, 1, sum)
  
  # Perform Bayes Linear adjustment to find Adjusted Expectation and Variance of df/dx(x)
  E_dfx <- c(rep(E_dfdx1, nP), rep(E_dfdx2, nP))
  Var_dfx <- c( rep(2*sigma^2/(theta[1]^2), nP ), rep( 2*sigma^2/(theta[2]^2), nP ))
  Cov_dfx_D <- rbind(
    cbind(Cov_dfdx1_fxdash(dist_matrix_xP_xD, dist_matrix_resc1_xP_xD), 
          Cov_dfdx1_df_dxdash1(dist_matrix_xP_xD, dist_matrix_resc1_xP_xD), 
          Cov_dfdx_df_dxdash(dist_matrix_xP_xD, dist_matrix_resc1_xP_xD, dist_matrix_resc2_xP_xD)),
    cbind(Cov_dfdx2_fxdash(dist_matrix_xP_xD, dist_matrix_resc2_xP_xD), 
          Cov_dfdx_df_dxdash(dist_matrix_xP_xD, dist_matrix_resc1_xP_xD, dist_matrix_resc2_xP_xD),
          Cov_dfdx2_df_dxdash2(dist_matrix_xP_xD, dist_matrix_resc2_xP_xD)))
  
  cov_dfx_D_Var_D_inv <- Cov_dfx_D %*% Var_D_inv
  
  ED_dfx <- E_dfx + cov_dfx_D_Var_D_inv %*% (D - E_D)
  VarD_dfx <- Var_dfx - apply(cov_dfx_D_Var_D_inv * Cov_dfx_D, 1, sum)
  
  # Return emulator expectation and variance
  return(as.matrix(data.frame("ExpD_f.x."=ED_fx, 
                              "VarD_f.x."=VarD_fx,
                              "ExpD_dfdx1.x."=ED_dfx[1:(0.5*length(xP))],
                              "VarD_dfdx1.x."=VarD_dfx[1:(0.5*length(xP))],
                              "ExpD_dfdx2.x."=ED_dfx[(0.5*length(xP)+1):length(xP)],
                              "VarD_dfdx2.x."=VarD_dfx[(0.5*length(xP)+1):length(xP)])))
}
