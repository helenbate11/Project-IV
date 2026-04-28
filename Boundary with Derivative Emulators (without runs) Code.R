# Code to implement the Boundary with Derivative Emulator on a 2d input space 
# with a single boundary or 2 perpendicular boundaries

# Note: calls on global functions of interest f() and gradf() (at the relevant boundary)

################################################################
# (Single-) Boundary with Derivative Emulator (without runs) on a 2d input space
################################################################

single_BwD_BL_emulator <- function(xP,                               # emulator prediction points
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
  
  # K projections of xP
  xK <- xP
  xK[,axis] <- rep(constant, nrow(xP))
  a <- xP[,axis] - rep(constant, nrow(xP))
  
  E_fx <- rep(E_f, nrow(xP))
  E_dfx <- rep(E_df, nrow(xP))
  Var_fx <- rep(sigma^2, nrow(xP))
  
  EK_fx <- E_fx + exp(-(a/theta[axis])^2) * (f(xK) - E_fx) +
    a*exp(-(a/theta[axis])^2) * (gradf(xK)[,axis] - E_dfx)
  VarK_fx <- Var_fx * ( 1 - exp(-2*(a/theta[axis])^2) - 2*(a/theta[axis])^2*exp(-2*(a/theta[axis])^2) )
  
  return(cbind("ExpD_f.x."=c(EK_fx), "VarD_f.x."=VarK_fx))
}

# boundary <- list("x1"=NULL, "x2"=0.4)
# em_out <- single_boundary_BL_emulator_fast(xP, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f)

# saving_boundary_with_runs_plots(em_out=em_out, boundary=boundary, x_grid=x_grid, xD=NULL,
#                                 filepath="Documents/Year 4/Diss/Code/Slow Boundary Wk_14/Boundary with Derivative/Images/", function_name="F1", 
#                                 EmulatorType="Boundary and Derivatives on the Boundary")

################################################################
# (Perpendicular-) Boundary with Derivative Emulator (without runs) on a 2d input space
################################################################

perp_BwD_BL_emulator <- function(xP,       # emulator prediction points
                                 boundary=list("x1"=0, "x2"=0), # boundary definitions: x1 = boundary$x1 AND x2 = boundary$x2
                                 theta=1,  # the correlation length vector
                                 sigma=1,  # the prior standard deviation
                                 E_f=0,    # the prior expectation of f(x)
                                 E_df=c(0,0,0)) { # the prior expectation of (df/dx_1(x), df/dx_2(x), d^2f/dx_1dx_2(x)) in the perpendicular directions to each boundary
  if(length(theta)==1){
    theta <- rep(theta, 2)
  }
  
  # K and L projections of xP
  xK <- xP
  xL <- xP
  xKL <- xP
  xK[,1] <- rep(boundary$x1, nrow(xP))  # (x1=constant[1])
  xL[,2] <- rep(boundary$x2, nrow(xP))  # (x1=constant[2])
  xKL <- matrix(c(rep(boundary$x1, nrow(xP)), rep(boundary$x2, nrow(xP))), ncol=2) # projection onto x1, x2 corner
  a <- xP[,1] - rep(boundary$x1, nrow(xP))
  b <- xP[,2] - rep(boundary$x2, nrow(xP))
  
  E_fx <- rep(E_f, nrow(xP))
  E_dfx <- matrix(rep(E_df, nrow(xP)), ncol=2, byrow=TRUE)
  Var_fx <- rep(sigma^2, nrow(xP))
  
  EKL_fx <- E_fx + exp(-(a/theta[1])^2) * (f(xK) - E_fx) + a*exp(-(a/theta[1])^2) * (gradf(xK)[,1] - E_dfx[1]) +
    exp(-(b/theta[2])^2) * (f(xL) - E_fx) + b*exp(-(b/theta[2])^2) * (gradf(xL)[,2] - E_dfx[2]) -
    exp(-(a/theta[1])^2) * exp(-(b/theta[2])^2) * (f(xKL) - E_fx) -
    a * exp(-(a/theta[1])^2) * exp(-(b/theta[2])^2) * (gradf(xKL)[,1] - E_dfx[1]) -
    b * exp(-(a/theta[1])^2) * exp(-(b/theta[2])^2) * (gradf(xKL)[,2] - E_dfx[2]) -
    a * b * exp(-(a/theta[1])^2) * exp(-(b/theta[2])^2) * (gradf(xKL)[,3] - E_dfx[3])
  
  VarKL_fx <- Var_fx * ( 1 - exp(-2*(a/theta[1])^2) - 2*(a/theta[1])^2*exp(-2*(a/theta[1])^2) ) *
    ( 1 - exp(-2*(b/theta[2])^2) - 2*(b/theta[2])^2*exp(-2*(b/theta[2])^2) )
  
  return(cbind("ExpD_f.x."=c(EKL_fx), "VarD_f.x."=VarKL_fx))
}


