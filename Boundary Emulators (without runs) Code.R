# (Separate) code to implement the Boundary Emulator on a 2d input space with a 
# single boundary, 2 perpendicular boundaries or 2 parallel boundaries and
# (Complete) function that will identify relevant emulator to use, 
# given the boundary definitions

# Note: calls on global function of interest f()

################################################################
# (Single-) Boundary Emulator (boundary only) on a 2d input space
################################################################

single_boundary_BL_emulator_fast <- function(xP,                               # emulator prediction points
                                             boundary=list("x1"=0, "x2"=NULL), # boundary definition: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL)
                                             theta=1,                          # the correlation length vector
                                             sigma=1,                          # the prior standard deviation
                                             E_f=0                             # the prior expectation of f(x)
                                             ){
  if(length(theta)==1){
    theta <- rep(theta, 2)
  }
  
  axis <- ifelse(is.null(boundary$x1), 2, 1)
  constant <- ifelse(is.null(boundary$x1), boundary$x2, boundary$x1)
  
  # K projections of xP
  xK <- xP
  xK[,axis] <- rep(constant, nrow(xP))
  a <- (xP[,axis] - rep(constant, nrow(xP))) / theta[axis]
  
  E_fx <- rep(E_f, nrow(xP))
  Var_fx <- rep(sigma^2, nrow(xP))
  
  ED_fx <- E_fx + exp(-a^2) * (f(xK) - E_fx)
  VarD_fx <- Var_fx - sigma^2 * exp(-2*a^2)
  
  return(cbind("ExpD_f.x."=c(ED_fx), "VarD_f.x."=VarD_fx))
}

################################################################
# (Perpendicular-) Boundary Emulator (boundary only) on a 2d input space
################################################################

perp_boundary_BL_emulator_fast <- function(xP,                            # emulator prediction points
                                           boundary=list("x1"=0, "x2"=1), # boundary definitions: x1 = boundary$x1 AND x2 = boundary$x2
                                           theta=1,                       # the correlation length vector
                                           sigma=1,                       # the prior standard deviation
                                           E_f=0) {                       # the prior expectation of f(x)
  if(length(theta)==1){
    theta <- rep(theta, 2)
  }
  
  # K projections of xP
  xK <- xP
  xL <- xP
  xKL <- xP
  xK[,2] <- rep(boundary$x2, nrow(xP))  # projection onto x1 boundary (x2=constant[2])
  xL[,1] <- rep(boundary$x1, nrow(xP))  # projection onto x2 boundary (x1=constant[1])
  xKL <- matrix(c(rep(boundary$x1, nrow(xP)), rep(boundary$x2, nrow(xP))), ncol=2) # projection onto x1, x2 corner
  a <- (xP[,2] - rep(boundary$x2, nrow(xP))) / theta[2]
  b <- (xP[,1] - rep(boundary$x1, nrow(xP))) / theta[1]
  
  E_fx <- rep(E_f, nrow(xP))
  Var_fx <- rep(sigma^2, nrow(xP))
  
  ED_fx <- E_fx + exp(-a^2)*(f(xK) - E_fx) + exp(-b^2)*(f(xL) - E_fx) - exp(-a^2)*exp(-b^2)*(f(xKL) - E_fx)
  VarD_fx <- Var_fx * (1 - exp(-2*a^2)) * (1 - exp(-2*b^2))
  
  return(cbind("ExpD_f.x."=c(ED_fx), "VarD_f.x."=VarD_fx))
}


################################################################
# (Parallel-) Boundary Emulator (boundary only) on a 2d input space
################################################################

parallel_boundary_BL_emulator_fast <- function(xP,                                    # emulator prediction points
                                               boundary=list("x1"=c(0,1), "x2"=NULL), # boundary definitions: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL)
                                               theta=1,                               # the correlation length vector
                                               sigma=1,                               # the prior standard deviation
                                               E_f=0) {                               # the prior expectation of f(x)
  if(length(theta)==1){
    theta <- rep(theta, 2)
  }
  double_axis <- ifelse(!is.null(boundary$x1), 1, 2)
  constant <- ifelse(c(!is.null(boundary$x1), !is.null(boundary$x1)), 
                     boundary$x1, boundary$x2)
  
  # K projections of xP
  axis_with_boundary <- c(2, 1)
  xK <- xP
  xL <- xP
  xK[,double_axis] <- rep(constant[1], nrow(xP))  # projection onto other boundary (double_axis=constant[1])
  xL[,double_axis] <- rep(constant[2], nrow(xP))  # projection onto other boundary (double_axis=constant[2])
  a <- (xP[,double_axis] - rep(constant[1], nrow(xP))) / theta[double_axis]
  b <- (xP[,double_axis] - rep(constant[2], nrow(xP))) / theta[double_axis]
  c <- (a + b) / theta[double_axis]
  
  coefficient <- function(y, z) {
    ( exp(-y^2) - exp(-z^2)*exp(-c^2) ) / ( 1 - exp(-2*c^2) )
  }
  
  E_fx <- rep(E_f, nrow(xP))
  Var_fx <- rep(sigma^2, nrow(xP))
  
  ED_fx <- E_fx + coefficient(a, b)*(f(xK) - E_fx) + coefficient(b, a)*(f(xL) - E_fx)
  VarD_fx <- Var_fx * 1/(1 - exp(-2*c^2)) * (1 - exp(-2*c^2) - exp(-2*a^2) - exp(-2*b^2) + 2*exp(-(c^2+a^2+b^2)))
  
  return(cbind("ExpD_f.x."=c(ED_fx), "VarD_f.x."=VarD_fx))
}


################################################################
# (Complete) Boundary Emulator (boundary only) on a 2d input space
################################################################

boundary_BL_emulator_fast <- function(xP,                            # emulator prediction points
                                      boundary=list("x1"=0, "x2"=0), # boundary definition: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL if single- or parallel-boundary case)
                                      theta=1,                       # the correlation length vector
                                      sigma=1,                       # the prior standard deviation
                                      E_f=0                          # the prior expectation of f(x)
                                      ){
  if(!is.null(boundary$x1) & !is.null(boundary$x2)){ # x1=c(), x2=c()
    return(perp_boundary_BL_emulator_fast(xP, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f))
  }
  if(!is.null(boundary$x1) & is.null(boundary$x2)){  # x1=c(), x2=NULL
    if(length(boundary$x1)==1){                      # x1=[], x2=NULL
      return(single_boundary_BL_emulator_fast(xP, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f))
    }else{
      return(parallel_boundary_BL_emulator_fast(xP, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f))
    }
  }
  if(is.null(boundary$x1) & !is.null(boundary$x2)){ # x1=NULL, x2=c()
    if(length(boundary$x2)==1){                     # x1=NULL, x2=[]
      return(single_boundary_BL_emulator_fast(xP, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f))
    }else{
      return(parallel_boundary_BL_emulator_fast(xP, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f))
    }
  }
}
