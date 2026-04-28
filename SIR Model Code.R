# Code to solve the SIR model over time for parameters beta and gamma

################################################################
# SIR Model
################################################################

library(deSolve)

# simple calculator of the derivatives at a given time point for given model 
# parameters and starting position (S_0, I_0, R_0)
# output = 3-element list
SIR_simple <- function(t,     # chosen time point
                       f,     # initial position (S_0, I_0, R_0)
                       parms  # values of model parameters: list(beta, gamma)
                       ){		
  with(as.list(parms),{
    
    # Now follows a list of differential equations for each of the 3 outputs, S, I and R:
    dS = -beta * f["S"] * f["I"]
    dI = beta * f["S"] * f["I"] - gamma * f["I"]
    dR = gamma * f["I"]
    
    f_result <- c(dS,dI,dR)
    list(f_result)			# this returns a list of the derivatives for use by the numerical solver below
  })
}



# simple calculator of the derivatives at a given time point for given model 
# parameters and starting position (S_0, I_0, R_0)
# output = 3-element data.frame
SIR_deriv <- function(t,     # chosen time point
                      f,     # initial position (S_0, I_0, R_0)
                      parms  # values of model parameters: list(beta, gamma)
                      ){		
  with(as.list(parms),{
    
    # Now follows a list of differential equations for each of the 3 outputs, S, I and R:
    dS = -beta * f[,"S"] * f[,"I"]
    dI = beta * f[,"S"] * f[,"I"] - gamma * f[,"I"]
    dR = gamma * f[,"I"]
    
    f_result <- cbind(dS,dI,dR)
    return(as.data.frame(f_result))
  })
}


# Calculates the number of people in each model compartment over a set of times
# for given model parameters and starting position (S_0, I_0, R_0)
# output = data.frame with columns: S(t), I(t), R(t), dS(t), dI(t), dR(t)
SIR_model <- function(beta,                          # value of parameter beta
                      gamma,                         # value of parameter gamma
                      f_start,                       # initial position (S_0, I_0, R_0)
                      times=seq(0, 2.5, length=200)  # series of time points at which to evaluate the model
                      ){
  
  x_inputs <- c(beta=beta, gamma=gamma)
  output <- as.matrix(lsoda(y=f_start, times = times, func= SIR_simple, parms=x_inputs))
  f_s <- output[,-1]
  derivs <- SIR_deriv(times, f=f_s, parms=x_inputs)
  
  return(cbind(as.data.frame(output), derivs))
}

# Calculates the number of people in each model compartment over a set of times
# for a set of different input parameters given starting position (S_0, I_0, R_0)
# output = list of data.frame for each time point having columns: S(t), I(t), R(t), dS(t), dI(t), dR(t) 
SIR_multi_inputs <- function(xP,                              # matrix of input values (beta, gamma)
                             timepoints,                      # sequence of time points at which to evaluate the model for each set of input parameters
                             f_0=c(S = 600, I = 100, R = 50)  # initial position (S_0, I_0, R_0)
                             ){
  # Arrange time points
  if(!(0 %in% timepoints)){
    timepoints <- c(0, timepoints)
  }
  timepoints <- sort(timepoints)
  
  # set beta-gamma vector p
  beta <- xP[,1]
  gamma <- xP[,2]
  
  if(length(timepoints)==1){ # only the case if timepoints=c(0)
    df <- as.data.frame(matrix(rep(NA, 7*nrow(xP)), ncol=7))
    colnames(df) <- c("time", "S", "I", "R", "dS", "dI", "dR")
    df[,"time"] <- rep(timepoints, nrow(df))
    df[,c("S", "I", "R")] <- rep(f_0, each=nrow(df))
    for(i in 1:nrow(df)){
      df[i,c("dS", "dI", "dR")] <- SIR_deriv(t=0, f=df[i,c("S", "I", "R")], parms=c(beta=beta[i], gamma=gamma[i]))
    }
    return(list("t=0"=df))
  }
  
  # define store for beta-gamma pair output in separate vectors for each time snapshot
  ts_names <- paste0("t=", timepoints)
  timepoint_store <- setNames(vector("list", length(timepoints)), ts_names)
  for(j in ts_names){
    timepoint_store[[j]] <- c()
  }
  
  # run model
  for(i in 1:nrow(xP)){
    out <- SIR_model(beta=beta[i], gamma=gamma[i], f_start=f_0, times=timepoints)
    for(j in 1:length(ts_names)){
      store <- ts_names[j]
      timepoint_store[[store]] <- rbind(timepoint_store[[store]], out[j,])
    }
  }
  return(timepoint_store)
}

# Calculates the number of people in chosen model compartments at a given time point
# for a set of different input parameters with starting position (S_0, I_0, R_0)
# output = vector of the number of people in the chosen compartment(s) at this time point
SIR_f <- function(xP,                             # matrix of input values (beta, gamma)
                  timepoint,                      # time point at which to calculate compartment sizes
                  variable,                       # the compartment(s) for which to calculate outputs
                  f_0=c(S = 600, I = 100, R = 50) # initial position (S_0, I_0, R_0)
                  ){
  out <- SIR_multi_inputs(xP=xP, timepoints=timepoint, f_0=f_0)
  label <- paste0("t=", timepoint)
  return(out[[label]][,variable])
}


################################################################
# SIR Model Derivative with respect to Beta (at beta=0)
################################################################

# S(t), I(t) and R(t) when beta=0
gradf_full <- function(xP,          # input matrix of (beta, gamma) on the boundary beta=0
                       f_0=f_start, # initial position (S_0, I_0, R_0)
                       variable     # the compartment(s) for which to calculate outputs
                       ){ 
  S0 <- f_0["S"]
  I0 <- f_0["I"]
  R0 <- f_0["R"]
  gammas <- xP[,2]
  
  dS_dbeta <- (1/gammas) * S0*I0 * exp(-gammas*timepoint) - (1/gammas) * S0*I0
  dI_dbeta <- S0*I0 * timepoint * exp(-gammas*timepoint)
  dR_dbeta <- - (timepoint + (1/gammas)) * S0*I0 * exp(-gammas*timepoint) + (1/gammas) * S0*I0
  
  # 
  gamma_0_points <- which(gammas==0)
  dS_dbeta[gamma_0_points] <- - S0*I0 * timepoint
  # dI_dbeta[gamma_0_points] <- S0*I0 * timepoint
  dR_dbeta[gamma_0_points] <- 0
  
  df <- list("S" = dS_dbeta,
             "I" = dI_dbeta,
             "R" = dR_dbeta)
  
  df_var <- df[[variable]]
  n <- length(df_var)
  return( cbind(df_var, rep(NA, n)) )
}
