# Code to generate all plots in the accompanying report depicting the SIR Model

source("./SIR Model Code.R")
source("./SIR Model over Time Plot Saving Function.R")
source("./SIR Emulator Output Plot Generators.R")
source("./SIR Emulator Output Plot Saving Functions.R")
source("./LHD Generator Functions.R")
source("./Simple Emulator Code.R")
source("./Derivative Emulator Code.R")
source("./Boundary Emulators (without runs) Code.R")
source("./Boundary with Runs Emulators Code.R")
source("./Boundary with Derivative Emulators (without runs) Code.R")
source("./Boundary with Derivative and Runs Emulators Code.R")

###########################################
# SIR change over time curves plot
###########################################

f_start <- c(S = 600, I = 100, R = 50)	
times <- seq(0, 2.5, length=200)

# beta = 0.01, gamma = 1 plot
beta <- 0.01
gamma <- 1

output <- SIR_model(beta=beta, gamma=gamma, f_start, times=times)

saving_SIR_plots(output=output, beta=beta, gamma=gamma, f_start=f_start, 
                 times=times, wd=800, ht=620, #ylim=ylim,
                 filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/",
                 deriv_plot=TRUE)

# beta = 0, gamma = 1 plot
beta <- 0
gamma <- 1
ylim <- c(0, 750)

output <- SIR_model(beta=beta, gamma=gamma, f_start, times=times)

saving_SIR_plots(output=output, beta=beta, gamma=gamma, f_start=f_start, 
                 times=times, wd=800, ht=620, ylim=ylim,
                 filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/",
                 deriv_plot=TRUE)



# showing our chosen time point (t=1) and input space over which to emulate 
# contains parameter combinations where the peak of infections has already passed, is
# currently happening and is still to come

f_start <- c(S = 600, I = 5, R = 10)	
times <- seq(0, 2.5, length=200)
ylim <- c(0, 750)
timepoint <- 1

# beta = 0.035, gamma = 3
beta <- 0.035
gamma <- 3

output <- SIR_model(beta=beta, gamma=gamma, f_start, times=times)

saving_SIR_plots(output=output, beta=beta, gamma=gamma, f_start=f_start, 
                 times=times, wd=800, ht=620, ylim=ylim,
                 filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/",
                 deriv_plot=TRUE, highlight_timepoint=TRUE, timepoint=timepoint)

# beta = 0.013, gamma = 1
beta <- 0.013
gamma <- 1

output <- SIR_model(beta=beta, gamma=gamma, f_start, times=times)

saving_SIR_plots(output=output, beta=beta, gamma=gamma, f_start=f_start, 
                 times=times, wd=800, ht=620, ylim=ylim,
                 filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/",
                 deriv_plot=TRUE, highlight_timepoint=TRUE, timepoint=timepoint)

# beta = 0.0055, gamma = 0.4
beta <- 0.0055
gamma <- 0.4 

output <- SIR_model(beta=beta, gamma=gamma, f_start, times=times)

saving_SIR_plots(output=output, beta=beta, gamma=gamma, f_start=f_start, 
                 times=times, wd=800, ht=620, ylim=ylim,
                 filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/",
                 deriv_plot=TRUE, highlight_timepoint=TRUE, timepoint=timepoint)

# following plot not used
# beta = 0.005, gamma = 1.5
beta <- 0.005
gamma <- 1.5

output <- SIR_model(beta=beta, gamma=gamma, f_start, times=times)

saving_SIR_plots(output=output, beta=beta, gamma=gamma, f_start=f_start, 
                 times=times, wd=800, ht=620, ylim=ylim,
                 filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/",
                 deriv_plot=TRUE, highlight_timepoint=TRUE, timepoint=timepoint)





########################################################
# SIR 2d Emulation - set up
########################################################

beta_lim <- c(0, 0.04)
gamma_lim <- c(0, 5)
beta_grid <- seq(beta_lim[1], beta_lim[2], len=80)
gamma_grid <- seq(gamma_lim[1], gamma_lim[2], len=80)
xP <- as.matrix(expand.grid("x1"=beta_grid, "x2"=gamma_grid))

set.seed(112)
n <- 20
xD <- best_lhd(n=n, M=50, x1_lim=beta_lim, x2_lim=gamma_lim, print_switch=FALSE)
nD <- nrow(xD)
plot(xD, xlim=beta_lim, ylim=gamma_lim, pch=16, xaxs="i", yaxs="i", col="blue", 
     xlab=expression(beta), ylab=expression(gamma), cex=1.4)
abline(v=beta_lim[1]+(beta_lim[2]-beta_lim[1])*(0:nD)/nD, col="grey60")
abline(h=gamma_lim[1]+(gamma_lim[2]-gamma_lim[1])*(0:nD)/nD, col="grey60")

# Emulation t=1
variable <- c("S", "I", "R")        # "S", "I" or "R"
timepoint <- 1
f_start <- c(S=600, I=5, R=10)

output <- SIR_f(xP=xD, timepoint=timepoint, variable=variable, f_0=f_start)

#########################################
# Emulator priors (and plot contour levels)
#########################################

sigma_S <- 110
E_f_S <- 250
cont_levs_exp_S <- seq(-50, 650, 50)
cont_levs_var_S <- seq(0, ceiling(sigma_S^2/1000)*1000, by=1000)

sigma_I <- 110
E_f_I <- 250
cont_levs_exp_I <- seq(-50, 650, 50)
cont_levs_var_I <- seq(0, ceiling(sigma_I^2/1000)*1000, by=1000)

sigma_R <- 110
E_f_R <- 300
cont_levs_exp_R <- seq(-50, 650, 50)
cont_levs_var_R <- seq(0, ceiling(sigma_R^2/1000)*1000, by=1000)

theta <- c(0.04/5, 1)*1
delta <- 1e-06

######################################
# Simple Emulation
######################################

# Emulation of S(t=1)
D <- output[,"S"]

sigma <- sigma_S
E_f <- E_f_S
cont_levs_exp <- cont_levs_exp_S
cont_levs_var <- cont_levs_var_S

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/Simple/")

simple_em_out <- simple_BL_emulator_2d_fast(xP=xP, xD=xD, D=D, 
                                            theta=theta, sigma=sigma, E_f=E_f)

saving_simple_em_plots(em_out=simple_em_out, xD=xD, 
                       beta_grid=beta_grid, gamma_grid=gamma_grid,
                       timepoint=timepoint, variable="S", 
                       cont_levs_exp=cont_levs_exp,
                       cont_levs_var=cont_levs_var,
                       wd=690, ht=620,
                       filepath=filepath,
                       EmulatorType="Simple",
                       diag_true_func_plot=TRUE)

# Emulation of I(t=1)
D <- output[,"I"]

sigma <- sigma_I
E_f <- E_f_I
cont_levs_exp <- cont_levs_exp_I
cont_levs_var <- cont_levs_var_I

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/Simple/")

simple_em_out <- simple_BL_emulator_2d_fast(xP=xP, xD=xD, D=D, 
                                            theta=theta, sigma=sigma, E_f=E_f)

saving_simple_em_plots(em_out=simple_em_out, xD=xD, 
                       beta_grid=beta_grid, gamma_grid=gamma_grid,
                       timepoint=timepoint, variable="I", 
                       cont_levs_exp=cont_levs_exp,
                       cont_levs_var=cont_levs_var,
                       wd=690, ht=620,
                       filepath=filepath,
                       EmulatorType="Simple",
                       diag_true_func_plot=TRUE)

# Emulation of R(t=1)
D <- output[,"R"]

sigma <- sigma_R
E_f <- E_f_R
cont_levs_exp <- cont_levs_exp_R

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/Simple/")

simple_em_out <- simple_BL_emulator_2d_fast(xP=xP, xD=xD, D=D, 
                                            theta=theta, sigma=sigma, E_f=E_f)

saving_simple_em_plots(em_out=simple_em_out, xD=xD, 
                       beta_grid=beta_grid, gamma_grid=gamma_grid,
                       timepoint=timepoint, variable="R", 
                       cont_levs_exp=cont_levs_exp,
                       cont_levs_var=cont_levs_var,
                       wd=690, ht=620,
                       filepath=filepath,
                       EmulatorType="Simple",
                       diag_true_func_plot=TRUE)

#########################################
# Boundary Emulation and BwR Emulation
#########################################

boundary <- list("x1"=0, "x2"=NULL)

# Emulation of S(t=1)
f <- function(xP, f_0=f_start){ # S(t), I(t) and R(t) when beta=0
  df <- list("S" = f_0["S"],
             "I" = f_0["I"]*exp(-xP[,2]*timepoint),
             "R" = f_0["I"]*(1 - exp(-xP[,2]*timepoint)) + f_0["R"])
  return( df[["S"]] )
}

D <- output[,"S"]

sigma <- sigma_S
E_f <- E_f_S
cont_levs_exp <- cont_levs_exp_S
cont_levs_var <- cont_levs_var_S

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/Boundary/")

boundary_em_out <- single_boundary_BL_emulator_fast(xP=xP, boundary=boundary,
                                                    theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=boundary_em_out, boundary=boundary, 
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="S", 
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary",
                                diag_true_func_plot=TRUE)

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/BwR/")

BwR_em_out <- single_boundary_with_runs_BL_emulator(xP=xP, xD=xD, D=D, 
                                                    boundary=boundary,
                                                    theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=BwR_em_out, boundary=boundary,
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="S", xD=xD,
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary with Runs",
                                diag_true_func_plot=TRUE)

# Emulation of I(t=1)
f <- function(xP, f_0=f_start){ # S(t), I(t) and R(t) when beta=0
  df <- list("S" = f_0["S"],
             "I" = f_0["I"]*exp(-xP[,2]*timepoint),
             "R" = f_0["I"]*(1 - exp(-xP[,2]*timepoint)) + f_0["R"])
  return( df[["I"]] )
}

D <- output[,"I"]

sigma <- sigma_I
E_f <- E_f_I
cont_levs_exp <- cont_levs_exp_I
cont_levs_var <- cont_levs_var_I

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/Boundary/")

boundary_em_out <- single_boundary_BL_emulator_fast(xP=xP, boundary=boundary,
                                                    theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=boundary_em_out, boundary=boundary, 
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="I", 
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary",
                                diag_true_func_plot=TRUE)

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/BwR/")

BwR_em_out <- single_boundary_with_runs_BL_emulator(xP=xP, xD=xD, D=D, 
                                                    boundary=boundary,
                                                    theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=BwR_em_out, boundary=boundary,
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="I", xD=xD,
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary with Runs",
                                diag_true_func_plot=TRUE)

# Emulation of R(t=1)
f <- function(xP, f_0=f_start){ # S(t), I(t) and R(t) when beta=0
  df <- list("S" = f_0["S"],
             "I" = f_0["I"]*exp(-xP[,2]*timepoint),
             "R" = f_0["I"]*(1 - exp(-xP[,2]*timepoint)) + f_0["R"])
  return( df[["R"]] )
}

D <- output[,"R"]

sigma <- sigma_R
E_f <- E_f_R
cont_levs_exp <- cont_levs_exp_R
cont_levs_var <- cont_levs_var_R

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/Boundary/")

boundary_em_out <- single_boundary_BL_emulator_fast(xP=xP, boundary=boundary,
                                                    theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=boundary_em_out, boundary=boundary, 
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="R", 
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary",
                                diag_true_func_plot=TRUE)

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/BwR/")

BwR_em_out <- single_boundary_with_runs_BL_emulator(xP=xP, xD=xD, D=D, 
                                                    boundary=boundary,
                                                    theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=BwR_em_out, boundary=boundary,
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="R", xD=xD,
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary with Runs",
                                diag_true_func_plot=TRUE)

##############################################################
# Boundary with Derivative Emulation and BwDR Emulation
##############################################################

boundary <- list("x1"=0, "x2"=NULL)

# Emulation of S(t=1)
f <- function(xP, f_0=f_start){ # S(t), I(t) and R(t) when beta=0
  df <- list("S" = f_0["S"],
             "I" = f_0["I"]*exp(-xP[,2]*timepoint),
             "R" = f_0["I"]*(1 - exp(-xP[,2]*timepoint)) + f_0["R"])
  return( df[["S"]] )
}

gradf <- function(xP){
  return(gradf_full(xP, f_0=f_start, variable="S"))
}

D <- output[,"S"]

sigma <- sigma_S
E_f <- E_f_S
E_df <- 0
cont_levs_exp <- cont_levs_exp_S
cont_levs_var <- cont_levs_var_S

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/BwD/")

BwD_em_out <- single_BwD_BL_emulator(xP=xP, boundary=boundary, 
                                     theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)

saving_boundary_with_runs_plots(em_out=BwD_em_out, boundary=boundary, 
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="S", 
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary and Derivative on the Boundary",
                                diag_true_func_plot=TRUE)

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/BwDR/")

BwDR_em_out <- single_BwD_with_runs_BL_emulator(xP=xP, xD=xD, D=D, boundary=boundary, 
                                                theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)

saving_boundary_with_runs_plots(em_out=BwDR_em_out, boundary=boundary, 
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="S", xD=xD,
                                cont_levs_exp=cont_levs_exp, # seq(-50, 650, 10)
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary and Derivative on the Boundary with Runs",
                                diag_true_func_plot=TRUE)

# Emulation of I(t=1)
f <- function(xP, f_0=f_start){ # S(t), I(t) and R(t) when beta=0
  df <- list("S" = f_0["S"],
             "I" = f_0["I"]*exp(-xP[,2]*timepoint),
             "R" = f_0["I"]*(1 - exp(-xP[,2]*timepoint)) + f_0["R"])
  return( df[["I"]] )
}

gradf <- function(xP){
  return(gradf_full(xP, f_0=f_start, variable="I"))
}

D <- output[,"I"]

sigma <- sigma_I
E_f <- E_f_I
E_df <- 0
cont_levs_exp <- cont_levs_exp_I
cont_levs_var <- cont_levs_var_I

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/BwD/")

BwD_em_out <- single_BwD_BL_emulator(xP=xP, boundary=boundary, 
                                     theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)

saving_boundary_with_runs_plots(em_out=BwD_em_out, boundary=boundary, 
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="I", 
                                cont_levs_exp=cont_levs_exp, # seq(-50, 500, 10)
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary and Derivative on the Boundary",
                                diag_true_func_plot=TRUE)

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/BwDR/")

BwDR_em_out <- single_BwD_with_runs_BL_emulator(xP=xP, xD=xD, D=D, boundary=boundary, 
                                                theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)

saving_boundary_with_runs_plots(em_out=BwDR_em_out, boundary=boundary, 
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="I", xD=xD,
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary and Derivative on the Boundary with Runs",
                                diag_true_func_plot=TRUE)

# Emulation of R(t=1)
f <- function(xP, f_0=f_start){ # S(t), I(t) and R(t) when beta=0
  df <- list("S" = f_0["S"],
             "I" = f_0["I"]*exp(-xP[,2]*timepoint),
             "R" = f_0["I"]*(1 - exp(-xP[,2]*timepoint)) + f_0["R"])
  return( df[["R"]] )
}

gradf <- function(xP){
  return(gradf_full(xP, f_0=f_start, variable="R"))
}

D <- output[,"R"]

sigma <- sigma_R
E_f <- E_f_R
E_df <- 0
cont_levs_exp <- cont_levs_exp_R
cont_levs_var <- cont_levs_var_R

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/BwD/")

BwD_em_out <- single_BwD_BL_emulator(xP=xP, boundary=boundary, 
                                     theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)

saving_boundary_with_runs_plots(em_out=BwD_em_out, boundary=boundary, 
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="R", 
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary and Derivative on the Boundary",
                                diag_true_func_plot=TRUE)

filepath <- paste0("/Users/helenbate/Documents/Year 4/Diss/Code/Images/SIR Model Example/DISS/t=1 Emulation/BwDR/")

BwDR_em_out <- single_BwD_with_runs_BL_emulator(xP=xP, xD=xD, D=D, boundary=boundary, 
                                                theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)

saving_boundary_with_runs_plots(em_out=BwDR_em_out, boundary=boundary, 
                                beta_grid=beta_grid, gamma_grid=gamma_grid,
                                timepoint=timepoint, variable="R", xD=xD,
                                cont_levs_exp=cont_levs_exp,
                                cont_levs_var=cont_levs_var,
                                wd=690, ht=620,
                                filepath=filepath,
                                EmulatorType="Boundary and Derivative on the Boundary with Runs",
                                diag_true_func_plot=TRUE)
