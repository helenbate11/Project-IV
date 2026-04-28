# Code to generate all plots in the accompanying report depicting the 1d or 
# 2d toy function examples

source("./LHD Generator Functions.R")
source("./Simple Emulator Code.R")
source("./Derivative Emulator Code.R")
source("./Boundary Emulators (without runs) Code.R")
source("./Boundary with Runs Emulators Code.R")
source("./Boundary with Derivative Emulators (without runs) Code.R")
source("./Boundary with Derivative and Runs Emulators Code.R")
source("./Emulator Output Plot Generators.R")
source("./Emulator Output Plot Saving Functions.R")
source("./Variance Comparators.R")
source("./Variance Comparison Plot Generators.R")
source("./Variance Comparison Plot Saving Function.R")

########################################################
# 1d Emulation - defining f
########################################################

# 1d toy function and its derivative
f <- function(x){
  cos(x) + sin(pi*x) + sin(x/10)
}

df <- function(x){
  -sin(x) + pi*cos(pi*x) + 1/10*cos(x/10)
}

# 1d run positions
xD <- seq(0.1, 4, by=0.8)
D <- f(xD)
xP <- seq(-0.001, 4.001, len=401)
nP <- length(xP)

########################################################
# 1d Simple Emulation (stand-alone)
########################################################

theta <- 0.7
sigma <- 0.6
delta <- 10^(-6)

E_f <- 0
simple_em_out <- simple_BL_emulator_1d(xP=xP, xD=xD, D=D, theta=theta, sigma=sigma, E_f=E_f)

# plot_BL_emulator_1d(simple_em_out, xP, xD, D, ylim=c(-2.5, 2.5))

saving_1d_BL_plot(em_out=simple_em_out, xP=xP, xD=xD, D=D, ylim=c(-2.5,2.5), 
                  EmulatorType="Simple", filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/1d Simple/",
                  function_name="F5", wd=800, ht=620)
saving_1d_BL_plot(em_out=simple_em_out, xP=xP, xD=xD, D=D, ylim=c(-2.5,2.5), 
                  EmulatorType="Simple", filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/1d Deriv/",
                  function_name="F5", wd=800, ht=620)

# following plot not used
E_dfdx <- 0
D <- c(f(xD), df(xD))
deriv_em_out <- deriv_BL_emulator_1d(xP=xP, xD=xD, D=D, theta=theta, sigma=sigma, E_f=E_f, E_dfdx=E_dfdx)

# plot_BL_emulator_1d(deriv_em_out, xP, xD, f(xD), ylim=c(-2.5, 2.5))

saving_1d_BL_plot(em_out=deriv_em_out, xP=xP, xD=xD, D=f(xD), ylim=c(-2.5,2.5), 
                  EmulatorType="Derivative", filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/1d Deriv/",
                  function_name="F5", wd=800, ht=620)

########################################################
# 1d Simple vs Derivative Emulator (comparison pair)
########################################################

theta <- 0.5
sigma <- 0.6
E_f <- 0
D <- f(xD)

simple_em_out <- simple_BL_emulator_1d(xP=xP, xD=xD, D=D, theta=theta, sigma=sigma, E_f=E_f)

# plot_BL_emulator_1d(simple_em_out, xP, xD, D, ylim=c(-2.5, 2.5))

saving_1d_BL_plot(em_out=simple_em_out, xP=xP, xD=xD, D=D, ylim=c(-3,3), 
                  EmulatorType="Simple", filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/1d Deriv/Joint/",
                  function_name="F5")


E_dfdx <- 0
D <- c(f(xD), df(xD))
deriv_em_out <- deriv_BL_emulator_1d(xP=xP, xD=xD, D=D, theta=theta, sigma=sigma, E_f=E_f, E_dfdx=E_dfdx)

# plot_BL_emulator_1d(deriv_em_out, xP, xD, f(xD), ylim=c(-2.5, 2.5))

saving_1d_BL_plot(em_out=deriv_em_out, xP=xP, xD=xD, D=f(xD), ylim=c(-3,3), 
                  EmulatorType="Derivative", filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/1d Deriv/Joint/",
                  function_name="F5")

saving_var_plots(simple_em_out=simple_em_out, deriv_em_out=deriv_em_out, sigma=sigma,
                 xD=xD, xP=xP, filepath="/Users/helenbate/Documents/Year 4/Diss/Code/Images/1d Deriv/Joint/", function_name="F5",
                 plot_xD=TRUE)

########################################################
# 2d Emulation - defining f
########################################################

# Functions
f <- function(x) {
  sin(pi*x[,1]) + cos(pi*x[,2]+0.5)*(x[,1]/2)/(x[,2]+1)
}

gradf <- function(x) {
  g1 <- pi*cos(pi*x[,1]) + (0.5/(x[,2]+1)) * cos(pi*x[,2]+0.5)
  g2 <- - (x[,1]/2)/(x[,2]+1) * ((1/(x[,2]+1)) * cos(pi*x[,2]+0.5) + pi*sin(pi*x[,2]+0.5) )
  # DECREASE g3 by 1/2 ----- DO THIS!!!!!
  g3 <- - (1/2)/(x[,2]+1) * ( 1/(x[,2]+1) * cos(pi*x[,2]+0.5) + pi*sin(pi*x[,2]+0.5) )
  cbind(g1, g2, g3)
}

# following three plots not used
fxP_mat <- matrix(f(xP), nrow=length(x_grid), ncol=length(x_grid))
gradfxP_mat_list <- list("x1"=matrix(gradf(xP)[,1], nrow=length(x_grid), ncol=length(x_grid)),
                         "x2"=matrix(gradf(xP)[,2], nrow=length(x_grid), ncol=length(x_grid)))
emul_fill_cont(cont_mat=fxP_mat, xD=xD, x_grid=x_grid, 
               color.palette=exp_cols, 
               main="True Function f(x)")
emul_fill_cont(cont_mat=gradfxP_mat_list[[1]], xD=xD, x_grid=x_grid, 
               color.palette=exp_cols, 
               main="True df/dx_1")
emul_fill_cont(cont_mat=gradfxP_mat_list[[2]], xD=xD, x_grid=x_grid, 
               color.palette=exp_cols, 
               main="True df/dx_2")

########################################################
# 2d Emulation - set up
########################################################

# Generating 2d plots for poster and presentation
n <- 12
# xD <- 2*best_lhd(n,M=100)
xD <- matrix(c(1.50000000, 1.5000000,
               1.91666667, 0.7500000,
               0.41666667, 0.4166667,
               1.20000000, 0.5833333,
               0.75000000, 0.7500000,
               0.25000000, 0.9166667,
               1.75000000, 1.0833333,
               1.25000000, 1.2500000,
               0.58333333, 1.4166667,
               0.08333333, 1.5833333,
               0.91666667, 1.9166667,
               1.90000000, 1.9166667), ncol=2, byrow=TRUE)
x1_lim <- c(0,2)
x2_lim <- c(0,2)

# following plot not used
nD <- nrow(xD)
plot(xD, xlim=x1_lim, ylim=x2_lim, pch=16, xaxs="i", yaxs="i", col="blue", 
     xlab="x1", ylab="x2", cex=1.4)
abline(v=x1_lim[1]+(x1_lim[2]-x1_lim[1])*(0:nD)/nD, col="grey60")
abline(h=x2_lim[1]+(x2_lim[2]-x2_lim[1])*(0:nD)/nD, col="grey60")


# Emulation space
x_grid <- seq(-0.001, 2.001, len=70)
xP <- as.matrix(expand.grid("x1"=x_grid, "x2"=x_grid))


# Initialisers
theta <- c(0.4, 0.5)
sigma <- 0.5
delta <- 1e-06
cont_levs_mat <- matrix(c(-1.5,1.5,0.2,
                          0,sigma^2,sigma^2/10,
                          -9,9,0.5,
                          0,13,1,
                          -5,5,0.5,
                          0,5,0.5), ncol=3, byrow=TRUE)

########################################################
# 2d Emulation
########################################################

# Simple Emulator
D <- f(xD)
E_f <- 0

em_out_simple <- simple_BL_emulator_2d_fast(xP=xP, xD=xD, D=D, 
                                            theta=theta, sigma=sigma, E_f=E_f)

# emul_fill_cont(cont_mat=matrix(em_out_simple[,1], nrow=length(x_grid), ncol=length(x_grid)),
#                cont_levs=seq(cont_levs_mat[1,1], cont_levs_mat[1,2], cont_levs_mat[1,3]), 
#                x_grid=x_grid, color.palette=exp_cols, plot_xD=TRUE, xD=xD)
# emul_fill_cont(cont_mat=matrix(em_out_simple[,2], nrow=length(x_grid), ncol=length(x_grid)),
#                cont_levs=seq(cont_levs_mat[2,1], cont_levs_mat[2,2], cont_levs_mat[2,3]), 
#                x_grid=x_grid, color.palette=var_cols, plot_xD=TRUE, xD=xD)

ED_fx_mat <- matrix(em_out_simple[,"ExpD_f.x."], nrow=length(x_grid), ncol=length(x_grid)) 
ED_fx_plot <- emul_fill_cont(cont_mat=ED_fx_mat, 
                             cont_levs=seq(cont_levs_mat[1,1],cont_levs_mat[1,2], cont_levs_mat[1,3]), 
                             xD=xD, x_grid=x_grid, color.palette=exp_cols, 
                             main=bquote("Emulator Adjusted Expectation " ~ E[D] * "[f(x)]" ~ "\n (" *
                                           theta == .(theta[1]) * "," ~ .(theta[2]) * "), " ~
                                           sigma == .(sigma)
                             ))

saving_fxD_plots(em_out=em_out_simple, xD=xD, x_grid=x_grid, 
                 cont_levs_mat=cont_levs_mat, wd=690,
                 filepath="Documents/Year 4/Diss/Code/Images/Simple/",
                 function_name="F", EmulatorType="Simple")

# Derivative Emulator
D <- c(f(xD), gradf(xD)[,1], gradf(xD)[,2])

E_f <- 0
E_dfdx1 <- 0
E_dfdx2 <- 0

em_out_deriv <- deriv_BL_emulator_2d_fast(xP=xP, xD=xD, D=D, 
                                          theta=theta, sigma=sigma, E_f=E_f, E_dfdx=c(E_dfdx1, E_dfdx2))

saving_fxD_plots(em_out=em_out_deriv, xD=xD, x_grid=x_grid, 
                 cont_levs_mat=cont_levs_mat, wd=690,
                 filepath="Documents/Year 4/Diss/Code/Images/Deriv/",
                 function_name="F", EmulatorType="Derivative")

# Boundary Emulator (single boundary)
E_f <- 0

boundary <- list("x1"=NULL, "x2"=0.1)
em_out_boundary <- boundary_BL_emulator_fast(xP=xP, 
                                             boundary=boundary,
                                             theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=em_out_boundary, x_grid=x_grid, wd=690,
                                cont_levs_mat=cont_levs_mat, boundary=boundary,
                                filepath="Documents/Year 4/Diss/Code/Images/Boundary/",
                                function_name="F", EmulatorType="Boundary")

# Boundary with Runs Emulator (single boundary)
D <- f(xD)
E_f <- 0

boundary <- boundary
em_out_boundary_wR <- single_boundary_with_runs_BL_emulator(xP=xP, xD=xD, D=D,
                                                            boundary=boundary,
                                                            theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=em_out_boundary_wR, x_grid=x_grid, xD=xD,
                                cont_levs_mat=cont_levs_mat, boundary=boundary, wd=690,
                                filepath="Documents/Year 4/Diss/Code/Images/Boundary with Runs/",
                                function_name="F", EmulatorType="Boundary with Runs")

# Boundary with Derivative on Boundary Emulator (single boundary)
E_f <- 0

boundary <- boundary
em_out_BwD <- single_BwD_BL_emulator(xP=xP, 
                                     boundary=boundary,
                                     theta=theta, sigma=sigma, E_f=E_f, E_df=0)

saving_boundary_with_runs_plots(em_out=em_out_BwD, x_grid=x_grid, wd=690,
                                cont_levs_mat=cont_levs_mat, boundary=boundary,
                                filepath="Documents/Year 4/Diss/Code/Images/BwD/",
                                function_name="F", EmulatorType="Boundary and Derivative on the Boundary")

# BwD and Runs Emulator (single boundary)
E_f <- 0

boundary <- boundary
em_out_BwDR <- single_BwD_with_runs_BL_emulator(xP=xP, xD=xD, D=D,
                                                boundary=boundary, 
                                                theta=theta, sigma=sigma, E_f=E_f, E_df=0)

saving_boundary_with_runs_plots(em_out=em_out_BwDR, xD=xD, x_grid=x_grid, wd=690,
                                cont_levs_mat=cont_levs_mat, boundary=boundary,
                                filepath="Documents/Year 4/Diss/Code/Images/BwD and Runs/",
                                function_name="F", EmulatorType="Boundary and Derivative on the Boundary with Runs")


########################################################
# 2d VARIANCE COMPARISON PLOTS
########################################################

################### simple vs deriv ################### 
number_to_av <- 500
M <- 35
n <- 12
n_lower <- 4
seed <- 10 # 1000 quite good for number_to_av=500

E_dfdx <- c(E_dfdx1, E_dfdx2)

# deriv_equilivator <- deriv_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
#                                               x_grid=x_grid, theta=theta, sigma=sigma, 
#                                               x1_lim=x1_lim, x2_lim=x2_lim, seed=seed, 
#                                               E_f=E_f, E_dfdx=E_dfdx,  
#                                               print_headline=FALSE, print_summary_table=FALSE,
#                                               return_summary_table=TRUE, return_results=TRUE)

# plotting_var_comp_boxplot(deriv_equilivator, type="Deriv")
# plotting_var_comp_mean_error_plot(deriv_equilivator, type="Deriv")

saving_var_comp_plots(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                      x_grid=x_grid, theta=c(0.4, 0.5), sigma=0.5, 
                      x1_lim=c(0,2), x2_lim=c(0,2), seed=seed, wd=690, ht=520,
                      filepath="Documents/Year 4/Diss/Code/Images/Var Comp Simple Deriv/",
                      function_name="F", type="Deriv")

################### simple vs BwR ################### 
n <- 12
boundary <- list("x1"=NULL, "x2"=0.1)

# BwR Emulator Variance examples for LHD on [0,2]x[0.5, 2]
# seeds = 14, 121, 124, 143
seed <- 14
set.seed(seed)
xD <- best_lhd(n, M=M, x1_lim=c(0,2), x2_lim=c(0.5,2), print_switch=FALSE)
em_out_boundary <- single_boundary_with_runs_BL_emulator(xP=xP, xD=xD, D=f(xD),
                                                         boundary=boundary,
                                                         theta=theta, sigma=sigma, E_f=E_f)
saving_boundary_with_runs_plots(em_out=em_out_boundary, xD=xD, x_grid=x_grid, 
                                boundary=list("x1"=NULL, "x2"=0.1),
                                wd=690, ht=620,
                                cont_levs_mat=cont_levs_mat,
                                filepath=paste0("Documents/Year 4/Diss/Code/Images/Var Comp Simple BwR/Variance Examples/DISS/seed=", seed, "/"),
                                function_name=paste0("F_seed_", seed),
                                EmulatorType="Boundary with Runs") 

# Boxplots and mean with error bars average variance plots
number_to_av <- 500
M <- 35
n <- 12
n_lower <- 4

# boundary_equilivator <- single_boundary_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
#                                                            x_grid=x_grid, theta=theta, sigma=sigma, 
#                                                            x1_lim=x1_lim, x2_lim=c(0,2), 
#                                                            x1_lim_2=x1_lim, x2_lim_2=c(0.5,2), seed=110, # Was 45 with single x2_lim
#                                                            E_f=E_f, boundary=list("x1"=NULL, "x2"=0.1),  
#                                                            print_headline=TRUE, print_summary_table=FALSE,
#                                                            return_summary_table=TRUE, return_results=TRUE)
# plotting_var_comp_boxplot(boundary_equilivator, type="BwR")
# plotting_var_comp_mean_error_plot(boundary_equilivator, type="BwR")

saving_var_comp_plots(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                      x_grid=x_grid, theta=c(0.4, 0.5), sigma=0.5, 
                      boundary=list("x1"=NULL, "x2"=0.1),
                      x1_lim=c(0,2), x2_lim=c(0,2), 
                      x1_lim_2=c(0,2), x2_lim_2=c(0.5,2), seed=45, wd=690, ht=520,
                      filepath="Documents/Year 4/Diss/Code/Images/Var Comp Simple BwR/",
                      function_name="F", type="BwR")


################### Simple / BwR vs BwDR ################### 
n <- 12
boundary <- list("x1"=NULL, "x2"=0.1)
x2_lim_2 <- c(0.6,2)

# BwDR Emulator Variance examples for LHD on [0,2]x[0.6, 2]
# seeds = 123, 124, 143, 444
seed <- 444
set.seed(seed)
xD <- best_lhd(n, M=M, x1_lim=c(0,2), x2_lim=x2_lim, print_switch=FALSE)
em_out_BwDR <- single_BwD_with_runs_BL_emulator(xP=xP, xD=xD, D=f(xD),
                                                boundary=boundary,
                                                theta=theta, sigma=sigma, E_f=E_f)
saving_boundary_with_runs_plots(em_out=em_out_BwDR, xD=xD, x_grid=x_grid, 
                                boundary=list("x1"=NULL, "x2"=0.1),
                                wd=690, ht=620,
                                cont_levs_mat=cont_levs_mat,
                                filepath=paste0("Documents/Year 4/Diss/Code/Images/Var Comp BwR BwDR/Variance Examples/seed=", seed, "/"),
                                function_name=paste0("F_seed_", seed),
                                EmulatorType="Boundary and Derivative on the Boundary with Runs") 

# Boxplots and mean with error bars average variance plots
number_to_av <- 500
M <- 35
n <- 12
n_lower <- 4
# BwDR_equilivator <- simple_BwDR_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
#                                                    x_grid=x_grid, theta=theta, sigma=sigma, 
#                                                    x1_lim=x1_lim, x2_lim=c(0.6,2), seed=45, 
#                                                    E_f=E_f, boundary=list("x1"=NULL, "x2"=0.1),  
#                                                    print_headline=FALSE, print_summary_table=FALSE,
#                                                    return_summary_table=TRUE, return_results=TRUE)
# plotting_var_comp_boxplot(BwDR_equilivator, type="Simple-BwDR")
# plotting_var_comp_mean_error_plot(BwDR_equilivator, type="BwDR")

# simple vs BwDR
saving_var_comp_plots(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                      x_grid=x_grid, theta=c(0.4, 0.5), sigma=0.5, 
                      boundary=list("x1"=NULL, "x2"=0.1),
                      x1_lim=c(0,2), x2_lim=c(0,2), 
                      x1_lim_2=c(0,2), x2_lim_2=c(0.6,2), seed=100, # 14 used with (100)
                      wd=690, ht=550,
                      filepath="Documents/Year 4/Diss/Code/Images/Var Comp Simple BwDR/",
                      function_name="F", type="Simple-BwDR")

# BwR vs BwDR
saving_var_comp_plots(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                      x_grid=x_grid, theta=c(0.4, 0.5), sigma=0.5, 
                      boundary=list("x1"=NULL, "x2"=0.1),
                      x1_lim=c(0,2), x2_lim=c(0.5,2), 
                      x1_lim_2=c(0,2), x2_lim_2=c(0.6,2), seed=1, wd=690, ht=550,
                      filepath="Documents/Year 4/Diss/Code/Images/Var Comp BwR BwDR/",
                      function_name="F", type="BwR-BwDR")

########################################################
# 2d BwD and BwDR Emulation with 2 Perpendicular Known Boundaries
########################################################

x_grid <- seq(0, 2, len=70)
xP <- as.matrix(expand.grid("x1"=x_grid, "x2"=x_grid))

E_f <- 0
E_df <- c(0, 0, 0)

boundary <- list("x1"=1, "x2"=0.1)

# Perpendicular 2 Boundary Emulator
em_out_boundary <- boundary_BL_emulator_fast(xP=xP, 
                                             boundary=boundary,
                                             theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=em_out_boundary, x_grid=x_grid, wd=690,
                                cont_levs_mat=cont_levs_mat, boundary=boundary,
                                filepath="Documents/Year 4/Diss/Code/Images/2 Boundary/",
                                function_name="F", EmulatorType="Boundary")

# Perpendicular 2 Boundary with Derivative Emulator
perpBwD_em_out <- perp_BwD_BL_emulator(xP=xP, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)

saving_boundary_with_runs_plots(em_out=perpBwD_em_out, x_grid=x_grid, wd=690, 
                                cont_levs_mat=cont_levs_mat, boundary=boundary, 
                                filepath="Documents/Year 4/Diss/Code/Images/2 BwD/",
                                function_name="F", EmulatorType="Boundary and Derivative on the Boundary")

# Perpendicular 2 Boundary with Runs Emulator
em_out_BwR <- perpendicular_boundary_with_runs_BL_emulator(xP=xP, xD=xD, D=D,
                                                           boundary=boundary,
                                                           theta=theta, sigma=sigma, E_f=E_f)

saving_boundary_with_runs_plots(em_out=em_out_BwR, x_grid=x_grid, xD=xD,
                                cont_levs_mat=cont_levs_mat, boundary=boundary, wd=690,
                                filepath="Documents/Year 4/Diss/Code/Images/2 BwR/",
                                function_name="F", EmulatorType="Boundary with Runs")

# Perpendicular Boundary with Derivative and Runs Emulator
perpBwDR_em_out <- perp_BwDR_BL_emulator(xP=xP, xD=xD, D=D, boundary=boundary, theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)

saving_boundary_with_runs_plots(em_out=perpBwDR_em_out, xD=xD, x_grid=x_grid, wd=690,
                                cont_levs_mat=cont_levs_mat, boundary=boundary, 
                                filepath="Documents/Year 4/Diss/Code/Images/2 BwDR/",
                                function_name="F", EmulatorType="Boundary and Derivative on the Boundary with Runs",)

