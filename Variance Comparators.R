# Code to compare the number of runs required to achieve a similar level of 
# average emulator uncertainty between:
# the Simple and Derivative Emulators, 
# the Simple and Boundary with Runs Emulators, 
# the Simple and Boundary with Derivative and Runs Emulators and 
# the Boundary with Runs and Boundary with Derivative and Runs Emulators

# Comparison Method:
# 1) Create for loop for simple emulator with n points
# 2) Calculate average variance for simple emulator with this n
# 3) Create for loop for derivative emulator with n, n-1, n-2, ..., n_lower points
# 4) Calculate average variance for derivative emulator with these m points
# 5) Compare which m has closest average variance to the simple emulator with n points

# Note: Must define model function f and model derivative function gradf (if 
# using derivatives) and x_grid (grid of the emulator evaluation space) before 
# running

# Possible outputs:
# print_headline = TRUE        =>  prints the number of runs for which the mean average alternative emulator variance is closest to that of the control with n runs and the absolute difference between these
# print_summary_table = TRUE   =>  prints a summary table for each emulator type and number of run combination calculated
# return_summary_table = TRUE  =>  returns the above summary table (usually set print_result=NULL in this case)
# return_results = TRUE        =>  returns a data.frame of the raw stored average variances calculated for each iteration of each emulator type and number of run combindation

# Example:
# f <- function(x) {
#   sin(2*pi*x[,1]) + cos(pi*x[,2]+0.5)*x[,1]/(x[,2]+1)
# }
# 
# gradf <- function(x) {
#   g1 <- 2*pi*cos(2*pi*x[,1]) + (1/(x[,2]+1)) * cos(pi*x[,2]+0.5)
#   g2 <- - x[,1]/(x[,2]+1) * ((1/(x[,2]+1)) * cos(pi*x[,2]+0.5) + pi*sin(pi*x[,2]+0.5) )
#   cbind(g1, g2)
# }
# 
# x_grid <- seq(-0.001, 2.001, len=70)
# deriv_av_var_equilivator(number_to_av=50, M=50, x_grid=x_grid, theta=c(0.2,0.4), sigma=0.5, E_f=0, E_dfdx=c(0,0),n=8, n_lower=4, x1_lim=c(0,2), x2_lim=c(0,2))
# Simple Emulator (n=8) and Derivative Emulator (n=5) Variances differ by only ~ 0.00336357

source("./LHD Generator Functions.R")
source("./Simple Emulator Code.R")
source("./Derivative Emulator Code.R")
source("./Boundary Emulators (without runs) Code.R")
source("./Boundary with Runs Emulators Code.R")
source("./Boundary with Derivative Emulators (without runs) Code.R")
source("./Boundary with Derivative and Runs Emulators Code.R")
source("./Variance Comparators.R")

################################################################
# Average Variance function
################################################################

av_var_func <- function(VarD_fx_mat){
  mean(diag(VarD_fx_mat))
}

################################################################
# Simple vs Derivative Variance Comparator
################################################################

deriv_av_var_equilivator <- function(number_to_av,                        # the number of LHDs for which to store the average variance of each emulator
                                     M=50,                                # the number of LHDs to check before choosing one to simulate emulator for (see best_lhd() from "LHD Generator Functions.R")
                                     x_grid,                              # the 1d grid of input points to expand for set of points for which we evaluate each emulator (expected square input space with equal scales)
                                     theta,                               # the correlation length vector
                                     sigma,                               # the prior standard deviation
                                     E_f,                                 # the prior expectation for f(x) 
                                     E_dfdx,                              # the prior expectation for df/dx(x), (df/dx_1(x), df/dx_2(x))
                                     n,                                   # the number of runs the control emulator uses (we place using LHD)
                                     n_lower=NULL,                        # the minimum number of runs for which summarise the average variance for the alternative emulator
                                     x1_lim=c(0,1), x2_lim=c(0,1),        # the x1- and x2-input space limits for LHD generation of the control emulator (see best_lhd() from "LHD Generator Functions.R")
                                     x1_lim_2=NULL, x2_lim_2=NULL,        # (optional) the x1- and x2-input space limits for LHD generation of the alternative emulator (if different from that of the control)
                                     seed=NULL,                           # (optional) run with set.seed() for LHD generation to ensure exact reproducibility
                                     print_headline=TRUE,                 # prints the number of runs for which the mean average alternative emulator variance is closest to that of the control with n runs and the absolute difference between these
                                     print_summary_table=FALSE,           # prints a summary table for each emulator type and number of run combination calculated
                                     return_summary_table=FALSE,          # returns the above summary table (usually set print_result=NULL in this case)
                                     return_results=FALSE                 # returns a data.frame of the raw stored average variances calculated for each iteration of each emulator type and number of run combindation
){
  # Simple emulator loop
  number_to_av <- number_to_av
  n <- n
  xP <- as.matrix(expand.grid("x1"=x_grid, "x2"=x_grid))
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(is.null(x1_lim_2)){
    x1_lim_2=x1_lim
  }
  if(is.null(x2_lim_2)){
    x2_lim_2=x2_lim
  }
  
  ############# Simple Loop #############
  simple_variance_storage <- rep(0,number_to_av)
  
  for(i in 1:number_to_av){
    xD <- best_lhd(n, M=M, x1_lim=x1_lim, x2_lim=x2_lim, print_switch=FALSE)
    D <- f(xD)
    em_out_simple <- simple_BL_emulator_2d_fast(xP=xP, xD=xD, D=D, 
                                                theta=theta, sigma=sigma, E_f=E_f)
    
    VarD_fx_mat <- matrix(em_out_simple[,"VarD_f.x."], nrow=length(x_grid), ncol=length(x_grid)) 
    simple_variance_storage[i] <- av_var_func(VarD_fx_mat)
  }
  
  av_var <- mean(simple_variance_storage)
  sd_var <- sd(simple_variance_storage)
  Q1_var <- quantile(simple_variance_storage, 0.25)
  Median_var <- median(simple_variance_storage)
  Q3_var <- quantile(simple_variance_storage, 0.75)
  simple_av_var <- list("n"=n,
                        "av_var"=av_var,
                        "sd_var"=sd_var,
                        "CI_lower_3sd_var"=max(0, av_var - 3*sd_var),
                        "CI_upper_3sd_var"=av_var + 3*sd_var,
                        "CI_lower_1sd_var"=max(0, av_var - sd_var),
                        "CI_upper_1sd_var"=av_var + sd_var,
                        "Q1"=Q1_var[[1]],
                        "Median"=Median_var,
                        "Q3"=Q3_var[[1]])
  simple_av_var
  
  if(is.null(n_lower)){
    n_seq <- seq(n,floor(0.5*n),by=-1)
  }else{
    n_seq <- seq(n, n_lower, by=-1)
  }
  
  ############# Derivative Loop #############
  # separate since for n_seq[i] < n LHD needs to be regnerated for fewer runs anyway
  
  deriv_variance_storage <- matrix(0, nrow=length(n_seq), ncol=number_to_av)
  # rows = sample size n, n-1, n-2, ..., 1
  # cols = iteration i=1, 2, 3, ..., number_to_av
  
  for(j in 1:length(n_seq)){
    for(i in 1:number_to_av){
      xD <- best_lhd(n_seq[j], M=M, x1_lim=x1_lim_2, x2_lim=x2_lim_2, print_switch=FALSE)
      D <- c(f(xD), gradf(xD)[,1], gradf(xD)[,2])
      em_out_deriv <- deriv_BL_emulator_2d_fast(xP=xP, xD=xD, D=D, theta=theta, 
                                                sigma=sigma, E_f=E_f, E_dfdx=E_dfdx)
      
      VarD_fx_mat <- matrix(em_out_deriv[,"VarD_f.x."], nrow=length(x_grid), ncol=length(x_grid)) 
      deriv_variance_storage[j,i] <- av_var_func(VarD_fx_mat)
    }
    print(j)
  }
  
  av_var <- rowMeans(deriv_variance_storage)
  sd_var <- apply(deriv_variance_storage, 1, sd)
  Q1_var <- apply(deriv_variance_storage, 1, function(x) quantile(x, 0.25)[[1]])
  Median_var <- apply(deriv_variance_storage, 1, median)
  Q3_var <- apply(deriv_variance_storage, 1, function(x) quantile(x, 0.75)[[1]])
  deriv_av_var <- list("n"=n_seq,
                       "av_var"=av_var,
                       "sd_var"=sd_var,
                       "CI_lower_3sd_var"=pmax(0, av_var - 3*sd_var),
                       "CI_upper_3sd_var"=av_var + 3*sd_var,
                       "CI_lower_1sd_var"=pmax(0, av_var - sd_var),
                       "CI_upper_1sd_var"=av_var + sd_var,
                       "Q1"=Q1_var,
                       "Median"=Median_var,
                       "Q3"=Q3_var)
  deriv_av_var
  
  ######### SUMMARY TABLE ##########
  simple_heading <- paste0("Simple_", n)
  deriv_headings <- rep(NA, length(n_seq))
  for(i in 1:length(n_seq)){
    deriv_headings[i] <- paste0("Deriv_", n_seq[i])
  }
  
  CI_3sd_contains_simple_mean <- (simple_av_var$av_var >= deriv_av_var$CI_lower_3sd_var & 
                                    simple_av_var$av_var <= deriv_av_var$CI_upper_3sd_var)
  CI_1sd_contains_simple_mean <- (simple_av_var$av_var >= deriv_av_var$CI_lower_1sd_var & 
                                    simple_av_var$av_var <= deriv_av_var$CI_upper_1sd_var)
  CI_1sd_overlaps_simple <- (simple_av_var$CI_lower_1sd_var >= deriv_av_var$CI_lower_1sd_var & 
                               simple_av_var$CI_lower_1sd_var <= deriv_av_var$CI_upper_1sd_var) |
    (simple_av_var$CI_upper_1sd_var >= deriv_av_var$CI_lower_1sd_var & 
       simple_av_var$CI_upper_1sd_var <= deriv_av_var$CI_upper_1sd_var) |
    (simple_av_var$CI_lower_1sd_var <= deriv_av_var$CI_lower_1sd_var &
       simple_av_var$CI_upper_1sd_var >= deriv_av_var$CI_upper_1sd_var)
  IQR_contains_simple_mean <- (simple_av_var$av_var >= deriv_av_var$Q1 & 
                                 simple_av_var$av_var <= deriv_av_var$Q3)
  IQR_contains_simple_median <- (simple_av_var$Median >= deriv_av_var$Q1 & 
                                   simple_av_var$Median <= deriv_av_var$Q3)
  IQR_overlaps_simple <- (simple_av_var$Q1 >= deriv_av_var$Q1 & 
                            simple_av_var$Q1 <= deriv_av_var$Q3) |
    (simple_av_var$Q3 >= deriv_av_var$Q1 & 
       simple_av_var$Q3 <= deriv_av_var$Q3) |
    (simple_av_var$Q1 <= deriv_av_var$Q1 &
       simple_av_var$Q3 >= deriv_av_var$Q3)
  
  summary_table <- data.frame(num_of_runs=c(simple_av_var$n, deriv_av_var$n), 
                              mean_average_var=c(simple_av_var$av_var, deriv_av_var$av_var),
                              sd_average_var=c(simple_av_var$sd_var, deriv_av_var$sd_var),
                              CI_3sd_lower=c(simple_av_var$CI_lower_3sd_var, deriv_av_var$CI_lower_3sd_var),
                              CI_3sd_upper=c(simple_av_var$CI_upper_3sd_var, deriv_av_var$CI_upper_3sd_var),
                              CI_3sd_contains_simple=c(NA, CI_3sd_contains_simple_mean),
                              CI_1sd_lower=c(simple_av_var$CI_lower_1sd_var, deriv_av_var$CI_lower_1sd_var),
                              CI_1sd_upper=c(simple_av_var$CI_upper_1sd_var, deriv_av_var$CI_upper_1sd_var),
                              CI_1sd_contains_simple=c(NA, CI_1sd_contains_simple_mean),
                              CI_1sd_overlaps_simple=c(NA, CI_1sd_overlaps_simple),
                              Median=c(simple_av_var$Median, deriv_av_var$Median),
                              Q1=c(simple_av_var$Q1, deriv_av_var$Q1),
                              Q3=c(simple_av_var$Q3, deriv_av_var$Q3),
                              IQR_contains_simple_mean=c(NA, IQR_contains_simple_mean),
                              IQR_contains_simple_median=c(NA, IQR_contains_simple_median),
                              IQR_overlaps_simple=c(NA, IQR_overlaps_simple),
                              row.names = c(simple_heading, deriv_headings))
  
  # colnames(summary_table) <- c("Number of Runs",
  #                   "Mean of Average Emulator Adjusted Variance",
  #                   "SD of Average Emulator Adjusted Variance",
  #                   "CI: lower bound",
  #                   "CI: upper bound",
  #                   "CI contains Simple Mean Average Adj Var: Yes/No",
  #                   "Median Average Emulator Adjusted Variance",
  #                   "Q1",
  #                   "Q3",
  #                   "IQR contains Simple Mean Average Adj Var: Yes/No",
  #                   "IQR contains Simple Median Average Adj Var: Yes/No",
  #                   "IQR overlaps with Simple Average Adj Var IQR: Yes/No")
  
  ######### RESULTS matrix ##########
  
  results_mat <- cbind(simple_variance_storage, t(deriv_variance_storage))
  results <- as.data.frame(results_mat)
  colnames(results) <- c(simple_heading, deriv_headings)
  
  ############# OUTPUT #############
  
  if(print_headline==TRUE){
    
    index <- which.min((simple_av_var$av_var - deriv_av_var$av_var)^2)
    cat("Simple Emulator Average Variance for LHD size", simple_av_var$n, ":", simple_av_var$av_var,
        "\nDerivative Emulator Average Variance for sizes: \n") 
    print(t(as.matrix(data.frame(av_var=deriv_av_var$av_var, row.names=as.character(deriv_av_var$n)))))
    cat("n which gives the closest Derivative Emulator average variance:", deriv_av_var$n[index],
        ", \nDifference between Simple Emulator Variance with n =", simple_av_var$n,
        " and Derivative Simple Emulator Variance with n =", deriv_av_var$n[index], 
        ":", abs(simple_av_var$av_var - deriv_av_var$av_var[index]))
    
  }
  if(print_summary_table==TRUE){
    
    # colnames(summary_table) <- c("Number of Runs",
    #                   "Mean of Average Emulator Adjusted Variance",
    #                   "SD of Average Emulator Adjusted Variance",
    #                   "CI: lower bound",
    #                   "CI: upper bound",
    #                   "CI contains Simple Mean Average Adj Var: Yes/No",
    #                   "Median Average Emulator Adjusted Variance",
    #                   "Q1",
    #                   "Q3",
    #                   "IQR contains Simple Mean Average Adj Var: Yes/No",
    #                   "IQR contains Simple Median Average Adj Var: Yes/No",
    #                   "IQR overlaps with Simple Average Adj Var IQR: Yes/No")
    
    print(summary_table)
    View(summary_table)
    # Median and IQR comparisons more helpful here as CIs too broad to be helpful
    
  }
  if(return_summary_table & return_results){
    return(list("summary_table"=summary_table,
                "results"=results))
  }else if(return_summary_table){
    return(summary_table)
  }else if(return_results){
    return(results)
  }
}





################################################################
# Simple vs Single Boundary with Runs Emulator Variance Comparator
################################################################

library(pdist)

single_boundary_av_var_equilivator <- function(number_to_av,                        # the number of LHDs for which to store the average variance of each emulator
                                               M=50,                                # the number of LHDs to check before choosing one to simulate emulator for (see best_lhd() from "LHD Generator Functions.R")
                                               x_grid,                              # the 1d grid of input points to expand for set of points for which we evaluate each emulator (expected square input space with equal scales)
                                               theta,                               # the correlation length vector
                                               sigma,                               # the prior standard deviation
                                               E_f,                                 # the prior expectation for f(x) 
                                               boundary=list("x1"=0, "x2"=NULL), # boundary definition: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL)
                                               n,                                   # the number of runs the control emulator uses (we place using LHD)
                                               n_lower=NULL,                        # the minimum number of runs for which summarise the average variance for the alternative emulator
                                               x1_lim=c(0,1), x2_lim=c(0,1),        # the x1- and x2-input space limits for LHD generation of the control emulator (see best_lhd() from "LHD Generator Functions.R")
                                               x1_lim_2=NULL, x2_lim_2=NULL,        # (optional) the x1- and x2-input space limits for LHD generation of the alternative emulator (if different from that of the control)
                                               seed=NULL,                           # (optional) run with set.seed() for LHD generation to ensure exact reproducibility
                                               print_headline=TRUE,                 # prints the number of runs for which the mean average alternative emulator variance is closest to that of the control with n runs and the absolute difference between these
                                               print_summary_table=FALSE,           # prints a summary table for each emulator type and number of run combination calculated
                                               return_summary_table=FALSE,          # returns the above summary table (usually set print_result=NULL in this case)
                                               return_results=FALSE                 # returns a data.frame of the raw stored average variances calculated for each iteration of each emulator type and number of run combindation
){
  # Simple emulator loop
  number_to_av <- number_to_av
  n <- n
  xP <- as.matrix(expand.grid("x1"=x_grid, "x2"=x_grid))
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(is.null(x1_lim_2)){
    x1_lim_2=x1_lim
  }
  if(is.null(x2_lim_2)){
    x2_lim_2=x2_lim
  }
  
  ############# Simple Loop #############
  simple_variance_storage <- rep(0, number_to_av)
  
  for(i in 1:number_to_av){
    xD <- best_lhd(n, M=M, x1_lim=x1_lim, x2_lim=x2_lim, print_switch=FALSE)
    D <- f(xD)
    em_out_simple <- simple_BL_emulator_2d_fast(xP=xP, xD=xD, D=D, 
                                                theta=theta, sigma=sigma, E_f=E_f)
    
    VarD_fx_mat <- matrix(em_out_simple[,"VarD_f.x."], nrow=length(x_grid), ncol=length(x_grid)) 
    simple_variance_storage[i] <- av_var_func(VarD_fx_mat)
  }
  
  av_var <- mean(simple_variance_storage)
  sd_var <- sd(simple_variance_storage)
  Q1_var <- quantile(simple_variance_storage, 0.25)
  Median_var <- median(simple_variance_storage)
  Q3_var <- quantile(simple_variance_storage, 0.75)
  simple_av_var <- list("n"=n,
                        "av_var"=av_var,
                        "sd_var"=sd_var,
                        "CI_lower_3sd_var"=max(0, av_var - 3*sd_var),
                        "CI_upper_3sd_var"=av_var + 3*sd_var,
                        "CI_lower_1sd_var"=max(0, av_var - sd_var),
                        "CI_upper_1sd_var"=av_var + sd_var,
                        "Q1"=Q1_var[[1]],
                        "Median"=Median_var,
                        "Q3"=Q3_var[[1]])
  simple_av_var
  
  if(is.null(n_lower)){
    n_seq <- seq(n,floor(0.5*n),by=-1)
  }else{
    n_seq <- seq(n, n_lower, by=-1)
  }
  
  ############# Boundary Loop #############
  boundary_variance_storage <- matrix(0, nrow=length(n_seq), ncol=number_to_av)
  # rows = sample size n, n-1, n-2, ..., 1
  # cols = iteration i=1, 2, 3, ..., number_to_av
  
  for(j in 1:length(n_seq)){
    for(i in 1:number_to_av){
      xD <- best_lhd(n_seq[j], M=M, x1_lim=x1_lim_2, x2_lim=x2_lim_2, print_switch=FALSE)
      D <- f(xD)
      em_out_boundary <- single_boundary_with_runs_BL_emulator(xP=xP, xD=xD, D=D,
                                                               boundary=boundary,
                                                               theta=theta, sigma=sigma, E_f=E_f)
      
      VarD_fx_mat <- matrix(em_out_boundary[,2], nrow=length(x_grid), ncol=length(x_grid)) 
      boundary_variance_storage[j,i] <- av_var_func(VarD_fx_mat)
    }
  }
  
  av_var <- rowMeans(boundary_variance_storage)
  sd_var <- apply(boundary_variance_storage, 1, sd)
  Q1_var <- apply(boundary_variance_storage, 1, function(x) quantile(x, 0.25)[[1]])
  Median_var <- apply(boundary_variance_storage, 1, median)
  Q3_var <- apply(boundary_variance_storage, 1, function(x) quantile(x, 0.75)[[1]])
  boundary_av_var <- list("n"=n_seq,
                          "av_var"=av_var,
                          "sd_var"=sd_var,
                          "CI_lower_3sd_var"=pmax(0, av_var - 3*sd_var),
                          "CI_upper_3sd_var"=av_var + 3*sd_var,
                          "CI_lower_1sd_var"=pmax(0, av_var - sd_var),
                          "CI_upper_1sd_var"=av_var + sd_var,
                          "Q1"=Q1_var,
                          "Median"=Median_var,
                          "Q3"=Q3_var)
  
  ######### SUMMARY TABLE ##########
  simple_heading <- paste0("Simple_", n)
  boundary_headings <- rep(NA, length(n_seq))
  for(i in 1:length(n_seq)){
    boundary_headings[i] <- paste0("Boundary_", n_seq[i])
  }
  
  CI_3sd_contains_simple_mean <- (simple_av_var$av_var >= boundary_av_var$CI_lower_3sd_var & 
                                    simple_av_var$av_var <= boundary_av_var$CI_upper_3sd_var)
  CI_1sd_contains_simple_mean <- (simple_av_var$av_var >= boundary_av_var$CI_lower_1sd_var & 
                                    simple_av_var$av_var <= boundary_av_var$CI_upper_1sd_var)
  CI_1sd_overlaps_simple <- (simple_av_var$CI_lower_1sd_var >= boundary_av_var$CI_lower_1sd_var & 
                               simple_av_var$CI_lower_1sd_var <= boundary_av_var$CI_upper_1sd_var) |
    (simple_av_var$CI_upper_1sd_var >= boundary_av_var$CI_lower_1sd_var & 
       simple_av_var$CI_upper_1sd_var <= boundary_av_var$CI_upper_1sd_var) |
    (simple_av_var$CI_lower_1sd_var <= boundary_av_var$CI_lower_1sd_var &
       simple_av_var$CI_upper_1sd_var >= boundary_av_var$CI_upper_1sd_var)
  IQR_contains_simple_mean <- (simple_av_var$av_var >= boundary_av_var$Q1 & 
                                 simple_av_var$av_var <= boundary_av_var$Q3)
  IQR_contains_simple_median <- (simple_av_var$Median >= boundary_av_var$Q1 & 
                                   simple_av_var$Median <= boundary_av_var$Q3)
  IQR_overlaps_simple <- (simple_av_var$Q1 >= boundary_av_var$Q1 & 
                            simple_av_var$Q1 <= boundary_av_var$Q3) |
    (simple_av_var$Q3 >= boundary_av_var$Q1 & 
       simple_av_var$Q3 <= boundary_av_var$Q3) |
    (simple_av_var$Q1 <= boundary_av_var$Q1 &
       simple_av_var$Q3 >= boundary_av_var$Q3)
  
  summary_table <- data.frame(num_of_runs=c(simple_av_var$n, boundary_av_var$n), 
                              mean_average_var=c(simple_av_var$av_var, boundary_av_var$av_var),
                              sd_average_var=c(simple_av_var$sd_var, boundary_av_var$sd_var),
                              CI_3sd_lower=c(simple_av_var$CI_lower_3sd_var, boundary_av_var$CI_lower_3sd_var),
                              CI_3sd_upper=c(simple_av_var$CI_upper_3sd_var, boundary_av_var$CI_upper_3sd_var),
                              CI_3sd_contains_simple=c(NA, CI_3sd_contains_simple_mean),
                              CI_1sd_lower=c(simple_av_var$CI_lower_1sd_var, boundary_av_var$CI_lower_1sd_var),
                              CI_1sd_upper=c(simple_av_var$CI_upper_1sd_var, boundary_av_var$CI_upper_1sd_var),
                              CI_1sd_contains_simple=c(NA, CI_1sd_contains_simple_mean),
                              CI_1sd_overlaps_simple=c(NA, CI_1sd_overlaps_simple),
                              Median=c(simple_av_var$Median, boundary_av_var$Median),
                              Q1=c(simple_av_var$Q1, boundary_av_var$Q1),
                              Q3=c(simple_av_var$Q3, boundary_av_var$Q3),
                              IQR_contains_simple_mean=c(NA, IQR_contains_simple_mean),
                              IQR_contains_simple_median=c(NA, IQR_contains_simple_median),
                              IQR_overlaps_simple=c(NA, IQR_overlaps_simple),
                              row.names = c(simple_heading, boundary_headings))
  
  # colnames(summary_table) <- c("Number of Runs",
  #                   "Mean of Average Emulator Adjusted Variance",
  #                   "SD of Average Emulator Adjusted Variance",
  #                   "CI: lower bound",
  #                   "CI: upper bound",
  #                   "CI contains Simple Mean Average Adj Var: Yes/No",
  #                   "Median Average Emulator Adjusted Variance",
  #                   "Q1",
  #                   "Q3",
  #                   "IQR contains Simple Mean Average Adj Var: Yes/No",
  #                   "IQR contains Simple Median Average Adj Var: Yes/No",
  #                   "IQR overlaps with Simple Average Adj Var IQR: Yes/No")
  
  ######### RESULTS matrix ##########
  
  results_mat <- cbind(simple_variance_storage, t(boundary_variance_storage))
  results <- as.data.frame(results_mat)
  colnames(results) <- c(simple_heading, boundary_headings)
  
  ############# OUTPUT #############
  
  if(print_headline==TRUE){
    
    index <- which.min((simple_av_var$av_var - boundary_av_var$av_var)^2)
    cat("Simple Emulator Average Variance for LHD size", simple_av_var$n, ":", simple_av_var$av_var,
        "\nBoundary with Runs Emulator Average Variance for sizes: \n") 
    print(t(as.matrix(data.frame(av_var=boundary_av_var$av_var, row.names=as.character(boundary_av_var$n)))))
    cat("n which gives the closest Boundary with Runs Emulator average variance:", boundary_av_var$n[index],
        ", \nDifference between Simple Emulator Variance with n =", simple_av_var$n,
        " and Boundary with Runs Simple Emulator Variance with n =", boundary_av_var$n[index], 
        ":", abs(simple_av_var$av_var - boundary_av_var$av_var[index]))
    
  }
  if(print_summary_table==TRUE){
    print(summary_table)
    View(summary_table)
  }
  if(return_summary_table & return_results){
    return(list("summary_table"=summary_table,
                "results"=results))
  }else if(return_summary_table){
    return(summary_table)
  }else if(return_results){
    return(results)
  }
}





################################################################
# Simple vs Single Boundary with Derivative and Runs Emulator Variance Comparator
################################################################

library(pdist)

simple_BwDR_av_var_equilivator <- function(number_to_av,                        # the number of LHDs for which to store the average variance of each emulator
                                           M=50,                                # the number of LHDs to check before choosing one to simulate emulator for (see best_lhd() from "LHD Generator Functions.R")
                                           x_grid,                              # the 1d grid of input points to expand for set of points for which we evaluate each emulator (expected square input space with equal scales)
                                           theta,                               # the correlation length vector
                                           sigma,                               # the prior standard deviation
                                           E_f,                                 # the prior expectation for f(x) 
                                           E_df,                                # the prior expectation for df/dx(x) in the direction perpendicularly to the boundary
                                           boundary=list("x1"=0, "x2"=NULL),    # boundary definition: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL)
                                           n,                                   # the number of runs the control emulator uses (we place using LHD)
                                           n_lower=NULL,                        # the minimum number of runs for which summarise the average variance for the alternative emulator
                                           x1_lim=c(0,1), x2_lim=c(0,1),        # the x1- and x2-input space limits for LHD generation of the control emulator (see best_lhd() from "LHD Generator Functions.R")
                                           x1_lim_2=NULL, x2_lim_2=NULL,        # (optional) the x1- and x2-input space limits for LHD generation of the alternative emulator (if different from that of the control)
                                           seed=NULL,                           # (optional) run with set.seed() for LHD generation to ensure exact reproducibility
                                           print_headline=TRUE,                 # prints the number of runs for which the mean average alternative emulator variance is closest to that of the control with n runs and the absolute difference between these
                                           print_summary_table=FALSE,           # prints a summary table for each emulator type and number of run combination calculated
                                           return_summary_table=FALSE,          # returns the above summary table (usually set print_result=NULL in this case)
                                           return_results=FALSE                 # returns a data.frame of the raw stored average variances calculated for each iteration of each emulator type and number of run combindation
){
  # Simple emulator loop
  number_to_av <- number_to_av
  n <- n
  xP <- as.matrix(expand.grid("x1"=x_grid, "x2"=x_grid))
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(is.null(x1_lim_2)){
    x1_lim_2=x1_lim
  }
  if(is.null(x2_lim_2)){
    x2_lim_2=x2_lim
  }
  
  ############# Simple Loop #############
  simple_variance_storage <- rep(0, number_to_av)
  
  for(i in 1:number_to_av){
    xD <- best_lhd(n, M=M, x1_lim=x1_lim, x2_lim=x2_lim, print_switch=FALSE)
    em_out_simple <- simple_BL_emulator_2d_fast(xP=xP, xD=xD, D=f(xD), 
                                                theta=theta, sigma=sigma, E_f=E_f)
    
    VarD_fx_mat <- matrix(em_out_simple[,"VarD_f.x."], nrow=length(x_grid), ncol=length(x_grid)) 
    simple_variance_storage[i] <- av_var_func(VarD_fx_mat)
  }
  
  av_var <- mean(simple_variance_storage)
  sd_var <- sd(simple_variance_storage)
  Q1_var <- quantile(simple_variance_storage, 0.25)
  Median_var <- median(simple_variance_storage)
  Q3_var <- quantile(simple_variance_storage, 0.75)
  simple_av_var <- list("n"=n,
                        "av_var"=av_var,
                        "sd_var"=sd_var,
                        "CI_lower_3sd_var"=max(0, av_var - 3*sd_var),
                        "CI_upper_3sd_var"=av_var + 3*sd_var,
                        "CI_lower_1sd_var"=max(0, av_var - sd_var),
                        "CI_upper_1sd_var"=av_var + sd_var,
                        "Q1"=Q1_var[[1]],
                        "Median"=Median_var,
                        "Q3"=Q3_var[[1]])
  
  if(is.null(n_lower)){
    n_seq <- seq(n,floor(0.5*n),by=-1)
  }else{
    n_seq <- seq(n, n_lower, by=-1)
  }
  
  ############# BwD Loop #############
  BwDR_variance_storage <- matrix(0, nrow=length(n_seq), ncol=number_to_av)
  # rows = sample size n, n-1, n-2, ..., 1
  # cols = iteration i=1, 2, 3, ..., number_to_av
  
  for(j in 1:length(n_seq)){
    for(i in 1:number_to_av){
      xD <- best_lhd(n_seq[j], M=M, x1_lim=x1_lim_2, x2_lim=x2_lim_2, print_switch=FALSE)
      em_out_BwDR <- single_BwD_with_runs_BL_emulator(xP=xP, xD=xD, D=f(xD),
                                                      boundary=boundary,
                                                      theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)
      
      VarD_fx_mat <- matrix(em_out_BwDR[,2], nrow=length(x_grid), ncol=length(x_grid)) 
      BwDR_variance_storage[j,i] <- av_var_func(VarD_fx_mat)
    }
  }
  
  av_var <- rowMeans(BwDR_variance_storage)
  sd_var <- apply(BwDR_variance_storage, 1, sd)
  Q1_var <- apply(BwDR_variance_storage, 1, function(x) quantile(x, 0.25)[[1]])
  Median_var <- apply(BwDR_variance_storage, 1, median)
  Q3_var <- apply(BwDR_variance_storage, 1, function(x) quantile(x, 0.75)[[1]])
  BwDR_av_var <- list(     "n"=n_seq,
                           "av_var"=av_var,
                           "sd_var"=sd_var,
                           "CI_lower_3sd_var"=pmax(0, av_var - 3*sd_var),
                           "CI_upper_3sd_var"=av_var + 3*sd_var,
                           "CI_lower_1sd_var"=pmax(0, av_var - sd_var),
                           "CI_upper_1sd_var"=av_var + sd_var,
                           "Q1"=Q1_var,
                           "Median"=Median_var,
                           "Q3"=Q3_var)
  
  ######### SUMMARY TABLE ##########
  simple_heading <- paste0("Simple_", n)
  BwDR_headings <- rep(NA, length(n_seq))
  for(i in 1:length(n_seq)){
    BwDR_headings[i] <- paste0("BwDR_", n_seq[i])
  }
  
  CI_3sd_contains_simple_mean <- (simple_av_var$av_var >= BwDR_av_var$CI_lower_3sd_var & 
                                    simple_av_var$av_var <= BwDR_av_var$CI_upper_3sd_var)
  CI_1sd_contains_simple_mean <- (simple_av_var$av_var >= BwDR_av_var$CI_lower_1sd_var & 
                                    simple_av_var$av_var <= BwDR_av_var$CI_upper_1sd_var)
  CI_1sd_overlaps_simple <- (simple_av_var$CI_lower_1sd_var >= BwDR_av_var$CI_lower_1sd_var & 
                               simple_av_var$CI_lower_1sd_var <= BwDR_av_var$CI_upper_1sd_var) |
    (simple_av_var$CI_upper_1sd_var >= BwDR_av_var$CI_lower_1sd_var & 
       simple_av_var$CI_upper_1sd_var <= BwDR_av_var$CI_upper_1sd_var) |
    (simple_av_var$CI_lower_1sd_var <= BwDR_av_var$CI_lower_1sd_var &
       simple_av_var$CI_upper_1sd_var >= BwDR_av_var$CI_upper_1sd_var)
  IQR_contains_simple_mean <- (simple_av_var$av_var >= BwDR_av_var$Q1 & 
                                 simple_av_var$av_var <= BwDR_av_var$Q3)
  IQR_contains_simple_median <- (simple_av_var$Median >= BwDR_av_var$Q1 & 
                                   simple_av_var$Median <= BwDR_av_var$Q3)
  IQR_overlaps_simple <- (simple_av_var$Q1 >= BwDR_av_var$Q1 & 
                            simple_av_var$Q1 <= BwDR_av_var$Q3) |
    (simple_av_var$Q3 >= BwDR_av_var$Q1 & 
       simple_av_var$Q3 <= BwDR_av_var$Q3) |
    (simple_av_var$Q1 <= BwDR_av_var$Q1 &
       simple_av_var$Q3 >= BwDR_av_var$Q3)
  
  summary_table <- data.frame(num_of_runs=c(simple_av_var$n, BwDR_av_var$n), 
                              mean_average_var=c(simple_av_var$av_var, BwDR_av_var$av_var),
                              sd_average_var=c(simple_av_var$sd_var, BwDR_av_var$sd_var),
                              CI_3sd_lower=c(simple_av_var$CI_lower_3sd_var, BwDR_av_var$CI_lower_3sd_var),
                              CI_3sd_upper=c(simple_av_var$CI_upper_3sd_var, BwDR_av_var$CI_upper_3sd_var),
                              CI_3sd_contains_simple=c(NA, CI_3sd_contains_simple_mean),
                              CI_1sd_lower=c(simple_av_var$CI_lower_1sd_var, BwDR_av_var$CI_lower_1sd_var),
                              CI_1sd_upper=c(simple_av_var$CI_upper_1sd_var, BwDR_av_var$CI_upper_1sd_var),
                              CI_1sd_contains_simple=c(NA, CI_1sd_contains_simple_mean),
                              CI_1sd_overlaps_simple=c(NA, CI_1sd_overlaps_simple),
                              Median=c(simple_av_var$Median, BwDR_av_var$Median),
                              Q1=c(simple_av_var$Q1, BwDR_av_var$Q1),
                              Q3=c(simple_av_var$Q3, BwDR_av_var$Q3),
                              IQR_contains_simple_mean=c(NA, IQR_contains_simple_mean),
                              IQR_contains_simple_median=c(NA, IQR_contains_simple_median),
                              IQR_overlaps_simple=c(NA, IQR_overlaps_simple),
                              row.names = c(simple_heading, BwDR_headings))
  
  # colnames(summary_table) <- c("Number of Runs",
  #                   "Mean of Average Emulator Adjusted Variance",
  #                   "SD of Average Emulator Adjusted Variance",
  #                   "CI: lower bound",
  #                   "CI: upper bound",
  #                   "CI contains Simple Mean Average Adj Var: Yes/No",
  #                   "Median Average Emulator Adjusted Variance",
  #                   "Q1",
  #                   "Q3",
  #                   "IQR contains Simple Mean Average Adj Var: Yes/No",
  #                   "IQR contains Simple Median Average Adj Var: Yes/No",
  #                   "IQR overlaps with Simple Average Adj Var IQR: Yes/No")
  
  ######### RESULTS matrix ##########
  
  results_mat <- cbind(simple_variance_storage, t(BwDR_variance_storage))
  results <- as.data.frame(results_mat)
  colnames(results) <- c(simple_heading, BwDR_headings)
  
  ############# OUTPUT #############
  
  if(print_headline==TRUE){
    
    index <- which.min((simple_av_var$av_var - BwDR_av_var$av_var)^2)
    cat("Simple Emulator Average Variance for LHD size", simple_av_var$n, ":", simple_av_var$av_var,
        "\nBoundary with Derivative and Runs Emulator Average Variance for sizes: \n") 
    print(t(as.matrix(data.frame(av_var=BwDR_av_var$av_var, row.names=as.character(BwDR_av_var$n)))))
    cat("n which gives the closest Boundary with Derivative and Runs Emulator average variance:", BwDR_av_var$n[index],
        ", \nDifference between Simple Emulator Variance with n =", simple_av_var$n,
        " and Boundary with Derivative and Runs Simple Emulator Variance with n =", BwDR_av_var$n[index], 
        ":", abs(simple_av_var$av_var - BwDR_av_var$av_var[index]))
    
  }
  if(print_summary_table==TRUE){
    print(summary_table)
    View(summary_table)
  }
  if(return_summary_table & return_results){
    return(list("summary_table"=summary_table,
                "results"=results))
  }else if(return_summary_table){
    return(summary_table)
  }else if(return_results){
    return(results)
  }
}

################################################################
# Single Boundary with Runs vs Single Boundary with Derivative and Runs Variance Comparator
################################################################

library(pdist)

BwR_BwDR_av_var_equilivator <- function(number_to_av,                        # the number of LHDs for which to store the average variance of each emulator
                                        M=50,                                # the number of LHDs to check before choosing one to simulate emulator for (see best_lhd() from "LHD Generator Functions.R")
                                        x_grid,                              # the 1d grid of input points to expand for set of points for which we evaluate each emulator (expected square input space with equal scales)
                                        theta,                               # the correlation length vector
                                        sigma,                               # the prior standard deviation
                                        E_f,                                 # the prior expectation for f(x) 
                                        E_df,                                # the prior expectation for df/dx(x) in the direction perpendicularly to the boundary
                                        boundary=list("x1"=0, "x2"=NULL),    # boundary definition: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL)
                                        n,                                   # the number of runs the control emulator uses (we place using LHD)
                                        n_lower=NULL,                        # the minimum number of runs for which summarise the average variance for the alternative emulator
                                        x1_lim=c(0,1), x2_lim=c(0,1),        # the x1- and x2-input space limits for LHD generation of the control emulator (see best_lhd() from "LHD Generator Functions.R")
                                        x1_lim_2=NULL, x2_lim_2=NULL,        # (optional) the x1- and x2-input space limits for LHD generation of the alternative emulator (if different from that of the control)
                                        seed=NULL,                           # (optional) run with set.seed() for LHD generation to ensure exact reproducibility
                                        print_headline=TRUE,                 # prints the number of runs for which the mean average alternative emulator variance is closest to that of the control with n runs and the absolute difference between these
                                        print_summary_table=FALSE,           # prints a summary table for each emulator type and number of run combination calculated
                                        return_summary_table=FALSE,          # returns the above summary table (usually set print_result=NULL in this case)
                                        return_results=FALSE                 # returns a data.frame of the raw stored average variances calculated for each iteration of each emulator type and number of run combindation
){
  # BwR emulator loop
  number_to_av <- number_to_av
  n <- n
  xP <- as.matrix(expand.grid("x1"=x_grid, "x2"=x_grid))
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if(is.null(x1_lim_2)){
    x1_lim_2=x1_lim
  }
  if(is.null(x2_lim_2)){
    x2_lim_2=x2_lim
  }
  
  ############# BwR Loop #############
  BwR_variance_storage <- rep(0, number_to_av)
  
  for(i in 1:number_to_av){
    xD <- best_lhd(n, M=M, x1_lim=x1_lim, x2_lim=x2_lim, print_switch=FALSE)
    em_out_BwR <- single_boundary_with_runs_BL_emulator(xP=xP, xD=xD, D=f(xD), 
                                                        theta=theta, sigma=sigma, E_f=E_f)
    
    VarD_fx_mat <- matrix(em_out_BwR[,"VarD_f.x."], nrow=length(x_grid), ncol=length(x_grid)) 
    BwR_variance_storage[i] <- av_var_func(VarD_fx_mat)
  }
  
  av_var <- mean(BwR_variance_storage)
  sd_var <- sd(BwR_variance_storage)
  Q1_var <- quantile(BwR_variance_storage, 0.25)
  Median_var <- median(BwR_variance_storage)
  Q3_var <- quantile(BwR_variance_storage, 0.75)
  BwR_av_var <- list(   "n"=n,
                        "av_var"=av_var,
                        "sd_var"=sd_var,
                        "CI_lower_3sd_var"=max(0, av_var - 3*sd_var),
                        "CI_upper_3sd_var"=av_var + 3*sd_var,
                        "CI_lower_1sd_var"=max(0, av_var - sd_var),
                        "CI_upper_1sd_var"=av_var + sd_var,
                        "Q1"=Q1_var[[1]],
                        "Median"=Median_var,
                        "Q3"=Q3_var[[1]])
  
  if(is.null(n_lower)){
    n_seq <- seq(n,floor(0.5*n),by=-1)
  }else{
    n_seq <- seq(n, n_lower, by=-1)
  }
  
  ############# BwD Loop #############
  BwDR_variance_storage <- matrix(0, nrow=length(n_seq), ncol=number_to_av)
  # rows = sample size n, n-1, n-2, ..., 1
  # cols = iteration i=1, 2, 3, ..., number_to_av
  
  for(j in 1:length(n_seq)){
    for(i in 1:number_to_av){
      xD <- best_lhd(n_seq[j], M=M, x1_lim=x1_lim_2, x2_lim=x2_lim_2, print_switch=FALSE)
      em_out_BwDR <- single_BwD_with_runs_BL_emulator(xP=xP, xD=xD, D=f(xD),
                                                      boundary=boundary,
                                                      theta=theta, sigma=sigma, E_f=E_f, E_df=E_df)
      
      VarD_fx_mat <- matrix(em_out_BwDR[,2], nrow=length(x_grid), ncol=length(x_grid)) 
      BwDR_variance_storage[j,i] <- av_var_func(VarD_fx_mat)
    }
  }
  
  av_var <- rowMeans(BwDR_variance_storage)
  sd_var <- apply(BwDR_variance_storage, 1, sd)
  Q1_var <- apply(BwDR_variance_storage, 1, function(x) quantile(x, 0.25)[[1]])
  Median_var <- apply(BwDR_variance_storage, 1, median)
  Q3_var <- apply(BwDR_variance_storage, 1, function(x) quantile(x, 0.75)[[1]])
  BwDR_av_var <- list(     "n"=n_seq,
                           "av_var"=av_var,
                           "sd_var"=sd_var,
                           "CI_lower_3sd_var"=pmax(0, av_var - 3*sd_var),
                           "CI_upper_3sd_var"=av_var + 3*sd_var,
                           "CI_lower_1sd_var"=pmax(0, av_var - sd_var),
                           "CI_upper_1sd_var"=av_var + sd_var,
                           "Q1"=Q1_var,
                           "Median"=Median_var,
                           "Q3"=Q3_var)
  
  ######### SUMMARY TABLE ##########
  BwR_heading <- paste0("BwR_", n)
  BwDR_headings <- rep(NA, length(n_seq))
  for(i in 1:length(n_seq)){
    BwDR_headings[i] <- paste0("BwDR_", n_seq[i])
  }
  
  CI_3sd_contains_BwR_mean <- (BwR_av_var$av_var >= BwDR_av_var$CI_lower_3sd_var & 
                                 BwR_av_var$av_var <= BwDR_av_var$CI_upper_3sd_var)
  CI_1sd_contains_BwR_mean <- (BwR_av_var$av_var >= BwDR_av_var$CI_lower_1sd_var & 
                                 BwR_av_var$av_var <= BwDR_av_var$CI_upper_1sd_var)
  CI_1sd_overlaps_BwR <- (BwR_av_var$CI_lower_1sd_var >= BwDR_av_var$CI_lower_1sd_var & 
                            BwR_av_var$CI_lower_1sd_var <= BwDR_av_var$CI_upper_1sd_var) |
    (BwR_av_var$CI_upper_1sd_var >= BwDR_av_var$CI_lower_1sd_var & 
       BwR_av_var$CI_upper_1sd_var <= BwDR_av_var$CI_upper_1sd_var) |
    (BwR_av_var$CI_lower_1sd_var <= BwDR_av_var$CI_lower_1sd_var &
       BwR_av_var$CI_upper_1sd_var >= BwDR_av_var$CI_upper_1sd_var)
  IQR_contains_BwR_mean <- (BwR_av_var$av_var >= BwDR_av_var$Q1 & 
                              BwR_av_var$av_var <= BwDR_av_var$Q3)
  IQR_contains_BwR_median <- (BwR_av_var$Median >= BwDR_av_var$Q1 & 
                                BwR_av_var$Median <= BwDR_av_var$Q3)
  IQR_overlaps_BwR <- (BwR_av_var$Q1 >= BwDR_av_var$Q1 & 
                         BwR_av_var$Q1 <= BwDR_av_var$Q3) |
    (BwR_av_var$Q3 >= BwDR_av_var$Q1 & 
       BwR_av_var$Q3 <= BwDR_av_var$Q3) |
    (BwR_av_var$Q1 <= BwDR_av_var$Q1 &
       BwR_av_var$Q3 >= BwDR_av_var$Q3)
  
  summary_table <- data.frame(num_of_runs=c(BwR_av_var$n, BwDR_av_var$n), 
                              mean_average_var=c(BwR_av_var$av_var, BwDR_av_var$av_var),
                              sd_average_var=c(BwR_av_var$sd_var, BwDR_av_var$sd_var),
                              CI_3sd_lower=c(BwR_av_var$CI_lower_3sd_var, BwDR_av_var$CI_lower_3sd_var),
                              CI_3sd_upper=c(BwR_av_var$CI_upper_3sd_var, BwDR_av_var$CI_upper_3sd_var),
                              CI_3sd_contains_BwR=c(NA, CI_3sd_contains_BwR_mean),
                              CI_1sd_lower=c(BwR_av_var$CI_lower_1sd_var, BwDR_av_var$CI_lower_1sd_var),
                              CI_1sd_upper=c(BwR_av_var$CI_upper_1sd_var, BwDR_av_var$CI_upper_1sd_var),
                              CI_1sd_contains_BwR=c(NA, CI_1sd_contains_BwR_mean),
                              CI_1sd_overlaps_BwR=c(NA, CI_1sd_overlaps_BwR),
                              Median=c(BwR_av_var$Median, BwDR_av_var$Median),
                              Q1=c(BwR_av_var$Q1, BwDR_av_var$Q1),
                              Q3=c(BwR_av_var$Q3, BwDR_av_var$Q3),
                              IQR_contains_BwR_mean=c(NA, IQR_contains_BwR_mean),
                              IQR_contains_BwR_median=c(NA, IQR_contains_BwR_median),
                              IQR_overlaps_BwR=c(NA, IQR_overlaps_BwR),
                              row.names = c(BwR_heading, BwDR_headings))
  
  # colnames(summary_table) <- c("Number of Runs",
  #                   "Mean of Average Emulator Adjusted Variance",
  #                   "SD of Average Emulator Adjusted Variance",
  #                   "CI: lower bound",
  #                   "CI: upper bound",
  #                   "CI contains BwR Mean Average Adj Var: Yes/No",
  #                   "Median Average Emulator Adjusted Variance",
  #                   "Q1",
  #                   "Q3",
  #                   "IQR contains BwR Mean Average Adj Var: Yes/No",
  #                   "IQR contains BwR Median Average Adj Var: Yes/No",
  #                   "IQR overlaps with BwR Average Adj Var IQR: Yes/No")
  
  ######### RESULTS matrix ##########
  
  results_mat <- cbind(BwR_variance_storage, t(BwDR_variance_storage))
  results <- as.data.frame(results_mat)
  colnames(results) <- c(BwR_heading, BwDR_headings)
  
  ############# OUTPUT #############
  
  if(print_headline==TRUE){
    
    index <- which.min((BwR_av_var$av_var - BwDR_av_var$av_var)^2)
    cat("BwR Emulator Average Variance for LHD size", BwR_av_var$n, ":", BwR_av_var$av_var,
        "\nBoundary with Derivative and Runs Emulator Average Variance for sizes: \n") 
    print(t(as.matrix(data.frame(av_var=BwDR_av_var$av_var, row.names=as.character(BwDR_av_var$n)))))
    cat("n which gives the closest Boundary with Derivative and Runs Emulator average variance:", BwDR_av_var$n[index],
        ", \nDifference between BwR Emulator Variance with n =", BwR_av_var$n,
        " and Boundary with Derivative and Runs Simple Emulator Variance with n =", BwDR_av_var$n[index], 
        ":", abs(BwR_av_var$av_var - BwDR_av_var$av_var[index]))
    
  }
  if(print_summary_table==TRUE){
    print(summary_table)
    View(summary_table)
  }
  if(return_summary_table & return_results){
    return(list("summary_table"=summary_table,
                "results"=results))
  }else if(return_summary_table){
    return(summary_table)
  }else if(return_results){
    return(results)
  }
}



