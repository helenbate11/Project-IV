
# The function to save average variance comparison plots into folder 
# with appropriate filenames
# Note: only valid for 2d input spaces and (if boundary known: 
# the single boundary case only)

# Saves with following file name:
# paste0(function_name, "_", control_type, "_", type, "_VC_mean_1sd_plot_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png")

source("./Variance Comparators.R")
source("./Variance Comparison Plot Generators.R")

################################################################
# Variance Comparison Plots Saving function
################################################################

saving_var_comp_plots <- function(number_to_av,                        # the number of LHDs for which to store the average variance of each emulator
                                  M,                                   # the number of LHDs to check before choosing one to simulate emulator for (see best_lhd() from "LHD Generator Functions.R")
                                  n,                                   # the number of runs the control emulator uses (we place using LHD)
                                  n_lower,                             # the minimum number of runs for which summarise the average variance for the alternative emulator
                                  x_grid,                              # the 1d grid of input points to expand for set of points for which we evaluate each emulator (expected square input space with equal scales)
                                  theta,                               # the correlation length vector
                                  sigma,                               # the prior standard deviation
                                  E_f=0,                               # the prior expectation for f(x)
                                  E_dfdx=c(0,0),                       # the prior expectation for df/dx(x)
                                  boundary=list("x1"=NULL, "x2"=NULL), # boundary definition: x1 = boundary$x1 OR x2 = boundary$x2 (other coordinate=NULL)
                                  x1_lim=c(0,1), x2_lim=c(0,1),        # the x1- and x2-input space limits for LHD generation of the control emulator (see best_lhd() from "LHD Generator Functions.R")
                                  x1_lim_2=NULL, x2_lim_2=NULL,        # (optional) the x1- and x2-input space limits for LHD generation of the alternative emulator (if different from that of the control)
                                  seed=NULL,                           # set the seed of random LHD generation
                                  wd=690, ht=620,                      # file width and height
                                  filepath="",                         # filepath for folder
                                  function_name="F",                   # function name (in file name)
                                  type=NULL                            # Comparison type: "Deriv", "BwR", "Simple-BwDR" or "BwR-BwDR" 
){
  if(type=="Deriv"){
    equilivator <- deriv_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                            x_grid=x_grid, theta=theta, sigma=sigma, 
                                            x1_lim=x1_lim, x2_lim=x2_lim, 
                                            x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed, 
                                            E_f=E_f, E_dfdx=E_dfdx,  
                                            print_headline=FALSE, print_summary_table=FALSE,
                                            return_summary_table=TRUE, return_results=TRUE)
    control_type <- "Simple"
  }else if(type=="BwR"){
    equilivator <- single_boundary_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                                      x_grid=x_grid, theta=theta, sigma=sigma, 
                                                      x1_lim=x1_lim, x2_lim=x2_lim, 
                                                      x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed, 
                                                      E_f=0, boundary=boundary, 
                                                      print_headline=FALSE, print_summary_table=FALSE,
                                                      return_summary_table=TRUE, return_results=TRUE)
    control_type <- "Simple"
  }else if(type=="Simple-BwDR"){
    equilivator <- simple_BwDR_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                                  x_grid=x_grid, theta=theta, sigma=sigma, 
                                                  x1_lim=x1_lim, x2_lim=x2_lim, 
                                                  x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed, 
                                                  E_f=0, boundary=boundary, 
                                                  print_headline=FALSE, print_summary_table=FALSE,
                                                  return_summary_table=TRUE, return_results=TRUE)
    control_type <- "Simple"
  }else if(type=="BwR-BwDR"){
    equilivator <- BwR_BwDR_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                               x_grid=x_grid, theta=theta, sigma=sigma, 
                                               x1_lim=x1_lim, x2_lim=x2_lim, 
                                               x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed, 
                                               E_f=0, boundary=boundary, 
                                               print_headline=FALSE, print_summary_table=FALSE,
                                               return_summary_table=TRUE, return_results=TRUE)
    control_type <- "BwR"
  }
  
  boxplot <- plotting_var_comp_boxplot(equilivator, type=type)
  mean_error_plot <- plotting_var_comp_mean_error_plot(equilivator, type=type)
  
  # boxplot
  png(filename = paste0(filepath, function_name, "_", control_type, "_", type, "_VC_boxplot_theta", theta[1], "_", theta[2], "_sigma", sigma, 
                        ".png"),
      width = wd, height = ht)
  print(boxplot)
  dev.off()
  
  # mean ± 1 sd error bar plot
  png(filename = paste0(filepath, function_name, "_", control_type, "_", type, "_VC_mean_1sd_plot_theta", theta[1], "_", theta[2], "_sigma", sigma, 
                        ".png"),
      width = wd, height = ht)
  print(mean_error_plot)
  dev.off()
}







