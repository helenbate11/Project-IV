# Code to visualise the number of runs required to achieve a similar level of 
# average emulator uncertainty between different emulators in
# quantile boxplot OR
# mean and standard deviation form

library(ggplot2)
library(tidyr)
library(tidyverse)

################################################################
# Quantile Boxplot of average emulator variance
################################################################

plotting_var_comp_boxplot <- function(equilivator,  # the variance equilivator output (from functions in: "Variance Comparators.R")
                                      type="Deriv", # the type of emulators being compared: "Deriv", "BwR", "Simple-BwDR" or "BwR-BwDR"
                                      
                                      # other inputs only if equilivator object not supplied
                                      number_to_av=NULL,   # see equilivator input variable definitions in "Variance Comparators.R"
                                      M=NULL, 
                                      n=NULL, 
                                      n_lower=NULL, 
                                      x_grid=NULL, 
                                      theta=NULL,
                                      sigma=NULL, 
                                      E_f=0, 
                                      E_dfdx=c(0,0), 
                                      boundary=NULL, 
                                      x1_lim=c(0,1), x2_lim=c(0,1),
                                      x1_lim_2=NULL, x2_lim_2=NULL,
                                      seed=NULL
                                      ){
  
  # generate equilivator object if not supplied
  if(is.null(equilivator) & type=="Deriv"){
    equilivator <- deriv_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                            x_grid=x_grid, theta=theta, sigma=sigma,
                                            x1_lim=x1_lim, x2_lim=x2_lim, 
                                            x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed,
                                            E_f=E_f, E_dfdx=E_dfdx,
                                            print_headline=FALSE, print_summary_table=FALSE,
                                            return_summary_table=TRUE, return_results=TRUE)
  }else if(is.null(equilivator) & type=="BwR"){
    equilivator <- boundary_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                               x_grid=x_grid, theta=theta, sigma=sigma,
                                               x1_lim=x1_lim, x2_lim=x2_lim, 
                                               x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed,
                                               E_f=E_f, boundary=boundary,
                                               print_headline=FALSE, print_summary_table=FALSE,
                                               return_summary_table=TRUE, return_results=TRUE)
  }else if(is.null(equilivator) & type=="Simple-BwDR"){
    equilivator <- simple_BwDR_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                                  x_grid=x_grid, theta=theta, sigma=sigma,
                                                  x1_lim=x1_lim, x2_lim=x2_lim, 
                                                  x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed,
                                                  E_f=E_f, boundary=boundary,
                                                  print_headline=FALSE, print_summary_table=FALSE,
                                                  return_summary_table=TRUE, return_results=TRUE)
  }else if(is.null(equilivator) & type=="BwR-BwDR"){
    equilivator <- BwR_BwDR_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                               x_grid=x_grid, theta=theta, sigma=sigma,
                                               x1_lim=x1_lim, x2_lim=x2_lim, 
                                               x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed,
                                               E_f=E_f, boundary=boundary,
                                               print_headline=FALSE, print_summary_table=FALSE,
                                               return_summary_table=TRUE, return_results=TRUE)
  }
  
  summary_table <- equilivator$summary_table
  results_by_n <- equilivator$results
  
  # boxplot labels
  set_levels <- colnames(results_by_n)
  if(type=="Deriv"){
    type_vec <- c("Simple", rep("Derivative", length(set_levels)-1))
    title <- "Comparison of the Average Variance for Simple and Derivative Emulators"
  }else if(type=="BwR"){
    type_vec <- c("Simple", rep("Boundary", length(set_levels)-1))
    title <- "Comparison of the Average Variance for Simple and Boundary with Runs Emulators"
  }else if(type=="Simple-BwDR"){
    type_vec <- c("Simple", rep("Boundary\nwith\nDerivative", length(set_levels)-1))
    title <- "Comparison of the Average Variance for\nSimple and Boundary with Derivative and Runs Emulators"
  }else if(type=="BwR-BwDR"){
    type_vec <- c("Boundary", rep("Boundary\nwith\nDerivative", length(set_levels)-1))
    title <- "Comparison of the Average Variance for\nBoundary with Runs and Boundary with Derivative and Runs Emulators"
  }
  factor_labels <- paste0(type_vec, "\n(n=", sub(".*_", "", set_levels), ")")
  
  # reshaping data
  results_long <- pivot_longer(results_by_n, cols = everything(),
                               names_to = "Set",
                               values_to = "Average_Variance") |>
    mutate(Set = factor(Set,
                        levels=set_levels,
                        labels=factor_labels))
  
  my_purple <- rgb(0.52, 0.15, 0.5)
  my_lightpurple <- rgb(0.82, 0.4, 0.8)
  light_purples <- colorRampPalette(c("white", my_lightpurple))(3)
  
  # generating plot
  boxplot <- ggplot(data=results_long, aes(x = Set, y = Average_Variance)) +
    geom_boxplot(fill="orange", outlier.size=1) +
    # light_purples[2]
    theme_bw() +
    geom_hline(yintercept=c(summary_table[1,"Q1"], summary_table[1,"Q3"]),
               colour="red", linetype="dashed", linewidth=0.7) +
    # my_purple
    geom_vline(xintercept=1.5, colour="darkgrey", linewidth=0.4) +
    labs(x="Emulator Type (n = number of runs)",
         y="Average Variance",
         title=title) +
    theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
          axis.title.x = element_text(colour="black", size = 14, margin=margin(t=10)),
          axis.title.y = element_text(colour="black", size = 14, margin=margin(r=5)),
          axis.text.x = element_text(colour="black", size=11),
          axis.text.y = element_text(colour="black", size=11)
    )
  return(boxplot)
}



################################################################
# Mean with  ± 1 standard deviation error bar plot for average emulator variance
################################################################

plotting_var_comp_mean_error_plot <- function(equilivator,  # the variance equilivator output (from functions in: "Variance Comparators.R")
                                              type="Deriv", # the type of emulators being compared: "Deriv", "BwR", "Simple-BwDR" or "BwR-BwDR"
                                              
                                              # other inputs only if equilivator object not supplied
                                              number_to_av=NULL,   # see equilivator input variable definitions in "Variance Comparators.R"
                                              M=NULL, 
                                              n=NULL, 
                                              n_lower=NULL, 
                                              x_grid=NULL, 
                                              theta=NULL,
                                              sigma=NULL, 
                                              E_f=0, 
                                              E_dfdx=c(0,0), 
                                              boundary=NULL, 
                                              x1_lim=c(0,1), x2_lim=c(0,1),
                                              x1_lim_2=NULL, x2_lim_2=NULL,
                                              seed=NULL
                                              ){
  
  if(is.null(equilivator) & type=="Deriv"){
    equilivator <- deriv_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                            x_grid=x_grid, theta=theta, sigma=sigma, 
                                            x1_lim=x1_lim, x2_lim=x2_lim, 
                                            x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed, 
                                            E_f=E_f, E_dfdx=E_dfdx,  
                                            print_headline=FALSE, print_summary_table=FALSE,
                                            return_summary_table=TRUE, return_results=TRUE)
  }else if(is.null(equilivator) & type=="BwR"){
    equilivator <- single_boundary_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                                      x_grid=x_grid, theta=theta, sigma=sigma, 
                                                      x1_lim=x1_lim, x2_lim=x2_lim, 
                                                      x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed, 
                                                      E_f=E_f, boundary=boundary, 
                                                      print_headline=FALSE, print_summary_table=FALSE,
                                                      return_summary_table=TRUE, return_results=TRUE)
  }else if(is.null(equilivator) & type=="Simple-BwDR"){
    equilivator <- simple_BwDR_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                                  x_grid=x_grid, theta=theta, sigma=sigma, 
                                                  x1_lim=x1_lim, x2_lim=x2_lim, 
                                                  x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed, 
                                                  E_f=E_f, boundary=boundary, 
                                                  print_headline=FALSE, print_summary_table=FALSE,
                                                  return_summary_table=TRUE, return_results=TRUE)
  }else if(is.null(equilivator) & type=="BwR-BwDR"){
    equilivator <- BwR_BwDR_av_var_equilivator(number_to_av=number_to_av, M=M, n=n, n_lower=n_lower,
                                               x_grid=x_grid, theta=theta, sigma=sigma, 
                                               x1_lim=x1_lim, x2_lim=x2_lim, 
                                               x1_lim_2=x1_lim_2, x2_lim_2=x2_lim_2, seed=seed, 
                                               E_f=E_f, boundary=boundary, 
                                               print_headline=FALSE, print_summary_table=FALSE,
                                               return_summary_table=TRUE, return_results=TRUE)
  }
  
  summary_table <- equilivator$summary_table
  # x-axis labels
  set_levels <- rownames(summary_table)
  if(type=="Deriv"){
    type_vec <- c("Simple", rep("Derivative", length(set_levels)-1))
    title <- "Comparison of the Average Variance for Simple and Derivative Emulators"
  }else if(type=="BwR"){
    type_vec <- c("Simple", rep("Boundary", length(set_levels)-1))
    title <- "Comparison of the Average Variance for Simple and Boundary with Runs Emulators"
  }else if(type=="Simple-BwDR"){
    type_vec <- c("Simple", rep("Boundary\nwith\nDerivative", length(set_levels)-1))
    title <- "Comparison of the Average Variance for\nSimple and Boundary with Derivative and Runs Emulators"
  }else if(type=="BwR-BwDR"){
    type_vec <- c("Boundary", rep("Boundary\nwith\nDerivative", length(set_levels)-1))
    title <- "Comparison of the Average Variance for\nBoundary with Runs and Boundary with Derivative and Runs Emulators"
  }
  factor_labels <- paste0(type_vec, "\n(n=", sub(".*_", "", set_levels), ")")
  
  my_purple <- rgb(0.52, 0.15, 0.5)
  my_lightpurple <- rgb(0.82, 0.4, 0.8)
  light_purples <- colorRampPalette(c("white", my_lightpurple))(3)
  
  # generating plot
  mean_error_plot <- ggplot(summary_table, aes(x = factor(1:nrow(summary_table), labels=factor_labels),
                                               y = mean_average_var)) +
    geom_errorbar(aes(ymin = mean_average_var - sd_average_var, 
                      ymax = mean_average_var + sd_average_var),
                  width = 0.2) +
    geom_point(shape=21, colour="black", fill="orange", size=3) +
    theme_bw() +
    geom_hline(yintercept=c(summary_table[1,"CI_1sd_lower"], summary_table[1,"CI_1sd_upper"]),
               colour="red", linetype="dashed", linewidth=0.7) +
    # my_lightpurple
    geom_vline(xintercept=1.5, colour="darkgrey", linewidth=0.4) +
    labs(x = "Emulator Type (n = number of runs)", 
         y = "Average Variance", 
         title = title,
         subtitle = "(with ± 1 standard deviation uncertainty bars)") +
    theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
          plot.subtitle = element_text(hjust = 0.5, face="bold", size=12),
          axis.title.x = element_text(colour="black", size = 14, margin=margin(t=10)),
          axis.title.y = element_text(colour="black", size = 14, margin=margin(t=5)),
          axis.text.x = element_text(colour="black", size=11),
          axis.text.y = element_text(colour="black", size=11)
    )
  
  return(mean_error_plot)
}

