# The function to save SIR emulator output plots into a folder 
# with appropriate filename, plot titles and colour scheme

# Note: file names often depend on global parameters of theta, sigma, E_f and E_df
# so save the relevant plots immediately after running that emulator

# Note: calls on global function of interest f()

################################################################
# 2d Emulator Output Plots Saving Function - SIR
# (for the 2d Simple Emulator)
################################################################

# Saves with following file names:
# paste0("SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, ".png")
# etc.

saving_simple_em_plots <- function(em_out,              # emulator output matrix
                                   xD,                  # the run input locations
                                   beta_grid,           # grid edge locations defining xP[,1], the beta input space
                                   gamma_grid,          # grid edge locations defining xP[,2], the gamma input space
                                   timepoint,           # the time at which SIR model evaulated
                                   variable="default",  # the model compartment being emulated: "S", "I" or "R"
                                   cont_levs_exp=NULL,  # the contour plot levels for the expectation E_D[f(x)] plot
                                   cont_levs_var=NULL,  # the contour plot levels for the variance Var_D[f(x)] plot
                                   EmulatorType="Simple",     # the type of emulator being plotted: "Simple"
                                   wd=690, ht=620,      # the height and width of the plots
                                   filepath="",         # the filepath to the folder in which to save plots (MUST SET)
                                   function_name="F",   # the model function name
                                   diag_true_func_plot=TRUE   # TRUE or FALSE: plot the diagnositcs and true function plots (saves time if FALSE)
){
  
  ED_fx_mat <- matrix(em_out[,"ExpD_f.x."], nrow=length(beta_grid), ncol=length(beta_grid))
  VarD_fx_mat <- matrix(em_out[,"VarD_f.x."], nrow=length(beta_grid), ncol=length(beta_grid))
  
  ED_main1 <- bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Expectation" ~ E[D] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
  ED_main2 <- bquote(bold("for SIR output" ~ .(variable) ~ "at timepoint" ~  t==.(as.character(timepoint))))
  VarD_main1 <- bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Variance" ~ Var[D] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
  VarD_main2 <- bquote(bold("for SIR output" ~ .(variable) ~ "at timepoint" ~  t==.(as.character(timepoint))))
  SD_main1 <- bquote(bold(.(EmulatorType) ~ "Emulator Diagnostics" ~ S[D] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
  SD_main2 <- bquote(bold("for SIR output" ~ .(variable) ~ "at timepoint" ~  t==.(as.character(timepoint))))
  true_main1 <- bquote(bold("True SIR Output" ~ .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")"))
  # true_main2 <- bquote(bold("for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))))

  # Emulator plots
  png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, ".png"),
      width = wd, height = ht)
  ED_fx_plot <- SIR_emul_fill_cont(cont_mat=ED_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, type=variable, xD=xD,
                                   cont_levs=cont_levs_exp,
                                   cex.main=2, cex.lab=1.7, cex.axis=1.4,
                                   key.axes = axis(4, cex.axis = 1.3))
  # title(xlab = bquote("Beta Input" ~ beta), line = 3, cex.lab=1.7)  
  title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
  title(main = ED_main1, cex.main = 2, outer=TRUE, line=-2.2)
  # title(main = ED_main1, line = -1.3, cex.main = 2, outer=TRUE)
  # title(main = ED_main2, line = -3.1, cex.main = 2, outer=TRUE)
  dev.off()
  
  png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType, "_VarDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, ".png"),
      width = wd, height = ht)
  VarD_fx_plot <- SIR_emul_fill_cont(cont_mat=VarD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, type="var", xD=xD,
                                     cont_levs=cont_levs_var,
                                     cex.main=2, cex.lab=1.7, cex.axis=1.4,
                                     key.axes = axis(4, cex.axis = 1.3))
  # title(xlab = bquote("Beta Input" ~ beta), line = 3, cex.lab=1.7)  
  title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
  title(main = VarD_main1, cex.main = 2, outer=TRUE, line=-2.2)
  # title(main = VarD_main1, line = -1.3, cex.main = 2, outer=TRUE)
  # title(main = VarD_main2, line = -3.1, cex.main = 2, outer=TRUE)
  dev.off()
  
  if(diag_true_func_plot){
    fxP <- SIR_f(xP=xP, timepoint=timepoint, variable=variable, f_0=f_start)
    fxP_mat <- matrix(fxP, nrow=length(beta_grid), ncol=length(beta_grid))
    
    SD_fx_mat <- (ED_fx_mat - fxP_mat) / sqrt(VarD_fx_mat)
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType, "_SDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, ".png"),
        width = wd, height = ht)
    SD_fx_plot <- SIR_emul_fill_cont(cont_mat=SD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid,
                                     type="diag", xD=xD, cont_levs=seq(-3,3,0.25), 
                                     cex.main=2, cex.lab=1.7, cex.axis=1.4,
                                     key.axes = axis(4, cex.axis = 1.3))
    # title(xlab = bquote("Beta Input" ~ beta), line = 3, cex.lab=1.7)  
    title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
    title(main = SD_main1, cex.main = 2, outer=TRUE, line=-2.2)
    # title(main = SD_main1, line = -1.3, cex.main = 2, outer=TRUE)
    # title(main = SD_main2, line = -3.1, cex.main = 2, outer=TRUE)
    dev.off()
    
    # True function plots
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_TrueFunction_fx.png"),
        width = wd, height = ht)
    
    SIR_emul_fill_cont(cont_mat=fxP_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, type=variable, xD=xD,
                       cont_levs=cont_levs_exp,
                       cex.main=2, cex.lab=1.7, cex.axis=1.4,
                       key.axes = axis(4, cex.axis = 1.3))
    # title(xlab = bquote("Beta Input" ~ beta), line = 3, cex.lab=1.7)  
    title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
    title(main = true_main1, cex.main = 2, outer=TRUE, line=-2.2)
    # title(main = true_main1, line = -1.3, cex.main = 2, outer=TRUE)
    # title(main = true_main2, line = -3.1, cex.main = 2, outer=TRUE)
    dev.off()
  }
}


################################################################
# 2d Emulator Output Plots Saving Function - SIR
# (for the 2d Boundary, Boundary with Runs, Boundary with Derivative
# and Boundary with Derivative and Runs Emulators)
################################################################

# Saves with following file names:
# paste0("SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, ".png")
# etc.

saving_boundary_with_runs_plots <- function(em_out,                    # the emulator output matrix
                                            boundary,                  # boundary definition: list(x1 = boundary$x1, x2 = boundary$x2)
                                            xD=NULL,                   # the design run positions
                                            beta_grid,                 # grid edge locations defining xP[,1], the beta input space
                                            gamma_grid,                # grid edge locations defining xP[,2], the gamma input space
                                            timepoint,                 # the time at which SIR model evaulated
                                            variable="default",        # the model compartment being emulated: "S", "I" or "R"
                                            cont_levs_exp=NULL,        # the contour plot levels for the expectation E_D[f(x)] plot
                                            cont_levs_var=NULL,        # the contour plot levels for the variance Var_D[f(x)] plot
                                            wd=800, ht=620,            # the width and height of the saved plots
                                            filepath="",               # the filepath to the folder in which to save plots (MUST SET)
                                            function_name="F",         # the model function name
                                            EmulatorType="Boundary",   # the emulator type: "Boundary", "Boundary with Runs", "Boundary with Derivative on the Boundary", "Boundary with Derivative on the Boundary and Runs")
                                            diag_true_func_plot=TRUE   # TRUE or FALSE: plot the diagnositcs and true function plots (saves time if FALSE)
){
  ED_fx_vec <- em_out[,1]
  VarD_fx_vec <- em_out[,2]
  
  neg_vars <- which(VarD_fx_vec<0)
  large_neg_var <- FALSE
  if( length(neg_vars)>=1 ){
    for(i in neg_vars){
      if( abs(VarD_fx_vec[i]) < 1e-06 ){
        VarD_fx_vec[i] <- 0
      }else{
        large_neg_var <- TRUE
      }
    }
    if(large_neg_var){
      print("Significantly large negative variance - needs investigating")
    }
  }
  
  ED_fx_mat <- matrix(ED_fx_vec, nrow=length(beta_grid), ncol=length(beta_grid)) 
  VarD_fx_mat <- matrix(VarD_fx_vec, nrow=length(beta_grid), ncol=length(beta_grid))
  
  EmulatorType_short <- ifelse(EmulatorType=="Boundary", "Boundary", 
                               ifelse(EmulatorType=="Boundary with Runs", "Boundary with Runs", 
                                      ifelse(EmulatorType=="Boundary with Derivative on the Boundary", "BwD",
                                             "BwD with Runs")))
  
  
  if(EmulatorType =="Boundary"){
    ED_main <- bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Expectation " ~ E[K] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
    VarD_main <- bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Variance " ~ Var[K] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
    SD_main <- bquote(bold(.(EmulatorType) ~ "Emulator Diagnostics " ~ S[K] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
  }else if(EmulatorType == "Boundary with Runs"){
    ED_main1 <- bquote(bold(.(EmulatorType) ~ "Emulator")) 
    ED_main2 <- bquote(bold("Adjusted Expectation " ~ E[D] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
    VarD_main1 <- bquote(bold(.(EmulatorType) ~ "Emulator"))
    VarD_main2 <- bquote(bold("Adjusted Variance " ~ Var[D] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
    SD_main1 <- bquote(bold(.(EmulatorType) ~ "Emulator"))
    SD_main2 <- bquote(bold("Diagnostics " ~ S[D] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
  }else if(EmulatorType == "Boundary with Derivative on the Boundary"){
    ED_main1 <- bquote(bold(.(EmulatorType) ~ "Emulator")) 
    ED_main2 <- bquote(bold("Adjusted Expectation " ~ E[K] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
    VarD_main1 <- bquote(bold(.(EmulatorType) ~ "Emulator"))
    VarD_main2 <- bquote(bold("Adjusted Variance " ~ Var[K] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
    SD_main1 <- bquote(bold(.(EmulatorType) ~ "Emulator"))
    SD_main2 <- bquote(bold("Diagnostics " ~ S[K] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
  }else if(EmulatorType == "Boundary with Derivative on the Boundary and Runs"){
    ED_main1 <- bquote(bold(.(EmulatorType))) 
    ED_main2 <- bquote(bold("Emulator Adjusted Expectation " ~ E[D] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
    VarD_main1 <- bquote(bold(.(EmulatorType)))
    VarD_main2 <- bquote(bold("Emulator Adjusted Variance " ~ Var[D] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
    SD_main1 <- bquote(bold(.(EmulatorType)))
    SD_main2 <- bquote(bold("Emulator Diagnostics " ~ S[D] * "[" * .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")]"))
  }
  
  true_main <- bquote(bold("True SIR Output" ~ .(variable) * "(" * beta * "," ~ gamma ~ "|" ~ t==.(as.character(timepoint)) *")"))

  if(EmulatorType=="Boundary"){
    # Emulator plots for E_D[f(x)]
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType_short, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, "_E_df", E_df[1], "_", E_df[2], "_", E_df[3],
                          "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
        width = wd, height = ht)
    ED_fx_plot <- SIR_emul_fill_cont(cont_mat=ED_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid,
                                     type=variable, xD=xD, cont_levs=cont_levs_exp,
                                     plot_boundary=TRUE, plot_boundary_line=boundary, 
                                     cex.main=2, cex.lab=1.7, cex.axis=1.4,
                                     key.axes = axis(4, cex.axis = 1.3))
    title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
    title(main = ED_main, cex.main = 2, outer=TRUE, line=-2.2)
    dev.off()
    
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType_short, "_VarDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, "_E_df", E_df[1], "_", E_df[2], "_", E_df[3],
                          "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
        width = wd, height = ht)
    VarD_fx_plot <- SIR_emul_fill_cont(cont_mat=VarD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, 
                                       type="var", xD=xD, cont_levs=cont_levs_var,
                                       plot_boundary=TRUE, plot_boundary_line=boundary, 
                                       cex.main=2, cex.lab=1.7, cex.axis=1.4,
                                       key.axes = axis(4, cex.axis = 1.3))
    title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
    title(main = VarD_main, cex.main = 2, outer=TRUE, line=-2.2)
    dev.off()
    
    if(diag_true_func_plot){
      fxP <- SIR_f(xP=xP, timepoint=timepoint, variable=variable, f_0=f_start)
      fxP_mat <- matrix(fxP, nrow=length(beta_grid), ncol=length(beta_grid))
      
      SD_fx_vec <- (ED_fx_vec - fxP) / sqrt(VarD_fx_vec)
      # Fixing diagnostic values on boundary (perfect knowledge => Var_D[f(x)]=0)
      beta_0_points <- which(xP[,1]==0)
      SD_fx_vec[beta_0_points] <- 0
      
      SD_fx_mat <- matrix(SD_fx_vec, nrow=length(beta_grid), ncol=length(beta_grid))
      
      png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType_short, "_SDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, "_E_df", E_df[1], "_", E_df[2], "_", E_df[3],
                            "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
          width = wd, height = ht)
      SD_fx_plot <- SIR_emul_fill_cont(cont_mat=SD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid,
                                       type="diag", xD=xD, cont_levs=seq(-3,3,0.25), 
                                       plot_boundary=TRUE, plot_boundary_line=boundary, 
                                       cex.main=2, cex.lab=1.7, cex.axis=1.4,
                                       key.axes = axis(4, cex.axis = 1.3))
      title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
      title(main = SD_main, cex.main = 2, outer=TRUE, line=-2.2)
      dev.off()
      
      # True function plots
      png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_TrueFunction_fx.png"),
          width = wd, height = ht)
      
      SIR_emul_fill_cont(cont_mat=fxP_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, 
                         type=variable, xD=xD,cont_levs=cont_levs_exp,
                         plot_boundary=TRUE, plot_boundary_line=boundary,
                         cex.main=2, cex.lab=1.7, cex.axis=1.4,
                         key.axes = axis(4, cex.axis = 1.3))
      title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
      title(main = true_main, cex.main = 2, outer=TRUE, line=-2.2)
      dev.off()
    }
  }else{
    # Emulator plots for E_D[f(x)]
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType_short, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, "_E_df", E_df[1], "_", E_df[2], "_", E_df[3],
                          "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
        width = wd, height = ht)
    ED_fx_plot <- SIR_emul_fill_cont(cont_mat=ED_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid,
                                     type=variable, xD=xD, cont_levs=cont_levs_exp,
                                     plot_boundary=TRUE, plot_boundary_line=boundary, 
                                     cex.main=2, cex.lab=1.7, cex.axis=1.4,
                                     key.axes = axis(4, cex.axis = 1.3))
    title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
    title(main = ED_main1, line = -1.3, cex.main = 2, outer=TRUE)
    title(main = ED_main2, line = -3.1, cex.main = 2, outer=TRUE)
    dev.off()
    
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType_short, "_VarDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, "_E_df", E_df[1], "_", E_df[2], "_", E_df[3],
                          "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
        width = wd, height = ht)
    VarD_fx_plot <- SIR_emul_fill_cont(cont_mat=VarD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, 
                                       type="var", xD=xD, cont_levs=cont_levs_var,
                                       plot_boundary=TRUE, plot_boundary_line=boundary, 
                                       cex.main=2, cex.lab=1.7, cex.axis=1.4,
                                       key.axes = axis(4, cex.axis = 1.3))
    title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
    title(main = VarD_main1, line = -1.3, cex.main = 2, outer=TRUE)
    title(main = VarD_main2, line = -3.1, cex.main = 2, outer=TRUE)
    dev.off()
    
    
    if(diag_true_func_plot){
        fxP <- SIR_f(xP=xP, timepoint=timepoint, variable=variable, f_0=f_start)
        fxP_mat <- matrix(fxP, nrow=length(beta_grid), ncol=length(beta_grid))
        
        SD_fx_vec <- (ED_fx_vec - fxP) / sqrt(VarD_fx_vec)
        # Fixing diagnostic values on boundary (perfect knowledge => Var_D[f(x)]=0)
        beta_0_points <- which(xP[,1]==0)
        SD_fx_vec[beta_0_points] <- 0
        
        SD_fx_mat <- matrix(SD_fx_vec, nrow=length(beta_grid), ncol=length(beta_grid))
        
        png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType_short, "_SDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, "_E_df", E_df[1], "_", E_df[2], "_", E_df[3],
                              "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
            width = wd, height = ht)
        SD_fx_plot <- SIR_emul_fill_cont(cont_mat=SD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid,
                                         type="diag", xD=xD, cont_levs=seq(-3,3,0.25), 
                                         plot_boundary=TRUE, plot_boundary_line=boundary, 
                                         cex.main=2, cex.lab=1.7, cex.axis=1.4,
                                         key.axes = axis(4, cex.axis = 1.3))
        title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
        title(main = SD_main1, line = -1.3, cex.main = 2, outer=TRUE)
        title(main = SD_main2, line = -3.1, cex.main = 2, outer=TRUE)
        dev.off()
        
        # True function plots
        png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_TrueFunction_fx.png"),
            width = wd, height = ht)
        
        SIR_emul_fill_cont(cont_mat=fxP_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, 
                           type=variable, xD=xD,cont_levs=cont_levs_exp,
                           plot_boundary=TRUE, plot_boundary_line=boundary,
                           cex.main=2, cex.lab=1.7, cex.axis=1.4,
                           key.axes = axis(4, cex.axis = 1.3))
        title(ylab = bquote("Gamma Input" ~ gamma), line = 2, cex.lab=1.7)
        title(main = true_main, cex.main = 2, outer=TRUE, line=-2.2)
        dev.off()
    }
  }
}

