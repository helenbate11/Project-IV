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
  
  # Emulator plots
  png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, ".png"),
      width = wd, height = ht)
  ED_fx_plot <- SIR_emul_fill_cont(cont_mat=ED_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, type=variable, xD=xD,
                                   cont_levs=cont_levs_exp,
                                   main=bquote(bold(atop(.(EmulatorType) ~ "Emulator Adjusted Expectation" ~ E[D] * "[f(x)]", "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))))),
                                   cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType, "_VarDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, ".png"),
      width = wd, height = ht)
  VarD_fx_plot <- SIR_emul_fill_cont(cont_mat=VarD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, type="var", xD=xD,
                                     cont_levs=cont_levs_var,
                                     main=bquote(bold(atop(.(EmulatorType) ~ "Emulator Adjusted Variance" ~ Var[D] * "[f(x)]", "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))))),
                                     cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  if(diag_true_func_plot){
    fxP <- SIR_f(xP=xP, timepoint=timepoint, variable=variable, f_0=f_start)
    fxP_mat <- matrix(fxP, nrow=length(beta_grid), ncol=length(beta_grid))
    
    SD_fx_mat <- (ED_fx_mat - fxP_mat) / sqrt(VarD_fx_mat)
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType, "_SDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, ".png"),
        width = wd, height = ht)
    SD_fx_plot <- SIR_emul_fill_cont(cont_mat=SD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid,
                                     type="diag", xD=xD, cont_levs=seq(-3,3,0.25), 
                                     main=bquote(bold(atop(.(EmulatorType) ~ "Emulator Diagnostics" ~ S[D] * "[f(x)]", "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))))),
                                     cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
    dev.off()
    
    # True function plots
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_TrueFunction_fx.png"),
        width = wd, height = ht)
    
    true_main <- bquote(bold(atop("True Function f(x)", 
                                  "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint)))))
    
    SIR_emul_fill_cont(cont_mat=fxP_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, type=variable, xD=xD,
                       cont_levs=cont_levs_exp,
                       main=true_main, 
                       cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
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
                                            EmulatorType="Boundary",   # the emulator type: "Boundary", "Boundary with Runs", "Boundary and Derivative on the Boundary", "Boundary and Derivative on the Boundary with Runs")
                                            diag_true_func_plot=TRUE   # TRUE or FALSE: plot the diagnositcs and true function plots (saves time if FALSE)
                                            ){
  ED_fx_vec <- em_out[,1]
  VarD_fx_vec <- em_out[,2]
  
  neg_vars <- which(VarD_fx_vec<0)
  if( length(neg_vars)>=1 ){
    for(i in neg_vars){
      large_neg_var <- FALSE
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
                                      ifelse(EmulatorType=="Boundary and Derivative on the Boundary", "BwD",
                                             "BwD with Runs")))
  
  
  if(EmulatorType =="Boundary"){
    ED_main <- bquote(bold(atop(.(EmulatorType) ~ "Emulator Adjusted Expectation" ~ E[K] * "[f(x)]", 
                                "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint)))))
    VarD_main <- bquote(bold(atop(.(EmulatorType) ~ "Emulator Adjusted Variance" ~ Var[K] * "[f(x)]", 
                                  "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint)))))
    SD_main <- bquote(bold(atop(.(EmulatorType) ~ "Emulator Diagnostics" ~ S[K] * "[f(x)]", 
                                "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint)))))
  }else if(EmulatorType == "Boundary with Runs"){
    ED_main <- bquote(bold(atop(.(EmulatorType) ~ "Emulator Adjusted Expectation" ~ E[D] * "[f(x)]", 
                                "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint)))))
    VarD_main <- bquote(bold(atop(.(EmulatorType) ~ "Emulator Adjusted Variance" ~ Var[D] * "[f(x)]", 
                                  "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint)))))
    SD_main <- bquote(bold(atop(.(EmulatorType) ~ "Emulator Diagnostics" ~ S[D] * "[f(x)]", 
                                "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint)))))
  }else if(EmulatorType == "Boundary and Derivative on the Boundary"){
    ED_main <- bquote(bold(
      atop(
        .(EmulatorType) ~ "Emulator",
        "Adjusted Expectation " ~ E[K] * "[f(x)]" ~ "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))
      )
    ))
    VarD_main <- bquote(bold(
      atop(
        .(EmulatorType) ~ "Emulator",
        "Adjusted Variance " ~ Var[K] * "[f(x)]" ~ "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))
      )
    ))
    SD_main <- bquote(bold(
      atop(
        .(EmulatorType),
        "Diagnostics " ~ S[K] * "[f(x)]" ~ "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))
      )
    ))
  }else if(EmulatorType == "Boundary and Derivative on the Boundary with Runs"){
    ED_main <- bquote(bold(
      atop(
        .(EmulatorType) ~ "Emulator",
        "Adjusted Expectation " ~ E[D] * "[f(x)]" ~ "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))
      )
    ))
    VarD_main <- bquote(bold(
      atop(
        .(EmulatorType) ~ "Emulator",
        "Adjusted Variance " ~ Var[D] * "[f(x)]" ~ "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))
      )
    ))
    SD_main <- bquote(bold(
      atop(
        .(EmulatorType),
        "Diagnostics " ~ S[D] * "[f(x)]" ~ "for SIR output" ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint))
      )
    ))
  }
  
  # Emulator plots for E_D[f(x)]
  png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType_short, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, "_E_df", E_df,
                        "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
      width = wd, height = ht)
  ED_fx_plot <- SIR_emul_fill_cont(cont_mat=ED_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid,
                                   type=variable, xD=xD, cont_levs=cont_levs_exp,
                                   plot_boundary=TRUE, plot_boundary_line=boundary, 
                                   main=ED_main,
                                   cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType_short, "_VarDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, "_E_df", E_df,
                        "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
      width = wd, height = ht)
  VarD_fx_plot <- SIR_emul_fill_cont(cont_mat=VarD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, 
                                     type="var", xD=xD, cont_levs=cont_levs_var,
                                     plot_boundary=TRUE, plot_boundary_line=boundary, 
                                     main=VarD_main,
                                     cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  if(diag_true_func_plot){
    fxP <- SIR_f(xP=xP, timepoint=timepoint, variable=variable, f_0=f_start)
    fxP_mat <- matrix(fxP, nrow=length(beta_grid), ncol=length(beta_grid))
    
    SD_fx_vec <- (ED_fx_vec - fxP) / sqrt(VarD_fx_vec)
    # Fixing diagnostic values on boundary (perfect knowledge => Var_D[f(x)]=0)
    beta_0_points <- which(xP[,1]==0)
    SD_fx_vec[beta_0_points] <- 0
    
    SD_fx_mat <- matrix(SD_fx_vec, nrow=length(beta_grid), ncol=length(beta_grid))
    
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_", EmulatorType_short, "_SDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, "_E_f", E_f, "_E_df", E_df,
                          "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
        width = wd, height = ht)
    SD_fx_plot <- SIR_emul_fill_cont(cont_mat=SD_fx_mat, beta_grid=beta_grid, gamma_grid=gamma_grid,
                                     type="diag", xD=xD, cont_levs=seq(-3,3,0.25), 
                                     plot_boundary=TRUE, plot_boundary_line=boundary, 
                                     main=SD_main,
                                     cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
    dev.off()
    
    # True function plots
    png(filename = paste0(filepath, "SIR_", variable, "_timepoint_", timepoint, "_TrueFunction_fx.png"),
        width = wd, height = ht)
    
    true_main <- bquote(bold(atop("True Function f(x)", 
                                  "for SIR output " ~ .(variable) ~ "at timepoint" ~ t==.(as.character(timepoint)))))
    
    SIR_emul_fill_cont(cont_mat=fxP_mat, beta_grid=beta_grid, gamma_grid=gamma_grid, 
                       type=variable, xD=xD,cont_levs=cont_levs_exp,
                       plot_boundary=TRUE, plot_boundary_line=boundary,
                       main=true_main, 
                       cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
    dev.off()
  }
}

