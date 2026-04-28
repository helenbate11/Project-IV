# The function to save emulator output plots into a folder 
# with appropriate filename

# Note: file names often depend on global parameters of theta, sigma, E_f and E_df
# so save the relevant plots immediately after running that emulator

# Note: calls on global function of interest f()

################################################################
# 1d Emulator Output Plots Saving Function
# (for the 1d Simple and Derivative Emulators)
################################################################

# Saves with following file names:
# paste0(function_name, "_Variances_psd_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png")

saving_1d_BL_plot <- function(em_out,            # emulator output matrix
                              xP,                # vector of inputs where the emulator was evaluated
                              xD,                # the run input locations
                              D,                 # the run outputs, D = (f(x^1),...,f(x^n))
                              ylim=c(0,1.1),     # the y-axis range
                              EmulatorType,      # the type of emulator being plotted: "Simple" or "Derivative"
                              wd=690, ht=620,    # the height and width of the plots
                              filepath="",       # the filepath to the folder in which to save plots (MUST SET)
                              function_name="F"  # the model function name
                              ){
  png(filename = paste0(filepath, function_name, "_1d_", EmulatorType, "_theta", theta, "_sigma", sigma, 
                        ".png"),
      width = wd, height = ht)
  plot_BL_emulator_1d(em_out=em_out, xP=xP, xD=xD, D=D, ylim=ylim,
                      maintitle=paste("1d", EmulatorType, "Emulator Output"))
  dev.off()
}


################################################################
# 1d Simple vs Derivative Emulator Variance Comparison Plot Saving Function
# (for the 1d Simple and Derivative Emulators)
################################################################

# Saves with following file name:
# paste0(function_name, "_1d_", EmulatorType, "_theta", theta, "_sigma", sigma, ".png")
# paste0(function_name, "_SDs_psd_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png")

saving_var_plots <- function(simple_em_out,    # the 1d Simple Emulator output matrix
                             deriv_em_out,     # the 1d Derivative Emulator output matrix
                             sigma,            # the prior standard deviation
                             plot_xD=FALSE,    # TRUE or FALSE: plot the design run positions 
                             xD,               # the design run positions
                             xP,               # the input parameters over which we emulated
                             col_vec=c("black", colorRampPalette(c("darkred", "red"))(3)[2], adjustcolor("red", alpha.f = 0.8)),
                             xD_col="green",   # the colour of the design runs if plotted
                             wd=800, ht=620,   # the width and height of saved plots
                             filepath="",      # the filepath to the folder in which to save plots (MUST SET)
                             function_name="F" # the model function name
                             ){
  ylim_var <- c(0, sigma^2+0.1)
  ylim_sd <- c(0, sigma+0.15)
  
  # Variance plot
  png(filename = paste0(filepath, function_name, "_Variances_psd_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png"),
      width = wd, height = ht)
  
  plot(xP, rep(sigma^2, length(xP)), ylim=ylim_var, type="l", col=col_vec[1], lwd=3,
       xlab="Input parameter x", ylab="Variance Var[f(x)]",
       main=bquote(bold("Variance of the Prior and the Simple and Derivative Emulators")),
       cex.main=1.5, cex.lab=1.4, mgp = c(2.5, 1, 0))
  lines(xP, simple_em_out[,"VarD_f(x)"], col=col_vec[2], lwd=3)
  lines(xP, deriv_em_out[,"VarD_f(x)"], col=col_vec[3], lwd=3)
  if(plot_xD){
    points(xD, rep(0, length(xD)), pch=21, col=1, bg=xD_col, cex=1.5)
    legend("topright", legend=c(bquote("Prior Variance" ~ sigma^2), 
                                bquote("Variance of Simple Emulator" ~ Var[D] * "[f(x)]"),
                                bquote("Variance of Derivative Emulator" ~ Var[D] * "[f(x)]"),
                                "Model Evaluations"),
           lty=c(1,1,1,NA), pch=c(NA, NA, NA, 16), col=c(col_vec, xD_col), 
           lwd=c(2.5, 2.5, 2.5, NA), pt.cex=1.3, cex=1.3)
  }else{
    legend("topright", legend=c(bquote("Prior Variance" ~ sigma^2), 
                                bquote("Variance of Simple Emulator" ~ Var[D] * "[f(x)]"),
                                bquote("Variance of Derivative Emulator" ~ Var[D] * "[f(x)]")),
           lty=c(1,1,1), pch=c(NA, NA, NA), col=col_vec,
           lwd=c(2.5, 2.5, 2.5), cex=1.3)
  }
  
  dev.off()
  
  # SD plot
  png(filename = paste0(filepath, function_name, "_SDs_psd_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png"),
      width = wd, height = ht)
  
  plot(xP, rep(sigma, length(xP)), ylim=ylim_sd, type="l", col=col_vec[1], lwd=3,
       xlab="Input parameter x", ylab="Standard Deviation SD[f(x)]",
       main=bquote(bold("Standard Deviation of the Prior and the Simple and Derivative Emulators")),
       cex.main=1.5, cex.lab=1.4, mgp = c(2.5, 1, 0))
  lines(xP, sqrt(simple_em_out[,"VarD_f(x)"]), col=col_vec[2], lwd=3)
  lines(xP, sqrt(deriv_em_out[,"VarD_f(x)"]), col=col_vec[3], lwd=3)
  if(plot_xD){
    points(xD, rep(0, length(xD)), pch=21, col=1, bg=xD_col, cex=1.5)
    legend("topright", legend=c(bquote("Prior Standard Deviation" ~ sigma), 
                                bquote("Standard Deviation of Simple Emulator" ~ SD[D] * "[f(x)]"),
                                bquote("Standard Deviation of Derivative Emulator" ~ SD[D] * "[f(x)]"),
                                "Model Evaluations"),
           lty=c(1,1,1,NA), pch=c(NA, NA, NA, 16), col=c(col_vec, xD_col), 
           lwd=c(2.5, 2.5, 2.5, NA), pt.cex=1.3, cex=1.3)
  }else{
    legend("topright", legend=c(bquote("Prior Standard Deviation" ~ sigma), 
                                bquote("Standard Deviation of Simple Emulator" ~ SD[D] * "[f(x)]"),
                                bquote("Standard Deviation of Derivative Emulator" ~ SD[D] * "[f(x)]")),
           lty=c(1,1,1), pch=c(NA, NA, NA), col=col_vec,
           lwd=c(2.5, 2.5, 2.5), cex=1.3)
  }
  
  dev.off()
}


################################################################
# 2d Emulator Output Plots Saving Function
# (for the 2d Simple and Derivative Emulators)
################################################################

# Saves with following file name:
# paste0(function_name, "_", EmulatorType, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png")
# etc.

saving_fxD_plots <- function(em_out,                    # the emulator output matrix
                             dfdx_plots=FALSE,          # TRUE of FALSE: save the plots from the emulation of the derivative of the function 
                             dfxP1_mat=NULL,            # (define only if dfdx_plots = TRUE) the emulator additional output for the emulation of df/dx_1
                             dfxP2_mat=NULL,            # (define only if dfdx_plots = TRUE) the emulator additional output for the emulation of df/dx_2
                             xD,                        # the design run positions
                             x_grid,                    # the coordinate grid from which xP defined
                             cont_levs_mat=matrix(c(-2,2,0.2,   # contour colour levels for the E_D[f(x)] plot
                                                    0,1,0.1,    # contour colour levels for the Var_D[f(x)] plot
                                                    -9,9,0.5,   # contour colour levels for the E_D[df/dx_1(x)] plot
                                                    0,13,1,     # contour colour levels for the Var_D[df/dx_1(x)] plot
                                                    -5,5,0.5,   # contour colour levels for the E_D[df/dx_2(x)] plot
                                                    0,5,0.5     # contour colour levels for the Var_D[df/dx_2(x)] plot
                                                    ), ncol=3, byrow=TRUE),
                             wd=800, ht=620,            # the width and height of the saved plots
                             filepath="",               # the filepath to the folder in which to save plots (MUST SET)
                             function_name="F",         # the model function name
                             EmulatorType="Derivative"  # the emulator type: "Simple" or "Derivative"
                             ){
  
  ED_fx_mat <- matrix(em_out[,"ExpD_f.x."], nrow=length(x_grid), ncol=length(x_grid)) 
  VarD_fx_mat <- matrix(em_out[,"VarD_f.x."], nrow=length(x_grid), ncol=length(x_grid)) 
  fxP_mat <- matrix(f(xP), nrow=length(x_grid), ncol=length(x_grid))
  
  # Emulator plots
  png(filename = paste0(filepath, function_name, "_", EmulatorType, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png"),
      width = wd, height = ht)
  ED_fx_plot <- emul_fill_cont(cont_mat=ED_fx_mat, 
                               cont_levs=seq(cont_levs_mat[1,1],cont_levs_mat[1,2], cont_levs_mat[1,3]), 
                               xD=xD, x_grid=x_grid, color.palette=exp_cols, 
                               main=bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Expectation " ~ E[D] * "[f(x)]")),
                               cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  png(filename = paste0(filepath, function_name, "_", EmulatorType, "_VarDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png"),
      width = wd, height = ht)
  VarD_fx_plot <- emul_fill_cont(cont_mat=VarD_fx_mat, 
                                 cont_levs=seq(cont_levs_mat[2,1],cont_levs_mat[2,2], cont_levs_mat[2,3]), 
                                 xD=xD, x_grid=x_grid, color.palette=var_cols, 
                                 main=bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Variance " ~ Var[D] * "[f(x)]")),
                                 cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  S_diag_mat <- (ED_fx_mat - fxP_mat) / sqrt(VarD_fx_mat)
  png(filename = paste0(filepath, function_name, "_", EmulatorType, "_SDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png"),
      width = wd, height = ht)
  SD_fx_plot <- emul_fill_cont(cont_mat=S_diag_mat, cont_levs=seq(-3,3,0.25), 
                               xD=xD, x_grid=x_grid, xD_col="purple", color.palette=diag_cols,
                               main=bquote(bold(.(EmulatorType) ~ "Emulator Diagnostics " ~ S[D] * "[f(x)]")),
                               cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  # True function plots
  png(filename = paste0(filepath, function_name, "_TrueFunction_fx.png"),
      width = wd, height = ht)
  emul_fill_cont(cont_mat=fxP_mat, 
                 cont_levs=seq(cont_levs_mat[1,1],cont_levs_mat[1,2], cont_levs_mat[1,3]), xD=xD, x_grid=x_grid,
                 color.palette=exp_cols, main="True Computer Model Function f(x)",
                 cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
}





################################################################
# 2d Emulator Output Plots Saving Function
# (for the 2d Boundary, Boundary with Runs, Boundary with Derivative
# and Boundary with Derivative and Runs Emulators)
################################################################

# Saves with following file name:
# paste0(function_name, "_", EmulatorType, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, ".png")
# etc.

saving_boundary_with_runs_plots <- function(em_out,                    # the emulator output matrix
                                            boundary,                  # boundary definition: list(x1 = boundary$x1, x2 = boundary$x2)
                                            xD=NULL,                   # the design run positions
                                            x_grid,                    # the coordinate grid from which xP defined
                                            cont_levs_mat=matrix(c(-2,2,0.2,   # contour colour levels for the E_D[f(x)] plot
                                                                   0,1,0.1,    # contour colour levels for the Var_D[f(x)] plot
                                            ), ncol=3, byrow=TRUE),
                                            wd=800, ht=620,            # the width and height of the saved plots
                                            filepath="",               # the filepath to the folder in which to save plots (MUST SET)
                                            function_name="F",         # the model function name
                                            EmulatorType="Boundary"  # the emulator type: "Boundary", "Boundary with Runs", "Boundary and Derivative on the Boundary", "Boundary and Derivative on the Boundary with Runs")
                                            ){

  ED_fx_vec <- em_out[,1]
  VarD_fx_vec <- em_out[,2]
  fxP_vec <- f(xP)
  
  ED_fx_mat <- matrix(ED_fx_vec, nrow=length(x_grid), ncol=length(x_grid)) 
  VarD_fx_mat <- matrix(VarD_fx_vec, nrow=length(x_grid), ncol=length(x_grid)) 
  fxP_mat <- matrix(fxP_vec, nrow=length(x_grid), ncol=length(x_grid))
  
  EmulatorType_short <- ifelse(EmulatorType=="Boundary", "Boundary", 
                               ifelse(EmulatorType=="Boundary with Runs", "Boundary with Runs", 
                                      ifelse(EmulatorType=="Boundary and Derivative on the Boundary", "BwD",
                                             "BwD with Runs")))
  
  if(EmulatorType =="Boundary"){
    ED_main <- bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Expectation " ~ E[K] * "[f(x)]"))
    VarD_main <- bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Variance " ~ Var[K] * "[f(x)]"))
    SD_main <- bquote(bold(.(EmulatorType) ~ "Diagnostics " ~ S[K] * "[f(x)]"))
  }else if(EmulatorType == "Boundary with Runs"){
    ED_main <- bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Expectation " ~ E[D] * "[f(x)]"))
    VarD_main <- bquote(bold(.(EmulatorType) ~ "Emulator Adjusted Variance " ~ Var[D] * "[f(x)]"))
    SD_main <- bquote(bold(.(EmulatorType) ~ "Diagnostics " ~ S[D] * "[f(x)]"))
  }else if(EmulatorType == "Boundary and Derivative on the Boundary"){
    ED_main <- bquote(bold(
      atop(
        .(EmulatorType) ~ "Emulator",
        "Adjusted Expectation " ~ E[K] * "[f(x)]"
      )
    ))
    VarD_main <- bquote(bold(
      atop(
        .(EmulatorType) ~ "Emulator",
        "Adjusted Variance " ~ Var[K] * "[f(x)]"
      )
    ))
    SD_main <- bquote(bold(
      atop(
        .(EmulatorType),
        "Diagnostics " ~ S[K] * "[f(x)]"
      )
    ))
  }else if(EmulatorType == "Boundary and Derivative on the Boundary with Runs"){
    ED_main <- bquote(bold(
      atop(
        .(EmulatorType) ~ "Emulator",
        "Adjusted Expectation " ~ E[D] * "[f(x)]"
      )
    ))
    VarD_main <- bquote(bold(
      atop(
        .(EmulatorType) ~ "Emulator",
        "Adjusted Variance " ~ Var[D] * "[f(x)]"
      )
    ))
    SD_main <- bquote(bold(
      atop(
        .(EmulatorType),
        "Diagnostics " ~ S[D] * "[f(x)]"
      )
    ))
  }
  
  # Emulator plots for E_D[f(x)]
  png(filename = paste0(filepath, function_name, "_", EmulatorType_short, "_EDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, 
                        "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
      width = wd, height = ht)
  ED_fx_plot <- emul_fill_cont(cont_mat=ED_fx_mat, 
                               cont_levs=seq(cont_levs_mat[1,1],cont_levs_mat[1,2], cont_levs_mat[1,3]), 
                               plot_boundary=TRUE, plot_boundary_line=boundary, xD=xD, x_grid=x_grid, color.palette=exp_cols, 
                               main=ED_main,
                               cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  png(filename = paste0(filepath, function_name, "_", EmulatorType_short, "_VarDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, 
                        "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
      width = wd, height = ht)
  VarD_fx_plot <- emul_fill_cont(cont_mat=VarD_fx_mat, 
                                 cont_levs=seq(cont_levs_mat[2,1],cont_levs_mat[2,2], cont_levs_mat[2,3]), 
                                 plot_boundary=TRUE, plot_boundary_line=boundary, xD=xD, x_grid=x_grid, color.palette=var_cols, 
                                 main=VarD_main,
                                 cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  S_diag_vec <- (ED_fx_vec - fxP_vec) / sqrt(VarD_fx_vec)
  # Fixing diagnostic values on boundary (perfect knowledge => Var_D[f(x)]=0)
  boundaryK_points <- which(xP[,1]==boundary$x1)
  boundaryL_points <- which(xP[,2]==boundary$x2)
  S_diag_vec[boundaryK_points] <- 0
  S_diag_vec[boundaryL_points] <- 0
  
  S_diag_mat <- matrix(S_diag_vec, nrow=length(x_grid), ncol=length(x_grid))
  
  png(filename = paste0(filepath, function_name, "_", EmulatorType_short, "_SDfx_theta", theta[1], "_", theta[2], "_sigma", sigma, 
                        "_boundary_x1_", boundary$x1[1], "_", boundary$x1[2], "_x2_", boundary$x2[1], "_", boundary$x2[2], ".png"),
      width = wd, height = ht)
  SD_fx_plot <- emul_fill_cont(cont_mat=S_diag_mat, cont_levs=seq(-3,3,0.25), 
                               plot_boundary=TRUE, plot_boundary_line=boundary, xD=xD, x_grid=x_grid, xD_col="purple", color.palette=diag_cols,
                               main=SD_main,
                               cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
  
  # True function plots
  png(filename = paste0(filepath, function_name, "_TrueFunction_fx.png"),
      width = wd, height = ht)
  emul_fill_cont(cont_mat=fxP_mat, 
                 cont_levs=seq(cont_levs_mat[1,1],cont_levs_mat[1,2], cont_levs_mat[1,3]), 
                 plot_boundary=TRUE, plot_boundary_line=boundary, xD=xD, x_grid=x_grid,
                 color.palette=exp_cols, main="True Computer Model Function f(x)",
                 cex.main=1.5, cex.lab=1.5, mgp = c(2.5, 1, 0))
  dev.off()
}

