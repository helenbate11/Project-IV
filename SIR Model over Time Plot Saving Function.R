# The function to save plots of the SIR model compartment number change over time
# into a folder with appropriate filename

################################################################
# SIR number of people in each compartment change over time 
# curves - plotting function
################################################################

# Saves with following file names:
# paste0(filepath, "SIR_model_", "S0_", f_start["S"], "_I0_", f_start["I"], "_R0_", f_start["R"], "_beta_", beta, "_gamma_", gamma, ".png")
# paste0(filepath, "SIR_model_derivs_", "S0_", f_start["S"], "_I0_", f_start["I"], "_R0_", f_start["R"], "_beta_", beta, "_gamma_", gamma, ".png")

SIR_cols <- c("red", "darkorange", "darkturquoise")

saving_SIR_plots <- function(output,       # SIR_model output from SIR_model() in "SIR Model Code.R"
                             beta,         # the beta parameter value
                             gamma,        # the gamma parameter value
                             f_start,      # the intial values (S_0, I_0, R_0)
                             times=seq(0, 2.5, length=200), # the sequence of timepoints at which to evaluate the SIR model
                             wd=800, ht=620,                # the height and width of the plots
                             ylim=NULL,    # the y-axis plot limits
                             filepath="",  # the filepath to the folder in which to save plots (MUST SET)
                             deriv_plot=FALSE,              # TRUE or FALSE: plot change over time over derivatives of model wrt time 
                             highlight_timepoint=FALSE,     # TRUE or FALSE add vertical line to plot at certain timepoint
                             timepoint=NULL                 # (defined only if highligh_timepoint = TRUE) the timepoint of interest
                             ){
  
  if(is.null(ylim)){
    ylim <- c(0, 1.15*max(output[,c("S", "I", "R")]))
  }
  png(filename = paste0(filepath, "SIR_model_", "S0_", f_start["S"], "_I0_", 
                        f_start["I"], "_R0_", f_start["R"], "_beta_", 
                        beta, "_gamma_", gamma, ".png"),
      width = wd, height = ht)
  plot(output[,"time"], output[,"S"], ty="n", 
       ylim=ylim,
       xlab="Time", ylab="Number in Each Compartment", 
       main=bquote(bold("SIR Curves (" *
                          S["0"]==bold(.(as.character(f_start["S"]))) * "," ~ 
                          I["0"]==bold(.(as.character(f_start["I"]))) * "," ~ 
                          R["0"]==bold(.(as.character(f_start["R"]))) * "," ~
                          beta==bold(.(as.character(beta))) * "," ~
                          gamma==bold(.(as.character(gamma))) * ")")),
       cex.lab=1.4, cex.main=1.5, mgp = c(2, 0.8, 0))
  for(i in c("S", "I", "R")){
    coli <- which(i==c("S", "I", "R")) 
    lines(output[,"time"], output[,i], col=SIR_cols[coli], lwd=3)
  }
  if(beta==0){
    abline(h=f_start["I"]+f_start["R"], col="black", lwd=1, lty=2)
    legend("topright", legend=c(paste("Susceptible", expression(S(t))), 
                                paste("Infected", expression(I(t))), 
                                paste("Recovered", expression(R(t))), 
                                paste0("Asymptote of ", expression(R(t)))), 
           col=c(SIR_cols, "black"), lty=c(1,1,1,2), lwd=c(2.5,2.5,2.5,1), cex=1.3)
  }else if(gamma==0){
    abline(h=f_start["S"]+f_start["I"], col="black", lwd=1, lty=2)
    legend("topright", legend=c(paste("Susceptible", expression(S(t))), 
                                paste("Infected", expression(I(t))), 
                                paste("Recovered", expression(R(t))), 
                                paste0("Asymptote of ", expression(I(t)))),
           col=c(SIR_cols, "black"), lty=c(1,1,1,2), lwd=c(2.5,2.5,2.5,1), cex=1.3)
  }else{
    legend("topright", legend=c(paste("Susceptible", expression(S(t))), 
                                paste("Infected", expression(I(t))), 
                                paste("Recovered", expression(R(t)))),  
           col=SIR_cols, lty=1, lwd=2.5, cex=1.3)
  }
  if(highlight_timepoint){
    abline(v=timepoint, col="darkgrey", lwd=2, lty=2)
  }
  dev.off()
  
  if(deriv_plot){
    png(filename = paste0(filepath, "SIR_model_derivs_", "S0_", f_start["S"], "_I0_", 
                          f_start["I"], "_R0_", f_start["R"], "_beta_", 
                          beta, "_gamma_", gamma, ".png"),
        width = wd, height = ht)
    plot(output[,"time"], output[,"dS"], ty="n", 
         ylim=c(1.1*min(output[,c("dS", "dI", "dR")]), 
                1.1*max(output[,c("dS", "dI", "dR")])),
         xlab="Time", ylab="Rate of Change of Each Compartment", 
         main=bquote(bold("Derivatives of the SIR Curves (" *
                            S["0"]==bold(.(as.character(f_start["S"]))) * "," ~ 
                            I["0"]==bold(.(as.character(f_start["I"]))) * "," ~ 
                            R["0"]==bold(.(as.character(f_start["R"]))) * "," ~
                            beta==bold(.(as.character(beta))) * "," ~
                            gamma==bold(.(as.character(gamma))) * ")")), 
         cex.lab=1.5, cex.main=1.5, mgp = c(2, 0.8, 0))
    for(i in c("dS", "dI", "dR")){
      coli <- which(i==c("dS", "dI", "dR")) 
      lines(output[,"time"], output[,i], col=SIR_cols[coli], lwd=3)
    }
    legend("topright", legend=c(expression(frac(dS, dt)(t)), expression(frac(dI, dt)(t)), expression(frac(dR, dt)(t))),
           col=SIR_cols, lty=1, lwd=2.5, y.intersp = 1.5, cex=1.3)
    if(highlight_timepoint){
      abline(v=timepoint, col="darkgrey", lwd=2, lty=2)
    }
    dev.off()
  }
}
