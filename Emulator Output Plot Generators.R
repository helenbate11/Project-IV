# Code to generate emulator output plots (as in the accompanying report)

################################################################
# Function to plot 1d emulator output 
# (for the 1d Simple and Derivative Emulators)
################################################################

my_purple <- rgb(0.52, 0.15, 0.5)

plot_BL_emulator_1d <- function(
    em_out,         # emulator output matrix
    xP,             # vector of inputs where the emulator was evaluated
    xD,             # the run input locations
    D,              # the run outputs, D = (f(x^1),...,f(x^n))
    ylim=c(0,1.1),  # the y-axis range
    maintitle=NULL, # the title of the plot (blank if NULL)
    ...             # other inputs for plot()
){
  
  plot(xP,em_out[,"ExpD_f(x)"],ylim=ylim,ty="l",col="hotpink",lwd=3,
       xlab="Input parameter x",
       ylab="Output f(x)",
       main=maintitle, 
       cex.main=1.5, cex.lab=1.4, mgp = c(2, 0.9, 0),
       ...) 
  lines(xP,em_out[,"ExpD_f(x)"]+3*sqrt(em_out[,"VarD_f(x)"]),col="darkorange",lwd=3)
  lines(xP,em_out[,"ExpD_f(x)"]-3*sqrt(em_out[,"VarD_f(x)"]),col="darkorange",lwd=3)
  
  # true function
  lines(xP,f(xP),col=my_purple, lwd=2.5, lty=1)
  
  # runs
  points(xD, D, pch=21, col=1, bg="green", cex=1.5)
  legend('topright', legend=c("Emulator Expectation",
                              "Emulator Prediction Interval",
                              "True function f(x)",
                              "Model Evaluations"),
         lty=c(1,1,1,NA), pch=c(NA,NA,NA,16), 
         col=c("hotpink", "darkorange", my_purple, "green"), 
         lwd=2.5, pt.cex=1.3, cex=1.3)
}




################################################################
# Function to plot 2d emulator outputs
################################################################

library(viridisLite)
exp_cols <- magma
var_cols <- function(n) {hcl.colors(n, "YlOrRd", rev=TRUE)}
diag_cols <- turbo

emul_fill_cont <- function(cont_mat,            # matrix of values over which we plot a contour plot
                           cont_levs=NULL,      # contour levels (NULL: auto selection)
                           nlev=20,             # the approx number of contour levels for auto selection
                           plot_xD=TRUE,        # TRUE or FALSE: plot the design run positions 
                           plot_boundary=FALSE, # TRUE or FALSE: plot the boundary line
                           plot_boundary_line=list("x1"=NULL, "x2"=NULL), # (only define if plot_boundary=TRUE) boundary definitions: x1 = boundary$x1, x2 = boundary$x2, (one coordinate=NULL)
                           xD=NULL,             # (only define if plot_xD=TRUE) the matrix of run positions
                           xD_col="green",      # colour of the design runs and boundary line if plotted
                           x_grid,              # the grid defining xP
                           ...                  # any extra arguments passed to filled.contour()
) {
  # Define contour levels if necessary
  if(is.null(cont_levs)){
    cont_levs <- pretty(cont_mat, n=nlev)
  }
  
  # create filled contour plot
  filled.contour(x_grid, x_grid, cont_mat, levels=cont_levs, 
                 xlab=expression(x[1]), ylab=expression(x[2]), ..., 
                 plot.axes={
                   axis(1)
                   axis(2)
                   contour(x_grid, x_grid, cont_mat, add=TRUE, levels=cont_levs, lwd=0.8)
                   if(plot_xD){
                     points(xD, pch=21, col=1, bg=xD_col, cex=1.5)
                   }
                   if(plot_boundary){
                     if(!is.null(plot_boundary_line$x1)){
                       abline(v=plot_boundary_line$x1, col=xD_col, lwd=3)
                     }
                     if(!is.null(plot_boundary_line$x2)){
                       abline(h=plot_boundary_line$x2, col=xD_col, lwd=3)
                     }
                   }
                 })
}




