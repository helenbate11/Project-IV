# Modification of emul_fill_cont() function from "Emulator Output Plot Generators.R"
# to plot with SIR colour schemes and plot and axis titles

################################################################
# Function to plot 2d emulator outputs with SIR colour schemes and titles
################################################################

library(viridisLite)
library(RColorBrewer)
exp_cols <- magma
var_cols <- function(n) {hcl.colors(n, "YlOrRd", rev=TRUE)}
diag_cols <- turbo

SIR_emul_fill_cont <- function(cont_mat,            # matrix of values over which we plot a contour plot
                               cont_levs=NULL,      # contour levels (NULL: auto selection)
                               nlev=20,             # the approx number of contour levels for auto selection
                               plot_xD=TRUE,        # TRUE or FALSE: plot the design run positions 
                               plot_boundary=FALSE, # TRUE or FALSE: plot the boundary line
                               plot_boundary_line=list("x1"=NULL, "x2"=NULL), # (only define if plot_boundary=TRUE) boundary definitions: x1 = boundary$x1, x2 = boundary$x2, (one coordinate=NULL)
                               xD=NULL,             # (only define if plot_xD=TRUE) the matrix of run positions
                               xD_col="green",      # colour of the design runs and boundary line if plotted
                               beta_grid,           # grid edge locations defining xP[,1], the beta input space
                               gamma_grid,          # grid edge locations defining xP[,2], the gamma input space
                               type="default",      # the plot type: "default"=exp_cols, 
                                                    #                "S"=reds (expectation plot for the SIR compartment S), 
                                                    #                "I"=orange (expectation plot for the SIR compartment I), 
                                                    #                "R"=turquoise (expectation plot for the SIR compartment R), 
                                                    #                "var"=var_cols, 
                                                    #                "diag"=diag_cols
                               ...                  # any extra arguments passed to filled.contour()
) {

  # Define contour levels if necessary
  if(is.null(cont_levs)){
    cont_levs <- pretty(cont_mat, n=nlev)
  }
  
  spectral <- colorRampPalette(brewer.pal(11, "Spectral"))(200)
  
  R_cols <- function(n) {
    pal <- hcl.colors(200, "Spectral", rev = TRUE)
    pal <- pal[1:100]                    # blue → yellow section
    colorRampPalette(pal)(n)
  }
  
  I_cols <- function(n) {
    # Base golden yellow tones
    pal <- hcl.colors(200, "YlOrBr", rev = FALSE)
    pal <- pal[40:180]  # remove very dark and extremely pale colors
    
    # Last base color for reference
    last_base <- pal[length(pal)]
    
    # Fade to pale yellow (monotonically increasing in lightness)
    fade <- colorRampPalette(
      c(last_base,
        hcl(h = 60, c = 30, l = 100),  # pale golden yellow
        hcl(h = 65, c = 10, l = 100)   # almost-white cream
      ),
      space = "Lab"
    )(40)
    
    # Combine base + fade
    full_pal <- c(pal, fade)
    
    # Interpolate to n colors
    colorRampPalette(full_pal, space = "Lab", bias = 1.2)(n)
  }
  
  S_cols <- function(n) {
    pal <- hcl.colors(200, "Peach", rev = FALSE)
    
    # create a long smooth fade to white
    fade <- colorRampPalette(c(pal[150], "#fff7f5"), space = "Lab")(80)
    
    full_pal <- c(pal[1:150], fade)
    
    colorRampPalette(full_pal, bias = 1.3)(n)
  }
  
  if(type=="default"){
    color.palette <- exp_cols
  }else if(type=="S"){
    color.palette <- S_cols
  }else if(type=="I"){
    color.palette <- I_cols
  }else if(type=="R"){
    color.palette <- R_cols
  }else if(type=="var"){
    color.palette <- var_cols
  }else if(type=="diag"){
    color.palette <- diag_cols
    xD_col <- "purple"
  }
  
  # create filled contour plot
  filled.contour(beta_grid, gamma_grid, cont_mat, levels=cont_levs, 
                 xlab=bquote("Beta Input" ~ beta), 
                 ylab=bquote("Gamma Input" ~ gamma), 
                 color.palette=color.palette,
                 xaxs="i",
                 yaxs="i",
                 ..., 
                 plot.axes={
                   axis(1)
                   axis(2)
                   contour(beta_grid, gamma_grid, cont_mat, add=TRUE, levels=cont_levs, lwd=0.8)
                   if(plot_xD){
                     points(xD, pch=21, col=1, bg=xD_col, cex=1.5)
                   }
                   if(plot_boundary){
                     if(!is.null(plot_boundary_line$x1)){
                       abline(v=plot_boundary_line$x1, col=xD_col, lwd=5)
                     }
                     if(!is.null(plot_boundary_line$x2)){
                       abline(h=plot_boundary_line$x2, col=xD_col, lwd=5)
                     }
                   }
                 })
}

