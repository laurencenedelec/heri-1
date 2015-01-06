
PlotHeatmapCorrelationMatrix <- function(X, 
                                         correlation = FALSE) {
  
  if(correlation == TRUE) {
    # Get the correlation matrix of X
    z <- cor(X)
  } else {
    z <- X
  }

  # Get colors between white and black
  col <- colorRampPalette(c("white", "black"), space="rgb")
  
  breaks <- seq(min(z), max(z), length.out=10)
  
  image(z, col=col(length(breaks)-1), breaks=breaks, xaxt="n", yaxt="n", ylab="", xlab="", useRaster=T)

}
