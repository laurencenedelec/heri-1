####**********************************************************************
####  Written and Developed by: 
####**********************************************************************
####
####    Julien Duvanel, Copyright 2014
####    email: duvanel@stanford.edu
####
####**********************************************************************
####**********************************************************************

####**********************************************************************
####
####  Diagnostic of methods
####
####**********************************************************************

library(ggplot2)
library(grid)

#' Compare methods via ScatterPlot
#' 
#' @title Scatter plot methods
#' @param full.res results obtained after estimating heritability
#' @return nothing but export a pdf
#' @author Julien Duvanel
PlotScatterMethods <- function(full.res, heridiff.file, editplotname.file, pdf.path = "") {
    
    # Used to give dynamic pdf name file
    datetime.stamp <- format(Sys.time(), "%d%m%Y_%H%M%S")
    
    
    # Give correct rownames (used to filter afterwards)
    for(i in 1:length(full.res)) {
        rownames(full.res[[i]]) <- as.character(full.res[[i]][,1])
    }
    
    # Load HeriDiff (from Chiara) and give correct rownames
    HeriDiff <- read.csv(heridiff.file)
    EDITED_PLOT_NAMES <- read.csv(editplotname.file, header = TRUE)
    
    # Select only the columns we are interested in
    HERIDIFF <- data.frame(Phenotypes = as.character(HeriDiff$FieldName), 
                           Value = as.numeric(HeriDiff$H2r))
    # give correct rownames
    rownames(HERIDIFF) <- as.character(HeriDiff$FieldName)
    
    # Filter data where we have data
    # (because we can only compare if we have data in both datasets)
    index <- intersect(HERIDIFF$Phenotypes, full.res[[1]][,1])
    HERIDIFF <- HERIDIFF[as.vector(index), ]
    HERIDIFF <- merge(HERIDIFF, EDITED_PLOT_NAMES, by.x = "Phenotypes", by.y = "FieldName", sort=F)
    
    # Required to have the correct "data format"
    for(i in 1:length(full.res)) {
        full.res[[i]] <- full.res[[i]][as.vector(index), ]
    }
    
    # Now we have 2 datasets with the same number of rows

    # For each methods stored in full.res
    # we go through an plot a scatterplot between
    # this method and heridiff
    for(i in 1:length(full.res)) {
        
        # Get data from raw results (thanks to the cluster)
        # Give correct name and transform data into numeric
        dat <- full.res[[i]][,-1]
        colnames(dat) <- c("LR", "DCOV", "DCOVPLUS")
        class(dat) <- "numeric"
        dat <- as.data.frame(dat)
        
        # For each combination
        p <- list()
        for(j in 1:ncol(dat)) {
            
            # Create the dataframe we want to plant
            datafrm <- data.frame(dat[, j], HeriDiff = HERIDIFF$Value, HERIDIFF$PlotGroup, HERIDIFF$PlotDomain)
            colnames(datafrm) <- c(colnames(dat)[j], "HeriDiff", "PlotGroup", "PlotDomain")
            
            # Store the plot (because we're gonna use it latter)
            p[[j]] <- ggplot(datafrm, 
                             aes_string(x = colnames(dat)[j],
                                        y = "HeriDiff",
                                        color = "PlotGroup")) +
                geom_point(aes_string(shape = "PlotGroup")) +    # Use hollow circles
                #geom_smooth(method=lm,   # Add linear regression line
                #            se=FALSE) +  # Don't add shaded confidence region
                GetCustomGgplotTheme()
        }
        
        # Plot results, 3 columns (= 3 methods)
        pdf(file = paste0(pdf.path, names(full.res)[i], "_", datetime.stamp, ".pdf"), width = 17, height = 7)
            
            # Create a grid with 2 rows, length(p) columns
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(2, length(p), heights = unit(c(0.5, 5), "null"))))   
            
            # First row contains the title
            grid.text(label = paste0("Comparison ", names(full.res)[i], " vs HeriDiff"),  
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:length(p)),
                gp = gpar(fontsize = 25))
            
            # Second row and length(p) columns contains the plot we saved earlier
            for(k in 1:length(p)) {
                print(p[[k]] + theme(legend.position = "none"), vp = vplayout(x=2,y=k))
            }
        
        dev.off()
        
    }
    
}
