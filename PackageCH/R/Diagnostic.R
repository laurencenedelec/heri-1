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
#' @param heridiff.file file containing estimates given by Chiara
#' @param editplotname.file file containing PlotGroup (needed to plot group)
#' @param pdf.path where the pdf are exported
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
            datafrm <- data.frame(dat[, j], HeriDiff = HERIDIFF$Value, 
                                  PlotGroup = HERIDIFF$PlotGroup, 
                                  PlotDomain = HERIDIFF$PlotDomain)
            colnames(datafrm) <- c(colnames(dat)[j], "HeriDiff", "PlotGroup", "PlotDomain")
            
            # Store the plot (because we're gonna use it latter)
            p[[j]] <- ggplot(datafrm, 
                             aes_string(x = colnames(dat)[j],
                                        y = "HeriDiff",
                                        color = "PlotGroup")) +
                geom_point(aes_string(shape = "PlotGroup")) +    # Use hollow circles
                #geom_smooth(method=lm,   # Add linear regression line
                #            se=FALSE) +  # Don't add shaded confidence region
                scale_shape_manual(values = 0:length(unique(datafrm$PlotGroup))) +
                guides(col = guide_legend(ncol = 1)) +
                GetCustomGgplotTheme()
        }
        
        # Plot results, 3 columns (= 3 methods)
        pdf(file = paste0(pdf.path, names(full.res)[i], "_", datetime.stamp, ".pdf"), width = 17, height = 7)
            
            # Create a grid with 2 rows, length(p) columns
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(2, length(p)+1, 
                                                       heights = unit(c(0.5, 5), "null"),
                                                       widths = unit(c(0.28, 0.28, 0.28, 0.16), "null"))))   
            
            # First row contains the title
            grid.text(label = paste0("Comparison ", names(full.res)[i], " vs HeriDiff"),  
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:(length(p)+1)),
                gp = gpar(fontsize = 25))
            
            # Second row and length(p) columns contains the plot we saved earlier
            for(k in 1:length(p)) {
                print(p[[k]] + theme(legend.position = "none"), vp = vplayout(x=2,y=k))
            }
        
            pushViewport(vp = vplayout(x=2, y = length(p)+1 ))
                grid.draw(GetLegendFromGgplot2(a.gplot = p[[1]]))    
            upViewport(0)

        dev.off()
        
    }
    
}


#' Compare method dcov+ via ScatterPlot
#' 
#' @title Scatter plot dcov+ using different matrices
#' @param full.res results obtained after estimating heritability
#' @param editplotname.file file containing PlotGroup (needed to plot group)
#' @param pdf.path where the pdf are exported
#' @return nothing but export a pdf
#' @author Julien Duvanel
PlotScatterDcovPlus <- function(full.res, editplotname.file, pdf.path = "") {
    
    # Used to give dynamic pdf name file
    datetime.stamp <- format(Sys.time(), "%d%m%Y_%H%M%S")
    
    # Give correct rownames (used to filter afterwards)
    for(i in 1:length(full.res)) {
        rownames(full.res[[i]]) <- as.character(full.res[[i]][,1])
    }
    
    # Load plot group
    EDITED_PLOT_NAMES <- read.csv(editplotname.file, header = TRUE)
    PlotNames <- data.frame(FieldName = as.character(EDITED_PLOT_NAMES$FieldName), 
                            PlotGroup = as.character(EDITED_PLOT_NAMES$PlotGroup))
    
    # Required to have the correct "data format"
    for(i in 1:length(full.res)) {
        full.res[[i]] <- full.res[[i]][as.vector(index), ]
    }
    
    # Now we have 2 datasets with the same number of rows
    
    # For each methods stored in full.res
    # we go through an plot a scatterplot between
    # 2 distance matrices
    combi <- combn(seq(length(full.res)), m=2)
    
    # For each combination
    p <- list()
    for(i in 1:ncol(combi)) {
        
        # Get data from raw results (thanks to the cluster)
        # Give correct name and transform data into numeric
        tmp.res <- merge(full.res[[combi[1,i]]][, c(1,4)], 
                         full.res[[combi[2,i]]][, c(1,4)],
                         by.x = "Phenotypes", by.y = "Phenotypes")
        colnames(tmp.res) <- c("Phenotypes", 
                               names(full.res)[combi[1,i]],
                               names(full.res)[combi[2,i]])
        
        dat <- merge(tmp.res, PlotNames, by.x = "Phenotypes", by.y = "FieldName")
        dat[,2] <- as.numeric(as.vector(dat[,2]))
        dat[,3] <- as.numeric(as.vector(dat[,3]))
        
        # Store the plot (because we're gonna use it latter)
        p[[i]] <- ggplot(dat, 
                         aes_string(x = colnames(dat)[2],
                                    y = colnames(dat)[3],
                                    color = "PlotGroup")) +
            geom_point(aes_string(shape = "PlotGroup")) +    # Use hollow circles
            #geom_smooth(method=lm,   # Add linear regression line
            #            se=FALSE) +  # Don't add shaded confidence region
            scale_shape_manual(values = 0:length(unique(dat$PlotGroup))) +
            guides(col = guide_legend(ncol = 1)) +
            GetCustomGgplotTheme()
        
    }
        
    # Plot results, 3 columns (= 3 methods)
    pdf(file = paste0(pdf.path, "DcovPlus_", datetime.stamp, ".pdf"), width = 17, height = 14)
    
        # Create a grid with 2 rows, length(p) columns
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(3, 4, 
                                                   heights = unit(c(0.5, 5, 5), "null"),
                                                   widths = unit(c(0.28, 0.28, 0.28, 0.16), "null"))))   
        
        # First row contains the title
        grid.text(label = paste0("Comparison dcovplus between methods"),  
                  vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4),
                  gp = gpar(fontsize = 25))
        
        # Second row and length(p) columns contains the plot we saved earlier
        for(k in 1:ncol(combi)) {
            print(p[[k]] + theme(legend.position = "none"), vp = vplayout(x = 2 + (k-1) %/% 3 ,y = (k-1) %% 3 + 1))
        }
        
        pushViewport(vp = vplayout(x = 2:3, y = 4 ))
        grid.draw(GetLegendFromGgplot2(a.gplot = p[[1]]))    
        upViewport(0)
        
    dev.off()
        
}
