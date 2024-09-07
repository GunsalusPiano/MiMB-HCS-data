########################
# R script is in support of the published MiMB book chapter (2024) and data analysis methods published in manuscript (https://doi.org/10.1038/s42003-022-04343-3) 
########################

########################
# Script is written by Yanthe E. Pearson
########################

########################
# R script for MiMB BookChapter - Step 3.1.3
########################

########################
# R script for figure  3: Sample of raw data as aggregate well medians   
########################

# clear variables
rm(list=ls())

########################
# libraries
########################

library(lattice)
library(RColorBrewer)
library(readr)


########################
# set working directory
########################

setwd("~/Location/of/Data")

########################
# (1) Import the annotated well median data file from panel A channels (Hoechst, syto14, mito, WGA)
########################

well_med_dat <- read_csv("data_step3.1.1_counts_and_medians_panelA.csv")

########################
# Check column names. Here column labels refer to metadata, cell counts and cellular features
# The data file includes median values for all 1848 wells (308 wells x 6 plates)
########################

names(well_med_dat)

########################
# identify 'empty' wells - these wells were not used for this HCS assay
########################

empty = which(grepl("empty", well_med_dat$Treatment)) 
# Check the empty wells
empty.check = well_med_dat[empty,]

########################
# identify low density wells - these wells have less than 40 cells as measured based on nucleus identification
########################

low.c = which(well_med_dat$cell_count < 40)

########################
# set the empty well data to NA - do not use these wells 
########################

well_med_dat[empty,13:dim(well_med_dat)[2]] = NA

########################
# numeric data
########################

well_med_dat.num = well_med_dat[,13:dim(well_med_dat)[2]]

########################
# identify each plate name
# "1A1" "1A2" "1A3" "2A1" "2A2" "2A3"
########################

plate.unique = sort(unique(well_med_dat$PlateNumber))

########################
# Initializing "lists" for heatmap plots of cell counts and well medians
# Well medians are formatted into data matrices according to well position and plate number 
# Matrices are saved to the lists for final heatmap plots
########################

datalist.plates = list()
datalist.plates.all = list()
datalist.plates.c = list()
datalist.plates.all.c = list()
merged.save = list()


########################
# j is a feature, here j = 8 is nucleus size. The user can pick any feature they are interested in plotting for visualizing well medians
# reasons to plot and visualize are to identify anomalies, examples include toxic treatments, plate effects, row effects and edge effects.
########################
 
  j = 8

########################
#  The index **i** is used to loop through each plate (here there are 6 plates)
#  plate.unique - unique plate names are "1A1" "1A2" "1A3" "2A1" "2A2" "2A3"
########################
  
  for (i in seq(1,length(plate.unique))) {     
    
    # identify wells within each plate
    plate.ind = which(well_med_dat$PlateNumber == plate.unique[i])
    plate.current = well_med_dat.num[plate.ind,j]
    plate.current.labels = well_med_dat[plate.ind,]
    plate.current.ord <- plate.current.labels[order(plate.current.labels$WellId),]
    
    # The data, and current feature
    plate.use <- plate.current.ord[,12+j]
    
    # The feature name
    feat.use <- names(plate.current.ord)[12+j]
    feat.name <- names(well_med_dat.num)[j]
    
    # The measured feature is formatted based on the 384-well plate layout
    # Raw data matrix 
    feature1_mat <- data.matrix(plate.use)
    feature1_mat1 <- matrix(feature1_mat, nrow = 14, ncol = 22, byrow = TRUE)
    rownames(feature1_mat1) <- LETTERS[2:15]
    colnames(feature1_mat1) <- 2:23
    
    # The measured feature is saved as a formatted matrix for plate **i**
    datalist.plates[[i]] <- feature1_mat1 
    
   # empty treatment wells set to NA
    ind.t = which(grepl("treated",plate.current.ord$Treatment)|grepl("empty",plate.current.ord$Treatment))
    ind.c = which(grepl("control",plate.current.ord$Treatment)|grepl("empty",plate.current.ord$Treatment))

    
    plate.use.c = plate.current.ord[,12+j] # setting all treatment wells to NA
    plate.use.c[ind.t,] = NA
    feature1_mat.ct <- data.matrix(plate.use.c)
    feature1_mat.ct1 <- matrix(feature1_mat.ct, nrow = 14, ncol = 22, byrow = TRUE) # raw data
    rownames(feature1_mat.ct1) <- LETTERS[2:15]
    colnames(feature1_mat.ct1) <- 2:23
    
    # save full matrix - control
    datalist.plates.c[[i]] <- feature1_mat.ct1 
    
  }
  
  
  
  ########################
  # Here each plate matrix is bound together
  ########################
  A = rbind(datalist.plates[[1]],datalist.plates[[2]], datalist.plates[[3]])
  B = rbind(datalist.plates[[4]],datalist.plates[[5]], datalist.plates[[6]])
  C = cbind(A,B)
  datalist.plates.all[[j]] <- C
  
  # Here, *full* use used to identify the full plate
  full = C
  
  ########################
  # Here each plate matrix is bound together, but the treatment values are muted. 
  # These plate matrices show only the control DMSO value per well.
  ########################
  A = rbind(datalist.plates.c[[1]],datalist.plates.c[[2]], datalist.plates.c[[3]])
  B = rbind(datalist.plates.c[[4]],datalist.plates.c[[5]], datalist.plates.c[[6]])
  C = cbind(A,B)
  datalist.plates.all.c[[j]] <- C
  
  # Here, *controls* is used used to identify the plate showing only the control wells
  controls = C

  # Here, the 12 plates are bound together into a larger matrix
  merged = cbind(controls, full)
  
  # merged.save is a list which saved the final data matrix used for the heatmap plot for feature **j**
  merged.save[[j]] = merged
  
  # all features
  feat.name <- names(well_med_dat.num) 

  # heatmap colors 
  #  https://cran.r-project.org/web/packages/colorBlindness/vignettes/colorBlindness.html
   new.palette <- colorRampPalette(brewer.pal(8, "RdBu"))(25)

  # Since we are only showing nucleus size here, the index i is equal to 8. For cell count, i = 1
  # If the user was interested in other features, i would be a different number
  
    
    i = 8
    
    ########################
    # the heatmap plot using the levelplot function and the lattice package
    ########################
  
    # Step 1: Call the pdf command to start the plot
    # pdf(file = "figure_step3.1.3_heatmap.pdf",   # The directory you want to save the file in
    #     width = 5, # The width of the plot in inches
    #     height = 8) # The height of the plot in inches
    # 
    levelplot(merged.save[[i]],panel = function(...){
    panel.levelplot(...)
    
    panel.abline(h = 22.5, col = "black", lwd = .5)
    panel.abline(h = 44.5, col = "black", lwd = .5)
    panel.abline(h = 67.5, col = "black", lwd = .5)
    
    panel.abline(v = 14.5, col = "black", lwd = .5)
    panel.abline(v = 28.5, col = "black", lwd = .5)
    },
    
    col.regions = new.palette,
    xlab =  list(label = "Nucleus size",cex = 1, rot = 180),
    ylab = NULL,
    #main = NULL,
    # sub="with log scales",
    # main=list(label="Nucleus size", cex=2),
    scale=list(x=list(at=c(1,14), rot=90, cex=1), y=list(at=c(1,22), rot=90, cex=1)), aspect = "fill")
  

  # dev.off()
  
  


