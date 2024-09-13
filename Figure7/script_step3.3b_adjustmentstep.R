
########################
# R script is in support of Methods in Molecular Biology book chapter and 
# data analysis methods published in manuscript (https://doi.org/10.1038/s42003-022-04343-3) 
# Script is written by Yanthe E. Pearson
########################

########################
# Quality control measures and statistical strategies to address the challenges of high-content phenotypic data
# by Yanthe E. Pearson and Kristin C. Gunsalus
########################

########################
# Per well adjustment of cells 
# 3.3 Plate normalization and cell standardization
########################

########################
# Description: 
# import adjustments file (adjustments for each well as calculated from the median polish analysis)
# import cell level data for panel A
# loop through features and plates and make adjustments to the single cell populations as needed.
# visualize the pre and post adjusted data (feature distributions)
########################

# clear variables
rm(list=ls())

########################
# libraries
########################

library(ggplot2)
library(RColorBrewer)
library(readr)
library(lattice)
library(gridExtra)

########################
# setwd("~/Directory.../.../")
########################

########################
# import data file data_step3.3b_adjustments.csv
########################

mp.adjust <- read_csv("data_step3.3b_adjustments.csv")
# mp.adjust = mp.adjust[,-2]
mp.adjust = data.frame(mp.adjust)

########################
# import panel A - cell level data
########################

cell.dat = read_csv("data_step3.3b_cell_panelA_subset.csv")
dat.sub = cell.dat[,1:32]

cell1 = cell.dat


    ########################  
    # Use *mp.adjust* and *data*
    ########################

    # Raw cell level data
    data = data.frame(cell1) 
  
    # identify number of features
    features.l = dim(data[,26:82])[2]
    
    # Panel A has 71 features
    features.l 
    
    # Create new data frame for saving newly adjusted cell features
    dat_modified = data.frame(data) 
    
    # 57 features in data frame *adjustment amount*
    feats = mp.adjust[,13:dim(mp.adjust)[2]] 
    
    # The list of features that have been checked for plate effects
    vars = names(feats)  
    vars
    
    
    #wells = mp.adjust$plate_well 
    #wells
  

    ########################
    # Description: 
    # For each feature, check each well in each plate (1848 wells across 6 plates) 
    # First feature is the area of the nucleus captured with the Hoechst stain in panel A
    # feature[1] = "ObjectArea_NUC_A"
    # First well is B02 plate 1 replicate 1 panel A
    # well[1] = "1A1_B02" 
    # try a subset of features
    ########################
    
    
    # Subset of features
    vars1 = vars[1:6]
    vars1
    
    j = vars[4]
    
    # loop through each feature
    for (j in vars) {   
  
      # loop through each well in each plate - 1848 wells total
     for (i in seq(1,length(mp.adjust$plate_well))) {         
    
      # ind will be the (number of) indices of cells that need adjusting
      ind = which(mp.adjust$plate_well[i] == data$plate_well) 
      
      # the adjustment amount
      adj_val = mp.adjust[i,j] 
      
      # each cell gets adjusted here
      the.change = adj_val + data[ind,j]  
      
      # dat_modified is the new dataframe with the updated *cell* values
      dat_modified[ind, j] = the.change 
    
     }
}

  
    
  # Save data output
  # setwd("~/Location/of/files")
  # write.csv(dat_modified, file = "name_of_your_file.csv")
  
  
    

  
  
 
