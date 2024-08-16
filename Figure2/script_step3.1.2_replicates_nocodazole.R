########################
# R script is in support of the published MiMB book chapter and data analysis methods published in manuscript (https://doi.org/10.1038/s42003-022-04343-3) 
########################

########################
# Script is written by Yanthe E. Pearson
########################

########################
# R script for MiMB BookChapter - Step 3.1.2: Cell feature distributions
########################

########################
# R script for figure  2: Visualizing raw data replicates from different plates as feature distributions   
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
library(dplyr)

########################
# setwd("~/Directory.../.../")
########################

setwd("~/Location/of/Data")

########################
# Import raw data - Nocodazole data
########################

cell.dat = read_csv("data_step3.1.2_density_replicates_nocodazole.csv")

########################
# Import raw data - DMSO cell data
########################

cell.dat.c = read_csv("data_step3.1.2_cell_data_control.csv")

########################
# Visualize the raw feature distributions - observe dose response and replicates
########################

########################
# Check annotations
########################

unique(cell.dat$plate_well)
unique(cell.dat$plate_well)
names(cell.dat)
unique(cell.dat$Nocodazole_uM)
unique(cell.dat$Metadata_dil_factor)
unique(cell.dat$PlateNumber)


######################
# convert the numeric column to a character column 
cell.dat$Nocodazole_uM = as.character(cell.dat$Nocodazole_uM)

unique(cell.dat$Nocodazole_uM)

# Check cell counts in each nocodazole treatment (with replicates merged) 
count(cell.dat, WellId)

# Check cell counts in each nocodazole treatment (without replicates merged) 
count(cell.dat, plate_well)
######################

# repeat for the control data
cell.dat.c$Nocodazole_uM = cell.dat.c$Drug_name
cell.dat.c$Nocodazole_uM = as.character(cell.dat.c$Nocodazole_uM)
# DMSO
unique(cell.dat.c$Nocodazole_uM)
######################
# Add annotations for plot labels
cell.dat$PlateNumber1 = cell.dat$PlateNumber
unique(cell.dat$PlateNumber)

ind1 = which(cell.dat$PlateNumber == "1A1")
ind2 = which(cell.dat$PlateNumber == "1A2")
ind3 = which(cell.dat$PlateNumber == "1A3")

cell.dat$PlateNumber1[ind1] = "Plate1 Rep1"
cell.dat$PlateNumber1[ind2] = "Plate1 Rep2"
cell.dat$PlateNumber1[ind3] = "Plate1 Rep3"
######################

######################
title4 = paste("Nocodazole", sep = "")
title4

title5 = paste("(\u03BC","M)", sep = "")
title5

title6 = paste(title4, title5, sep = " ")
title6
######################


########################
# Setup ggplot parameters - Y axis numbers, X axis numbers, X and Y axis labels, Title text.
########################

Y.text <- element_text( color = "black", size = 10, family = "Helvetica")
X.text <- element_text( color = "black",  size = 10, family = "Helvetica") 
a.text = element_text(color = "black", size= 12, family = "Helvetica")
t.text = element_text(color = "black", size= 10, family = "Helvetica")  

########################
# Object size - Nucleus
########################


# Step 1: Call the pdf command to start the plot
# pdf(file = "figure_step3.1.2_density.pdf",   # The directory you want to save the file in
#     width = 10, # The width of the plot in inches
#     height = 3) # The height of the plot in inches



  ggplot(cell.dat, aes(ObjectSizeCh1, color = Nocodazole_uM)) +
  geom_density(size = 0.8) + 
  # Add control distribution
  geom_density(data = cell.dat.c, aes(ObjectSizeCh1), color = "darkgrey", size = 0.9, linetype = 2) + 
  # color of curves
  scale_color_brewer(palette = "Reds", name = title6) +
  # set themes
  theme_classic() +
  theme(legend.position ='none') + 
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text,
        title = t.text, text = t.text ) +
  theme(strip.text.x = X.text) +
  labs(x = "Nucleus Size", title = NULL, y = 'Density') + 
  # separate data by plate number/replicate
  facet_wrap(~PlateNumber1, nrow = 1) 
  
  # dev.off()
  
  
  
  
  #######################
  #######################
  # Nocodazole - Legend
  ######################
  #######################
  
  mypalette <- brewer.pal(7,"Reds")
  
  col.tox = mypalette

  labels1 = c("0.3125", "0.625", "1.25", "2.5", "5","10", "20")
  
  # Step 1: Call the pdf command to start the plot
  pdf(file = "figure_step3.1.2_density_legend.pdf",   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  
  #######################
  #######################
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

  legend("topleft", legend = labels1, pch=15, pt.cex=3, cex=1.5, bty='n',
         col = col.tox)
  
  mtext("Nocodazole (uM)", at = 0.08, cex = 2)
  #######################
  #######################
  
  dev.off()
  
  
