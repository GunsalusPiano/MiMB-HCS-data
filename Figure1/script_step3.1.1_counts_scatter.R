########################
# R script is in support of the published MiMB book chapter and data analysis methods published in manuscript (https://doi.org/10.1038/s42003-022-04343-3) 
########################

########################
# Script is written by Yanthe E. Pearson
########################

########################
# R script for MiMB BookChapter - Step 3.1.1
########################

########################
# R script for figure  1: Raw cell counts   
########################

# clear variables
rm(list=ls())

########################
# libraries
########################

library(ggplot2)
library(RColorBrewer)
library(readr)

########################
# set working directory
########################

setwd("~/Documents/NYU_work/Sept_17_2018/R_code_b/MiMB_BookChapter_2024")

########################
# (1) Import the annotated well median data file from panel A channels (Hoechst, syto14, mito, WGA)
########################

well_med_dat <- read_csv("data_step3.1.1_counts_and_medians_panelA.csv")
  
########################
# identify control wells
ind_control = grep('control', well_med_dat$Treatment)

########################
# dataframe of control wells
annot.count.c = well_med_dat[ind_control,]
  
########################
# summary of control well counts, including range, median
########################

control_range = range(annot.count.c$cell_count)
# The range is 419 to 990 cells
control_range

control_med = median(annot.count.c$cell_count)
# The median is 812 cells
control_med 
########################

########################
# Check column names and number of plates
names(well_med_dat)
unique(well_med_dat$PlateNumber)

########################
# Remove empty treatments
annot.count = well_med_dat[-grep("empty",well_med_dat$Treatment),]
  
########################
# grab the control wells to create a second layer of control (purple) points
ind.c = which(annot.count$Drug_name == "DMSO")
df.c = annot.count[ind.c,]
######################## 

########################
# grab the some treatment wells (Nocodazole)
ind.t = which(annot.count$Drug_name == "nocodazole")
df.t = annot.count[ind.t,]
########################

########################
# Scatterplot paramters
########################
     
Y.text <- element_text( color = "black", size = 12, family = "Helvetica")
X.text <- element_text( color = "black", size = 12, family = "Helvetica")
a.text = element_text( color = "black", size = 14, family = "Helvetica")
strip.text.a = element_text(color = "black", size = 14, family = "Helvetica")
t.text = element_text(color = "black", size= 12, family = "Helvetica") 

########################
# Scatter plot updates - grey for treatments
# purple for control
########################

########################
# Scatter plot implementation using ggplot
########################


# Step 1: Call the pdf command to start the plot
pdf(file = "figure_step3.1.1_scatter.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches


  ggplot(data = annot.count, aes(Position_index, cell_count, color = Treatment), shape = 20, size = 2, alpha = .8) +
  geom_point( ) + 
  # add layer of controls
  geom_point(data = df.c, aes(Position_index, cell_count), shape = 20, size = 2, alpha = .8 ) +
  # add layer for plotting nocodazole
  geom_point(data = df.t, aes(Position_index, cell_count), shape = 20, size = 2, alpha = .8, color = "red" ) +
  
  ###
  scale_color_manual(values=c("purple", "lightgrey")) +
  ###
  geom_hline(yintercept= control_med, linetype ="dashed", color = "darkgrey", size = .5) +
  geom_hline(yintercept = control_range[1], linetype="dashed", color = "darkgrey") +
  geom_hline(yintercept = control_range[2], linetype="dashed", color = "darkgrey") +
  ###
  ylim(-10, 1110) +
  
  theme_classic() +
  theme(strip.text = strip.text.a, axis.text.x = X.text, 
        axis.text.y = Y.text, axis.title = a.text, title = t.text, text = t.text) + 
  labs(x = "Plate # rep #", y = "Cell count", title = NULL) +
  scale_x_continuous(breaks = c(1, 309, 617, 925, 1233, 1541), 
                     labels=c("p1.r1", "p1.r2", "p1.r3", "p2.r1", "p2.r2", "p2.r3")) +
  theme(legend.position = "none") 
  
  dev.off()

########################
# Legend
########################
  
  # Step 1: Call the pdf command to start the plot
  pdf(file = "figure_step3.1.1_scatter_legend.pdf",   # The directory you want to save the file in
      width = 5, # The width of the plot in inches
      height = 4) # The height of the plot in inches
  
  
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

legend("topleft", legend =c("DMSO control", "Treatments", "Nocodazole"),
                             
       pch = 19, pt.cex = 3, cex = 1.5, bty='n',
       
       col = c("purple","lightgrey", "red"))

mtext("Treatment", at = 0.2, cex = 2)

dev.off()
     
 
     
