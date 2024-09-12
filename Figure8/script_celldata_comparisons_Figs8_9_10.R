########################
# R script is in support Methods in Molecular Biology book chapter and 
# data analysis methods published in manuscript (https://doi.org/10.1038/s42003-022-04343-3) 
# Script is written by Yanthe E. Pearson
########################

########################
# Quality control measures and statistical strategies to address the challenges of high-content phenotypic data
# by Yanthe E. Pearson and Kristin C. Gunsalus
########################

########################
# Import cell level data to plot pre and post normalized feature distributions
# Figures 8, 9 and 10
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

setwd("~/Documents/NYU_work/Sept_17_2018/R_code_b/MiMB_BookChapter_2024")

########################
# import **data_step4_adjustments.csv**
########################

df <- read_csv("cell_data_allplates.csv")
df = data.frame(df)

df.sub <- read_csv("cell_data_plate1rep1.csv")
df.sub = data.frame(df.sub)



          