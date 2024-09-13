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

df <- read_csv("data_step3.3d_cell_data_allplates.csv")
df = data.frame(df)
class(df)

########################
# Identify binning parameters in preparation for plot
########################

bin.old = diff(range(df$FeatureX_old))

bin.sta = diff(range(df$FeatureX_sta))


unique(df$PlateNumber)
ind1 = which(df$PlateNumber == "1A1")
ind2 = which(df$PlateNumber == "1A2")
ind3 = which(df$PlateNumber == "1A3")

ind4 = which(df$PlateNumber == "2A1")
ind5 = which(df$PlateNumber == "2A2")
ind6 = which(df$PlateNumber == "2A3")

df$plate_rep = df$PlateNumber
df$plate_rep[ind1] = "Plate1 rep1"
df$plate_rep[ind2] = "Plate1 rep2"
df$plate_rep[ind3] = "Plate1 rep3"
df$plate_rep[ind4] = "Plate2 rep1"
df$plate_rep[ind5] = "Plate2 rep2"
df$plate_rep[ind6] = "Plate2 rep3"

unique(df$plate_rep)




Y.text <- element_text( size =10, color = "black", family = "Helvetica") 
X.text <- element_text(size = 10, color = "black", family = "Helvetica")
a.text = element_text(size = 12, color = "black", family = "Helvetica")
t.text.o = element_text(size = 8 ,color = "black", family = "Helvetica")
t.text.n = element_text(size = 12 ,color = "black", family = "Helvetica")
strip.text.a = element_text(size = 10, color = "black", family = "Helvetica")


########################
# Raw data one histogram showing all 6 plates worth of control wells
# Figure 8
########################

o1 = ggplot(df, aes(FeatureX_old, colour = "blue")) + 
  geom_histogram(binwidth = bin.old/130, alpha = 0.7, show.legend = FALSE, 
                 color="darkgrey", fill="darkgrey") +
  theme_bw() +
  labs(title = "Raw data distribution of cells from all plates", y = "Distribution", x = "Total nucleus intensity\n(DMSO control)") +
  # ylim(0, 4000) +
  theme(axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text.o) 


########################
# Raw data - per plate histograms
# Figure 9
########################

Y.text <- element_text( size = 10, color = "black", family = "Helvetica") 
X.text <- element_text(size = 10, color = "black", family = "Helvetica")


o2 = ggplot(df, aes(FeatureX_old, colour = "blue")) + 
  geom_histogram(binwidth = bin.old/130, alpha = 0.7, show.legend = FALSE, 
                 color="darkgrey", fill="darkgrey") +
  labs(title = "Raw data distribution of cells per plate", y = "Distribution", x = "Total nucleus intensity\n(DMSO control)") +
  theme_bw() +
  ylim(0, 4000) +
  theme(strip.background = element_rect(fill="white"), strip.text = strip.text.a, 
        axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text.o) +
  facet_wrap(~ plate_rep,  nrow = 2) 


########################
# Plate -- 1A1 -- 
########################

names(df)
unique(df$plate_well)
unique(df$Plate)
unique(df$Plate_number)
unique(df$PlateNumber)

ind = which(df$PlateNumber == "1A1")

df.sub = df[ind,]

########################
# Prepare plot parameters 
########################

Y.text <- element_text( size =10, color = "black", family = "Helvetica") 
X.text <- element_text(size = 10, color = "black", family = "Helvetica")
a.text = element_text(size = 12, color = "black", family = "Helvetica")
t.text.o = element_text(size = 10 ,color = "black", family = "Helvetica")
t.text.n = element_text(size = 12 ,color = "black", family = "Helvetica")
strip.text.a = element_text(size = 12, color = "black", family = "Helvetica")


########################
# Raw data one histogram showing all 6 plates worth of control wells
########################

o1_sub =  ggplot(df.sub, aes(FeatureX_old, colour = "blue")) + 
  geom_histogram(binwidth = bin.old/130, alpha = 0.7, show.legend = FALSE, 
                 color="darkgrey", fill="darkgrey") +
  ylim(0,3500) +  
  theme_bw() +
  labs(title = "Raw data distribution: Plate1 rep1", y = "Distribution", x = "Total nucleus intensity\n(DMSO control cells)") + 
  theme(axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text.o) 



########################
# Raw data - per well and per plate density
########################

o4_sub = ggplot(df.sub, aes(x=FeatureX_old, color=plate_well)) +
  geom_density(show.legend = FALSE) + theme(legend.position='none') + 
  scale_color_grey() + 
  theme_bw() +       
  labs(title = "Per well density of raw data: Plate1 rep1", y = "Density", x = "Total nucleus intensity\n(DMSO control)") +
  theme(strip.text = strip.text.a, axis.text.x = X.text, 
        axis.text.y = Y.text, axis.title = a.text, title = t.text.o) 
# +      
# facet_wrap(~ PlateNumber,  nrow = 2) + 



##################
# Standardized data one histogram showing all 6 plates worth of control wells
##################

s1_sub = ggplot(df.sub, aes(FeatureX_sta, colour = "blue")) + 
  geom_histogram(binwidth = bin.sta/130, alpha = 0.7, show.legend = FALSE, 
                 color="darkgrey", fill="darkgrey") +
  ylim(0,3500) +     
  theme_bw() +
  labs(title = "Standardized data: Plate1 rep1", y = "Distribution", x = "Total nucleus intensity\n(DMSO control)") + 
  theme(axis.text.x = X.text, axis.text.y = Y.text, axis.title = a.text, title = t.text.o) 



##################
# Standardized data - per well and per plate density
##################

s4_sub = ggplot(df.sub, aes(x = FeatureX_sta, color = plate_well)) +
    geom_density(show.legend = FALSE) + theme(legend.position='none') + 
    scale_color_grey() + 
    theme_bw() +
    labs(title = "Per well density of standardized data: Plate1 rep1", y = "Density", x = "Total nucleus intensity\n(DMSO control)") +
    theme(strip.text = strip.text.a, axis.text.x = X.text, 
        axis.text.y = Y.text, axis.title = a.text, title = t.text.o) 



# Fig 8 plot to pdf
pdf(file = "Figure8.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4) # The height of the plot in inches

o1

dev.off()

# Fig 9 plot to pdf
pdf(file = "Figure9.pdf",   # The directory you want to save the file in
    width = 9, # The width of the plot in inches
    height = 4.5) # The height of the plot in inches

o2

dev.off()


pdf(file = "Figure10a_left.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4) # The height of the plot in inches

o1_sub

dev.off()

pdf(file = "Figure10b_left.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4) # The height of the plot in inches

s1_sub

dev.off()

pdf(file = "Figure10a_right.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4) # The height of the plot in inches

o4_sub

dev.off()


pdf(file = "Figure10b_right.pdf",   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4) # The height of the plot in inches

s4_sub

dev.off()


o1
o2


o1_sub
s1_sub

o4_sub
s4_sub



          
