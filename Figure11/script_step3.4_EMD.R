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
# Methods section 3.4: Earth mover's distance
# Figure 11: Simulated cell features to illustrate the EMD calculation.  
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
library(reshape2)
library(plyr)

########################
# setwd("~/Directory.../.../")
########################



########################
# Case 1
########################
x1 = rnorm(n = 900, mean = 0, sd = .5)
y1 = rnorm(n = 900, mean = 0, sd = 1.5)
########################
# Case 2
########################
x2 = rnorm(n= 900, mean = 0, sd = .5)
y2 = rnorm(n= 900, mean = 1, sd = .5)

########################
# Pick one of the five cases above, or define a new case by changing the parameters below:
# Data distribution: number of data points (n), mean, standard deviation.
########################

x = x2
y = y2

########################
# Calculate the EMD - area between the two CDF curves
########################

xy = c(x,y)
ord = order(xy);
v = c( rep(1/length(x), length(x) ), rep(-1/length(y), length(y) ) )
emd1 = sum( abs( diff(xy[ord]) * ( cumsum(v[ord]) [-length(v)] ) ),  na.rm = T )
emd1


########################
# Melt the data
########################
df1 = data.frame(x,y)
df2 = melt(df1)

########################
# Calculate median of each sample
########################
med <- ddply(df2, "variable", summarise, grp.median = median(value))
head(med)

mu <- ddply(df2, "variable", summarise, grp.mean = mean(value))
head(mu)

########################
# Setup ggplot parameters
########################
# Yaxis numbers
Y.text <- element_text( color = "black", size = 10, family = "Helvetica")
# X axis numbers
X.text <- element_text( color = "black",  size = 10, family = "Helvetica") 
# X and Y axis labels
a.text = element_text(color = "black", size= 10, family = "Helvetica")
# Title text
t.text = element_text(color = "black", size= 8, family = "Helvetica")  

########################
# Fig 1
# (a) plot density curves, (b) plot CDF curves
########################

a = ggplot(df2, aes(x=value, color=variable)) +
    geom_density(size=1.2,) + 
    
    #geom_vline(data=mu, aes(xintercept=grp.mean, color=variable), linetype="dashed") +
    geom_vline(data=med, aes(xintercept=grp.median, color=variable), linetype="dashed") +
  
    scale_color_manual(values=c("black", "steelblue2")) +
    # # set up the y-axis labels
    scale_y_continuous(limits=c(0,.9), breaks=c(0, 0.3,
                                              0.6, 0.9), expand = c(0,0)) +
  scale_x_continuous(limits=c(-3,6), breaks=c(-3, -2,
                                              -1, 0, 1, 2, 3, 4, 5, 6), expand = c(0,0)) +
  
    theme_bw() + 
    theme(legend.position = 'none') +
    theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, 
        title = t.text, text = t.text) +
    labs(x = NULL, y = "Density curves", title = "N = 900, D1: (mean = 0, sd = .5), D2: (mean = 2, sd = 0.5)")
    # labs(x = NULL, y = "Density curves", title = "N = 300, D1: (mean = 0, sd = 1), D2: (mean = 0, sd = 2)")


b = ggplot(df2, aes(x=value, color=variable)) + stat_ecdf(size=1.2)+
    theme_bw() + 
    #scale_color_manual(values=c("black", "#E69F00")) +
    scale_color_manual(values=c("black", "steelblue2")) +
    theme(legend.position = 'none') +
    #scale_y_continuous(limits=c(0,1), breaks=c(0, 0.2, 0.4,
    #                                           0.6, 0.8, 1), expand = c(0,0)) +
    scale_x_continuous(limits=c(-3,5), breaks=c(-3, -2,
                                              -1, 0, 1, 2, 3, 4, 5), expand = c(0,0)) +
  
    theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, 
        title = t.text, text = t.text) +
    labs(title=NULL, y = "Cumulative distribution", x = "Feature") 

pdf(file = "Figure11_b.pdf",   # The directory you want to save the file in
    width = 4.5, # The width of the plot in inches
    height = 6) # The height of the plot in inches
########################
# Plot the two figures together
grid.arrange(a,b,nrow = 2)
########################

dev.off()








########################
# Pick one of the five cases above, or define a new case by changing any of the following:
# Data distribution, number of data points (n), mean, standard deviation.
########################
x = x1
y = y1

########################
# Calculate the EMD - area between the two CDF curves
########################

xy = c(x,y)
ord = order(xy);
v = c( rep(1/length(x), length(x) ), rep(-1/length(y), length(y) ) )
emd2 = sum( abs( diff(xy[ord]) * ( cumsum(v[ord]) [-length(v)] ) ),  na.rm = T )
emd2


########################
# Melt the data
########################
df1 = data.frame(x,y)
df2 = melt(df1)

########################
# Calculate median of each sample
########################
med1 <- ddply(df2, "variable", summarise, grp.median = median(value))
head(med1)

mu1 <- ddply(df2, "variable", summarise, grp.mean = mean(value))
head(mu1)


########################
# Fig 1
# (a) plot density curves, (b) plot CDF curves
########################

a1 = ggplot(df2, aes(x=value, color=variable)) +
  geom_density(size=1.2,) + 
  
  geom_vline(data=med1, aes(xintercept=grp.median, color=variable), linetype="dashed") +
  
  scale_color_manual(values=c("black", "steelblue2")) +
  # # set up the y-axis labels
  scale_y_continuous(limits=c(0,.9), breaks=c(0, 0.3,
                                              0.6, 0.9), expand = c(0,0)) +
  scale_x_continuous(limits=c(-3,6), breaks=c(-3, -2,
                                              -1, 0, 1, 2, 3, 4, 5, 6), expand = c(0,0)) +
  
  theme_bw() + 
  theme(legend.position = 'none') +
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, 
        title = t.text, text = t.text) +
  labs(x = NULL, y = "Density curves", title = "N = 900, D1: (mean = 0, sd = 0.5), D2: (mean = 0, sd = 1.5)")


b1 = ggplot(df2, aes(x=value, color=variable)) + stat_ecdf(size=1.2)+
  theme_bw() + 
  #scale_color_manual(values=c("black", "#E69F00")) +
  scale_color_manual(values=c("black", "steelblue2")) +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits=c(0,1), breaks=c(0, 0.2, 0.4,
  #                                           0.6, 0.8, 1), expand = c(0,0)) +
  scale_x_continuous(limits=c(-3,5), breaks=c(-3, -2,
                                              -1, 0, 1, 2, 3, 4, 5), expand = c(0,0)) +
  
  theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, 
        title = t.text, text = t.text) +
  labs(title=NULL, y = "Cumulative distribution", x = "Feature") 

# Fig 11 plot to pdf
pdf(file = "Figure11_a.pdf",   # The directory you want to save the file in
    width = 4.5, # The width of the plot in inches
    height = 6) # The height of the plot in inches

########################
# Plot the two figures together
grid.arrange(a1,b1,nrow = 2)
########################

dev.off()

emd1
emd2
