########################
# R script is in support of the published manuscript (https://doi.org/10.1038/s42003-022-04343-3) 
# The book chapter for MiMb - May 2024 
########################

########################
# Script is written by Yanthe E. Pearson
########################
# May 2nd 2024

########################
# R script for Fig. 4: Statistical metric comparison and feature reproducibility.
# (b) Hypothetical probability density (PDF) and cumulative density (CDF) curves for 
# two random samples of the same feature are illustrated to show how the Kolmogorov–Smirnov (KS) distance and Wasserstein metric (EMD) are estimated.
########################

# Generate different probability distributions

# https://www.stat.umn.edu/geyer/old/5101/rlook.html

# runif, rexp, rnorm, rbinom, rgeom, 



########################
# Description: 
#
#
#
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
# Case 1
########################
x1 = rnorm(n = 300, mean = 0, sd = 1)
y1 = rnorm(n = 300, mean = 1.5, sd = 1.5)
########################
# Case 2
########################
x2 = rnorm(n= 300, mean = 0, sd = 1)
y2 = rnorm(n= 300, mean = 0, sd = 2)

########################
# Case 3
########################
x3 = rnorm(n = 300, mean = 0, sd = 1)
y3 = rnorm(n = 300, mean = 2, sd = 1)

########################
# Case 4
########################
x4 = rnorm(n = 300, mean = 0, sd = 1)
y4 = rnorm(n = 300, mean = 2.5, sd = 1.75)

########################
# Case 5
########################
 x5 = rnorm(n = 300, mean = 0, sd = 1)
 y5 = rnorm(n = 300, mean = 3, sd = 2)

########################
# Case 6
########################
# rlnorm(n, meanlog = 0, sdlog = 1)
# rgamma(n, shape, rate = 1, scale = 1/rate)
# rbeta(n, shape1, shape2, ncp = 0)
# rweibull(n, shape, scale = 1)


########################
# Pick one of the five cases above, or define a new case by changing any of the following:
# Data distribution, number of data points (n), mean, standard deviation.
########################
x = x3
y = y3

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
Y.text <- element_text( color = "black", size = 10, angle = 0, family = "Helvetica")
# X axis numbers
X.text <- element_text( color = "black",  size = 10, angle = 0, family = "Helvetica") 
# X and Y axis labels
a.text = element_text(size= 10, color = "black", angle = 0, family = "Helvetica")
# Title text
t.text = element_text(size= 8, color = "black", family = "Helvetica")  




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
    # scale_y_continuous(limits=c(0,.6), breaks=c(0, 0.1, 0.2,
    #                                                0.3, 0.4, 0.5, 0.6), expand = c(0,0)) +
    theme_bw() + 
    theme(legend.position = 'none') +
    theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, 
        title = t.text, text = t.text) +
    labs(x = NULL, y = "Density curves", title = "N = 300, D1: (mean = 0, sd = 1), D2: (mean = 0, sd = 0.5)")
     #labs(x = NULL, y = "Density curves", title = "N = 300, D1: (mean = 0, sd = 1), D2: (mean = 0, sd = 2)")


b = ggplot(df2, aes(x=value, color=variable)) + stat_ecdf(size=1.2)+
    theme_bw() + 
    #scale_color_manual(values=c("black", "#E69F00")) +
    scale_color_manual(values=c("black", "steelblue2")) +
    theme(legend.position = 'none') +
    #scale_y_continuous(limits=c(0,1), breaks=c(0, 0.2, 0.4,
    #                                           0.6, 0.8, 1), expand = c(0,0)) +
    theme(axis.text.y = Y.text, axis.text.x = X.text, axis.title = a.text, 
        title = t.text, text = t.text) +
    labs(title=NULL, y = "Cumulative distribution", x = "Normal distribution") 


########################
# Plot the two figures together
grid.arrange(a,b,nrow = 2)
########################

######################
# Fig 2
######################

sample1 = x
sample2 = y


######################
# create ECDF of data
cdf1 <- ecdf(sample1) 
cdf2 <- ecdf(sample2) 

# Using line segments to shade in the area between two ECDFs

ecdf0 <- cdf1 
ecdf1 <- cdf2

plot(ecdf0, xlim=range(sample1, sample2), verticals=TRUE, do.points=FALSE, col="black", lwd = 3, main = NULL,
     col.main="blue", col.lab="black", col.sub="black") 
plot(ecdf1, verticals=TRUE, do.points=FALSE, col="steelblue2", lwd = 3, add=TRUE) 

# ----- this will plot lots of segments (2001 of them?)
tmpx <- seq(min(sample1, sample2), max(sample1, sample2), len=2001) 
p0 <- ecdf0(tmpx) # ----- y upper and y lower
p1 <- ecdf1(tmpx) # ----- y upper and y lower
segments(x0=tmpx, y0=p0, y1=p1, col=adjustcolor(ifelse(p0<p1, "steelblue", "grey"), alpha.f=0.1)) 

# Show KS distance with vertical red line
# points(c(x0, x0), c(y0, y1), pch=16, col="red") 
# segments(x0, y0, x0, y1, col="red", lty="dotted") 


