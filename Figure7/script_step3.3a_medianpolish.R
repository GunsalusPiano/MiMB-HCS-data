########################
# R script is in support of the published MiMB book chapter (2024) and data analysis methods published in manuscript (https://doi.org/10.1038/s42003-022-04343-3) 
########################

########################
# Script is written by Yanthe E. Pearson
########################

########################
# R script for MiMB BookChapter - Methods section 3.3: Plate normalization and cell standardization
########################


########################
# Figure 7: Scatter plot of well-level data after plate correction and data standardization for six plates.
########################


########################
# Description: 
# import well medians
# for each feature, apply the two-way ANOVA
# if significant row or column effects are detected (p value < 0.0006)
# apply median polish to whole plate, and calculate B scores
# if no plate effects detected, simply calculate the the per well Z scores.
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
# setwd("~/Location/of/your/files")
########################

setwd("~/Documents/NYU_work/Sept_17_2018/R_code_b/MiMB_BookChapter_2024/fig3C")

########################
# Import the updated medians file 
########################

df.features <- read_csv("data_step3.3a_well_medians.csv")

df.features.num = df.features[,12:dim(df.features)[2]]

plate.unique = sort(unique(df.features$PlateNumber))

# initialize
mp = 0

# initialize
z = 0

########################
# FULL FEATURE SET 
# (1) For each feature 
# (2) Apply the two way anova model
# (3) Apply median polish + trend/adjust/residuals and B scores
########################

########################
# Create lists for saving raw well medians
datalist.plates = list()
datalist.plates.all = list()

########################
# Create lists for saving post normalized medians
plates.polish = list()
plates.polish.all = list()
plates.mpmedian.save.all = list()

########################
# Stores observed trend
plates.trend = list()
plates.trend.all = list()

########################
# Stores adjustment amount
plates.adjusted = list()
plates.adjusted.all = list()
plates.adjusted.save.all = list()

########################
# Stores b scores
plates.Bc = list()
plates.Bc.all = list()
bscores.save = list()
bscores.save.all = list()

########################
# j loops through each feature, i loops through each plate
########################

for (j in seq(1,length(df.features.num))) { 
  
  feat.current = names(df.features.num)[j]
  feat.current
  
  # i loops through each plate
  for (i in seq(1,length(plate.unique))) {   
    

    plate.ind = which(df.features$PlateNumber == plate.unique[i])
    plate.current = df.features.num[plate.ind,j]
    plate.current.labels = df.features[plate.ind,]
    plate.current.ord <- plate.current.labels[order(plate.current.labels$WellId),]
    
    # The data
    plate.use <- plate.current.ord[,11+j]
    # The feature name
    feat.use <- names(plate.current.ord)[11+j]
    
    ########################
    # Raw data matrix 
    feature1_mat <- data.matrix(plate.use)
    feature1_mat1 <- matrix(feature1_mat, nrow = 14, ncol = 22, byrow = TRUE)
    rownames(feature1_mat1) <- LETTERS[2:15]
    colnames(feature1_mat1) <- 2:23
    
    ########################
    # Saving raw data matrix 
    datalist.plates[[i]] <- feature1_mat1 
    
    ########################
    # Apply the two - factor ANOVA to determine if there is positional dependency within each plate
    # This is applied to *control* wells only
    ########################
    #
    # The feature ~ factors (row and column position) - The data (cell feature) is treated as the dependent variable
    ########################
    
    # "Plate"  "Treatment" "PlateNumber"  "Row"  "Column"  "feature"  "feature"
    plate.use1 <- plate.current.ord[,c(3,6,9,10,11,11+j,11+j)]
    
    # re-name replicate column
    names(plate.use1)[6] = "feature_x"
    
    # Prepare control well data 
    plate.c = plate.current.ord[which('control' == plate.current.ord$Treatment),c(3,6,9, 10,11,11+j,11+j)]
    
    # re-name replicate column
    names(plate.c)[6] = "feature_x"
    
    ########################
    # A two-way ANOVA assesses whether the factors (categories) - categories are row and columnn number (plate position)
    # have a significant effect on the outcome (measurement/feature). 
    # Its implementation can be conceptualized as a regression analysis where each row and column 
    # level is treated as a dummy variable. 
    ########################

    plate.c$Row = as.factor(plate.c$Row)
    plate.c$Column = as.factor(plate.c$Column)
    
    ########################
    # Implementation of the two factor anova
    ########################
    
    M <- lm(feature_x ~ Row + Column, data=plate.c)
    # summary(M)
    
    MA = Anova(M)
    
    row_F = MA$`F value`[1]
    col_F = MA$`F value`[2]
    
    row_p = MA$`Pr(>F)`[1]
    col_p = MA$`Pr(>F)`[2]
    
    print(names(plate.use))
    print(row_p)
    print(col_p)
    
    
    
    
    ########################
    # check for significance row or columns effects via p-value
    ########################
    
    if (row_p < 0.00006 || col_p < 0.00006)  {
    
      
     mp = mp+1
      
 
      ########################
      # Performs median polish on whole plate
      # (1) overall/grand median --	the fitted constant term.
      # (2) row -- the fitted row effects.
      # (3) col -- the fitted column effects.
      # (4) residuals	the residuals.
      ########################
     
      print("Median Polishing")
      
      plate.mp = medpolish(feature1_mat1, na.rm = T) 

      # plot(plate.mp)
      # TREND = RAW DATA - REDISUALS

      trend = feature1_mat1 - plate.mp$residuals
      normed = plate.mp$overall + plate.mp$residuals   # MEDIAN POLISH VALUES (NEW MEDIANS)
      raw = feature1_mat1
      adjustment = normed - raw
      res = plate.mp$residuals

      ########################
      # Next, use residual values to calculate the B score using control wells
      # (1) create dataframe of residuals and identify treatment/control wells
      # (2) use the control wells, find MED (control) and MAD (control wells)
      ########################
      
      # transform matrix to a vector
      df.Bc = data.frame(as.vector(t(res)))      
      df.Bc.ann = cbind(plate.current.labels[,1:11],df.Bc)
      names(df.Bc.ann)[dim(df.Bc.ann)[2]] = "mp.residuals"

      # identify controls
      ind_control <- which(grepl("control", df.Bc.ann$Treatment)) 
      control_wells <- df.Bc.ann[ind_control,]
      MED.RES.C = median(control_wells$mp.residuals)
      MAD.RES.C = mad(control_wells$mp.residuals)

      
      if (MAD.RES.C > 0) {
        df.Bc.ann$BSCORE.C = (df.Bc.ann$mp.residuals - MED.RES.C)/MAD.RES.C
        BC.save = (df.Bc.ann$mp.residuals - MED.RES.C)/MAD.RES.C
      } else {
        MED.RES.C = mean(control_wells$mp.residuals)  # Can't use median if too many zeros and denominator is 0 --- use MEAN instead --- RNA only
        MAD.RES.C = sd(control_wells$mp.residuals)
        df.Bc.ann$BSCORE.C = (df.Bc.ann$mp.residuals - MED.RES.C)/MAD.RES.C
        BC.save = (df.Bc.ann$mp.residuals - MED.RES.C)/MAD.RES.C
      }
      

      ########################
      # SAVING THE B SCORES -- VECTOR OF DATA HERE
      ########################
      
      bscores.save[[i]] = BC.save  

      ########################
      # CONVERT BACK TO MATRIX
      ########################
      
      Bcontrol <- data.matrix(df.Bc.ann$BSCORE.C)
      Bcontrol <- matrix(Bcontrol, nrow = 14, ncol = 22, byrow = TRUE)
      rownames(Bcontrol) <- LETTERS[2:15]
      colnames(Bcontrol) <- 2:23 
      
      ########################
      # For each feature, store the "trend", "normed", "adjustment" and "Bcontrol" matrices
      # All six plates 
      ########################
      
      plates.trend[[i]] <- trend 
      plates.polish[[i]] <- normed
      plates.adjusted[[i]] <- adjustment
      plates.Bc[[i]] <- Bcontrol
      
   
      
  } 
  
    ########################
    # If no significant plate effects are detected 
    # if p >= 0.00006
    # Do not apply median polish to the well medians. 
    # The plate is simply normalized to the DMSO-control wells (Z scores, see paper), 
    # this will correct for plate to plate variation
    ########################
    
  
  
  else {
    
    z = z+1
    
    print("Z scoring")
    
    ########################
    # Set entries to NA so as to preserve the dimensionality of the data
    ########################
    
    trend = matrix(NA, nrow = 14, ncol = 22)
    normed = matrix(NA, nrow = 14, ncol = 22)   # Since we aren't doing median-polish --- set to NA
    raw = feature1_mat1
    adjustment = matrix(0, nrow = 14, ncol = 22)
    res = matrix(NA, nrow = 14, ncol = 22)
    
    ########################
    # Set entries to NA so as to preserve the dimensionality of the data
    ########################
    
    
    
    
    
    ########################
    # Next, Calculate the robust Z scores using control wells
    # (1) create data frame of medians and identify treatment/control wells
    # (2) use the control wells, find MED (control) and MAD (control wells)
    ########################
    
    ########################
    # the raw data
    ########################
    df.Bc.ann = plate.use1 
    
    ########################
    # Identify controls
    ########################
    ind_control <- which(grepl("control", df.Bc.ann$Treatment)) 
    control_wells <- df.Bc.ann[ind_control,]
    MED.RES.C = median(control_wells$feature_x)
    MAD.RES.C = mad(control_wells$feature_x)
    
    ########################
    # Check if the MAD of a plate/feature is zero!
    ########################
    
    if (MAD.RES.C > 0) {
      df.Bc.ann$BSCORE.C = (df.Bc.ann$feature_x - MED.RES.C)/MAD.RES.C
      BC.save = (df.Bc.ann$feature_x - MED.RES.C)/MAD.RES.C
    } else {
      MED.RES.C = mean(control_wells$feature_x)  # Can't use median if too many zeros and denominator is 0 --- use MEAN instead --- RNA only
      MAD.RES.C = sd(control_wells$feature_x)
      df.Bc.ann$BSCORE.C = (df.Bc.ann$feature_x - MED.RES.C)/MAD.RES.C
      BC.save = (df.Bc.ann$feature_x - MED.RES.C)/MAD.RES.C
    }
    
    
    ########################
    # Saving the B (or Z) scores 
    ########################
    
    bscores.save[[i]] = BC.save 
    
    ########################
    # Convert back to matrix
    ########################
    
    Bcontrol <- data.matrix(df.Bc.ann$BSCORE.C)
    Bcontrol <- matrix(Bcontrol, nrow = 14, ncol = 22, byrow = TRUE)
    rownames(Bcontrol) <- LETTERS[2:15]
    colnames(Bcontrol) <- 2:23 
    
    ########################
    # For each feature, store all six plates
    ########################
    
    plates.trend[[i]] <- trend 
    plates.polish[[i]] <- normed
    plates.adjusted[[i]] <- adjustment
    
    plates.Bc[[i]] <- Bcontrol
    
    
  } # ----- closing *if* statement
  
     ############################################# 
      
  } # ----- loop for plates
  

  
  
  ###################
  # Store plate matrices here
  # Per well medians
  A = rbind(datalist.plates[[1]],datalist.plates[[2]], datalist.plates[[3]])
  B = rbind(datalist.plates[[4]],datalist.plates[[5]], datalist.plates[[6]])
  C = cbind(A,B)
  datalist.plates.all[[j]] <- C # store the six plates
  
  ###################
  # Store plate matrices here
  # Per well medians after median polish normalization
  A.mp = rbind(plates.polish[[1]], plates.polish[[2]], plates.polish[[3]])
  B.mp = rbind(plates.polish[[4]], plates.polish[[5]], plates.polish[[6]])    #, plates.polish[[8]])
  C.mp = cbind(A.mp,B.mp)
  plates.polish.all[[j]] <- C.mp # store the six plates
  
  ###################
  # Store plate matrices here
  # Observed trens
  A.t = rbind(plates.trend[[1]], plates.trend[[2]], plates.trend[[3]])
  B.t = rbind(plates.trend[[4]], plates.trend[[5]], plates.trend[[6]])
  C.t = cbind(A.t,B.t)
  plates.trend.all[[j]] <- C.t # store the six plates

  ###################
  # Store plate matrices here
  # Adjustment amount
  A.a = rbind(plates.adjusted[[1]], plates.adjusted[[2]], plates.adjusted[[3]])
  B.a = rbind(plates.adjusted[[4]], plates.adjusted[[5]], plates.adjusted[[6]])
  C.a = cbind(A.a,B.a)
  plates.adjusted.all[[j]] <- C.a # store the six plates
  
  ###################
  # Store plate matrices here
  # B scores
  A.bc = rbind(plates.Bc[[1]], plates.Bc[[2]], plates.Bc[[3]])
  B.bc = rbind(plates.Bc[[4]], plates.Bc[[5]], plates.Bc[[6]])
  C.bc = cbind(A.bc,B.bc)
  plates.Bc.all[[j]] <- C.bc # store the six plates
  

  

  ########################
  # Save B scores to a dataframe 
  ########################
  
  one = as.vector(t(plates.Bc[[1]]))
  two = as.vector(t(plates.Bc[[2]]))
  three = as.vector(t(plates.Bc[[3]]))
  four = as.vector(t(plates.Bc[[4]]))
  five = as.vector(t(plates.Bc[[5]]))
  six = as.vector(t(plates.Bc[[6]]))
  

  hi2 = c(one, two, three, four, five, six)  
  bscores.save.all[[j]] = hi2   

  
  

  
  ########################
  # Save well adjustment data to a dataframe 
  ########################
  
  one.a = as.vector(t(plates.adjusted[[1]]))
  two.a = as.vector(t(plates.adjusted[[2]]))
  three.a = as.vector(t(plates.adjusted[[3]]))
  four.a = as.vector(t(plates.adjusted[[4]]))
  five.a = as.vector(t(plates.adjusted[[5]]))
  six.a = as.vector(t(plates.adjusted[[6]]))
  
  
  hi2.a = c(one.a,two.a,three.a,four.a,five.a,six.a) #,seven.a,eight.a)
  plates.adjusted.save.all[[j]] = hi2.a
  
  
  
  
  ##########################
  # loop for features closes here
  ##########################
  
} 




##########################
# prep data frame of B scores
##########################

hi = do.call(cbind,bscores.save.all)
colnames(hi) = names(df.features.num)
hi1 = cbind(df.features[,1:11],hi)

setwd("~/Documents/NYU_work/Sept_17_2018/R_code_b/MiMB_BookChapter_2024/fig3C")
write.csv(hi1, file = "data_output_Bscores_panelA_March14_2024.csv")  

##########################
# prep data frame of plate adjustments
##########################

hi2 = do.call(cbind,plates.adjusted.save.all)
colnames(hi2) = names(df.features.num)
hi3 = cbind(df.features[,1:11],hi2)

setwd("~/Documents/NYU_work/Sept_17_2018/R_code_b/MiMB_BookChapter_2024/fig3C")
write.csv(hi3, file = "data_output_adjustments_panelA_March14_2024.csv") 



##########################
# "hi1" is the B-scores
##########################

names1 = names(hi1)[12:dim(hi1)[2]]

##########################
# Prepare plot parameters
##########################

Y.text <- element_text( color = "black", size =10, family = "Helvetica")
X.text <- element_text( color = "black", size = 10, family = "Helvetica")
a.text = element_text(color = "black", size =10, family = "Helvetica")

##########################
# Index "j" calls feature # 4
# "ObjectTotalInten_NUC_A" or Total nucleus intensity
##########################

j = 4

# for (j in  seq(1,length(names1))) {
  
  feat = names1[j]
  feat.df.all = hi1[,c(1:11,j+11,j+11)]
  colnames(feat.df.all)[13] <- "feature_x"
  
  ind.c = which(grepl('control',feat.df.all$Treatment))
  ctrl = feat.df.all[ind.c,]
  
  r.ctrl = range(ctrl$feature_x)
  
  ind.t = which(grepl('treated',feat.df.all$Treatment))
  trt = feat.df.all[ind.t,]
  trt$X1 = 1:1491
  r.trt = range(trt$feature_x)
  feat.df.all$X1 = 1:1848
 
  
  ##########################
  # Plot off B scores for feature j = 4, Total intensity of the nucleaus (Hoechst channel)
  ##########################
    
  
   # pdf(file = "figure_step3.3a_Bscores_totalinten.pdf",   # The directory you want to save the file in
   #     width = 7, # The width of the plot in inches
   #     height = 5) # The height of the plot in inches
  
  
    ggplot(feat.df.all, aes(X1, feature_x)) + geom_point(aes(colour = feature_x)) + 
    theme_bw()+
    geom_hline(yintercept=r.ctrl[1], linetype="dashed", color = "black") +
    geom_hline(yintercept=r.ctrl[2], linetype="dashed", color = "black") +
    scale_colour_gradient2(low = "blue", high = "red") + 
    # labs(title = feat, y = "B scores" , x = "plate and replicate number") +
    labs(title = "Feature: Total nucleus intensity", y = "B scores" , x = "plate # replicate #") +
    scale_x_continuous(breaks = c(1, 309, 617, 925, 1233, 1541), 
                         labels=c("p1.r1", "p1.r2", "p1.r3", "p2.r1", "p2.r2", "p2.r3")) +
    ylim(-36.97533, 32.45456) +
    theme(plot.title = element_text(color="black", size=13),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    axis.text.x = X.text,
    axis.text.y = X.text,
    legend.title = element_blank()
    
   )
      
     
#  dev.off() 

#}

# dev.off()




  ##########################
  # Index "j" calls feature # 4
  # "ObjectSize_NUC_A"
  ##########################
  
  j = 7
  
  # for (j in  seq(1,length(names1))) {
  
  feat = names1[j]
  feat
  feat.df.all = hi1[,c(1:11,j+11,j+11)]
  colnames(feat.df.all)[13] <- "feature_x"
  
  ind.c = which(grepl('control',feat.df.all$Treatment))
  ctrl = feat.df.all[ind.c,]
  
  r.ctrl = range(ctrl$feature_x)
  
  ind.t = which(grepl('treated',feat.df.all$Treatment))
  trt = feat.df.all[ind.t,]
  trt$X1 = 1:1491
  r.trt = range(trt$feature_x)
  feat.df.all$X1 = 1:1848
  
  
  ##########################
  # Plot off B scores for feature j = 4, Total intensity of the nucleaus (Hoechst channel)
  ##########################
  
  
  # pdf(file = "figure_step3.3a_Zscores_nucsize.pdf",   # The directory you want to save the file in
  #     width = 7, # The width of the plot in inches
  #     height = 5) # The height of the plot in inches
  
  
  ggplot(feat.df.all, aes(X1, feature_x)) + geom_point(aes(colour = feature_x)) + 
    theme_bw()+
    geom_hline(yintercept=r.ctrl[1], linetype="dashed", color = "black") +
    geom_hline(yintercept=r.ctrl[2], linetype="dashed", color = "black") +
    scale_colour_gradient2(low = "blue", high = "red") + 
    # labs(title = feat, y = "Z scores" , x = "plate and replicate number") +
    labs(title = "Feature: Nucleus size", y = "Z scores" , x = "plate # replicate #") +
    scale_x_continuous(breaks = c(1, 309, 617, 925, 1233, 1541), 
                       labels=c("p1.r1", "p1.r2", "p1.r3", "p2.r1", "p2.r2", "p2.r3")) +
    ylim(-36.97533, 32.45456) +
    theme(plot.title = element_text(color="black", size=13),
          axis.title.x = element_text(color="black", size=14),
          axis.title.y = element_text(color="black", size=14),
          axis.text.x = X.text,
          axis.text.y = X.text,
          legend.title = element_blank()
          
    )
  
  
#  dev.off() 
  



