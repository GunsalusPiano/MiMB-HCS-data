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
# R script for MiMB BookChapter - Methods section 3.2: Two-way ANOVA
########################

########################
# R script for figure 4: Two-way ANOVA test for detection of plate effects
# R script for figure 5: Summary of row effect detection for RNA channel (Syto14) features.   
########################

########################
# Description
# Step 1: Import panel A well level median aggregates measured from six plates
# Step 2: Set up a data frame for saving the two-way ANOVA output: test statistic (F) and corresponding p-values.
# Step 3: Run the analysis one plate and one feature at a time
# Step 4: Add annotations to the two-way ANOVA analysis results file
# Step 5: Add additional metadata for distinguishing intensity features from non intensity features
# Step 6: Observe which p-value or a range of p-values distinguishes plates and features with and without detectable technical variation among the control wells.
# Step 7: Add metadata for distinguishing features measured in different fluorescent channels
########################

 # clear variables
 rm(list=ls())

 ########################
 # libraries
 ########################

 library(ggplot2)
 library(RColorBrewer)
 library(readr)
 library(car)

########################
# set working directory
########################

setwd("~/Location/of/Data")

  ########################
  # (1) Import the annotated well median data file from panel A channels (Hoechst, syto14, mito, WGA)
  ########################

  well_med_dat <- read_csv("data_step3.1.1_counts_and_medians_panelA.csv")
  names(well_med_dat)[1:20]

  # Identify empty wells 
  empty = which(grepl("empty", well_med_dat$Treatment)) 
  well_med_dat[empty,13:dim(well_med_dat)[2]] = NA

  #######################
  # (2) Set up a data frame for saving the ANOVA test statistic (F statistic) and corresponding p values.
  # Column names are 'feature', 'plate', 'F_r', 'P_r', 'F_c', 'P_c', the subscripts 'r' and 'c' stand for row and column.
  #######################
    
    feature <- c("feature","feature","feature")
    plate <- c("plate1","plate1","plate1")
    
    F_r <- c(0.11, 0.11, 0.11)
    P_r <- c(0.11, 0.11, 0.11)
    F_c <- c(0.11, 0.11, 0.11)
    P_c <- c(0.11, 0.11, 0.11)
    
    ANOVA.results <- data.frame(feature, plate, F_r, P_r, F_c, P_c, stringsAsFactors=FALSE)
    ANOVA.results

    #######################
    # (3) Run the ANOVA analysis one plate and one feature at a time
    # The ANOVA model is analyzing only the DMSO control wells
    #######################
    
    # Numerical data (without metadata columns)
    well_med_dat.num = well_med_dat[,13:dim(well_med_dat)[2]]
    names(well_med_dat.num)
    
    # Plates are names based on plate # and replicate #, 1A1 is plate 1 replicate 1: "1A1" "1A2" "1A3" "2A1" "2A2" "2A3", 
    plate.unique = sort(unique(well_med_dat$PlateNumber))
    
    # initialize ss, ss is used for saving model outputs into "Anova.results" dataframe.
    ss = 0
    
    #######################
    # Parse through each column (feature j) and each plate (i)
    #######################
    
    for (j in seq(1,length(well_med_dat.num))) { 
  
            ff = names(well_med_dat.num)[j]
          
    for (i in seq(1,length(plate.unique))) {   
    
    pp = plate.unique[i]
            
    plate.ind = which(well_med_dat$PlateNumber == plate.unique[i])
    plate.current = well_med_dat.num[plate.ind,j]
    plate.current.labels = well_med_dat[plate.ind,]
    plate.current.ord <- plate.current.labels[order(plate.current.labels$WellId),]
    
    # The data (feature only -- dependent variable)
    plate.use <- plate.current.ord[,12+j]
       
    # The data (feature only -- dependent variable ~ factors (row and column position))
    plate.use1 <- plate.current.ord[,c(4,7,11,12,12+j,12+j)]
    names(plate.use1)[6] = "feature_x"
    
    # Model the controls only
    plate.c = plate.current.ord[which('control' == plate.current.ord$Treatment),c(4,7,11,12,12+j,12+j)]
    names(plate.c)[6] = "feature_x"
    
    ###########################################################
    # (2) Set up "control" data for Two way anova Row + Col --- 
    # A two-way ANOVA assesses whether the factors (categories) 
    # have a significant effect on the outcome (measurement/feature). 
    # Its implementation can be conceptualized as a regression analysis where each row and column 
    # level is treated as a dummy variable. 
    ###########################################################
    
     plate.c$Row = as.factor(plate.c$Row)
     plate.c$Column = as.factor(plate.c$Column)
  
    ###########################################################
    # Implementation of the two factor anova
    ###########################################################
   
    M <- lm(feature_x ~ Row + Column, data=plate.c)
    # summary(M)

    MA = Anova(M)
    row_F = MA$`F value`[1]
    col_F = MA$`F value`[2]
    
    row_p = MA$`Pr(>F)`[1]
    col_p = MA$`Pr(>F)`[2]
    
    ss = ss + 1
    
    ANOVA.results[ss,1:6] <- data.frame(ff, pp, row_F, row_p, col_F, col_p, stringsAsFactors=FALSE)
    
    print(names(plate.use))

          }
          }
    
    
    #######################
    # After all plates and features have been analysed based on their control well values we add metadata to this *ANOVA.results* data frame.
    # (4) Add annotations to the ANOVA results file *ANOVA.results*
    # The columns *Channel* and *Channel_color* 
    #######################
    
    df.new = ANOVA.results
    
    
    #######################
    # (5) Add more metadata (column: feature_type) for distinguishing intensity features from non intensity features
    #######################
    
    # Intensity features
    df.new$feature_type = "non-int"
    ind.int = which(grepl("Inten", df.new$feature))
    df.new$feature_type[ind.int] = "int"

    # Log transform the p-values
    df.new$P_r_log = log(df.new$P_r)
    df.new$P_c_log = log(df.new$P_c)
   
    # add columns for Negative (Log (p) )
    df.new$P_r_log_neg = -(df.new$P_r_log)
    df.new$P_c_log_neg = -(df.new$P_c_log)
    
    # add column for plate labels
    df.new$plate_rep = df.new$plate
    in1 = which(df.new$plate == "1A1")
    df.new$plate_rep[in1] = "p1.r1" 
    in2 = which(df.new$plate == "1A2")
    df.new$plate_rep[in2] = "p1.r2" 
    in3 = which(df.new$plate == "1A3")
    df.new$plate_rep[in3] = "p1.r3" 
    
    in4 = which(df.new$plate == "2A1")
    df.new$plate_rep[in4] = "p2.r1" 
    in5 = which(df.new$plate == "2A2")
    df.new$plate_rep[in5] = "p2.r2" 
    in6 = which(df.new$plate == "2A3")
    df.new$plate_rep[in6] = "p2.r3" 
    
    
    #######################
    # (6) COLOR CODE THE CURVES BY CHANNEL
    #######################
    
       # df.save = df.new
      
       
      for (i in seq(1,length(df.new$feature))) {

        feat = df.new$feature[i]

        print(feat)

        if (grepl("cell_count",feat)) {
          print("cell_count")
          df.new$Channel[i] = "cell_count"
          df.new$Channel_color[i] = "purple"
        }

        if (grepl("NUC_A",feat)) {
          print("NUC_A")
          df.new$Channel[i] = "NUC_A"
          df.new$Channel_color[i] = "blue"

        }

        if (grepl("RNA",feat)) {
          print("RNA")
          df.new$Channel[i] = "RNA"
          df.new$Channel_color[i] = "yellow"

        }

        if (grepl("PMG",feat)) {
          print("PMG")
          df.new$Channel[i] = "PMG"
          df.new$Channel_color[i] = "yellow"

        }

        if (grepl("MITO",feat)) {
          print("MITO")
          df.new$Channel[i] = "MITO"
          df.new$Channel_color[i] = "red"


        }

        if (grepl("Tex",feat)) {
          print("Tex")
          df.new$Channel[i] = "NUC_A"
          df.new$Channel_color[i] = "lightblue"

        }

      }



        ##############################
        # Add columns for fluorescent channel annotations
        ##############################
        
        ind1 = which(df.new$Channel == "cell_count")
        ind2 = which(df.new$Channel == "NUC_A")
        ind6 = which(df.new$Channel == "NUC_Tex")
        ind3 = which(df.new$Channel == "RNA")
        ind4 = which(df.new$Channel == "PMG")
        ind5 = which(df.new$Channel == "MITO")
        
        df.new$Channel_new = df.new$Channel

        df.new$Channel_new[ind1] = "Cell_count_panelA"
        df.new$Channel_new[ind2] = "Hoechst_panelA"
        df.new$Channel_new[ind6] = "Hoechst_panelA"
        df.new$Channel_new[ind3] = "SYTO14"
        df.new$Channel_new[ind4] = "WGA-AlexaFluor555"
        df.new$Channel_new[ind5] = "MitoTrackerDR"
        
        
        ########################
        # Add new column
        ########################
        
        df.new$cell_struc = df.new$Channel
        
        ########################
        # editing labels 
        ########################
        
        df.new$cell_struc[ind1] = "Cell count"
        df.new$cell_struc[ind2] = "DNA (Hoechst)"
        df.new$cell_struc[ind6] = "DNA (Hoechst)"
        df.new$cell_struc[ind3] = "RNA"
        df.new$cell_struc[ind4] = "Golgi and membranes"
        df.new$cell_struc[ind5] = "Mitochondria"
        

        ########################
        # Prepare ggplot plot parameters
        ########################
       
        Y.text <- element_text(color = "black", size =10, family = "Helvetica")
        X.text <- element_text(color = "black", size = 10, family = "Helvetica")
        a.text = element_text(color = "black", size = 10, family = "Helvetica")
        t.text = element_text(color = "black", size = 10, family = "Helvetica")
        strip.text.a = element_text(color = "black", size = 12, family = "Helvetica")
          
          ################################
          # Column effects: Plot of all negative p-values (P_c) after testing for possible row effects (row dependencies or patterns)
          ################################
          
          # Remove cell counts
          ind = which(df.new$feature == "cell_count")
          df.new1 = df.new[-ind,]
          
          
          ################################
          # Row effects: Plot of all negative p-values (P_r) after testing for possible row effects (row dependencies or patterns)
          ################################
          
          pdf(file = "figure_step3.2_row.pdf",   # The directory you want to save the file in
              width = 11, # The width of the plot in inches
              height = 3) # The height of the plot in inches
          
          ggplot(df.new1) + 
          aes(x = plate_rep, y = P_r_log_neg, group = feature, colour = cell_struc) +
          geom_line(linetype = 1, size = .6, show.legend = FALSE) + 
          ylim(min(df.new$P_r_log_neg),max(df.new$P_r_log_neg)) +
          theme_bw() +
          # adding horizontal lines to the graphs here 
          geom_hline(yintercept=-log(0.00001), linetype = "dashed", color = "black", size = 0.3) +
          geom_hline(yintercept=-log(0.001), linetype = "dotdash", color = "black", size = 0.3) +
          # graph labels are specified here
          labs(x = "Plate # replicate #", y = "- Log(P-value)", title = "Two-way ANOVA test for row effects") + 
          theme(strip.background = element_rect(fill="white"), strip.text = strip.text.a, axis.text.x = X.text, 
                  text = t.text, axis.text.y = Y.text, axis.title = a.text, title = t.text) +
          facet_wrap(~cell_struc, nrow=1)
          
          dev.off()
          
          ################################
          # Column effects: Plot of all negative p-values (P_c) after testing for possible row effects (row dependencies or patterns)
          ################################
          
          
          pdf(file = "figure_step3.2_col.pdf",   # The directory you want to save the file in
              width = 11, # The width of the plot in inches
              height = 3) # The height of the plot in inches
          
          ggplot(df.new1) + 
            
          aes(x = plate_rep, y = P_c_log_neg, group = feature, colour = cell_struc) +
          geom_line(linetype = 1, size = .6, show.legend = FALSE) + 
            
          ylim(min(df.new$P_c_log_neg),max(df.new$P_r_log_neg)) +
            
          theme_bw() +
            
          geom_hline(yintercept=-log(0.00001), linetype = "dashed", color = "black", size = 0.3) +
          geom_hline(yintercept=-log(0.001), linetype = "dotdash", color = "black", size = 0.3) +
            
          labs(x = "Plate # replicate #", y = "- Log(P-value)", title = "Two-way ANOVA test for column effects") + 
            
          theme(strip.background = element_rect(fill="white"), strip.text = strip.text.a, axis.text.x = X.text, 
                  axis.text.y = Y.text, axis.title = a.text, title = t.text, text = t.text) +
            
          facet_wrap(~cell_struc, nrow = 1)
          
          dev.off()
          
          ################################
          # Check RNA - Row effect
          ################################
          
          
           # Create new dataframe for RNA features only - 16 unique features
          df.new2 = df.new1[grep("RNA", df.new1$feature),]
          # 16 unique features
          length(unique(df.new2$feature))
          
          
          ################################
          # Create color wheel of 16 colors
          ################################
          colors25 <- c(
            "dodgerblue2", "#E31A1C", # red
            "green4",
            "#6A3D9A", # purple
            "#FF7F00", # orange
            "black", "gold1",
            "skyblue2", "#FB9A99", # lt pink
            "palegreen2",
            "#CAB2D6", # lt purple
            "#FDBF6F", # lt orange
            "gray70", "khaki2",
            "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
            "darkturquoise", "green1", "yellow4", "yellow3",
            "darkorange4", "brown"
          )
          pie(rep(1, 25), col = colors25)
          
          
          ################################
          # The plot of RNA feature performance 
          ################################
          
           pdf(file = "figure_step3.2_RNA_all.pdf",   # The directory you want to save the file in
               width = 7, # The width of the plot in inches
               height = 5) # The height of the plot in inches
          
          
            ggplot(df.new2, aes(x = plate_rep, y = P_r_log_neg, group = feature)) + 
            geom_line(linetype = 1, size = .8, show.legend = TRUE, aes(color = feature)) + 
            
            # color of curves
            scale_color_manual(values = colors25, name = "Features") +
            theme_bw() +
            
            geom_hline(yintercept=-log(0.00001), linetype = "dashed", color = "black", size = 0.3) +
            geom_hline(yintercept=-log(0.001), linetype = "dotdash", color = "black", size = 0.3) +
            
            labs(x = "Plate # replicate #", y = "- Log(P-value)", title = "RNA channel features") + 
            theme(strip.background = element_rect(fill="white"), strip.text = strip.text.a, axis.text.x = X.text, 
                  axis.text.y = Y.text, axis.title = a.text, title = t.text, text = t.text) 
          
           dev.off()
          
          ################################
          # The plot of RNA feature performance - color features based on intensity or non- intensity related features
          ################################
          
          
          # Step 1: Call the pdf command to start the plot
          pdf(file = "figure_step3.2_RNA_intensity.pdf",   # The directory you want to save the file in
              width = 5.5, # The width of the plot in inches
              height = 5) # The height of the plot in inches

          ggplot(df.new2, aes(x = plate_rep, y = P_r_log_neg, group = feature)) + 
          geom_line(linetype = 1, size = .8, aes(color = feature_type)) + 
              
          # color of curves
          # scale_color_manual(values = c(colors25[18],colors25[8]), name = "Feature type") +
          
          scale_color_manual(values = c(colors25[4],colors25[16]), name = "Feature type") +
        
          theme_bw() +
          geom_hline(yintercept=-log(0.00001), linetype = "dashed", color = "black", size = 0.3) +
          geom_hline(yintercept=-log(0.001), linetype = "dotdash", color = "black", size = 0.3) +
              
          theme(legend.position="top") +
          labs(x = "Plate # replicate #", y = "- Log(P-value)", title = NULL) + 
          theme(strip.background = element_rect(fill="white"), strip.text = strip.text.a, axis.text.x = X.text, 
                    axis.text.y = Y.text, axis.title = a.text, title = t.text, text = t.text,
                legend.text = element_text(colour="black", size=11, family = "Helvetica"))
                                            
              
          dev.off()
          
          
