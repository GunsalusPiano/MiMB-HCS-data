# Outline

Quality control measures and statistical strategies to address the challenges of high-content phenotypic data

Yanthe E. Pearson and Kristin C. Gunsalus


Outline:  <br />

Introduction  <br />
Material <br />
2.1 Data format <br />
2.2 Software <br />
2.3 R packages <br />
Methods <br />
3.1 Raw data visualizations <br />
3.1.1 Cell counts scatterplot <br />
3.1.2 Cell feature distributions <br />
3.1.3 Well aggregate heatmaps <br />
3.2 Two-way ANOVA <br />
3.3 Plate normalization and cell standardization <br />
3.4: Earth mover’s distance <br />

# R scripts:
script_step3.1.1_counts_scatter.R <br />
script_step3.1.2_replicates_nocodazole.R <br />
script_step3.1.3_nucsize_heatmap.R <br />
script_step3.2_anova.R <br />
script_step3.3a_medianpolish.R <br />
script_step3.3b_adjustmentstep.R <br />
script_step3.3c_BZscores.R <br />
script_step3.3d_compare_data.R<br /> 

# Data files:
data_step3.1.1_counts_and_medians_panelA.csv* <br />
data_step3.1.2_density_replicates_nocodazole.csv <br />
data_step3.1.2_cell_data_control.csv <br />
data_step3.3a_well_medians.csv <br />
data_step3.3b_adjustments.csv (maybe change this to adjustment_amount) <br />
data_step3.3b_cell_panelA_subset.csv <br />
data_step3.3c_cells_adjusted.csv <br />


# Figure files:
Figure_step3.1.1_scatter.pdf (Fig 1) <br />
Figure_step3.1.1_scatter_legend.pdf (Fig 1 legend) <br />
Figure_step3.1.2_density.pdf (Fig 2) <br />
Figure_step3.1.2_density_legend (Fig 2 legend) <br />
Figure_step3.1.3_heatmap.pdf (Fig 3)  <br />
Figure_step3.2_row.pdf (Fig 4 top) <br />
Figure_step3.2_col.pdf (Fig 4 bottom) <br />
Figure_step3.2_RNA_all.pdf (Fig 5 left) <br />
Figure_step3.2_RNA_intensity.pdf (Fig 5 right) <br />
Figure_step3.3_flowchart.pdf (Fig 6: Flowchart) <br />
Figure_step3.3a_Bscores_totalinten.pdf (Fig 7: B scores) <br />
Figure_step3.3a_Zscores_nucsize.pdf (Fig 7: Z scores) <br />

*data file to be used three times, for steps 3.1.1, 3.1.3 and 3.2.
