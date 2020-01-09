library(RColorBrewer)
library(ggplot2)
library(GGally)
library(latex2exp)
library(grid)
library(gridExtra)
library(tidyr)
library(dplyr)

source("plottingCode/plottingFunctions.R")

methsAll <- c("Lasso", "AdaLasso", "LENet", "HENet", "Ridge", "SCAD", "Stability", "Dantzig")

###############################
# Main Text Figures
#####

# Fig 1
methsIdx <- 1:8
r_lineplot(meths = methsIdx, corDesign = "Independence", SNR = 2,
           printPlot = TRUE, savePlot = TRUE, 
           plotWidth = 6, plotHeight = 5,
           plotPath = "manuscriptFigures", plotFileName = "Fig1.pdf")

# Fig 2
methsIdx <- c(1, 2, 4, 5, 6, 7)
lineplot_separate_factors(metric = "pAUC", meths = methsIdx, corDesign = "Independence", 
                          facetGridFormula = SNR ~ s0 + n,
                          s0 = c(10, 40), SNR = c(2, 0.5), n = c(300, 100), p = NULL,
                          x = "p", rotate_x_tickLabels = TRUE, errorbars = FALSE, 
                          printPlot = TRUE, savePlot = TRUE, 
                          plotWidth = 6, plotHeight = 4.66,
                          plotPath = "manuscriptFigures", plotFileName = "Fig2.pdf")

# Fig 3
methsIdx <- c(1, 2, 4, 5, 6, 7)
lineplot_separate_factors(metric = "pAUC", meths = methsIdx, corDesign = "high", 
                          facetGridFormula = SNR ~ s0 + n,
                          s0 = c(10, 40), SNR = c(2, 0.5), n = c(300, 100), p = 2000,
                          s0B = NULL,
                          includeIndependence = TRUE,
                          x = "s0B", rotate_x_tickLabels = FALSE, errorbars = FALSE, 
                          printPlot = TRUE, savePlot = TRUE, 
                          plotWidth = 6, plotHeight = 5,
                          plotPath = "manuscriptFigures", plotFileName = "Fig3.pdf")

# Fig 4 (to be included)


# Fig 5
methsIdx <- c(1, 2, 4, 5, 6)
lineplot_separate_factors(metric = "RMSE", meths = methsIdx, corDesign = "Independence", 
                          facetGridFormula = SNR ~ s0 + n,
                          s0 = c(10, 40), SNR = c(2, 0.5), n = c(300, 100), p = NULL,
                          x = "p", rotate_x_tickLabels = TRUE, errorbars = FALSE, 
                          scales = "free_y", 
                          printPlot = TRUE, savePlot = TRUE, 
                          plotWidth = 6, plotHeight = 4.66,
                          plotPath = "manuscriptFigures", plotFileName = "Fig5.pdf")

# Fig 6
methsIdx <- c(1, 2, 4, 5, 6)
lineplot_separate_factors(metric = "RMSE", meths = methsIdx, corDesign = "high", 
                          facetGridFormula = SNR ~ s0 + n,
                          s0 = c(10, 40), SNR = c(2, 0.5), n = c(300, 100), p = 2000,
                          s0B = NULL,
                          includeIndependence = TRUE,
                          x = "s0B", rotate_x_tickLabels = FALSE, errorbars = FALSE, 
                          scales = "free_y",
                          printPlot = TRUE, savePlot = TRUE, 
                          plotWidth = 6, plotHeight = 5,
                          plotPath = "manuscriptFigures", plotFileName = "Fig6.pdf")

# Fig 7 (to be included)


# Fig 8
methsIdx <- c(1, 2, 4, 6, 7)
lineplot_separate_factors(metric = c("TPR", "PPV"), meths = methsIdx, corDesign = "Independence", 
                          facetGridFormula = SNR ~ s0 + n,
                          s0 = c(10, 40), SNR = c(2, 0.5), n = c(300, 100), p = NULL,
                          x = "p", rotate_x_tickLabels = TRUE, errorbars = FALSE,
                          ylim = c(0, 1),
                          printPlot = TRUE, savePlot = TRUE, 
                          plotWidth = 6, plotHeight = 7,
                          plotPath = "manuscriptFigures", plotFileName = "Fig8.pdf")

# Fig 9
methsIdx <- c(1, 2, 4, 6, 7)
lineplot_separate_factors(metric = c("TPR", "PPV"), meths = methsIdx, corDesign = "high", 
                          facetGridFormula = SNR ~ s0 + n,
                          s0 = c(10, 40), SNR = c(2, 0.5), n = c(300, 100), p = 2000,
                          s0B = NULL,
                          includeIndependence = TRUE,
                          x = "s0B", rotate_x_tickLabels = FALSE, errorbars = FALSE,
                          ylim = c(0, 1),
                          printPlot = TRUE, savePlot = TRUE, 
                          plotWidth = 6, plotHeight = 7,
                          plotPath = "manuscriptFigures", plotFileName = "Fig9.pdf")

# Fig 10 (to be included)

# Figs 11-13 (not included)


###############################
# Supplementary Figures
#####

# Fig S1
methsIdx <- 1:7
r_lineplot(meths = methsIdx, corDesign = "Independence", SNR = 0.5,
           printPlot = TRUE, savePlot = TRUE, 
           plotWidth = 6, plotHeight = 5,
           plotPath = "manuscriptFigures", plotFileName = "FigS1.pdf")

# Fig S2
methsIdx <- c(1, 2, 4, 5, 6, 7)
r_lineplot(meths = methsIdx, corDesign = "high", SNR = 2, s0B = 5,
           printPlot = TRUE, savePlot = TRUE, 
           plotWidth = 6, plotHeight = 5,
           plotPath = "manuscriptFigures", plotFileName = "FigS2.pdf")

# Fig S3 (not included)

# Fig S4
methsIdx <- c(1, 2, 4, 5, 6, 7)
methods_vs_methods("pAUC", methsIdx, corDesign = "Independence",
                   colMap = "SNR", shapeMap = "r",
                   xylim = NULL, log_xy = FALSE, 
                   rotate_x_tickLabels = TRUE,
                   printPlot = TRUE, savePlot = TRUE, 
                   plotWidth = 6.5, plotHeight = 7,
                   plotPath = "manuscriptFigures", plotFileName = "FigS4.pdf")

# Fig S5 (not included)

# Fig S6
methsIdx <- c(1, 2, 4, 5, 6, 7)
methods_vs_methods("pAUC", methsIdx, corDesign = "low",
                   colMap = "SNR", shapeMap = "r",
                   xylim = NULL, log_xy = FALSE, 
                   rotate_x_tickLabels = TRUE,
                   printPlot = TRUE, savePlot = TRUE, 
                   plotWidth = 6.5, plotHeight = 7,
                   plotPath = "manuscriptFigures", plotFileName = "FigS6.pdf")

# Fig S7
methsIdx <- c(1, 2, 4, 5, 6, 7)
lineplot_separate_factors(metric = "pAUC", meths = methsIdx, corDesign = "high", 
                          facetGridFormula = SNR ~ s0 + n,
                          s0 = c(10, 40), SNR = c(4, 2, 1, 0.5), n = c(300, 100), p = 500,
                          s0B = NULL,
                          includeIndependence = TRUE,
                          x = "s0B", rotate_x_tickLabels = FALSE, errorbars = FALSE, 
                          printPlot = TRUE, savePlot = TRUE, 
                          plotWidth = 6, plotHeight = 8,
                          plotPath = "manuscriptFigures", plotFileName = "FigS7.pdf")


# Fig S8 (to be included)


# Fig S9
methsIdx <- c(1, 2, 4, 5, 6)
methods_vs_methods("RMSE", methsIdx, corDesign = "Independence",
                   colMap = "SNR", shapeMap = "r",
                   xylim = NULL, log_xy = FALSE, 
                   rotate_x_tickLabels = TRUE,
                   printPlot = TRUE, savePlot = TRUE, 
                   plotWidth = 6.5, plotHeight = 7,
                   plotPath = "manuscriptFigures", plotFileName = "FigS9.pdf")
# Fig S10
methods_vs_methods("RMSE", methsIdx, corDesign = "low",
                   colMap = "SNR", shapeMap = "r",
                   xylim = NULL, log_xy = FALSE, 
                   rotate_x_tickLabels = TRUE,
                   printPlot = TRUE, savePlot = TRUE, 
                   plotWidth = 6.5, plotHeight = 7,
                   plotPath = "manuscriptFigures", plotFileName = "FigS10.pdf")

# Fig S11
methsIdx <- c(1, 2, 4, 5, 6)
lineplot_separate_factors(metric = "RMSE", meths = methsIdx, corDesign = "high", 
                          facetGridFormula = SNR ~ s0 + n,
                          s0 = c(10, 40), SNR = c(4, 2, 1, 0.5), n = c(300, 100), p = 500,
                          s0B = NULL,
                          includeIndependence = TRUE,
                          x = "s0B", rotate_x_tickLabels = FALSE, errorbars = FALSE, 
                          scales = "free_y",
                          printPlot = TRUE, savePlot = TRUE, 
                          plotWidth = 6, plotHeight = 8,
                          plotPath = "manuscriptFigures", plotFileName = "FigS11.pdf")

# Fig S12 (to be included)


# Fig S13
methsIdx <- c(1, 2, 4, 6, 7)
methods_vs_methods("TPR", methsIdx, corDesign = "Independence",
                   colMap = "SNR", shapeMap = "r",
                   xylim = NULL, log_xy = FALSE, 
                   rotate_x_tickLabels = TRUE,
                   printPlot = TRUE, savePlot = TRUE, 
                   plotWidth = 6.5, plotHeight = 7,
                   plotPath = "manuscriptFigures", plotFileName = "FigS13.pdf")
# Fig S14
methods_vs_methods("PPV", methsIdx, corDesign = "Independence",
                   colMap = "SNR", shapeMap = "r",
                   xylim = NULL, log_xy = FALSE, 
                   rotate_x_tickLabels = TRUE,
                   printPlot = TRUE, savePlot = TRUE, 
                   plotWidth = 6.5, plotHeight = 7,
                   plotPath = "manuscriptFigures", plotFileName = "FigS14.pdf")
# Fig S15
methods_vs_methods("TPR", methsIdx, corDesign = "low",
                   colMap = "SNR", shapeMap = "r",
                   xylim = NULL, log_xy = FALSE, 
                   rotate_x_tickLabels = TRUE,
                   printPlot = TRUE, savePlot = TRUE, 
                   plotWidth = 6.5, plotHeight = 7,
                   plotPath = "manuscriptFigures", plotFileName = "FigS15.pdf")
# Fig S16
methods_vs_methods("PPV", methsIdx, corDesign = "low",
                   colMap = "SNR", shapeMap = "r",
                   xylim = NULL, log_xy = FALSE, 
                   rotate_x_tickLabels = TRUE,
                   printPlot = TRUE, savePlot = TRUE, 
                   plotWidth = 6.5, plotHeight = 7,
                   plotPath = "manuscriptFigures", plotFileName = "FigS16.pdf")

# Fig S17
methsIdx <- c(1, 2, 4, 6, 7)
lineplot_separate_factors(metric = c("TPR", "PPV"), meths = methsIdx, corDesign = "high", 
                          facetGridFormula = SNR ~ s0 + n,
                          s0 = c(10, 40), SNR = c(4, 2, 1, 0.5), n = c(300, 100), p = 500,
                          s0B = NULL,
                          includeIndependence = TRUE,
                          x = "s0B", rotate_x_tickLabels = FALSE, errorbars = FALSE,
                          ylim = c(0, 1),
                          printPlot = TRUE, savePlot = TRUE, 
                          plotWidth = 6, plotHeight = 9,
                          plotPath = "manuscriptFigures", plotFileName = "FigS17.pdf")

# Fig S18 (to be included)


# Figs S19-S21 (not included)