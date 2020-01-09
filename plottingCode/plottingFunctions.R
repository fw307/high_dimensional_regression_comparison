
#############################
# ggplot theme
#####
ggplotThemePub <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(text = element_text(color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(face = "bold", size = rel(1.2)),
          axis.text = element_text(size = rel(1), colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.title = element_text(face = "bold", size = rel(1)),
          legend.text = element_text(size = rel(1)),
          legend.background = element_rect(fill = "transparent"),
          legend.key.size = unit(1, 'lines'),
          legend.position = "bottom",
          legend.direction = "horizontal",
          panel.border = element_rect(colour = "black"),
          panel.grid.minor = element_line(colour = NA),
          strip.text = element_text(face = "bold", size = rel(1), colour = "black",
                                    margin = margin(1.5, 1.5, 1.5, 1.5)),
          plot.tag = element_text(face = "bold", size = rel(1.4), hjust = 0, vjust = 1), 
          plot.tag.position = c(0, 1)
    )
}


#############################
# helper functions
#####
getAxisLegendLabel <- function(param) {
  switch(param,
         param,
         rho = TeX("$\\textbf{\\rho}$"),
         pB = TeX("$\\textbf{p^B}$"),
         s0B = TeX("$\\textbf{s_0^B}$"),
         s0 = TeX("$\\textbf{s_0}$"),
         rDisc = "r",
         corDesign = "Design")
}

getFactorLevelLabels <- function(param, lvls = NULL) {
  switch(param,
         rDisc = unname(TeX(c("$0 < r \\leq 1$", "$1 < r \\leq 2$", "$2 < r \\leq 3$", "$r > 3$"))),
         rho = unname(TeX(paste0("$\\textbf{\\rho = ", lvls, "}$"))),
         pB = unname(TeX(paste0("$\\textbf{p^B = ", lvls, "}$"))),
         s0B = unname(TeX(paste0("$\\textbf{s_0^B = ", lvls, "}$"))),
         s0 = unname(TeX(paste0("\\textbf{$s_0 = ", lvls, "}$"))),
         n = unname(TeX(paste0("\\textbf{n = ", lvls, "}"))),
         p = unname(TeX(paste0("\\textbf{p = ", lvls, "}"))),
         SNR = unname(TeX(paste0("\\textbf{SNR = ", lvls, "}"))))
}

setNullParams <- function(df, corDesign = NA, rho = NA, pB = NA, s0B = NA,
                          n = NA, p = NA, s0 = NA, SNR = NA) {
  out <- list()
  
  if(is.null(corDesign)) {
    corDesign <- sort(unique(df$corDesign))
    out <- c(out, list(corDesign = corDesign))
  }
  if(is.null(pB)) {
    pB <- setdiff(sort(unique(df$pB)), 0)
    out <- c(out, list(pB = pB))
  }
  if(is.null(rho)) {
    rho <- setdiff(sort(unique(df$rho)), 0)
    out <- c(out, list(rho = rho))
  }
  if(is.null(s0B)) {
    s0B <- setdiff(sort(unique(df$s0B)), 0)
    out <- c(out, list(s0B = s0B))
  }
  if(is.null(n)) {
    n <- sort(unique(df$n))
    out <- c(out, list(n = n))
  }
  if(is.null(p)) {
    p <- sort(unique(df$p))
    out <- c(out, list(p = p))
  }
  if(is.null(s0)) {
    s0 <- sort(unique(df$s0))
    out <- c(out, list(s0 = s0))
  }
  if(is.null(SNR)) {
    SNR <- sort(unique(df$SNR))
    out <- c(out, list(SNR = SNR))
  }
  return(out)
}

lineplot_combine_2metrics <- function(gt1, gt2) {
  
  nvalsRow <- length(grep("axis-r", gt1$layout$name))
  
  ylab1grob <- gtable::gtable_filter(gt1, "ylab-l")$grobs[[1]]
  ylab2grob <- gtable::gtable_filter(gt2, "ylab-l")$grobs[[1]]
  hposPanels <- unique(gt1$layout$t[grep("panel", gt1$layout$name)])
  gt1 <- gtable::gtable_add_grob(gt1, rep(list(ylab1grob), length(hposPanels)), t = hposPanels, l = 3)
  gt2 <- gtable::gtable_add_grob(gt2, rep(list(ylab2grob), length(hposPanels)), t = hposPanels, l = 3)
  
  gt <- gt1[1:(hposPanels[1]-1), ]
  for(i in 1:(length(hposPanels) - 1)) {
    gt <- gridExtra::gtable_rbind(gt, 
                                  gt1[hposPanels[i]:(hposPanels[i] + 1), ],
                                  gt2[hposPanels[i]:(hposPanels[i] + 1), ], 
                                  size = "first")
  }
  gt <- gridExtra::gtable_rbind(gt, 
                                gt1[hposPanels[i + 1], ],
                                gt2[(hposPanels[i + 1]-1):nrow(gt2), ], size = "first")
  
  hposPanels <- unique(gt$layout$t[grep("panel", gt$layout$name)])
  hposSmallGaps <- hposPanels[seq(1, length(hposPanels), 2)] + 1
  gt$heights[hposSmallGaps] <- 0 * gt$heights[hposSmallGaps] # smaller gaps betwen TPR/PPV rows
  hposLargeGaps <- hposPanels[seq(2, length(hposPanels) - 1, 2)] + 1
  gt$heights[hposLargeGaps] <- 1.2 * gt$heights[hposLargeGaps]
  
  
  # right facet strips span two rows
  # code below adapted from https://stackoverflow.com/questions/40316169/nested-facets-in-ggplot2-spanning-groups
  stripLocations <- grep("strip-r", gt$layout$name)
  strip <- gtable::gtable_filter(gt, "strip-r", trim = FALSE)
  nStrip <- nrow(strip$layout)
  nval <- length(grep("axis-r", gt1$layout$name))
  
  r <- strip$layout$r[1]
  top <- strip$layout$t[nStrip/nval * (0:(nval - 1)) + 1]
  b   <- strip$layout$b[nStrip/nval * (1:nval)]
  
  vec   <- vector("list", length = 3)
  vec[] <- list(zeroGrob())
  multiRowStrip <- gtable::gtable_col("strip-r", vec, width = unit(1, "null"), 
                                      heights = unit(c(1, 0, 1), "null"))
  
  for(i in 1:nval) {
    combinedStrip <- gtable::gtable_add_grob(multiRowStrip, 
                                             gt$grobs[[stripLocations[2 * (i-1) + 1]]]$grobs[[1]], 
                                             t = 1, l = 1, b = length(multiRowStrip))
    gt <- gtable::gtable_add_grob(gt, combinedStrip, 
                                  t = top[i], b = b[i],  l = r, 
                                  name = paste0("strip-r-", 2 * (i-1) + 1, "-", 2 * i))
  }
  
  invisible(gt)
  
}

#############################
# rescaled sample size, r, line plots (used for Figs 1, S1, S2)
#####
r_lineplot <- function(meths, corDesign, SNR, rho = NULL, pB = NULL, s0B = NULL,
                       printPlot = TRUE,
                       savePlot = FALSE, plotWidth = 7, plotHeight = 7,
                       plotUnits = "in",
                       plotPath = ".", plotFileName = "r_lineplot.pdf") {
  #####
  # meths: subset of c("Lasso", "AdaLasso", "LENet", "HENet", "Ridge", "SCAD", "Stability", "Dantzig")
  #        or numeric vector giving corresponding indexes
  # corDesign: one of "Independence", "Pairwise" and "Toeplitz" (for synthetic data), or,
  #            "low" and "high" (for semi-synthetic data)
  # SNR: one of c(0.5, 1, 2, 4)
  # rho: one of c(0.5, 0.7, 0.9); applies for Pairwise corDesign only (Toeplitz has rho fixed to 0.7)
  # pB: subset of c(10, 100); applies for Pairwise corDesign only (Toeplitz and high have pB fixed to 100 and 10 respectively)
  # s0B: subset of c(1, 2, 5); applies for Pairwise and high corDesigns only (Toeplitz has s0B fixed to 2)
  # printPlot: logical; print plot?
  # savePlot: logical: save plot to pdf?
  # plotWidth, plotHeight, plotUnits: dimensions for pdf
  # plotPath, plotFileName: path and filename to save pdf to; defaults to working directory with filename r_lineplot.pdf
  # 
  
  #####  
  methsAll <- c("Lasso", "AdaLasso", "LENet", "HENet", "Ridge", "SCAD", "Stability", "Dantzig")
  if((is.character(meths) && !all(tolower(meths) %in% tolower(methsAll))) ||
     (is.numeric(meths) && !all(meths %in% 1:8))) {
    stop("meths should be a subset of c('Lasso', 'AdaLasso', 'LENet', 'HENet', 'Ridge', 'SCAD', 'Stability', 'Dantzig')
          or a subset of 1:8")
  } 
  
  if(is.numeric(meths)) {
    meths <- methsAll[meths]
  }
  
  colours <- (brewer.pal(9, "Set1"))
  colours <- c(colours[c(1, 2, 9, 3)], "#EEEE00", "#000000", colours[c(4, 7)])
  names(colours) <- tolower(methsAll)
  colours <- colours[methsAll %in% meths]
  
  metricNames <- c("pAUC", "RMSE", "TPR", "PPV")
  
  if(corDesign %in% c("Independence", "Pairwise", "Toeplitz")) {
    load("results/synthetic/allSettingsAveraged.RData")
  } else if(corDesign %in% c("low", "high")) {
    load("results/semi-synthetic/allSettingsAveraged.RData")  
  } else {
    stop("corDesign should be one of 'Independence', 'Pairwise', 'Toeplitz', 'low' or 'high'")
  }
  
  if(corDesign %in% c("Independence", "low", "Toeplitz")) {
    if(!is.null(pB)) {
      warning("Ignoring argument: pB")
    }
    if(!is.null(s0B)) {
      warning("Ignoring argument: s0B")
    }
    if(!is.null(rho)) {
      warning("Ignoring argument: rho")
    }
    plotData <- lapply(scores, function(df) {
      tmp <- df[df$SNR == SNR & df$corDesign == corDesign & 
                  df$method %in% tolower(meths), ]
      tmp$method <- factor(tmp$method, tolower(meths))
      return(tmp)
    })
  } else if(corDesign == "high") {
    if(!is.null(pB)) {
      warning("Ignoring argument: pB")
    }
    if(!is.null(rho)) {
      warning("Ignoring argument: rho")
    }
    if(is.null(s0B) || is.na(s0B)) {
      stop("Argument s0B is required")
    }
    plotData <- lapply(scores, function(df) {
      tmp <- df[df$SNR == SNR & df$corDesign == corDesign & 
                  df$s0B == s0B &
                  df$method %in% tolower(meths), ]
      tmp$method <- factor(tmp$method, tolower(meths))
      return(tmp)
    })
  } else {
    if(is.null(pB) || is.na(pB)) {
      stop("Argument pB is required")
    }
    if(is.null(rho) || is.na(rho)) {
      stop("Argument rho is required")
    }
    if(is.null(s0B) || is.na(s0B)) {
      stop("Argument s0B is required")
    }
    plotData <- lapply(scores, function(df) {
      tmp <- df[df$SNR == SNR & df$corDesign == corDesign & 
                  df$rho == rho & df$pB == pB & df$s0B == s0B &
                  df$method %in% tolower(meths), ]
      tmp$method <- factor(tmp$method, tolower(meths))
      return(tmp)
    })
  }
  
  mrgn <- list(c(5.5, 10.5, -2,  5.5),
               c(5.5, 5.5,  -2,  5.5),
               c(5.5, 10.5, 5.5, 5.5),
               c(5.5, 5.5,  5.5, 5.5))
  
  ggp <- mapply(function(m, sl, mrgn) {
    ggp <- ggplot(plotData[[tolower(m)]], aes(x = r, y = mean, colour = method)) + 
      geom_line(size = 0.5, show.legend = sl) +
      geom_point(size = 0.9, show.legend = sl) +
      guides(colour = guide_legend(nrow = 1)) +
      scale_x_log10(breaks = c(seq(0.3, 0.9, 0.1), 1:5), 
                    labels = c(cbind(
                      rbind("", seq(0.4, 1, 0.2)),
                      rbind(seq(2, 4, 2), "")) )) +
      scale_colour_manual(name = "Method", values = colours, labels = meths) +
      ylab(m) +
      ggplotThemePub() +
      theme(plot.margin = unit(mrgn, "pt"))
    
    if(tolower(m) != "rmse") {
      ggp <- ggp + ylim(0,1)
    }
    return(ggp)
  }, 
  metricNames, c(TRUE, FALSE, FALSE, FALSE), mrgn,
  SIMPLIFY = FALSE)
  
  ggp[[1]] <- ggp[[1]] + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank(), 
                               legend.title = element_blank())
  ggp[[2]] <- ggp[[2]] + theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank())
  ggp[[3]] <- ggp[[3]] + theme(axis.title.x = element_blank())
  ggp[[4]] <- ggp[[4]] + xlab(TeX("$\\textbf{r = n / s_0 \\log(p-s_0)}$"))
  
  
  row1 <- do.call(gtable_cbind, c(lapply(list(ggp[[1]] + theme(legend.position = "none"), ggp[[2]]), ggplotGrob), size = "last"))
  row2 <- do.call(gtable_cbind, c(lapply(list(ggp[[3]], ggp[[4]]), ggplotGrob), size = "last"))
  gt <- do.call(gtable_rbind, c(list(row1, row2), size  = "last"))
  
  # move x-axis label to centre of bottom plots
  xlabGrobIdx <- which(gt$layout$name == "xlab-b" & !grepl("zeroGrob", sapply(gt$grobs, "[[", "name")))[1]
  lExtent <-  min(gt$layout$l[which(gt$layout$name == "panel")])
  gt$layout[xlabGrobIdx, "l"] <- lExtent
  
  leg <- ggpubr::get_legend(ggp[[1]])
  h <- grobHeight(leg) + unit(0.5, "line")
  gt <- gtable::gtable_add_rows(gt, heights = h, -1)
  gt <- gtable::gtable_add_grob(gt, leg, t = nrow(gt), l = 1, r = ncol(gt))
  
  if(printPlot) {
    grid.newpage()
    grid.draw(gt)
  }
  
  if(savePlot) {
    ggsave(filename = plotFileName, path = plotPath, 
           plot = gt, width = plotWidth, height = plotHeight, units = plotUnits)
  } 
} 



#############################
# methods vs methods scatter plot (used for Figs S4, S6, S9, S10, S13-S16)
#####
methods_vs_methods <- function(metric, meths,
                               corDesign = NULL, 
                               s0 = NULL, SNR = NULL, 
                               n = NULL, p = NULL, 
                               rho = NULL, pB = NULL, s0B = NULL, 
                               colMap = "SNR", shapeMap = "r",
                               xylim = NULL, log_xy = FALSE, rotate_x_tickLabels = TRUE,
                               printPlot = TRUE,
                               savePlot = FALSE, plotWidth = 7, plotHeight = 7,
                               plotUnits = "in",
                               plotPath = ".", plotFileName = paste0(metric, "_methods_vs_methods.pdf")) {
  #####
  # metric: one of "pAUC", "RMSE", "TPR", "PPV" 
  # meths: subset of c("Lasso", "AdaLasso", "LENet", "HENet", "Ridge", "SCAD", "Stability", "Dantzig")
  #        or numeric vector giving corresponding indexes
  # corDesign: subset of "Independence", "Pairwise" and "Toeplitz" (for synthetic data), and,
  #            "low" and "high" (for semi-synthetic data), or,
  #            can be set to NULL to get all values
  # s0: subset of c(10, 20, 40) or set to NULL for all values
  # SNR: subset of c(0.5, 1, 2, 4) or set to NULL for all values
  # n: subset of c(100, 200 ,300) or set to Null for all values
  # p: subset of c(500, 1000, 2000, 4000) or set to NULL for all values
  # rho: subset of c(0.5, 0.7, 0.9) or set to NULL for all values; applies for Pairwise corDesign only (Toeplitz has rho fixed to 0.7)
  # pB: subset of c(10, 100) or set to NULL for all values; applies for Pairwise corDesign only (Toeplitz and high have pB fixed to 100 and 10 respectively)
  # s0B: subset of c(1, 2, 5) or set to NULL for all values; applies for Pairwise and high corDesigns only (Toeplitz has s0B fixed to 2)
  # colMap: choose one of the factors (or rescaled sample size r) to map to colour aesethetic (Default is SNR)
  # shapeMap: choose one of the factors (or rescaled sample sieze r) to map to shape aesethetic (Default is r)
  # xylim: numeric vector of length 2 giving axes limits, or set to Null to use automatic limits
  # log_xy: logical; log-scale for axes?
  # rotate_x_tickLabels: logical: rotate x-axis tick labels by 90 degrees?
  # printPlot: logical; print plot?
  # savePlot: logical: save plot to pdf?
  # plotWidth, plotHeight, plotUnits: dimensions for pdf
  # plotPath, plotFileName: path and filename to save pdf to; defaults to working directory with filename <metric>_methods_vs_methods.pdf
  # 
  
  #####
  if(length(metric) != 1) {
    stop("metric should be character vector of length = 1")
  } else if(!(tolower(metric) %in% c("pauc", "rmse", "tpr", "ppv"))) {
    stop("metric should be one of 'pAUC', 'RMSE', 'TPR' or 'PPV'")
  }
  
  methsAll <- c("Lasso", "AdaLasso", "LENet", "HENet", "Ridge", "SCAD", "Stability", "Dantzig")
  if((is.character(meths) && !all(tolower(meths) %in% tolower(methsAll))) ||
            (is.numeric(meths) && !all(meths %in% 1:8))) {
    stop("meths should be a subset of c('Lasso', 'AdaLasso', 'LENet', 'HENet', 'Ridge', 'SCAD', 'Stability', 'Dantzig')
          or a subset of 1:8")
  } 
  
  if(is.numeric(meths)) {
    meths <- methsAll[meths]
  }
  

  xylabel <- metric
  metric <- tolower(metric)
  
  load("results/synthetic/allSettingsAveraged.RData")
  dfSyn <- scores[[metric]]
  load("results/semi-synthetic/allSettingsAveraged.RData")
  dfReal <- scores[[metric]]
  
  df <- bind_rows(dfSyn, dfReal)
  
  list2env(setNullParams(df, corDesign = corDesign, rho = rho, pB = pB, s0B = s0B,
                n = n, p = p, s0 = s0, SNR = SNR), envir = environment())
  
 
  if(any(c("Independence", "low", "Toeplitz") %in% corDesign)) {
    plotData <- df[df$SNR %in% SNR & df$s0 %in% s0 & df$n %in% n & df$p %in% p &
                     df$corDesign %in% intersect(corDesign, c("Independence", "low", "Toeplitz")) &
                     df$method %in% tolower(meths), ]
  } else {
    plotData <- df[NULL, ]
  }
  if("high" %in% corDesign) {
    plotData <- rbind(plotData,
                      df[df$SNR %in% SNR & df$s0 %in% s0 & df$n %in% n & df$p %in% p &
                           df$corDesign == "high" &
                           df$s0B %in% s0B & 
                           df$method %in% tolower(meths), ])
  }
  if("Pairwise" %in% corDesign) {
    plotData <- rbind(plotData,
                      df[df$SNR %in% SNR & df$s0 %in% s0 & df$n %in% n & df$p %in% p &
                           df$corDesign == "Pairwise" &
                           df$rho %in% rho & df$pB %in% pB & df$s0B %in% s0B & 
                           df$method %in% tolower(meths), ])
  }
  
  # add in discretised r values and dummy variable to use to create facet strip
  plotData <- cbind(plotData,
                    rDisc = factor((plotData$r > 0) + (plotData$r > 1) +
                                     (plotData$r > 2) + (plotData$r > 3)))
  
  if(colMap == "r") {
    colMap <- "rDisc"
    colLabels <- getFactorLevelLabels("rDisc")
    
  } else {
    colLabels <- waiver()
  }
  if(shapeMap == "r") {
    shapeMap <- "rDisc"
    shapeLabels <- getFactorLevelLabels("rDisc")
  } else {
    shapeLabels <- waiver()
  }
  
  plotData[[colMap]] <- factor(plotData[[colMap]])
  plotData[[shapeMap]] <- factor(plotData[[shapeMap]])
  
  if(log_xy) {
    trans <- "log10"
  } else {
    trans <- "identity"
  }
  
  plotData <- plotData[, names(plotData) != "sd"]
  plotData <- spread(plotData, method, mean)
  
  lowerfun <- function(data, mapping, limits = NULL, breaks = waiver()) {
    ggplot(data = data, mapping = mapping) +
      geom_point(size = 1.2) +
      geom_abline(slope = 1, intercept = 0) +
      scale_x_continuous(limits = limits, breaks = breaks) +
      scale_y_continuous(limits = limits, breaks = breaks) +
      scale_colour_manual(name = getAxisLegendLabel(colMap), values = brewer.pal(6, "Dark2"), labels = colLabels) +
      scale_shape_discrete(name = getAxisLegendLabel(shapeMap), solid = FALSE, 
                           labels = shapeLabels) +
      guides(colour = guide_legend(title.position = "top", title.hjust = 0.5),
             shape = guide_legend(title.position = "top", title.hjust = 0.5))
  }  
  
  if(metric == "rmse") {
    limits = NULL
    breaks = waiver()
  } else {
    limits = c(0, 1)
    breaks = c(0, 0.5, 1)
  }
  
  ggp <- ggpairs(plotData, columns = intersect(tolower(meths), colnames(plotData)), 
                 lower = list(continuous = wrap(lowerfun, limits = limits,
                                                breaks = breaks)),
                 upper = list(continuous = "blank"),
                 diag = list(continuous = "blankDiag"),
                 mapping = aes_string(colour = colMap, shape = shapeMap),
                 legend = c(2,1),
                 xlab = xylabel, ylab = xylabel,
                 switch = "both",
                 columnLabels = meths) +
    ggplotThemePub() %+replace%
    theme(panel.grid.minor = theme_bw()$panel.grid.minor) +
    theme(legend.position = "bottom")
  
  if(rotate_x_tickLabels) {
    ggp <- ggp + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  }
  
  gpairs_lower <- function(g){
    g$plots <- g$plots[-(1:g$nrow)]
    g$yAxisLabels <- g$yAxisLabels[-1]
    g$nrow <- g$nrow -1
    
    g$plots <- g$plots[-(seq(g$ncol, length(g$plots), by = g$ncol))]
    g$xAxisLabels <- g$xAxisLabels[-g$ncol]
    g$ncol <- g$ncol - 1
    
    g
  }
  
  ggp <- gpairs_lower(ggp)
  
  if(printPlot) {
    print(ggp)
  }
  
  if(savePlot) {
    ggsave(filename = plotFileName, path = plotPath, 
           plot = ggp, width = plotWidth, height = plotHeight, units = plotUnits)
  }  
  
  invisible(ggp)
}




#############################
# lineplots where up to five factors are varied (five of n, p, s0, SNR, rho, pB, s0B) (used for Figs 2, 3, 5, 6, 8, 9)
#####
lineplot_separate_factors <- function(metric, meths,
                                      corDesign, 
                                      facetGridFormula,
                                      s0 = c(10, 40), SNR = c(4, 1), 
                                      n = c(300, 200, 100), p = NULL, 
                                      rho = 0, pB = 0, s0B = 0, 
                                      x = "p",
                                      includeIndependence = FALSE,
                                      ylim = NULL, rotate_x_tickLabels = FALSE,
                                      scales = "fixed", 
                                      errorbars = FALSE, printPlot = TRUE,
                                      savePlot = FALSE, plotWidth = 7, plotHeight = 7,
                                      plotUnits = "in",
                                      plotPath = ".", 
                                      plotFileName = paste0(metric, "_lineplot_", corDesign, "_", 
                                                            paste0(as.character.default(facetGridFormula)[c(2,1,3)], collapse = " "),
                                                            ".pdf")) {
  #####
  # metric: one of "pAUC", "RMSE", "TPR", "PPV" or c("TPR", "PPV") 
  # meths: subset of c("Lasso", "AdaLasso", "LENet", "HENet", "Ridge", "SCAD", "Stability", "Dantzig")
  #        or numeric vector giving corresponding indexes
  # corDesign: one of "Independence", "Pairwise" or "Toeplitz" (for synthetic data), or,
  #            "low" and "high" (for semi-synthetic data)
  # facetGridFormula: formula for ggplot2 facet_grid function. 
  #                   Specifies which factors are plotted as rows and columns.
  #                   Allows one factor for rows and two factors for columns. e.g. SNR ~ s0 + n
  # s0: subset of c(10, 20, 40) or set to NULL for all values
  # SNR: subset of c(0.5, 1, 2, 4) or set to NULL for all values
  # n: subset of c(100, 200 ,300) or set to Null for all values
  # p: subset of c(500, 1000, 2000, 4000) or set to NULL for all values
  # rho: subset of c(0.5, 0.7, 0.9) or set to NULL for all values; applies for Pairwise corDesign only (Toeplitz has rho fixed to 0.7)
  # pB: subset of c(10, 100) or set to NULL for all values; applies for Pairwise corDesign only (Toeplitz and high have pB fixed to 100 and 10 respectively)
  # s0B: subset of c(1, 2, 5) or set to NULL for all values; applies for Pairwise and high corDesigns only (Toeplitz has s0B fixed to 2)
  # x: factor to plot on x-axis
  # includeIndependence: when (i) corDesign != "Independence"; (ii) x is one of "s0B", "rho" or "pB";
  #                      and (iii) rows/columns in facetGridFormula are a subset of c("n", "p", "s0", "SNR"),
  #                      include "Independence" results in each panel
  # ylim: y-axis limits. Set automatically by default.
  # rotate_x_tickLabels: logical: rotate x-axis tick labels by 45 degrees?
  # scales: argument for ggplot2 facet_grid function. 
  #         "fixed" sets all panels to have the same y-axis limits. 
  #         "free_y" allows each row in grid to have its own scale.
  # errorbars: logical: include errorbars?
  # printPlot: logical; print plot?
  # savePlot: logical: save plot to pdf?
  # plotWidth, plotHeight, plotUnits: dimensions for pdf
  # plotPath, plotFileName: path and filename to save pdf to; 
  #                         defaults to working directory with filename <metric>_lineplot_<corDesign>_<facetGridFormula>.pdf
  #####
  if(length(metric) == 2) {
    gt1 <- lineplot_separate_factors(metric = metric[1], meths = meths,
                                     corDesign = corDesign,
                                   facetGridFormula = facetGridFormula,
                                   s0 = s0, SNR = SNR, n = n, p = p, 
                                   s0B = s0B, pB = pB, rho = rho,
                                   x = x,
                                   includeIndependence = includeIndependence,
                                   ylim = ylim, rotate_x_tickLabels = rotate_x_tickLabels,
                                   scales = scales, 
                                   errorbars = errorbars, printPlot = FALSE, 
                                   savePlot = FALSE)
    
    gt2 <- lineplot_separate_factors(metric = metric[2], meths = meths, 
                                     corDesign = corDesign,
                                   facetGridFormula = facetGridFormula,
                                   s0 = s0, SNR = SNR, n = n, p = p, 
                                   s0B = s0B, pB = pB, rho = rho,
                                   x = x,
                                   includeIndependence = includeIndependence,
                                   ylim = ylim, rotate_x_tickLabels = rotate_x_tickLabels,
                                   scales = scales, 
                                   errorbars = errorbars, printPlot = FALSE, 
                                   savePlot = FALSE)
    
    gt <- lineplot_combine_2metrics(gt1, gt2)
    
  } else {
  
    rows <- all.vars(facetGridFormula[[2]])
    cols <- all.vars(facetGridFormula[[3]])
    
    if(!all(rows %in% c("n", "p", "s0", "SNR", "rho", "pB", "s0B"))) {
      stop("facetGridFormula elements should be a subset of c('n', 'p', 's0', 'SNR', 'rho', 'pB', 's0B')")
    }
    if(!all(cols %in% c("n", "p", "s0", "SNR", "rho", "pB", "s0B"))) {
      stop("facetGridFormula elements should be a subset of c('n', 'p', 's0', 'SNR', 'rho', 'pB', 's0B')")
    }
    if(length(cols) > 2) {
      stop("Up to two factors can be varied by column in the facet_grid")
    }
    if(length(rows) > 2) {
      stop("Up to two factors can be varied by row in the facet_grid")
    }
    if(length(intersect(cols, rows)) > 1) {
      stop("RHS and LHS of facetGridFormula should have no common elements")
    }
    if(corDesign == "Independence" && includeIndependence) {
      warning("Ignoring argument: includeIndependence = TRUE")
    }
    if(includeIndependence && x %in% c("p", "n", "SNR", "s0")) {
      stop("includeIndependence can only be TRUE when the x-axis variable is s0B, pB or rho")
    }
    
    methsAll <- c("Lasso", "AdaLasso", "LENet", "HENet", "Ridge", "SCAD", "Stability", "Dantzig")
    if((is.character(meths) && !all(tolower(meths) %in% tolower(methsAll))) ||
       (is.numeric(meths) && !all(meths %in% 1:8))) {
      stop("meths should be a subset of c('Lasso', 'AdaLasso', 'LENet', 'HENet', 'Ridge', 'SCAD', 'Stability', 'Dantzig')
          or a subset of 1:8")
    } 
    
    if(is.numeric(meths)) {
      meths <- methsAll[meths]
    }
    
    colours <- (brewer.pal(9, "Set1"))
    colours <- c(colours[c(1, 2, 9, 3)], "#EEEE00", "#000000", colours[c(4, 7)])
    names(colours) <- tolower(methsAll)
    colours <- colours[methsAll %in% meths]
    
    ylabel <- metric
    metric <- tolower(metric)
    
    xlabel <- getAxisLegendLabel(x)
    
    if(corDesign %in% c("Independence", "Pairwise", "Toeplitz")) {
      load("results/synthetic/allSettingsAveraged.RData")
      df <- scores[[metric]]
      
      list2env(setNullParams(df, corDesign = corDesign, rho = rho, pB = pB, s0B = s0B,
                             n = n, p = p, s0 = s0, SNR = SNR), envir = environment())
      
      if(corDesign %in% "Pairwise") {
        fixedVar <- setdiff(c("n", "p", "s0", "SNR", "rho", "pB", "s0B"), c(rows, cols, x))
        
        plotData <- df[df$SNR %in% SNR & df$s0 %in% s0 & df$n %in% n & df$p %in% p &
                         df$corDesign == corDesign &
                         df$rho %in% rho & df$pB %in% pB & df$s0B %in% s0B & 
                         df$method %in% tolower(meths), ]
      } else {
        fixedVar <- setdiff(c("n", "p", "s0", "SNR"), c(rows, cols, x))
        
        plotData <- df[df$SNR %in% SNR & df$s0 %in% s0 & df$n %in% n & df$p %in% p &
                         df$corDesign == corDesign &
                         df$method %in% tolower(meths), ]
      }
      
    } else {
      load("results/semi-synthetic/allSettingsAveraged.RData")
      df <- scores[[metric]]
      
      list2env(setNullParams(df, corDesign = corDesign, rho = rho, pB = pB, s0B = s0B,
                             n = n, p = p, s0 = s0, SNR = SNR), envir = environment())
      
      if(corDesign %in% "low") {
        fixedVar <- setdiff(c("n", "p", "s0", "SNR"), c(rows, cols, x))
        
        plotData <- df[df$SNR %in% SNR & df$s0 %in% s0 & df$n %in% n & df$p %in% p &
                         df$corDesign == corDesign &
                         df$method %in% tolower(meths), ]
      } else {
        fixedVar <- setdiff(c("n", "p", "s0", "SNR", "s0B"), c(rows, cols, x))
        
        plotData <- df[df$SNR %in% SNR & df$s0 %in% s0 & df$n %in% n & df$p %in% p &
                         df$corDesign == corDesign &
                         df$s0B %in% s0B &
                         df$method %in% tolower(meths), ]
      }
    }
    plotData$method <- factor(plotData$method, tolower(meths))
    
    if(includeIndependence) {
      load("results/synthetic/allSettingsAveraged.RData")
      dfInd <- scores[[metric]]
      dfInd <- dfInd[dfInd$corDesign == "Independence", ]
      plotDataInd <- dfInd[dfInd$n %in% n & dfInd$p %in% p & 
                             dfInd$SNR %in% SNR & dfInd$s0 %in% s0 &
                             dfInd$corDesign == "Independence" &
                             dfInd$method %in% tolower(meths), ]
      plotDataInd$method <- factor(plotDataInd$method, tolower(meths))
      
      if(corDesign == "Pairwise") {
        legendLabelCorDesign <- "Synthetic Correlation Design" 
      } else {
        legendLabelCorDesign <- "Semisynthetic Correlation Design" 
      }
    }
    
    for(prm in c("n", "p", "s0", "SNR", "rho", "pB", "s0B")) {
      if(x != prm && !is.null(plotData[[prm]])) {
        plotData[[prm]] <- factor(plotData[[prm]], levels = get(prm))
        levels(plotData[[prm]]) <- getFactorLevelLabels(prm, levels(plotData[[prm]]))
        
        if(includeIndependence){
          plotDataInd[[prm]] <- factor(plotDataInd[[prm]], levels = get(prm))
          levels(plotDataInd[[prm]]) <- getFactorLevelLabels(prm, levels(plotDataInd[[prm]]))
        }
      }
    }
    
    if(x %in% c("n", "rho")) {
      log_x <- FALSE
    } else {
      log_x <- TRUE
    }
    
    if(log_x) {
      if(includeIndependence) {
        indX <- get(x)[1] / 2
        plotDataInd[[x]] <- indX
      }
      
      trans <- "log10"
    } else {
      if(includeIndependence) {
        indX <- get(x)[1] - (get(x)[2] - get(x)[1])
        plotDataInd[[x]] <- indX
      }
      trans <- "identity"
    }
    
    if(metric=="rmse") {
      yBreaks <- waiver()
    } else {
      yBreaks <- c(0, 0.25, 0.5, 0.75, 1)
    }
    
    if(x == "p") {
      pd <- position_dodge(width = 0.02)
    } else {
      pd <- position_dodge(width = 0.05)
    }
    
    ggp <- ggplot(plotData, aes_string(x = x, y = "mean", colour = "method", shape = "corDesign")) +
      ggplotThemePub() +
      geom_line(size = 0.5, position = pd) +
      geom_point(size = 1, position = pd) +
      facet_grid(facetGridFormula, 
                 labeller = label_parsed,
                 scales = scales)
    
    if(includeIndependence) {
      ggp <- ggp + geom_point(data = plotDataInd,
                              size = 1, position = position_dodge(width = 0.15)) +
        geom_vline(xintercept = (get(x)[1] + indX) / 2, colour = "black", linetype = 2) +
        scale_x_continuous(trans = trans, breaks = c(indX, get(x)), 
                           labels = c("Ind", get(x))) +
        scale_shape_manual(name = "Design", values = c(1, 4), 
                           breaks = c("Independence", corDesign),
                           labels = c("Synthetic Independence Design", 
                                      legendLabelCorDesign))
      
    } else {
      ggp <- ggp +  scale_x_continuous(trans = trans, breaks = unique(plotData[[x]]))
    }
    
    ggp <- ggp + scale_y_continuous(breaks = yBreaks) +
      scale_colour_manual(name = "Method", values = colours,
                          labels = meths) +
      ylab(ylabel) +
      xlab(getAxisLegendLabel(x)) +
      coord_cartesian(ylim = ylim) +
      theme(legend.title = element_blank(),
            legend.box = "vertical")
    
    if(includeIndependence) {
      ggp <- ggp + guides(colour = guide_legend(nrow = 1, order = 2),
                          shape = guide_legend(order = 1))
    } else {
      ggp <- ggp + guides(colour = guide_legend(nrow = 1),
                          shape = FALSE)
    }
    
    if(rotate_x_tickLabels) {
      ggp <- ggp + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1))
    }
    if(errorbars) {
      ggp <- ggp + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                                 size = 0.5, position = pd)
      if(includeIndependence) {
        ggp <- ggp + geom_errorbar(data = plotDataInd,
                                   aes(ymin = mean - sd, ymax = mean + sd),
                                   size = 0.5, position = position_dodge(width = 0.15),
                                   width = 1/4) 
      }
    }
    # browser()
    
    gt <- ggplotGrob(ggp)
    
    # top facet strips span two columns
    # code below adapted from https://stackoverflow.com/questions/40316169/nested-facets-in-ggplot2-spanning-groups
    stripLocations <- grep("strip-t", gt$layout$name)
    strip <- gtable::gtable_filter(gt, "strip-t", trim = FALSE)
    nStripCol <- nrow(strip$layout)
    nCol1 <- length(get(cols[1]))
    nCol2 <- length(get(cols[2]))
    
    top <- strip$layout$t[1]
    l   <- strip$layout$l[nStripCol/nCol1 * (1:nCol1) - nCol2 +1]
    r   <- strip$layout$r[nStripCol/nCol1 * (1:nCol1)]
    
    mat   <- matrix(vector("list", length = 2 * (2 * nCol2 - 1)), nrow = 2)
    mat[] <- list(zeroGrob())
    multiColStrip <- gtable::gtable_matrix("strip-t", mat, widths = unit(c(1, rep(c(0, 1), nCol2 - 1)), "null"), 
                                           heights = unit(rep(1, 2), "null"))
    
    for(i in 1:nCol1) {
      combinedStrip <- gtable::gtable_add_grob(multiColStrip, 
                                               gt$grobs[[stripLocations[nCol2 * (i-1) + 1]]]$grobs[[1]], 
                                               t = 1, l = 1, r = ncol(multiColStrip))
      gt <- gtable::gtable_add_grob(gt, combinedStrip, 
                                    t = top,  l = l[i], r = r[i], 
                                    name = paste0("strip-t-", nCol2 * (i-1) + 1, "-", nCol2 * i))
    }
    
    # larger gap between top level facets
    gt$widths[4 + (2 * nCol2) * 1:(nCol2 - 1)] <- 1.5 * gt$widths[4 + (2 * nCol2) * 1:(nCol2 - 1)]
    
  }
  
  if(printPlot) {
    grid.newpage()
    grid.draw(gt)
  }
  
  if(savePlot) {
    ggsave(filename = plotFileName, path = plotPath, 
           plot = gt, width = plotWidth, height = plotHeight, units = plotUnits)
  }  
  
  invisible(gt)
}