violinPlot <- function(ctGenes, byFactor, factorOrder, groupLabel, 
                       extraLabel, dotSize, dotAlpha = 0.5){
  
  ## --- normalize grouping column (patient / probe / etc.) ---
  grp <- trimws(as.character(ctGenes[, byFactor]))
  grp[is.na(grp) | grp == ""] <- "Unknown"
  ctGenes[, byFactor] <- factor(grp)
  
  ## factorOrder is mainly for kmeans.cluster / infl-uninfl-blood logic
  if (missing(factorOrder) || is.null(factorOrder)) {
    factorOrder <- levels(ctGenes[, byFactor])
  } else {
    factorOrder <- unique(as.character(factorOrder))
  }
  
  if (nrow(ctGenes) <= 4) return(invisible(NULL))
  
  ## account for kmeans.cluster column when in use
  if (byFactor != "kmeans.cluster") {
    idCols <- grep("kmeans.cluster", colnames(ctGenes))
  } else {
    idCols <- grep("kmeans.cluster", colnames(ctGenes))
    ctGenes <- subset(ctGenes, kmeans.cluster %in% factorOrder)
  }
  
  write.csv(ctGenes, paste0(baseDir, "results/ctGenes_Print1.csv"))
  
  ## -------- get p-values for differential expression --------
  if (length(idCols) == 0) {
    stop("kmeans.cluster column not found in ctGenes.")
  }
  
  ## patientColor may or may not exist (only if you made patient TSNE earlier)
  EndColumn <- grep("^patientColor$", colnames(ctGenes))
  if (length(EndColumn) == 0) {
    ## no patientColor -> treat all columns after kmeans.cluster as genes
    EndColumn <- ncol(ctGenes) + 1
  }
  
  ## gene columns = everything between kmeans.cluster and patientColor (or end)
  geneCols <- (idCols + 1):(EndColumn - 1)
  if (length(geneCols) == 0) {
    stop("No gene columns found between kmeans.cluster and patientColor / end of ctGenes.")
  }
  
  ## force numeric expression values
  for (g in geneCols) {
    ctGenes[, g] <- suppressWarnings(as.numeric(ctGenes[, g]))
  }
  
  ## --- SAFE Kruskal–Wallis loop (skip genes with only 1 group) ---
  groups <- ctGenes[, byFactor]
  if (!is.factor(groups)) groups <- factor(groups)
  
  pvals <- sapply(geneCols, function(col) {
    vals <- ctGenes[, col]
    ok   <- !is.na(vals) & !is.na(groups)
    v    <- vals[ok]
    g    <- droplevels(groups[ok])
    
    ## if all obs in one group (or no data), return NA
    if (length(v) == 0 || length(unique(g)) < 2) {
      return(NA_real_)
    }
    
    suppressWarnings(
      tryCatch(kruskal.test(v, g)$p.value,
               error = function(e) NA_real_)
    )
  })
  
  ## name pvals by gene names & BH adjust
  names(pvals) <- colnames(ctGenes)[geneCols]
  pvals <- p.adjust(pvals, method = "BH")
  
  ## order genes by significance, dropping all-NA genes
  ord <- order(pvals, na.last = NA)
  orderedGenes <- geneCols[ord]
  
  if (byFactor == "IBD_Status") {
    AddedColumn <- grep("^IBD_Status$", colnames(ctGenes))
    ctOrderGenes <- ctGenes[, c(1:idCols, AddedColumn, orderedGenes)]
  } else {
    ctOrderGenes <- ctGenes[, c(1:idCols, orderedGenes)]
  }
  
  sigGeneNames <- names(pvals)[ord]
  sigGeneVals  <- signif(pvals[ord], digits = 4)
  sigLength    <- sum(sigGeneVals < 0.05 & !is.na(sigGeneVals))
  
  cat(paste("Differentially expressed genes between",
            groupLabel, extraLabel, ":\n"))
  if (sigLength > 0) {
    print(paste(sigGeneNames[1:sigLength],
                sigGeneVals[1:sigLength], sep=": "),
          quote = FALSE)
  }
  
  ## -------- melt for plotting --------
  if (byFactor != "IBD_Status") {
    ctMelt <- reshape2::melt(ctOrderGenes,
                             id.vars = 1:idCols,
                             variable.name = "gene")
  } else {
    AddedColumn <- grep("^IBD_Status$", colnames(ctOrderGenes))
    ctMelt <- reshape2::melt(ctOrderGenes,
                             id.vars = c(1:AddedColumn),
                             variable.name = "gene")
  }
  names(ctMelt)[names(ctMelt) == "value"] <- "log2Ex"
  ctMelt$log2Ex <- suppressWarnings(as.numeric(ctMelt$log2Ex))
  
  genesToPlotGG <- list()
  ctMeltSubGG   <- list()
  genesGG       <- list()
  fillLabel     <- list()
  
  for (i in 1:24) {
    genesToPlotGG[[i]] <- sigGeneNames[(i*6 - 5):(i*6)]
    genesToPlotGG[[i]] <- genesToPlotGG[[i]][!is.na(genesToPlotGG[[i]])]
    
    if (length(genesToPlotGG[[i]]) == 0) next
    
    dat <- ctMelt[ctMelt$gene %in% genesToPlotGG[[i]], ]
    if (nrow(dat) == 0) next
    
    genesGG[[i]]   <- factor(dat$gene)
    fillLabel[[i]] <- droplevels(factor(dat[, byFactor]))
    
    ctMeltSubGG[[i]] <- dat
  }
  
  ## -------- plotting --------
  if (identical(factorOrder, c("infl","uninfl","blood"))) {
    
    for (i in seq(1, 24, 4)) {
      if (is.null(ctMeltSubGG[[i]])) next
      
      vplot1 <- ggplot(ctMeltSubGG[[i]],
                       aes(genesGG[[i]], log2Ex, fill = fillLabel[[i]])) +
        geom_violin(scale = "width",
                    position = position_dodge(width = 0.75),
                    bw = "nrd0",
                    adjust = 0.8) +
        geom_point(size = dotSize,
                   alpha = dotAlpha,
                   shape = 19,
                   position = position_jitterdodge(jitter.width = 0.3,
                                                   jitter.height = 0,
                                                   dodge.width = 0.75)) +
        scale_y_continuous() +
        guides(fill = guide_legend(title = groupLabel)) +
        ggtitle(paste("most significant expression differences between ",
                      groupLabel, " ", extraLabel, " (plot #", i, ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text = element_text(size = 25),
              panel.grid.minor = element_line(color = "gray90"),
              panel.grid.major = element_line(color = "gray75", size = 0.4),
              panel.grid.major.x = element_blank(),
              axis.title.x = element_text(vjust = -0.5),
              plot.margin = unit(c(1,0,1,0),"cm")) +
        scale_fill_manual(values = c("blood"  = "#f25e3b",
                                     "infl"   = "#3939FF",
                                     "uninfl" = "deepskyblue2"))
      
      print(vplot1)
    }
    
  } else {
    
    for (i in seq(1, 24, 4)) {
      if (is.null(ctMeltSubGG[[i]])) next
      
      vplot1 <- ggplot(ctMeltSubGG[[i]],
                       aes(genesGG[[i]], log2Ex, fill = fillLabel[[i]])) +
        geom_violin(scale = "width",
                    position = position_dodge(width = 0.75),
                    bw = "nrd0",
                    adjust = 0.8) +
        geom_point(size = dotSize,
                   alpha = dotAlpha,
                   shape = 19,
                   position = position_jitterdodge(jitter.width = 0.3,
                                                   jitter.height = 0,
                                                   dodge.width = 0.75)) +
        scale_fill_brewer(palette = "Set1") +
        scale_y_continuous() +
        guides(fill = guide_legend(title = groupLabel)) +
        ggtitle(paste("most significant expression differences between ",
                      groupLabel, " ", extraLabel, " (plot #", i, ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text = element_text(size = 25),
              panel.grid.minor = element_line(color = "gray90"),
              panel.grid.major = element_line(color = "gray75", size = 0.4),
              panel.grid.major.x = element_blank(),
              axis.title.x = element_text(vjust = -0.5),
              plot.margin = unit(c(1,0,1,0),"cm"))
      
      print(vplot1)
    }
  }
}
