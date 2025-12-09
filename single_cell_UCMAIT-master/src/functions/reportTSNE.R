##### function to generate variable number tsne plots for cluster reports ####
reportTSNE <- function(ctPlot, tsne_y, deGenes, clusters = clusters, cluster = cluster){
  
  pointSize <- 2.9
  myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  
  if(length(clusters) > 8){
    clusterPalette <- "Set3"
  }else{
    clusterPalette <- "Accent"
  }
  
  shapeVals <- c(19, 17, 15, 18)
  numTissues <- length(unique(tsne_y$tissue))
  shapeVals <- shapeVals[1:numTissues]
  
  ## plot tsne colored by cluster
  tsne <- ggplot(tsne_y, aes(y1, y2)) +
    geom_point(aes(color=kmeans.cluster, shape=tissue), size = pointSize, alpha = 0.8) +
    scale_colour_brewer(palette = clusterPalette) +
    scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                       minor_breaks = NULL) +
    scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                       minor_breaks = NULL) +
    scale_shape_manual(values = shapeVals) +
    guides(color=guide_legend(title="cluster", order = 1)) +
    guides(shape = guide_legend(order = 2)) +
    ggtitle("t-SNE between tissues (colored by kmeans.cluster)") +
    theme_minimal() +
    theme(#axis.line=element_blank(),
          panel.border=element_rect(fill=NA, color="gray75", size=0.4),
          panel.grid.minor=element_line(color="gray90"),
          panel.grid.major=element_line(color="gray85", size=0.3),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          text=element_text(size=28),
          plot.margin=unit(c(14,9,14,9),"cm"))
  print(tsne)
  
  ##create gene tsne ggplot function
  ggtsne <- function(tsne_y, gene, tissue, pointSize, myPalette){
    scaleDist <- tsne_y[,gene] - min(tsne_y[,gene])
    scalePerc <- as.vector(quantile(scaleDist, probs = seq(0, 1, by= 0.01)))
    scalePercTest <- c(scalePerc[-1], scalePerc[length(scalePerc)])
    scalePercDiff <- scalePercTest - scalePerc
    scalePercDiff <- scalePercDiff[-c(1:20,82:101)]
    scalePercDiffMax <- which(scalePercDiff == max(scalePercDiff)) + 20
    
    colorBreakLow <- scalePerc[scalePercDiffMax - 10]/max(scaleDist)
    colorBreakHi <- scalePerc[scalePercDiffMax + 10]/max(scaleDist)
    
    colorValsLow <- seq(0, colorBreakLow, length.out = 30)
    colorValsMid <- seq(colorBreakLow, colorBreakHi, length.out = 50)
    colorValsMid <- colorValsMid[-1]
    colorValsHi <- seq(colorBreakHi, 1, length.out=30)
    colorValsHi <- colorValsHi[-1]
    colorVals <- c(colorValsLow, colorValsMid, colorValsHi)
    
    tsnePlot <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=tsne_y[,gene], shape=tissue), size = pointSize, alpha = 0.8) +
      scale_colour_gradientn(colours = myPalette(108), values=colorVals) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_colorbar(title=gene, order = 1)) +
      guides(shape = guide_legend(order = 2)) +
      ggtitle(paste("cluster ", cluster, ": ", gene, sep = "")) +
      theme_minimal() +
      theme(#axis.line=element_blank(),
            panel.border=element_rect(fill=NA, color="gray75", size=0.4),
            panel.grid.minor=element_line(color="gray90"),
            panel.grid.major=element_line(color="gray85", size=0.3),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            text=element_text(size=18),
            plot.margin=unit(c(2,1,2,1),"cm"))
    
    return(tsnePlot)
  }
  
  myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  
  geneNumber <- length(unlist(deGenes))
  plotNumber <- geneNumber - geneNumber%%4
  residualPlots <- geneNumber - plotNumber
  
  if(plotNumber >= 4){
    for(i in seq(1, plotNumber, by=4)){
      gene <- unlist(deGenes)[i]
      tsne1 <- ggtsne(tsne_y, gene, tissue, pointSize, myPalette)
      
      gene2 <- unlist(deGenes)[i+1]
      tsne2 <- ggtsne(tsne_y, gene2, tissue, pointSize, myPalette)
      
      gene3 <- unlist(deGenes)[i+2]
      tsne3 <- ggtsne(tsne_y, gene3, tissue, pointSize, myPalette)
      
      gene4 <- unlist(deGenes)[i + 3]
      tsne4 <- ggtsne(tsne_y, gene4, tissue, pointSize, myPalette)
      
      grid.arrange(tsne1, tsne2, tsne3, tsne4, ncol=2, nrow=2)
      
    }
  }
  
  ## create full residual plot lists
  if(residualPlots > 0){
    i <- plotNumber + 1
    gene <- unlist(deGenes)[i]
    tsne1 <- ggtsne(tsne_y, gene, tissue, pointSize, myPalette)
    
    if((i+1) <= length(unlist(deGenes))){
      gene2 <- unlist(deGenes)[i+1]
      tsne2 <- ggtsne(tsne_y, gene2, tissue, pointSize, myPalette)
    }else{
      tsne2 <- nullGrob()
    }
    
    if((i+2) <= length(unlist(deGenes))){
      gene3 <- unlist(deGenes)[i+2]
      tsne3 <- ggtsne(tsne_y, gene3, tissue, pointSize, myPalette)
    }else{
      tsne3 <- nullGrob()
    }
    
    if((i+3) <= length(unlist(deGenes))){
      gene4 <- unlist(deGenes)[i + 3]
      tsne4 <- ggtsne(tsne_y, gene4, tissue, pointSize, myPalette)
    }else{
      tsne4 <- nullGrob()
    }
    
    grid.arrange(tsne1, tsne2, tsne3, tsne4, ncol=2, nrow=2)
  }
}