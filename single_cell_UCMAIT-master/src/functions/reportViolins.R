##### function to generate variable number violin plots for cluster reports ####
reportViolins <- function(ctPlot, deGenes, clusters = clusters, cluster = cluster){
  
  ## melt dataframe for plotting
  ctPlotMelt <- melt(ctPlot, id.vars = 1:9, variable.name = "gene")
  names(ctPlotMelt)[which(names(ctPlotMelt)=="value")] <- "log2Ex"
  
  geneNumber <- length(unlist(deGenes))
  plotNumber <- geneNumber - geneNumber%%8
  residualPlots <- geneNumber - plotNumber
  
  ## create empty plotting lists
  genesToPlotGG <- list()
  ctPlotMeltSubGG <- list()
  genesGG <- list()
  fillLabel <- list()
  
  ## plot full grid arrange frames
  if(plotNumber >= 8){
    for(i in 1:(plotNumber/2)){
      genesToPlotGG[[i]] <- unlist(deGenes)[(i*2 - 1):(i*2)]
      ctPlotMeltSubGG[[i]] <- ctPlotMelt[which(ctPlotMelt$gene %in% genesToPlotGG[[i]]),]
      
      ctPlotMeltSubGG[[i]] <- ctPlotMeltSubGG[[i]][order(factor(ctPlotMeltSubGG[[i]][, "kmeans.cluster"],
                                                              levels=clusters)),]
      
      genesGG[[i]] <- factor(ctPlotMeltSubGG[[i]]$gene, levels=genesToPlotGG[[i]])
      fillLabel[[i]] <- factor(ctPlotMeltSubGG[[i]][, "kmeans.cluster"], levels=clusters)
    }
    
    for(i in seq(1, (plotNumber/2), 4)){
      vplot1 = ggplot(ctPlotMeltSubGG[[i]], aes(genesGG[[i]], log2Ex, fill=fillLabel[[i]])) +
        geom_violin(scale = "count", 
                    position=position_dodge(width = 0.75), 
                    # adjust=3) +
                    bw = "SJ",
                    adjust=2) +
        geom_point(size = 1.5, 
                   alpha = 0.3, 
                   shape = 19, 
                   # stroke=1.35,   # stroke is for shape=1 circles
                   show.legend = F,
                   position=position_jitterdodge(jitter.width = 0.3, 
                                                 jitter.height = 0, 
                                                 dodge.width = 0.75)) +
        scale_fill_brewer(palette = "Set3") +
        scale_y_continuous() +
        guides(fill=guide_legend(title="cluster")) +
        ggtitle(paste("most significant expression differences for cluster ",
                      cluster, " (plot #", i, ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text=element_text(size=25),
              panel.grid.minor=element_line(color="gray90"),
              panel.grid.major=element_line(color="gray75", size=0.4),
              panel.grid.major.x=element_blank(),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(1,0,1,0),"cm"))
      
      vplot2 = ggplot(ctPlotMeltSubGG[[i+1]], aes(genesGG[[i+1]], log2Ex, fill=fillLabel[[i+1]])) +
        geom_violin(scale = "count", 
                    position=position_dodge(width = 0.75), 
                    # adjust=3) +
                    bw = "SJ",
                    adjust=2) +
        geom_point(size = 1.5, 
                   alpha = 0.3, 
                   shape = 19, 
                   # stroke=1.35,   # stroke is for shape=1 circles
                   show.legend = F,
                   position=position_jitterdodge(jitter.width = 0.3, 
                                                 jitter.height = 0, 
                                                 dodge.width = 0.75)) +
        scale_fill_brewer(palette = "Set3") +
        scale_y_continuous() +
        guides(fill=guide_legend(title="cluster")) +
        ggtitle(paste("most significant expression differences for cluster ",
                      cluster, " (plot #", (i+1), ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text=element_text(size=25),
              panel.grid.minor=element_line(color="gray90"),
              panel.grid.major=element_line(color="gray75", size=0.4),
              panel.grid.major.x=element_blank(),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(1,0,1,0),"cm"))
      
      vplot3 = ggplot(ctPlotMeltSubGG[[i+2]], aes(genesGG[[i+2]], log2Ex, fill=fillLabel[[i+2]])) +
        geom_violin(scale = "count", 
                    position=position_dodge(width = 0.75), 
                    # adjust=3) +
                    bw = "SJ",
                    adjust=2) +
        geom_point(size = 1.5, 
                   alpha = 0.3, 
                   shape = 19, 
                   # stroke=1.35,   # stroke is for shape=1 circles
                   show.legend = F,
                   position=position_jitterdodge(jitter.width = 0.3, 
                                                 jitter.height = 0, 
                                                 dodge.width = 0.75)) +
        scale_fill_brewer(palette = "Set3") +
        scale_y_continuous() +
        guides(fill=guide_legend(title="cluster")) +
        ggtitle(paste("most significant expression differences for cluster ",
                      cluster, " (plot #", (i+2), ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text=element_text(size=25),
              panel.grid.minor=element_line(color="gray90"),
              panel.grid.major=element_line(color="gray75", size=0.4),
              panel.grid.major.x=element_blank(),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(1,0,1,0),"cm"))
      
      vplot4 = ggplot(ctPlotMeltSubGG[[i+3]], aes(genesGG[[i+3]], log2Ex, fill=fillLabel[[i+3]])) +
        geom_violin(scale = "count", 
                    position=position_dodge(width = 0.75), 
                    # adjust=3) +
                    bw = "SJ",
                    adjust=2) +
        geom_point(size = 1.5, 
                   alpha = 0.3, 
                   shape = 19, 
                   # stroke=1.35,   # stroke is for shape=1 circles
                   show.legend = F,
                   position=position_jitterdodge(jitter.width = 0.3, 
                                                 jitter.height = 0, 
                                                 dodge.width = 0.75)) +
        scale_fill_brewer(palette = "Set3") +
        scale_y_continuous() +
        guides(fill=guide_legend(title="cluster")) +
        ggtitle(paste("most significant expression differences for cluster ",
                      cluster, " (plot #", (i+3), ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text=element_text(size=25),
              panel.grid.minor=element_line(color="gray90"),
              panel.grid.major=element_line(color="gray75", size=0.4),
              panel.grid.major.x=element_blank(),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(1,0,1,0),"cm"))
      
      grid.arrange(vplot1, vplot2, vplot3, vplot4, ncol=1, nrow=4)
    }
  }
  ## create full residual plot lists
  blankGeneNum <- 0
  
  if(residualPlots > 0){
    ## set blank plot variables
    ngNum <- floor((8-residualPlots)/2)
    blankGeneNum <- (2 - residualPlots)%%2
    
    if(residualPlots>=2){
      for(i in ((plotNumber/2) + 1):((plotNumber/2) + floor(residualPlots/2))){
        genesToPlotGG[[i]] <- unlist(deGenes)[(i*2 - 1):(i*2)]
        ctPlotMeltSubGG[[i]] <- ctPlotMelt[which(ctPlotMelt$gene %in% genesToPlotGG[[i]]),]
        
        ctPlotMeltSubGG[[i]] <- ctPlotMeltSubGG[[i]][order(factor(ctPlotMeltSubGG[[i]][, "kmeans.cluster"],
                                                                  levels=clusters)),]
        
        genesGG[[i]] <- factor(ctPlotMeltSubGG[[i]]$gene, levels=genesToPlotGG[[i]])
        fillLabel[[i]] <- factor(ctPlotMeltSubGG[[i]][, "kmeans.cluster"], levels=clusters)
      }
    }
    
    ## create partial residual plot entry
    if(blankGeneNum > 0){
      i <- ((plotNumber/2) + floor(residualPlots/2)) + 1

      genesToPlotGG[[i]] <- unlist(deGenes)[(i*2 - 1):(i*2 - blankGeneNum)]
      ctPlotMeltSubGG[[i]] <- ctPlotMelt[which(ctPlotMelt$gene %in% genesToPlotGG[[i]]),]
      
      ctPlotMeltSubGG[[i]] <- ctPlotMeltSubGG[[i]][order(factor(ctPlotMeltSubGG[[i]][, "kmeans.cluster"],
                                                                levels=clusters)),]
      
      genesGG[[i]] <- factor(ctPlotMeltSubGG[[i]]$gene, levels=genesToPlotGG[[i]])
      fillLabel[[i]] <- factor(ctPlotMeltSubGG[[i]][, "kmeans.cluster"], levels=clusters)
    }
  
    i <- (plotNumber/2) + 1
    
    rightMargin <- 0
    rightMarginAdjust <- 27
    if(length(genesToPlotGG[[i]]) < 2){
      rightMargin <- rightMarginAdjust
    }
      
    vplot1 <- ggplot(ctPlotMeltSubGG[[i]], 
                    aes(genesGG[[i]], log2Ex, fill=fillLabel[[i]])) +
                geom_violin(scale = "count", 
                            position=position_dodge(width = 0.75), 
                            # adjust=3) +
                            bw = "SJ",
                            adjust=2) +
                geom_point(size = 1.5, 
                           alpha = 0.3, 
                           shape = 19, 
                           # stroke=1.35,   # stroke is for shape=1 circles
                           show.legend = F,
                           position=position_jitterdodge(jitter.width = 0.3, 
                                                         jitter.height = 0, 
                                                         dodge.width = 0.75)) +
                scale_fill_brewer(palette = "Set3") +
                scale_y_continuous() +
                guides(fill=guide_legend(title="cluster")) +
                ggtitle(paste("most significant expression differences for cluster ",
                              cluster, " (plot #", i, ")", sep="")) +
                xlab("\ngenes") +
                ylab("log2(gene exp. / Gapdh exp.)\n") +
                theme_minimal() +
                theme(text=element_text(size=25),
                      panel.grid.minor=element_line(color="gray90"),
                      panel.grid.major=element_line(color="gray75", size=0.4),
                      panel.grid.major.x=element_blank(),
                      axis.title.x=element_text(vjust=-0.5),
                      plot.margin=unit(c(1,rightMargin,1,0),"cm"))
    
    if((i+1) <= length(genesToPlotGG)){
      if(length(genesToPlotGG[[i+1]]) < 2){
        rightMargin <- rightMarginAdjust
      }
      
      vplot2 <- ggplot(ctPlotMeltSubGG[[i+1]], aes(genesGG[[i+1]], log2Ex, fill=fillLabel[[i+1]])) +
        geom_violin(scale = "count", 
                    position=position_dodge(width = 0.75), 
                    # adjust=3) +
                    bw = "SJ",
                    adjust=2) +
        geom_point(size = 1.5, 
                   alpha = 0.3, 
                   shape = 19, 
                   # stroke=1.35,   # stroke is for shape=1 circles
                   show.legend = F,
                   position=position_jitterdodge(jitter.width = 0.3, 
                                                 jitter.height = 0, 
                                                 dodge.width = 0.75)) +
        scale_fill_brewer(palette = "Set3") +
        scale_y_continuous() +
        guides(fill=guide_legend(title="cluster")) +
        ggtitle(paste("most significant expression differences for cluster ",
                      cluster, " (plot #", (i+1), ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text=element_text(size=25),
              panel.grid.minor=element_line(color="gray90"),
              panel.grid.major=element_line(color="gray75", size=0.4),
              panel.grid.major.x=element_blank(),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(1,rightMargin,1,0),"cm"))
    }else{
      vplot2 <- nullGrob()
    }
    
    if((i+2) <= length(genesToPlotGG)){
      if(length(genesToPlotGG[[i+2]]) < 2){
        rightMargin <- rightMarginAdjust
      }
      
      vplot3 <- ggplot(ctPlotMeltSubGG[[i+2]], aes(genesGG[[i+2]], log2Ex, fill=fillLabel[[i+2]])) +
        geom_violin(scale = "count", 
                    position=position_dodge(width = 0.75), 
                    # adjust=3) +
                    bw = "SJ",
                    adjust=2) +
        geom_point(size = 1.5, 
                   alpha = 0.3, 
                   shape = 19, 
                   # stroke=1.35,   # stroke is for shape=1 circles
                   show.legend = F,
                   position=position_jitterdodge(jitter.width = 0.3, 
                                                 jitter.height = 0, 
                                                 dodge.width = 0.75)) +
        scale_fill_brewer(palette = "Set3") +
        scale_y_continuous() +
        guides(fill=guide_legend(title="cluster")) +
        ggtitle(paste("most significant expression differences for cluster ",
                      cluster, " (plot #", (i+2), ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text=element_text(size=25),
              panel.grid.minor=element_line(color="gray90"),
              panel.grid.major=element_line(color="gray75", size=0.4),
              panel.grid.major.x=element_blank(),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(1,rightMargin,1,0),"cm"))
    }else{
      vplot3 <- nullGrob()
    }
    
    if((i+3) <= length(genesToPlotGG)){
      if(length(genesToPlotGG[[i+3]]) < 2){
        rightMargin <- rightMarginAdjust
      }
      
      vplot4 <- ggplot(ctPlotMeltSubGG[[i+3]], aes(genesGG[[i+3]], log2Ex, fill=fillLabel[[i+3]])) +
        geom_violin(scale = "count", 
                    position=position_dodge(width = 0.75), 
                    # adjust=3) +
                    bw = "SJ",
                    adjust=2) +
        geom_point(size = 1.5, 
                   alpha = 0.3, 
                   shape = 19, 
                   # stroke=1.35,   # stroke is for shape=1 circles
                   show.legend = F,
                   position=position_jitterdodge(jitter.width = 0.3, 
                                                 jitter.height = 0, 
                                                 dodge.width = 0.75)) +
        scale_fill_brewer(palette = "Set3") +
        scale_y_continuous() +
        guides(fill=guide_legend(title="cluster")) +
        ggtitle(paste("most significant expression differences for cluster ",
                      cluster, " (plot #", (i+3), ")", sep="")) +
        xlab("\ngenes") +
        ylab("log2(gene exp. / Gapdh exp.)\n") +
        theme_minimal() +
        theme(text=element_text(size=25),
              panel.grid.minor=element_line(color="gray90"),
              panel.grid.major=element_line(color="gray75", size=0.4),
              panel.grid.major.x=element_blank(),
              axis.title.x=element_text(vjust=-0.5),
              plot.margin=unit(c(1,rightMargin,1,0),"cm"))
    }else{
      vplot4 <- nullGrob()
    }
  
    grid.arrange(vplot1, vplot2, vplot3, vplot4, ncol=1, nrow=4)
  }
}



