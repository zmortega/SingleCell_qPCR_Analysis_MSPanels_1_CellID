#### generate cluster reports ####

clusterReport <- function(ctRep = ctClust, clusters = c(2,3), pValThreshold = 5e-15, plotType = "violins"){
  
  ## gene sets for grouping reported genes
  geneGroups <- list(Chemokines=c("Ccr1", "Ccr7", "Ccr2", "Ccr3", "Ccr4", "Ccr5", "Ccr6", "Cxcr3", "Cxcr4", "Cxcl10"), 
                     Cell_surface_receptors=c("Tnfrsf1a", "Tnfrsf1b", "Il.7r", "Il.5ra", "Il.4ra", "Il.2ra", "Ly6e", "Ifngr1", "Il.12rb", "Il.18r1", "Il.27r", "Icam1", "Cd44", "Ceacam1", "Il.1r2", "Cd28", "Cd3e", "Cd4", "Cd40", "Cd80", "Cd86", "Cd8a", "Ctla4", "Pdl.1", "Pd1", "Icos", "Tgfbr2"), 
                     Signaling_molecules=c("Stat1", "Stat3", "Stat4", "Stat5", "Map2k6", "Mapk8", "Fyn", "Jak1", "Jak2", "Pten", "Socs3", "Tnfaip3", "Traf2", "Vav1", "Zap70"), 
                     Cytokines=c("Ifng", "Il.10", "Il.12b", "Il.17A", "Il.2", "IL.21", "Il.25", "Il.27", "Il.3", "Il.4", "Il.5", "Il.6", "Il.7", "Tnf"), 
                     Transcription_factors=c("Bcl6", "Foxp3", "Gata4", "Ppara", "Pparg", "Ppargc1a", "Tbx21", "Bcl2", "Nfkb1", "Nur77", "Zeb2"), 
                     Metabolism=c("Gsk3a", "Gsk3b", "Hprt"), 
                     Interferon_response=c("Ifi44", "Ifi44l", "Ifit1", "Ifit3", "Irf1", "Irf2", "Irf4", "Irf7", "Isg15", "Mx1", "Aim2", "Oas1b", "Oas2", "Oasl1", "Rsad2"))
  
  ## create differential expression dataframe
  ctDE <- list()
  
  ## remove clusters that will not be reported
  ctRep <- subset(ctRep, kmeans.cluster %in% clusters)
  
  if(plotType == "tsne"){
    ## remove duplicated rows
    if(anyDuplicated(ctRep[,10:ncol(ctRep)]) > 0){
      ctRep <- ctRep[-which(duplicated(ctRep[,10:ncol(ctRep)])),]
    }
    
    ## remove genes without variance
    vars <- NULL
    for(i in 10:ncol(ctRep)){
      vars <- c(vars, var(ctRep[,i], na.rm=T))
    }
    ctGenesNoVar <- ctRep[which(vars == 0 | is.na(vars))+9]
    ctRep <- ctRep[,c(1:9, (which(!is.na(vars) & vars!=0)+9))]
    
    ## run t-SNE
    set.seed(1)
    tsne_out <- Rtsne(as.matrix(ctRep[, 10:ncol(ctRep)]), perplexity = 30)
    
    tsne_y <- as.data.frame(cbind(tsne_out$Y, 
                                  ctRep$cellSource,
                                  ctRep$kmeans.cluster,
                                  ctRep[, 10:ncol(ctRep)]))
    
    names(tsne_y)[1:4] <- c("y1", "y2", "tissue", "kmeans.cluster")
    tsne_y$kmeans.cluster <- factor(tsne_y$kmeans.cluster, levels=clusters)
    relTissues <- unique(tsne_y$tissue)
    relTissues <- relTissues[match(c("I", "LN", "S", "PB"), relTissues, nomatch=F)]
    tsne_y$tissue <- factor(tsne_y$tissue, levels=relTissues)

    for(i in c(1,2,5:ncol(tsne_y))){
      tsne_y[, i] <- as.numeric(tsne_y[, i])
    }
  }
  
  ## iterate through clusters
  for(cluster in clusters){
    
    ## create cluster specific dataframe
    ctCluster <- ctRep
    ctCluster[which(ctCluster$kmeans.cluster != cluster), "kmeans.cluster"] <- "others"
    
    #### Tissue Counts ####
    ## calculate and print tissue proportions for the current cluster
    cellSources <- unique(ctCluster$cellSource)
    cellSources <- cellSources[match(c("I", "LN", "S", "PB"), cellSources, nomatch=F)]
    
    tissueCount <- NULL
    for(i in 1:length(cellSources)){
      tissueCount[i] <- nrow(subset(ctClust, cellSource==cellSources[i] & kmeans.cluster==cluster))
    }
    totalClusterCellCount <- sum(tissueCount)
    tissueCountCluster <- tissueCount/sum(tissueCount)*100
    tissueCountCluster <- round(tissueCountCluster, digits = 1)
    tissueCountCluster <- paste(tissueCountCluster, "%", sep = "")
    
    print(paste("Cluster ", cluster, " Report:", sep = ""), row.names=F, quote = F)
    print(paste("Tissue proportions within cluster ", cluster, ":", sep=""), row.names=F, quote = F)
    print(paste(cellSources, tissueCountCluster, sep = "=", collapse = "; "), row.names=F, quote = F)
    print("     ", row.names=F, quote=F)
    
    ## calculate and print tissue proportions for the current cluster
    totalTissueCount <- NULL
    for(i in 1:length(cellSources)){
      totalTissueCount[i] <- nrow(subset(ctClust, cellSource==cellSources[i]))
    }
    tissueCountTotal <- tissueCount/totalTissueCount*100
    tissueCountTotal <- round(tissueCountTotal, digits = 1)
    tissueCountTotal <- paste(tissueCountTotal, "%", sep = "")
    
    print(paste("Tissue proportions across clusters (% of total tissue) for cluster ", cluster, ":", sep=""), row.names=F, quote = F)
    print(paste(tissueCountTotal, cellSources, sep = " of ", collapse = "; "), row.names=F, quote = F)
    print("     ", row.names=F, quote=F)
    
    #### Differential Expression ####
    ## calculate differential expression p-values
    pvals <- NULL
    for(gene in 10:ncol(ctCluster)){
      pvals <- c(pvals, kruskal.test(ctCluster[, gene], factor(ctCluster[, "kmeans.cluster"]))$p.value)
    }
    
    pvals <- p.adjust(pvals, method = "BH")
    
    ## calculate differences between means 
    ctCluster$kmeans.cluster <- factor(ctCluster$kmeans.cluster, levels = c(cluster, "others"))
    diffMean <- NULL
    for(gene in names(ctCluster)[(10):ncol(ctCluster)]){
      diffMean <- c(diffMean, t.test(ctCluster[, gene] ~ ctCluster$kmeans.cluster)$statistic)
    }
    names(diffMean) <- NULL #remove vector names
    upGenes <- names(ctCluster[,(which(diffMean > 0)+9)])

    genes <- names(ctCluster)[10:ncol(ctCluster)]
    
    geneDEtable <- as.data.frame(cbind(genes, pvals, diffMean), stringsAsFactors = F)
    geneDEtable$pvals <- as.numeric(geneDEtable$pvals)
    geneDEtable$diffMean <- as.numeric(geneDEtable$diffMean)
    
    geneDEup <- subset(geneDEtable, pvals < pValThreshold & diffMean > 0)
    geneDEup <- geneDEup[order(geneDEup$pvals),]
    
    geneDEup$pvals <- signif(geneDEup$pvals, digits = 3)
    
    
    if(nrow(geneDEup) > 0){
      row.names(geneDEup) <- NULL
      row.names(geneDEup) <- 1:nrow(geneDEup)
      
      groupedDEgenes <- sapply(geneGroups, function(x) x[match(geneDEup$genes, x, nomatch = F)])
      groupedDEgenes <- groupedDEgenes[lapply(groupedDEgenes, length) > 0]
      
      groupOrder <- sapply(geneGroups, function(x) match(x, geneDEup$genes, nomatch = F))
      groupOrder <- sapply(groupOrder, function(x) x[x>0])
      groupOrder <- groupOrder[lapply(groupOrder, length) > 0]
      groupOrder <- sapply(groupOrder, min)
      groupOrder <- groupOrder[order(unlist(groupOrder))]
      
      groupedDEgenes <- groupedDEgenes[names(groupOrder)]
      
      ctDE[[cluster]] <- groupedDEgenes
      
      maxListLength <- max(sapply(groupedDEgenes, length))
      deGeneFilledList <- lapply(groupedDEgenes, function(x) x <- c(x, rep(NA, maxListLength - length(x))))
      deGeneTable <- as.data.frame(do.call(cbind, deGeneFilledList), stringsAsFactors = F)
      
      
      ## print gene information
      print(paste("Differentially expressed genes (p-value < ", pValThreshold, ") for cluster ", cluster, ":", sep = ""), row.names=F, quote=F)
      print(geneDEup[, 1:2])
      cat("\n")
      
      print(paste("Cluster ", cluster, " categorized DE genes:", sep = ""), row.names=F, quote=F)
      print(deGeneTable, na.print = "", justify = "left", print.gap = 3, row.names = F)
      
      if(plotType == "violins"){
        ## plot violins for differentially expressed genes
        reportViolins(ctPlot = ctRep, deGenes = groupedDEgenes, clusters, cluster)
      }
      if(plotType == "tsne"){
        ## plot tSNE for differentially expressed genes
        reportTSNE(ctPlot = ctRep, tsne_y = tsne_y, deGenes = groupedDEgenes, clusters, cluster)
      }

    } else{
      print(paste("Differentially expressed genes (p-value < ", pValThreshold, ") for cluster ", cluster, ": NONE", sep = ""), row.names=F, quote=F)
      cat("\n")
      cat("\n")
      cat("\n")
      cat("\n")
      cat("\n")
    }
  }
  
  return(ctDE)
  
}