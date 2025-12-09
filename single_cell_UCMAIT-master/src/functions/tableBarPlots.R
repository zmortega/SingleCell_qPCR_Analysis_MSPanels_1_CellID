tableBarPlots <- function(plotFactor = "tissues", plotFactorVals = cellSources, clusterVals, factorTable = sourceTable, colorKey = sourceColorKey){

  rawFactorTable <- factorTable
  
  ## set color key
  #print(as.vector(colorKey[,2])) ##SMS Edit, because Louis sucks.
  factorColorKeyVec <- as.vector(colorKey[,2])
  #print(factorColorKeyVec) ##SMS Edit, because Louis sucks.
  #print(plotFactorVals) ##SMS Edit, because Louis sucks.
  names(factorColorKeyVec) <- plotFactorVals
  #print(factorColorKeyVec) ##SMS Edit, because Louis sucks.
  
  ## create factor table data frame
  factorTable <- as.data.frame(factorTable)
  names(factorTable) <- plotFactorVals
  
  factorTable$cluster <- as.factor(clusterVals)
  factorTableMelt <- melt(factorTable)
  
  ## calculate positions for text labels
  factorTableMelt <- ddply(factorTableMelt, .(variable), transform, pos = sum(value) - (cumsum(value) - (0.5 * value)))
  factorTableMelt$count <- factorTableMelt$value
  factorTableMelt[which(factorTableMelt$value < 5), "count"] <- ""
  
  ##SMS edit, to see what is in the factorTableMelt variable
  write.csv(factorTableMelt, paste0(baseDir, "results/tableBarPlots-001.csv"))
  
  factorTablePlot <- ggplot(factorTableMelt, aes(x = factorTableMelt$variable, y = value)) +
    geom_bar(aes(fill = factorTableMelt$cluster), stat = "identity", position = "stack") +
    geom_text(aes(label = count, y = pos), size = 7) +
    scale_fill_brewer(palette = "Set3") +
    labs(title="") +
    labs(x=plotFactor) +
    labs(y=paste(plotFactor, " counts for each cluster", sep = "")) +
    guides(fill=guide_legend(title="cluster", keyheight = 2)) +
    theme_minimal() +
    scale_y_continuous() +
    theme(text=element_text(size=27),
          panel.grid.minor=element_line(color="gray90"),
          panel.grid.major=element_line(color="gray75", size=0.4),
          panel.grid.major.x=element_blank(),
          axis.title.x=element_text(vjust=0),
          axis.text.x=element_text(angle = 45, vjust = 0.5),
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
  #print(factorTablePlot)
  
  factorTableMelt$cluster2 <- factor(factorTableMelt$cluster, levels = rev(clusterVals))
  
  factorTableHeatPlot <- ggplot(factorTableMelt, aes(x = factorTableMelt$variable, y = factorTableMelt$cluster2)) +
    geom_raster(aes(fill = value)) +
    geom_text(aes(label = value), size = 7) +
    scale_fill_gradient(low = "white", high = "#50A6C2") +
    labs(title="") +
    labs(x=plotFactor) +
    labs(y="cluster") +
    theme_minimal() +
    theme(text=element_text(size=27),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.title.x=element_text(vjust=0),
          axis.text.x=element_text(angle = 45, vjust = 0.5),
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
  
  
  ## compute column proportions
  factorTableCol <- as.data.frame(t(t(rawFactorTable)/colSums(rawFactorTable)))
  
  ## melt data frame
  factorTableCol$cluster <- as.factor(clusterVals)
  factorTableColMelt <- melt(factorTableCol)
  
  ## calculate positions for text labels
  factorTableColMelt <- ddply(factorTableColMelt, .(variable), transform, pos = cumsum(value) - (0.5 * value))
  factorTableColMelt$pos <- 1 - factorTableColMelt$pos
  
  ## create rounded frequency column and remove small values
  factorTableColMelt$freq <- round(factorTableColMelt$value, digits = 2)
  factorTableColMelt$freq2 <- factorTableColMelt$freq
  factorTableColMelt[which(factorTableColMelt$freq < 0.05), "freq2"] <- ""
  
  ## create plot
  factorTableColPlot <- ggplot(factorTableColMelt, aes(x = factorTableColMelt$variable, y = value)) +
    geom_bar(aes(fill = factorTableColMelt$cluster), stat = "identity", position = "fill") +
    geom_text(aes(label = factorTableColMelt$freq2, y = pos), size = 7) +
    scale_fill_brewer(palette = "Set3") +
    labs(title="") +
    labs(x=plotFactor) +
    labs(y=paste(plotFactor, " fractions for each cluster", sep = "")) +
    guides(fill=guide_legend(title="cluster", keyheight = 2)) +
    theme_minimal() +
    scale_y_continuous() +
    theme(text=element_text(size=27),
          panel.grid.minor=element_line(color="gray90"),
          panel.grid.major=element_line(color="gray75", size=0.4),
          panel.grid.major.x=element_blank(),
          axis.title.x=element_text(vjust=0),
          axis.text.x=element_text(angle = 45, vjust = 0.5),
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
  
  #print(factorTableColPlot)
  
  ## row proportions
  factorTableRow <- as.data.frame(rawFactorTable/rowSums(rawFactorTable))
  
  names(factorTableRow) <- plotFactorVals
  factorTableRow$cluster <- row.names(factorTableRow)
  factorTableRowMelt <- melt(factorTableRow)
  
  factorTableRowMelt <- ddply(factorTableRowMelt, .(cluster), transform, pos = cumsum(value) - (0.5 * value))
  factorTableRowMelt$pos <- factorTableRowMelt$pos
  factorTableRowMelt$freq <- round(factorTableRowMelt$value, digits = 2)
  factorTableRowMelt$freq2 <- factorTableRowMelt$freq
  factorTableRowMelt[which(factorTableRowMelt$freq < 0.1), "freq2"] <- ""
  
  ## reverse order of axis values for coord flip
  factorTableRowMelt$cluster <- factor(factorTableRowMelt$cluster, levels = rev(unique(factorTableRowMelt$cluster)))
  factorTableRowMelt$variable <- factor(factorTableRowMelt$variable, levels = rev(unique(factorTableRowMelt$variable)))
  
  ##SMS edit
  write.csv(factorTableRowMelt, paste0(baseDir, "results/tableBarPlots-002.csv"))
  
  factorTableRowPlot <- ggplot(factorTableRowMelt, aes(x = factorTableRowMelt$cluster, y = value)) +
    geom_bar(aes(fill = factorTableRowMelt$variable), stat = "identity", position = "fill") +
    geom_text(aes(label = factorTableRowMelt$freq2, y = pos), size = 7) +
    scale_fill_manual(values = contrastColorList, drop = TRUE) +
    labs(title="") +
    labs(x="") +
    labs(y=paste("cluster fractions for each ", plotFactor, sep="")) +
    guides(fill=guide_legend(title=plotFactor, keyheight = 2, reverse = T)) +
    theme_minimal() +
    coord_flip() +
    scale_y_continuous() +
    theme(text=element_text(size=27),
          panel.grid.minor=element_line(color="gray90"),
          panel.grid.major=element_line(color="gray75", size=0.4),
          panel.grid.major.y=element_blank(),
          axis.title.x=element_text(vjust=0),
          plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
  
  return(list(factorTableColPlot, factorTableRowPlot, factorTablePlot, factorTableHeatPlot))
  
}