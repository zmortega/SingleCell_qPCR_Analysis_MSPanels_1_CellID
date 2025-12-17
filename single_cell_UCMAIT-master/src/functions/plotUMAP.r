plotUMAP <- function(ctClust,
                     colorby = c("kmeans.cluster", "probe", "patient", "age",
                                 "FRorFZ", "tissue", "IBD_Status", "Gene_List"),
                     Genes) {
  
  ctRep <- ctClust
  write.csv(ctRep, paste0(baseDir, "results/plotUMAP-001.csv"))
  
  write.csv(ctRep, paste0(baseDir, "results/plotUMAP-002.csv"))
  
  ## Identify gene columns (everything after kmeans.cluster)
  KmeansCol <- grep("kmeans.cluster", colnames(ctRep))
  FirstGene <- (KmeansCol + 1)
  LastGene  <- ncol(ctRep)
  
  ## remove duplicated rows based on gene-expression columns
  if (anyDuplicated(ctRep[, (FirstGene:LastGene)]) > 0) {
    ctRep <- ctRep[-which(duplicated(ctRep[, (FirstGene:LastGene)])), ]
  }
  
  ## remove genes without variance
  vars <- NULL
  for (i in FirstGene:LastGene) {
    vars <- c(vars, var(ctRep[, i], na.rm = TRUE))
  }
  
  ## preserve your original column math (+8) as-written (so behavior matches)
  ctGenesNoVar <- ctRep[which(vars == 0 | is.na(vars)) + 8]
  ctRep <- ctRep[, c(1:KmeansCol, (which(!is.na(vars) & vars != 0) + 8))]
  write.csv(ctRep, paste0(baseDir, "results/plotUMAP-003.csv"))
  
  ## recompute gene boundaries after filtering
  FirstGene <- (grep("kmeans.cluster", colnames(ctRep)) + 1)
  FirstGeneName <- colnames(ctRep)[FirstGene]
  LastGene <- ncol(ctRep)
  LastGeneName <- colnames(ctRep)[LastGene]
  
  ## ---- run UMAP (replaces t-SNE) ----
  set.seed(100)
  X <- as.matrix(ctRep[, c(FirstGene:LastGene)])
  
  umap_Y <- uwot::umap(
    X,
    n_components = 2,
    n_neighbors  = 30,
    min_dist     = 0.3,
    metric       = "euclidean",
    verbose      = FALSE
  )
  
  tsne_y <- as.data.frame(cbind(
    umap_Y,
    ctRep$cellSource,
    ctRep$kmeans.cluster,
    ctRep$age,
    ctRep$patient,
    ctRep$probe,
    ctRep$cellType,
    ctRep[, FirstGene:LastGene]
  ))
  write.csv(tsne_y, paste0(baseDir, "results/plotUMAP-006.csv"))
  
  ## --- Metadata parsing block (unchanged, still keys off Species) ---
  if (Species == "HU") {
    
    tsne_y$FRorFZ <- sapply(1:nrow(tsne_y), function(x) {
      str_split(tsne_y[x, "ctRep$cellType"], "_")[[1]][1]
    })
    
    CohortVector <- sapply(1:nrow(tsne_y), function(x) {
      str_split(tsne_y[x, "ctRep$cellType"], "_")[[1]][10]
    })
    
    for (i in 1:length(CohortVector)) {
      if (grepl("RUS*", CohortVector[i])) {
        CohortVector[i] <- "RUS"
      } else if (grepl("UCM*", CohortVector[i])) {
        CohortVector[i] <- "UCM"
      } else if (grepl("OXF*", CohortVector[i])) {
        CohortVector[i] <- "OXF"
      } else if (grepl("NBD*", CohortVector[i])) {
        CohortVector[i] <- "NBD"
      }
    }
    tsne_y$cohort <- CohortVector
    
    OrganPattern <- c("-CO-", "-SI-", "_BL_")
    Organ <- c("colon", "small_intestine", "blood")
    for (i in 1:length(OrganPattern)) {
      tsne_y[grepl(OrganPattern[i], tsne_y$`ctRep$cellType`, fixed = TRUE), "organ"] <- Organ[i]
    }
    
    OrganLocationPattern <- c("-CECU-","-ASCN-","-TRVZ-","-DESC-","-SIGM-","-RECT-","-DUO","_BL_")
    OrganLocation <- c("cecum","ascending_colon","transverse_colon","descending_colon","sigmoid_colon","rectum","duodenum","blood")
    for (i in 1:length(OrganLocationPattern)) {
      tsne_y[grepl(OrganLocationPattern[i], tsne_y$`ctRep$cellType`, fixed = TRUE), "organlocation"] <- OrganLocation[i]
    }
    
    TissueLocationPattern <- c("-IEL_","-LPL_","-ALL_","_BL_")
    TissueLocation <- c("intra-epithelial","lamina-propria","total_digestion","blood")
    for (i in 1:length(TissueLocationPattern)) {
      tsne_y[grepl(TissueLocationPattern[i], tsne_y$`ctRep$cellType`, fixed = TRUE), "tissuelocation"] <- TissueLocation[i]
    }
    
  } else if (Species == "MS") {
    
    tsne_y$FRorFZ <- sapply(1:nrow(tsne_y), function(x) {
      str_split(tsne_y[x, "ctRep$cellType"], "_")[[1]][1]
    })
    tsne_y$cohort <- sapply(1:nrow(tsne_y), function(x) {
      str_split(tsne_y[x, "ctRep$cellType"], "_")[[1]][8]
    })
    tsne_y$organ <- sapply(1:nrow(tsne_y), function(x) {
      str_split(tsne_y[x, "ctRep$cellType"], "_")[[1]][9]
    })
    tsne_y$organlocation <- sapply(1:nrow(tsne_y), function(x) {
      paste(str_split(tsne_y[x, "ctRep$cellType"], "_")[[1]][13],
            str_split(tsne_y[x, "ctRep$cellType"], "_")[[1]][8],
            sep = "_")
    })
    tsne_y$tissuelocation <- sapply(1:nrow(tsne_y), function(x) {
      str_split(tsne_y[x, "ctRep$cellType"], "_")[[1]][13]
    })
  }
  
  ## column reorder + rename (unchanged)
  FirstGenePosition <- grep(FirstGeneName, colnames(tsne_y))
  LastGenePosition  <- grep(paste0(LastGeneName, "$"), colnames(tsne_y))
  
  ColNumOrder <- c()
  for (i in c("^1$", "^2$", "^ctRep\\$cellSource$", "^ctRep\\$kmeans.cluster$", "^ctRep\\$age$",
              "^ctRep\\$patient$", "^ctRep\\$probe$", "^cohort$", "^FRorFZ$", "^organ$",
              "^organlocation$", "^tissuelocation$", "^ctRep\\$cellType$")) {
    ColNumOrder <- c(ColNumOrder, grep(i, colnames(tsne_y)))
  }
  
  tsne_y <- tsne_y[, c(ColNumOrder, c(FirstGenePosition:LastGenePosition))]
  
  names(tsne_y)[1:13] <- c(
    "y1","y2","cellSource","kmeans.cluster","age","patient","probe",
    "cohort","FRorFZ","organ","organlocation","tissuelocation","cellType"
  )
  
  write.csv(tsne_y, paste0(baseDir, "results/plotUMAP-007.csv"))
  
  ## numeric conversion (unchanged)
  FirstGene <- (grep(pattern = "cellType", x = colnames(tsne_y)) + 1)
  for (i in c(1, 2, (FirstGene):(ncol(tsne_y)))) {
    tsne_y[, i] <- as.numeric(tsne_y[, i])
  }
  
  ## plotting settings (unchanged except titles)
  pointSize <- 4
  myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  clusterPalette <- "Set1"
  shapeVals <- c(19, 17, 15, 18)
  
  ## NOTE: your original uses tsne_y$tissue here; leaving as-is to match behavior
  numTissues <- length(unique(tsne_y$tissue))
  shapeVals <- shapeVals[1:numTissues]
  
  if ("kmeans.cluster" %in% colorby) {
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$kmeans.cluster)), size = pointSize, alpha = 1) +
      scale_colour_brewer(palette = clusterPalette) +
      scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color = guide_legend(title = "cluster", order = 1)) +
      ggtitle("UMAP between tissues (colored by kmeans.cluster)") +
      theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 28),
        plot.margin = unit(c(14,9,14,9), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  if ("probe" %in% colorby) {
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$probe)), size = pointSize, alpha = 1) +
      scale_color_manual(values = c(contrastColorList[1:length(unique(tsne_y$probe))])) +
      scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color = guide_legend(title = "probe", order = 1)) +
      ggtitle("UMAP between tissues (colored by probe)") +
      theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 28),
        plot.margin = unit(c(14,9,14,9), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  if ("organ" %in% colorby) {
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$organ)), size = pointSize, alpha = 1) +
      scale_color_manual(values = c(contrastColorList[1:length(unique(tsne_y$organ))])) +
      scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color = guide_legend(title = "organ", order = 1)) +
      ggtitle("UMAP between tissues (colored by organ)") +
      theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 28),
        plot.margin = unit(c(14,9,14,9), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  if ("patient" %in% colorby) {
    if (!"patientColor" %in% colnames(tsne_y)) tsne_y$patientColor <- NA
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$patient)), size = pointSize, alpha = 1) +
      scale_color_manual(values = contrastColorList[1:length(unique(tsne_y$patient))]) +
      scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color = guide_legend(title = "patient", order = 1)) +
      ggtitle("UMAP colored by patient") +
      theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 28),
        plot.margin = unit(c(14,9,14,9), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  if ("age" %in% colorby) {
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$age)), size = pointSize, alpha = 1) +
      scale_colour_manual(values = c("orangered", "navy", "deepskyblue2")) +
      scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color = guide_legend(title = "inflamation", order = 1)) +
      ggtitle("UMAP colored by tissue source") +
      theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 28),
        plot.margin = unit(c(14,9,14,9), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  if ("FRorFZ" %in% colorby) {
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$FRorFZ)), size = pointSize, alpha = 1) +
      scale_colour_manual(values = c("deepskyblue2", "orangered")) +
      scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color = guide_legend(title = "Frozen or Fresh", order = 1)) +
      ggtitle("UMAP colored by fresh or frozen status") +
      theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 28),
        plot.margin = unit(c(14,9,14,9), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  if ("tissue" %in% colorby) {
    write.csv(tsne_y, paste0(baseDir, "results/plotUMAP-007.csv"))
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$organlocation)), size = pointSize, alpha = 1) +
      scale_colour_manual(values = c("deepskyblue2","firebrick1","springgreen3","gold","purple2",
                                     "darkorange2","turquoise3","orchid1","chartreuse1","dodgerblue4")) +
      scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color = guide_legend(title = "Tissue Location", order = 1)) +
      ggtitle("UMAP colored by organ location") +
      theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 28),
        plot.margin = unit(c(14,9,14,9), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  if ("tissue" %in% colorby) {
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$tissuelocation)), size = pointSize, alpha = 1) +
      scale_colour_manual(values = c("deepskyblue2","firebrick1","springgreen3","gold","purple2",
                                     "darkorange2","turquoise3","orchid1","chartreuse1","dodgerblue4")) +
      scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color = guide_legend(title = "Tissue Location", order = 1)) +
      ggtitle("UMAP colored by tissue location") +
      theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 28),
        plot.margin = unit(c(14,9,14,9), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  UC <- c("RUS001","RUS003","RUS005","RUS009","RUS010","RUS012","RUS013","RUS014","RUS015","RUS016","RUS017","RUS018","RUS020","RUS022",
          "UCM10","UCM12","UCM13","UCM14","UCM15","UCM16")
  CD <- c("RUS002","RUS004","RUS006","RUS007","RUS008","RUS011","RUS019","RUS021")
  tsne_y[str_detect(pattern = paste(UC, collapse="|"), tsne_y$cellType), "IBD_Status"] <- "UC"
  tsne_y[str_detect(pattern = paste(CD, collapse="|"), tsne_y$cellType), "IBD_Status"] <- "CD"
  
  if ("IBD_Status" %in% colorby) {
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$IBD_Status)), size = pointSize, alpha = 1) +
      scale_colour_manual(values = c("#B2DF8A","thistle3")) +
      scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
      scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color = guide_legend(title = "IBD Status", order = 1)) +
      ggtitle("UMAP colored by IBD Status") +
      theme_minimal() +
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size = 28),
        plot.margin = unit(c(14,9,14,9), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  if ("Gene_List" %in% colorby & any("None" != Genes) == TRUE) {
    GenestoPull <- grep(pattern = paste(Genes, collapse="|"), x = colnames(tsne_y))
    
    for (i in GenestoPull) {
      plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
        geom_point(aes(color = tsne_y[, i]), size = pointSize, alpha = 1) +
        scale_colour_gradient(low = "grey75", high = "red", guide = "colourbar") +
        scale_x_continuous(breaks = seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10), minor_breaks = NULL) +
        scale_y_continuous(breaks = seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10), minor_breaks = NULL) +
        scale_shape_manual(values = shapeVals) +
        guides(color = guide_legend(title = i, order = 1)) +
        ggtitle(paste("UMAP colored by", colnames(tsne_y)[i], "expression")) +
        theme_minimal() +
        theme(
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          text = element_text(size = 28),
          plot.margin = unit(c(14,9,14,9), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      print(plTSNE)
    }
  }
  
  ## fuse embedding into ctClust for downstream plots
  ctClust <- tsne_y
  
  CellTypeColumn  <- grep("cellType",       colnames(ctClust))
  KmeansColumn    <- grep("kmeans.cluster", colnames(ctClust))
  PatientColorCol <- grep("patientColor",   colnames(ctClust))
  EndColumn       <- ncol(ctClust)
  
  if (length(PatientColorCol) == 0) {
    ctClust <- ctClust[c(
      1:(KmeansColumn - 1),
      (KmeansColumn + 1):CellTypeColumn,
      KmeansColumn,
      (CellTypeColumn + 1):EndColumn
    )]
  } else {
    ctClust <- ctClust[c(
      1:(KmeansColumn - 1),
      (KmeansColumn + 1):CellTypeColumn,
      KmeansColumn,
      (CellTypeColumn + 1):(PatientColorCol - 1),
      PatientColorCol:EndColumn
    )]
  }
  
  write.csv(ctClust, paste0(baseDir, "results/plotUMAP-008.csv"))
  return(ctClust)
}
