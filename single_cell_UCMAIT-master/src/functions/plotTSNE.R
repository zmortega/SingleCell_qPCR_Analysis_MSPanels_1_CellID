plotTSNE <- function(ctClust, colorby = c("kmeans.cluster", "probe", "patient", "age", "FRorFZ", "tissue", "IBD_Status", "Gene_List"), Genes) {
  ctRep <- ctClust
  #print(colnames(ctRep))
  write.csv(ctRep, paste0(baseDir, "results/plotTSNE-001.csv"))
  
  ## Decide which metadata column to color by
  #if (colorby == "kmeans.cluster") {
  #  color_var <- factor(ctRep$kmeans.cluster)
    
  #} else if (colorby == "probe") {
  #  color_var <- factor(ctRep$probe)
    
  #} else if (colorby == "patient") {
  #  color_var <- factor(ctRep$patient)
    
  #} else if (colorby == "age") {
  #  color_var <- factor(ctRep$age)
    
  #} else if (colorby == "tissue") {
  #  color_var <- factor(ctRep$cellSource)
    
  #} else {
  #  stop(paste("Unknown colorby:", colorby))
  #}
  
  #ctRep$tsne_color <- color_var
  
  
  # ##SMS_edit: Trying to break out Fresh vs Frozen status from "SPAM" column
  # ctRep[, "FRorFZ"] <- sapply(1:nrow(ctRep), function (x) {Status<- str_split(ctRep[x,"SPAM"], "\\.")[[1]][1]})
  # FRorFZ<- ctRep[,"FRorFZ"]
  write.csv(ctRep, paste0(baseDir, "results/plotTSNE-002.csv"))
  
  ###SMS edit:
  KmeansCol<- grep("kmeans.cluster", colnames(ctRep))
  FirstGene<- (KmeansCol+1)
  LastGene<- ncol(ctRep)
  
  ## remove duplicated rows
  #if(anyDuplicated(ctRep[,10:(ncol(ctRep))]) > 0){
  if(anyDuplicated(ctRep[,(FirstGene:LastGene)]) > 0){
    #ctRep <- ctRep[-which(duplicated(ctRep[,10:ncol(ctRep)])),]
    ctRep <- ctRep[-which(duplicated(ctRep[,(FirstGene:LastGene)])),]
  }
  
  ## remove genes without variance
  #print(colnames(ctRep))
  #LastColumn<- grep("cohort", colnames(ctRep)) ### SMS edit- Trying to find the last column to use for the for loop below. It is using a character column and then replacing that column with "NA".
  vars <- NULL
  #for(i in 10:(LastColumn-1)){
  for(i in FirstGene:LastGene){
    vars <- c(vars, var(ctRep[,i], na.rm=T))
  }
  ctGenesNoVar <- ctRep[which(vars == 0 | is.na(vars))+8] ###SMS: switched it to +8 from +9
  ctRep <- ctRep[,c(1:KmeansCol, (which(!is.na(vars) & vars!=0)+8))] ###Switched 8 here to KmeansCol and +8 from +9
  #ctRep[, "FRorFZ"] <- FRorFZ
  write.csv(ctRep, paste0(baseDir, "results/plotTSNE-003.csv"))
  
  #write.csv(ctRep, paste0(baseDir, "results/plotTSNE-004.csv"))

  #write.csv(ctRep, paste0(baseDir, "results/plotTSNE-005.csv"))
  
  ## run t-SNE
  set.seed(100)
  #tsne_out <- Rtsne(as.matrix(ctRep[, 10:ncol(ctRep)]), perplexity = 30)
  #tsne_out <- Rtsne(as.matrix(ctRep[, c(8,10:ncol(ctRep))]), perplexity = 35)
  #tsne_out <- Rtsne(as.matrix(ctRep[, c(8,10:(ncol(ctRep)-2))]), perplexity = 35)
  
  FirstGene<- (grep("kmeans.cluster", colnames(ctRep))+1)
  FirstGeneName<- colnames(ctRep)[FirstGene]
  LastGene<- ncol(ctRep)
  LastGeneName<- colnames(ctRep)[LastGene]
  
  tsne_out <- Rtsne(as.matrix(ctRep[,c(FirstGene:LastGene)]), perplexity = 30) #35) ###SMS: Switching to First Gene and LastGene.
  
  tsne_y <- as.data.frame(cbind(tsne_out$Y, 
                                ctRep$cellSource,
                                ctRep$kmeans.cluster,
                                ctRep$age,
                                ctRep$patient,
                                ctRep$probe,
                                #ctRep$FRorFZ,
                                #ctRep$cohort,
                                ctRep$cellType,
                                #ctRep[, 10:(ncol(ctRep)-2)]))
                                ctRep[,FirstGene:LastGene]))
  write.csv(tsne_y, paste0(baseDir, "results/plotTSNE-006.csv"))
  
  ###SMS: Break out ALL metadata post tSNE analysis here. This should eliminate the need to break out meta-data just before each tSNE plot and standardize the dataframe within a species
  if (Species == "HU") {
    tsne_y$FRorFZ<- sapply(1:nrow(tsne_y), function (x) {str_split(tsne_y[x,"ctRep$cellType"], "_")[[1]][1] })
    CohortVector<- sapply(1:nrow(tsne_y), function (x) {str_split(tsne_y[x,"ctRep$cellType"], "_")[[1]][10] })
    for (i in 1:length(CohortVector)) {
      RUS_Test<- grepl("RUS*", CohortVector[i])
      #print(RUS_Test)
      UCM_Test<- grepl("UCM*", CohortVector[i])
      #print(UCM_Test)
      OXF_Test<- grepl("OXF*", CohortVector[i])
      #print(OXF_Test)
      NBD_Test<- grepl("NBD*", CohortVector[i])
      #print(NBD_Test)
      if (RUS_Test == TRUE) {
        CohortVector[i] <- "RUS"
      } else if (UCM_Test == TRUE) {
        CohortVector[i] <- "UCM"
      } else if (OXF_Test == TRUE) {
        CohortVector[i] <- "OXF"
      } else if (NBD_Test == TRUE) {
        CohortVector[i] <- "NBD"
      }
    }
    tsne_y$cohort <- CohortVector

    OrganPattern<- c("-CO-","-SI-","_BL_")
    Organ<- c("colon","small_intestine","blood")
    for (i in 1:length(OrganPattern)) {
      tsne_y[grepl(OrganPattern[i], tsne_y$'ctRep$cellType', fixed=TRUE), "organ"] <- Organ[i]
    }

    OrganLocationPattern <- c("-CECU-","-ASCN-","-TRVZ-","-DESC-","-SIGM-","-RECT-","-DUO","_BL_")
    OrganLocation<- c("cecum","ascending_colon","transverse_colon","descending_colon","sigmoid_colon","rectum","duodenum","blood")
    for (i in 1:length(OrganLocationPattern)) {
      tsne_y[grepl(OrganLocationPattern[i], tsne_y$'ctRep$cellType', fixed=TRUE), "organlocation"] <- OrganLocation[i]
    }

    TissueLocationPattern <- c("-IEL_","-LPL_","-ALL_","_BL_")
    TissueLocation<- c("intra-epithelial","lamina-propria","total_digestion","blood")
    for (i in 1:length(TissueLocationPattern)) {
      tsne_y[grepl(TissueLocationPattern[i], tsne_y$'ctRep$cellType', fixed=TRUE), "tissuelocation"] <- TissueLocation[i]
    }
  } else if (Species == "MS") {
    tsne_y$FRorFZ<- sapply(1:nrow(tsne_y), function (x) {str_split(tsne_y[x,"ctRep$cellType"], "_")[[1]][1] })
    tsne_y$cohort<- sapply(1:nrow(tsne_y), function (x) {str_split(tsne_y[x,"ctRep$cellType"], "_")[[1]][8] })
    tsne_y$organ<- sapply(1:nrow(tsne_y), function (x) {str_split(tsne_y[x,"ctRep$cellType"], "_")[[1]][9] })
    tsne_y$organlocation <- sapply(1:nrow(tsne_y), function (x) {paste(str_split(tsne_y[x,"ctRep$cellType"], "_")[[1]][13], str_split(tsne_y[x,"ctRep$cellType"], "_")[[1]][8], sep = "_") })
    tsne_y$tissuelocation<- sapply(1:nrow(tsne_y), function (x) {str_split(tsne_y[x,"ctRep$cellType"], "_")[[1]][13] })
    
  }
  
  FirstGenePosition<- grep(FirstGeneName, colnames(tsne_y))
  LastGenePosition<- grep(paste0(LastGeneName,"$"), colnames(tsne_y))
  
  ColNumOrder<- c()
  
  for (i in c("^1$", "^2$", "^ctRep\\$cellSource$", "^ctRep\\$kmeans.cluster$", "^ctRep\\$age$", "^ctRep\\$patient$", "^ctRep\\$probe$", "^cohort$", "^FRorFZ$","^organ$","^organlocation$","^tissuelocation$", "^ctRep\\$cellType$")) {
    ColNumOrder<- c(ColNumOrder, grep(i, colnames(tsne_y)))
  }
  
  tsne_y<- tsne_y[,c(ColNumOrder,c(FirstGenePosition:LastGenePosition))]
  
  names(tsne_y)[1:13] <- c(
    "y1",
    "y2",
    "cellSource",
    "kmeans.cluster",
    "age",
    "patient",
    "probe",
    "cohort",
    "FRorFZ",
    "organ",
    "organlocation",
    "tissuelocation",
    "cellType"
  )
  
  
  write.csv(tsne_y, paste0(baseDir, "results/plotTSNE-007.csv"))
  
  #tsne_y$kmeans.cluster <- factor(tsne_y$kmeans.cluster, levels=clusters)
  #relTissues <- unique(tsne_y$tissue)
  #relTissues <- relTissues[match(c("I", "LN", "S", "PB"), relTissues, nomatch=F)]
  #tsne_y$tissue <- factor(tsne_y$tissue, levels=relTissues)
  
  FirstGene<- (grep(pattern = "cellType", x = colnames(tsne_y))+1)
  for(i in c(1,2,(FirstGene):(ncol(tsne_y)))){
    tsne_y[, i] <- as.numeric(tsne_y[, i])
  }

  ## set color for t-SNE plot
  pointSize <- 4
  myPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  clusterPalette <- "Set1" ##SMS Edit: LG original setting "Set3"
  shapeVals <- c(19, 17, 15, 18)
  numTissues <- length(unique(tsne_y$tissue))
  shapeVals <- shapeVals[1:numTissues]
  
  ## label cells by kmeans.cluster
  if("kmeans.cluster" %in% colorby) {
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$kmeans.cluster)), size = pointSize, alpha = 1) +
      scale_colour_brewer(palette = clusterPalette) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="cluster", order = 1)) +
      ggtitle("t-SNE between tissues (colored by kmeans.cluster)") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        # panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        # panel.grid.minor=element_line(color="gray90"),
        # panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(plTSNE)
  }
  
  
  ## label cells by probe
  if("probe" %in% colorby) {
    # tsne_y$probeColor <- NA
    # tsne_y[which(tsne_y$probe=="UCM5"), "probeColor"] <- "#A6CEE3"
    # tsne_y[which(tsne_y$probe=="UCM6"), "probeColor"] <- "#1F78B4"
    # tsne_y[which(tsne_y$probe=="UCM9"), "probeColor"] <- "blue3"
    # tsne_y[which(tsne_y$probe=="UCM10"), "probeColor"] <- "turquoise"
    # tsne_y[which(tsne_y$probe=="UCM12"), "probeColor"] <- "navy"
    # tsne_y[which(tsne_y$probe=="UCM13"), "probeColor"] <- "maroon"
    # tsne_y[which(tsne_y$probe=="UCM14"), "probeColor"] <- "darkmagenta"
    # tsne_y[which(tsne_y$probe=="UCM15"), "probeColor"] <- "thistle3"
    # tsne_y[which(tsne_y$probe=="UCM16"), "probeColor"] <- "peru"
    # tsne_y[which(tsne_y$probe=="UCM17"), "probeColor"] <- "purple4"
    # tsne_y[which(tsne_y$probe=="NBD1"), "probeColor"] <- "#B2DF8A"
    # tsne_y[which(tsne_y$probe=="NBD3"), "probeColor"] <- "#33A02C"
    # tsne_y[which(tsne_y$probe=="NBD4"), "probeColor"] <- "yellow2"
    # tsne_y[which(tsne_y$probe=="RUS002"), "probeColor"] <- "#003333"
    # tsne_y[which(tsne_y$probe=="RUS008"), "probeColor"] <- "#9933FF"
    # tsne_y[which(tsne_y$probe=="OXFGI4330"), "probeColor"] <- "orangered"
    # tsne_y[which(tsne_y$probe=="OXFGI4870"), "probeColor"] <- "plum4"
    # tsne_y[which(tsne_y$probe=="OXFIBD1322"), "probeColor"] <- "violetred4"
    # tsne_y[which(tsne_y$probe=="RUS001"), "probeColor"] <- "indianred1"
    # tsne_y[which(tsne_y$probe=="RUS011"), "probeColor"] <- "#3333FF"
    # tsne_y[which(tsne_y$probe=="RUS007"), "probeColor"] <- "#003300"
    #Probes <- tsne_y$probeColor
    #tsne_y <- subset(tsne_y, select=-c(probeColor))
    #AnnoColors <- cbind(AnnoColors, Probes)
    
    #plot(tsne_y$y1, tsne_y$y2, main="t-SNE between tissues (colored by probe)", col=tsne_y$probeColor, asp = 1, pch = 20)
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$probe)), size = pointSize, alpha = 1) +
      scale_color_manual(values=c(contrastColorList[1:length(unique(tsne_y$probe))])) + 
      #SMS edit out: used Karli's list of colors instead of assigning individuals manually; LG code: scale_color_manual(values=c("orangered", "purple4","#A6CEE3", "#1F78B4", "blue3", "turquoise", "navy", "maroon", "darkmagenta", "thistle3", "peru", "#B2DF8A", "#33A02C", "yellow2", "#003333", "#9933FF", "plum4", "violetred4", "indianred1", "#3333FF", "#003300","#AA4371","#c00000")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="probe", order = 1)) +
      ggtitle("t-SNE between tissues (colored by probe)") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        # panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        # panel.grid.minor=element_line(color="gray90"),
        # panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(plTSNE)
  }
  
  ## label cells by probe
  if("organ" %in% colorby) {
    # tsne_y$probeColor <- NA
    # tsne_y[which(tsne_y$probe=="UCM5"), "probeColor"] <- "#A6CEE3"
    # tsne_y[which(tsne_y$probe=="UCM6"), "probeColor"] <- "#1F78B4"
    # tsne_y[which(tsne_y$probe=="UCM9"), "probeColor"] <- "blue3"
    # tsne_y[which(tsne_y$probe=="UCM10"), "probeColor"] <- "turquoise"
    # tsne_y[which(tsne_y$probe=="UCM12"), "probeColor"] <- "navy"
    # tsne_y[which(tsne_y$probe=="UCM13"), "probeColor"] <- "maroon"
    # tsne_y[which(tsne_y$probe=="UCM14"), "probeColor"] <- "darkmagenta"
    # tsne_y[which(tsne_y$probe=="UCM15"), "probeColor"] <- "thistle3"
    # tsne_y[which(tsne_y$probe=="UCM16"), "probeColor"] <- "peru"
    # tsne_y[which(tsne_y$probe=="UCM17"), "probeColor"] <- "purple4"
    # tsne_y[which(tsne_y$probe=="NBD1"), "probeColor"] <- "#B2DF8A"
    # tsne_y[which(tsne_y$probe=="NBD3"), "probeColor"] <- "#33A02C"
    # tsne_y[which(tsne_y$probe=="NBD4"), "probeColor"] <- "yellow2"
    # tsne_y[which(tsne_y$probe=="RUS002"), "probeColor"] <- "#003333"
    # tsne_y[which(tsne_y$probe=="RUS008"), "probeColor"] <- "#9933FF"
    # tsne_y[which(tsne_y$probe=="OXFGI4330"), "probeColor"] <- "orangered"
    # tsne_y[which(tsne_y$probe=="OXFGI4870"), "probeColor"] <- "plum4"
    # tsne_y[which(tsne_y$probe=="OXFIBD1322"), "probeColor"] <- "violetred4"
    # tsne_y[which(tsne_y$probe=="RUS001"), "probeColor"] <- "indianred1"
    # tsne_y[which(tsne_y$probe=="RUS011"), "probeColor"] <- "#3333FF"
    # tsne_y[which(tsne_y$probe=="RUS007"), "probeColor"] <- "#003300"
    #Probes <- tsne_y$probeColor
    #tsne_y <- subset(tsne_y, select=-c(probeColor))
    #AnnoColors <- cbind(AnnoColors, Probes)
    
    #plot(tsne_y$y1, tsne_y$y2, main="t-SNE between tissues (colored by probe)", col=tsne_y$probeColor, asp = 1, pch = 20)
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$organ)), size = pointSize, alpha = 1) +
      scale_color_manual(values=c(contrastColorList[1:length(unique(tsne_y$organ))])) + 
      #SMS edit out: used Karli's list of colors instead of assigning individuals manually; LG code: scale_color_manual(values=c("orangered", "purple4","#A6CEE3", "#1F78B4", "blue3", "turquoise", "navy", "maroon", "darkmagenta", "thistle3", "peru", "#B2DF8A", "#33A02C", "yellow2", "#003333", "#9933FF", "plum4", "violetred4", "indianred1", "#3333FF", "#003300","#AA4371","#c00000")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="organ", order = 1)) +
      ggtitle("t-SNE between tissues (colored by organ)") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        # panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        # panel.grid.minor=element_line(color="gray90"),
        # panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(plTSNE)
  }
  
  ## label cells by patient
  if("patient" %in% colorby) {
    # keep a placeholder patientColor column so the bottom grep("patientColor") still works
    if(!"patientColor" %in% colnames(tsne_y)) {
      tsne_y$patientColor <- NA
    }
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color = factor(tsne_y$patient)), size = pointSize, alpha = 1) +
      scale_color_manual(values = contrastColorList[1:length(unique(tsne_y$patient))]) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="patient", order = 1)) +
      ggtitle("t-SNE colored by patient") +
      theme_minimal() +
      theme(
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
      )
    print(plTSNE)
  }
  
  ## label cells by age
  if("age" %in% colorby) {
    tsne_y$ageColor <- NA
    tsne_y[which(tsne_y$age=="infl"), "ageColor"] <- "deepskyblue2"
    tsne_y[which(tsne_y$age=="uninfl"), "ageColor"] <- "navy"
    tsne_y[which(tsne_y$age=="blood"), "ageColor"] <- "orangered"
    #Ages <- tsne_y$ageColor
    #tsne_y <- subset(tsne_y, select=-c(ageColor))
    #AnnoColors <- cbind(AnnoColors, Ages)
    
    #plot(tsne_y$y1, tsne_y$y2, main="t-SNE colored by inflamed site", col=tsne_y$ageColor, asp = 1, pch = 20)
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$age)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("orangered", "navy", "deepskyblue2")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="inflamation", order = 1)) +
      ggtitle("t-SNE colored by tissue source") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        # panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        # panel.grid.minor=element_line(color="gray90"),
        # panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(plTSNE)
  }
  
  ## label cells by Fresh or Frozen - SMS edit
  if("FRorFZ" %in% colorby) {
    tsne_y$ageColor <- NA
    tsne_y[which(tsne_y$FRorFZ=="FR"), "ageColor"] <- "deepskyblue2"
    tsne_y[which(tsne_y$FRorFZ=="FZ"), "ageColor"] <- "orangered"
    #tsne_y[which(tsne_y$age=="infl"), "ageColor"] <- "deepskyblue2"
    #tsne_y[which(tsne_y$age=="uninfl"), "ageColor"] <- "navy"
    #tsne_y[which(tsne_y$age=="blood"), "ageColor"] <- "orangered"
    #Ages <- tsne_y$ageColor
    #tsne_y <- subset(tsne_y, select=-c(ageColor))
    #AnnoColors <- cbind(AnnoColors, Ages)
    
    #plot(tsne_y$y1, tsne_y$y2, main="t-SNE colored by inflamed site", col=tsne_y$ageColor, asp = 1, pch = 20)
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$FRorFZ)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("deepskyblue2", "orangered")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="Frozen or Fresh", order = 1)) +
      ggtitle("t-SNE colored by fresh or frozen status") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        # panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        # panel.grid.minor=element_line(color="gray90"),
        # panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(plTSNE)
  }
  
  ## label cells by tissue (location) - SMS edit
  if("tissue" %in% colorby) {
    #tsne_y$ageColor <- NA
    # tsne_y[grepl("-CO-", tsne_y$cellType, fixed=TRUE), "organ"] <- "colon"
    # tsne_y[grepl("-SI-", tsne_y$cellType, fixed=TRUE), "organ"] <- "small_intestine"
    # tsne_y[grepl("_BL_", tsne_y$cellType, fixed=TRUE), "organ"] <- "blood"
    # tsne_y[grepl("-CECU-", tsne_y$cellType, fixed=TRUE), "organlocation"] <- "cecum"
    # tsne_y[grepl("-ASCN-", tsne_y$cellType, fixed=TRUE), "organlocation"] <- "ascending_colon"
    # tsne_y[grepl("-TRVZ-", tsne_y$cellType, fixed=TRUE), "organlocation"] <- "transverse_colon"
    # tsne_y[grepl("-DESC-", tsne_y$cellType, fixed=TRUE), "organlocation"] <- "descending_colon"
    # tsne_y[grepl("-SIGM-", tsne_y$cellType, fixed=TRUE), "organlocation"] <- "sigmoid_colon"
    # tsne_y[grepl("-RECT-", tsne_y$cellType, fixed=TRUE), "organlocation"] <- "rectum"
    # tsne_y[grepl("-DUO", tsne_y$cellType, fixed=TRUE), "organlocation"] <- "duodenum"
    # tsne_y[grepl("_BL_", tsne_y$cellType, fixed=TRUE), "organlocation"] <- "blood"
    # tsne_y[grepl("-ALL", tsne_y$cellType, fixed=TRUE), "organlocation"] <- "total_digestion"
    # tsne_y[grepl("-IEL_", tsne_y$cellType, fixed=TRUE), "tissuelocation"] <- "intra-epithelial"
    # tsne_y[grepl("-LPL_", tsne_y$cellType, fixed=TRUE), "tissuelocation"] <- "lamina-propria"
    # tsne_y[grepl("-ALL_", tsne_y$cellType, fixed=TRUE), "tissuelocation"] <- "total_digestion"
    # tsne_y[grepl("_BL_", tsne_y$cellType, fixed=TRUE), "tissuelocation"] <- "blood"
    # 
    write.csv(tsne_y, paste0(baseDir, "results/plotTSNE-007.csv"))
 
  
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$organlocation)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("deepskyblue2","firebrick1","springgreen3","gold","purple2","darkorange2","turquoise3","orchid1","chartreuse1","dodgerblue4")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="Tissue Location", order = 1)) +
      ggtitle("t-SNE colored by organ location") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        # panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        # panel.grid.minor=element_line(color="gray90"),
        # panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(plTSNE)
  }
  
  ## label cells by tissue (location) - SMS edit
  if("tissue" %in% colorby) {

    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$tissuelocation)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("deepskyblue2","firebrick1","springgreen3","gold","purple2","darkorange2","turquoise3","orchid1","chartreuse1","dodgerblue4")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="Tissue Location", order = 1)) +
      ggtitle("t-SNE colored by tissue location") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        # panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        # panel.grid.minor=element_line(color="gray90"),
        # panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(plTSNE)
  }
  
  ## label cells by IBD status - SMS edit
  UC<- c("RUS001","RUS003","RUS005","RUS009","RUS010","RUS012","RUS013","RUS014","RUS015","RUS016", "RUS017", "RUS018", "RUS020", "RUS022",
         "UCM10", "UCM12", "UCM13", "UCM14", "UCM15", "UCM16")
  CD<- c("RUS002","RUS004","RUS006","RUS007","RUS008","RUS011","RUS019", "RUS021")
  tsne_y[str_detect(pattern = paste(UC, collapse ="|"), tsne_y$cellType), "IBD_Status"] <- "UC"
  tsne_y[str_detect(pattern = paste(CD, collapse ="|"), tsne_y$cellType), "IBD_Status"] <- "CD"
  
  if("IBD_Status" %in% colorby) {
    
    plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
      geom_point(aes(color=factor(tsne_y$IBD_Status)), size = pointSize, alpha = 1) +
      scale_colour_manual(values=c("#B2DF8A","thistle3")) +
      scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                         minor_breaks = NULL) +
      scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                         minor_breaks = NULL) +
      scale_shape_manual(values = shapeVals) +
      guides(color=guide_legend(title="IBD Status", order = 1)) +
      ggtitle("t-SNE colored by IBD Status") +
      theme_minimal() +
      theme(#axis.line=element_blank(),
        # panel.border=element_rect(fill=NA, color="gray75", size=0.4),
        # panel.grid.minor=element_line(color="gray90"),
        # panel.grid.major=element_line(color="gray85", size=0.3),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        text=element_text(size=28),
        plot.margin=unit(c(14,9,14,9),"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(plTSNE)
  }
  
  if("Gene_List" %in% colorby & any("None" != Genes) == TRUE) {
    print(Genes)
    #FirstNumber<- grep(x = colnames(tsne_y), pattern = "cellType")
    #LastNumber<-  grep(x = colnames(tsne_y), pattern = "probeColor")
    GenestoPull<- grep(pattern = paste(Genes, collapse = "|"), x = colnames(tsne_y))
    
    for (i in GenestoPull) {
      plTSNE <- ggplot(tsne_y, aes(y1, y2)) +
        geom_point(aes(color=tsne_y[,i]), size = pointSize, alpha = 1) +
        scale_colour_gradient(low = "grey75", high = "red", guide = "colourbar") +
        scale_x_continuous(breaks=seq(min(tsne_y$y1), max(tsne_y$y1), length.out = 10),
                           minor_breaks = NULL) +
        scale_y_continuous(breaks=seq(min(tsne_y$y2), max(tsne_y$y2), length.out = 10),
                           minor_breaks = NULL) +
        scale_shape_manual(values = shapeVals) +
        guides(color=guide_legend(title=i, order = 1)) +
        ggtitle(paste("t-SNE colored by ", colnames(tsne_y)[i], " expression")) +
        theme_minimal() +
        theme(#axis.line=element_blank(),
          # panel.border=element_rect(fill=NA, color="gray75", size=0.4),
          # panel.grid.minor=element_line(color="gray90"),
          # panel.grid.major=element_line(color="gray85", size=0.3),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          text=element_text(size=28),
          plot.margin=unit(c(14,9,14,9),"cm"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      print(plTSNE)
    }
  }
  
  ctClust <- tsne_y  ### fuse tsne_y into ctClust for downstream plots
  
  CellTypeColumn  <- grep("cellType",       colnames(ctClust))
  KmeansColumn    <- grep("kmeans.cluster", colnames(ctClust))
  PatientColorCol <- grep("patientColor",   colnames(ctClust))
  EndColumn       <- ncol(ctClust)
  
  if (length(PatientColorCol) == 0) {
    ## No patientColor column present (e.g. colorby doesn't include "patient"):
    ## just keep everything after cellType in its existing order
    ctClust <- ctClust[ c(
      1:(KmeansColumn - 1),
      (KmeansColumn + 1):CellTypeColumn,
      KmeansColumn,
      (CellTypeColumn + 1):EndColumn
    ) ]
  } else {
    ## Original behavior when patientColor exists
    ctClust <- ctClust[ c(
      1:(KmeansColumn - 1),
      (KmeansColumn + 1):CellTypeColumn,
      KmeansColumn,
      (CellTypeColumn + 1):(PatientColorCol - 1),
      PatientColorCol:EndColumn
    ) ]
  }
  
  write.csv(ctClust, paste0(baseDir, "results/plotTSNE-008.csv"))
  return(ctClust)
  #invisible(ctClust)
}