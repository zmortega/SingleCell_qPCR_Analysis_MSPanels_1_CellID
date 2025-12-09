#### cluster and filter to remove outliers

clusterFilter <- function(ctInput, testK = F, numCenters = 2, plotHeatmap=F, plotClustOnly=F, 
                          heatmapFactor, heatmapColorBy, heatmapTissueLabel,
                          # fisherTests = c("probe", "tissue", "age", "mouse"),
                          fisherTests = NULL,
                          cumulativeExpHist = T, filterClusters = F, clustersToRemove = NULL,
                          removalText = T){
  
  #### kmeans analysis ####
  ## order by ctSum
  
  ###SMS edit: This function had extensive additions made to it at the top of the code. Because it is used twice, successively, a number of 'if/else' statements added to identify the first 
  ###use of the code from the second time it is applied to the same dataframe.
  
  ###SMS edit: This generates a csv file for the input dataframe so that we can figure out where an issues may lie.
  write.csv(ctInput, paste0(baseDir, "results/clusterFilter-001-ctInput-InputDataFrame.csv"))
  
  assign("ctInput", ctInput, envir = .GlobalEnv)
  print(ctInput)
  
  if (Troubleshoot == TRUE) {
    print(paste0("Column Names are: "))
    print(colnames(ctInput))
  }
  

  ###SMS edit: This creates the panel numbers that we will then use below. "If-else" statement to figure out what gene names we are going to use for each panel.
  if (Species == "HU") {
    if (PanelNumberDetected == "1" || PanelNumberDetected == "3") {
        Panel1<- c("AIM2","B2M","BCL2","BCL6","C3","C5","CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CD274","CD28","CD3E","CD4","CD40","CD44","CD52","CD53",
                   "CD8A","CDC42","CTLA4","CXCL10","CXCL9","CXCR3","CXCR4","FOXP3","FYN","GAPDH","GATA3","GSK3A","GSK3B","HLA-DRA","ICAM1","ICOS","IFI44","IFI44L",  
                   "IFIT1","IFIT3","IFNAR1", "IFNAR2","IFNG","IFNGR1","IL10","IL12RB1","IL17A","IL2","IL21","IL25","IL2RA","IL3","IL4","IL4R","IL5","IL5RA","IL7",
                   "IL7R","IRF1","IRF2","IRF4","IRF7","IRF9","ISG15","JAK1","JAK2","LY6E","NFKB1","NLRP3","NR4A1","PDCD1","PPARA","PPARG","PPARGC1A","PTEN",
                   "RAC1","RAC2","RORC","RPL13A","SOCS3","STAT1","STAT3","STAT4","STAT5B","TBX21","TGFBR2","TNF","TNFAIP3","TNFRSF1A","TNFRSF1B",
                   "TRAF2","TYK2","VAV1","ZAP70","ZEB2")
        Panel3<- c("AKT1","AKT1S1","ARPC2","BIRC5","BTG1","CALM1","CBLB","CCL5","CD200","CD200R1","CD27","CD40LG","CD48","CD69","CRIP1","CXCR6","DNM1L","DPP4","EVL",
                   "FAIM2","FAS","FOS","FOXO1","FOXP3","GATA3","GRB2","HAVCR2","HIF1A","HLA-DQB1","HRAS","ID2","IKZF2","IL6ST","IL7R","ITGA1","ITGA4","ITGB1","ITK",
                   "IZUMO1R","JUN","KLF2","KRAS","LAG3","LAT","LCK","LEF1","LGALS1","LYPD6","MAF","MALAT1","MAPK1","MAPK3","MAPK8","MTOR","MYB","MYC","NFATC1","NKG7",
                   "NRAS","ORAI1","PI3","PIK3CA","PLCG1","PML","PPP3CC","PRDM1","PRKAA2","PRKCQ","PTMA","PTPN6","PTPRC","RAP1A","RGS1","RICTOR","RPTOR","RUNX1","S100A4",
                   "S100A6","S1PR1","SELL","SLC2A1","SLC2A3","SREBF1","STMN1","TCF7","TIGIT","TLN1","TMSB10","TNFRSF18","TNFRSF4","TNFRSF9","TNFSF8","TRAF3IP2","TRAF5",
                   "TXK","ZBTB16")
        Panel1and3<- sort(c(Panel1, Panel3))
        print("Using gene names WITHOUT the underscores") 
        } else if (CombineCommonGenes == FALSE && PanelNumberDetected == "1 and 3") {
        ###But, there's an issue (there's always an issue), if you run this function with Panels "1 and 3" then you need to include underscores in the gene names to differentiate one set of genes from the other. 
        Panel1<- c("AIM2","B2M","BCL2","BCL6","C3_1","C5_1","CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CD274","CD28","CD3E","CD4","CD40","CD44","CD52","CD53",
                   "CD8A","CDC42_1","CTLA4","CXCL10_1","CXCL9","CXCR3","CXCR4","FOXP3","FYN","GAPDH","GATA3","GSK3A","GSK3B","HLA-DRA","ICAM1","ICOS","IFI44","IFI44L",  
                   "IFIT1","IFIT3","IFNAR1", "IFNAR2","IFNG_1","IFNGR1_1","IL10","IL12RB1_1","IL17A_1","IL2","IL21","IL25","IL2RA","IL3","IL4","IL4R","IL5","IL5RA","IL7",
                   "IL7R","IRF1","IRF2","IRF4","IRF7","IRF9","ISG15","JAK1","JAK2","LY6E","NFKB1","NLRP3","NR4A1","PDCD1_1","PPARA","PPARG","PPARGC1A","PTEN",
                   "RAC1_1","RAC2_1","RORC","RPL13A","SOCS3","STAT1","STAT3","STAT4","STAT5B","TBX21","TGFBR2","TNF_1","TNFAIP3","TNFRSF1A_1","TNFRSF1B",
                   "TRAF2","TYK2","VAV1","ZAP70","ZEB2")
        Panel3<- c("AKT1","AKT1S1","ARPC2","BIRC5","BTG1","CALM1","CBLB","CCL5","CD200","CD200R1","CD27","CD40LG","CD48","CD69","CRIP1","CXCR6","DNM1L","DPP4","EVL",
                   "FAIM2","FAS","FOS","FOXO1","FOXP3","GATA3","GRB2","HAVCR2","HIF1A","HLA-DQB1","HRAS","ID2","IKZF2","IL6ST","IL7R","ITGA1","ITGA4","ITGB1","ITK",
                   "IZUMO1R","JUN","KLF2","KRAS","LAG3","LAT","LCK","LEF1","LGALS1","LYPD6","MAF","MALAT1","MAPK1","MAPK3","MAPK8","MTOR","MYB","MYC","NFATC1","NKG7",
                   "NRAS","ORAI1","PI3","PIK3CA","PLCG1","PML","PPP3CC","PRDM1","PRKAA2","PRKCQ","PTMA","PTPN6","PTPRC","RAP1A","RGS1","RICTOR","RPTOR","RUNX1","S100A4",
                   "S100A6","S1PR1","SELL","SLC2A1","SLC2A3","SREBF1","STMN1","TCF7","TIGIT","TLN1","TMSB10","TNFRSF18","TNFRSF4","TNFRSF9","TNFSF8","TRAF3IP2","TRAF5",
                   "TXK","ZBTB16")
        Panel1and3<- sort(c(Panel1, Panel3))
        print("Using gene names WITH the underscores") 
        } else if (CombineCommonGenes == TRUE && PanelNumberDetected == "1 and 3") {
        Panel1<- c("AIM2","B2M","BCL2","BCL6","C3","C5","CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CD274","CD28","CD3E","CD4","CD40","CD44","CD52","CD53",
                   "CD8A","CDC42","CTLA4","CXCL10","CXCL9","CXCR3","CXCR4","FOXP3","FYN","GAPDH","GATA3","GSK3A","GSK3B","HLA-DRA","ICAM1","ICOS","IFI44","IFI44L",  
                   "IFIT1","IFIT3","IFNAR1", "IFNAR2","IFNG","IFNGR1","IL10","IL12RB1","IL17A","IL2","IL21","IL25","IL2RA","IL3","IL4","IL4R","IL5","IL5RA","IL7",
                   "IL7R","IRF1","IRF2","IRF4","IRF7","IRF9","ISG15","JAK1","JAK2","LY6E","NFKB1","NLRP3","NR4A1","PDCD1","PPARA","PPARG","PPARGC1A","PTEN",
                   "RAC1","RAC2","RORC","RPL13A","SOCS3","STAT1","STAT3","STAT4","STAT5B","TBX21","TGFBR2","TNF","TNFAIP3","TNFRSF1A","TNFRSF1B",
                   "TRAF2","TYK2","VAV1","ZAP70","ZEB2")
        Panel3<- c("AKT1","AKT1S1","ARPC2","BIRC5","BTG1","CALM1","CBLB","CCL5","CD200","CD200R1","CD27","CD40LG","CD48","CD69","CRIP1","CXCR6","DNM1L","DPP4","EVL",
                   "FAIM2","FAS","FOS","FOXO1","FOXP3","GATA3","GRB2","HAVCR2","HIF1A","HLA-DQB1","HRAS","ID2","IKZF2","IL6ST","IL7R","ITGA1","ITGA4","ITGB1","ITK",
                   "IZUMO1R","JUN","KLF2","KRAS","LAG3","LAT","LCK","LEF1","LGALS1","LYPD6","MAF","MALAT1","MAPK1","MAPK3","MAPK8","MTOR","MYB","MYC","NFATC1","NKG7",
                   "NRAS","ORAI1","PI3","PIK3CA","PLCG1","PML","PPP3CC","PRDM1","PRKAA2","PRKCQ","PTMA","PTPN6","PTPRC","RAP1A","RGS1","RICTOR","RPTOR","RUNX1","S100A4",
                   "S100A6","S1PR1","SELL","SLC2A1","SLC2A3","SREBF1","STMN1","TCF7","TIGIT","TLN1","TMSB10","TNFRSF18","TNFRSF4","TNFRSF9","TNFSF8","TRAF3IP2","TRAF5",
                   "TXK","ZBTB16")
        Panel3<- Panel3[-which(Panel3 %in% CommonGenes)]
        Panel1and3<- sort(c(Panel1, Panel3))
        }
  } else if (Species == "MS") {
    Panel1<- c("Aim2","Bcl2","Bcl6","Ccr1","ccr2","ccr3","ccr4","ccr5","ccr6","Ccr7","cd28","cd3e","cd4","cd40","Cd44","cd80","cd86","cd8a","Ceacam1","ctla4",
               "Cxcl10","Cxcr3","Cxcr4","foxp3","Fyn","gapdh","gata4","gsk3a","gsk3b","Hprt","icam1","Icos","Ifi44","Ifi44l","Ifit1","Ifit3","ifng","Ifngr1",
               "IL-21","il10","il12b","Il12rb","il17A","Il18r1","il1r2","il2","il25","Il27","Il27r","il2ra","il3","il4","il4ra","il5","il5ra","il6","il7","il7r",
               "Irf1","Irf2","Irf4","Irf7","Isg15","Jak1","Jak2","Ly6e","Map2k6","Mapk8","Mx1","nfkb1","Nur77","Oas1b","Oas2","Oasl1","Pd1","Pdl-1","ppara","pparg",
               "ppargc1a","pten","Rsad2","Socs3","Stat1","Stat3","Stat4","Stat5","Tbx21","Tgfbr2","tnf","Tnfaip3","tnfrsf1a","tnfrsf1b","Traf2","Vav1","Zap70","Zeb2")
    #Panel2<- c()
    Panel3<- c("AKT1","AKT1S1","ARPC2","BIRC5","BTG1","CALM1","CBLB","CCL5","CD200","CD200R1","CD27","CD40LG","CD48","CD69","CRIP1","CXCR6","DNM1L","DPP4","EVL",
               "FAIM2","FAS","FOS","FOXO1","FOXP3","GATA3","GRB2","HAVCR2","HIF1A","HLA-DQB1","HRAS","ID2","IKZF2","IL6ST","IL7R","ITGA1","ITGA4","ITGB1","ITK",
               "IZUMO1R","JUN","KLF2","KRAS","LAG3","LAT","LCK","LEF1","LGALS1","LYPD6","MAF","MALAT1","MAPK1","MAPK3","MAPK8","MTOR","MYB","MYC","NFATC1","NKG7",
               "NRAS","ORAI1","PI3","PIK3CA","PLCG1","PML","PPP3CC","PRDM1","PRKAA2","PRKCQ","PTMA","PTPN6","PTPRC","RAP1A","RGS1","RICTOR","RPTOR","RUNX1","S100A4",
               "S100A6","S1PR1","SELL","SLC2A1","SLC2A3","SREBF1","STMN1","TCF7","TIGIT","TLN1","TMSB10","TNFRSF18","TNFRSF4","TNFRSF9","TNFSF8","TRAF3IP2","TRAF5",
               "TXK","ZBTB16")
    Panel1and3<- sort(c(Panel1, Panel3))
  }

  ## KA edit: list of colors with some contrast to use for various graphs/plots
  contrastColorList <- c("deepskyblue2", "firebrick", "springgreen4", "gold", "royalblue",
                         "magenta4", "cyan", "darkseagreen1", "blue3",
                         "orchid", "khaki1", "deepskyblue", "tomato", "lightseagreen",
                         "violetred", "palegreen1", "slateblue", "coral", "darkolivegreen4",
                         "yellow", "mediumorchid", "lightskyblue1", "sienna", "turquoise",
                         "maroon", "paleturquoise1", "mediumpurple4", "darkgoldenrod1", "plum",
                         "seagreen", "lightsteelblue", "deeppink4", "palevioletred", "green4",
                         "navajowhite", "darkorchid", "lightsalmon3", "steelblue", "chocolate",
                         "cyan4", "tan", "hotpink", "darkolivegreen3", "mistyrose",
                         "dodgerblue", "darkorange1", "mediumpurple", "springgreen4", "lightpink3",
                         "olivedrab", "orangered", "cadetblue1", "green",
                         "lightpink", "paleturquoise", "brown", "lightskyblue2", "indianred1",
                         "lightgoldenrod", "yellow4", "darkseagreen", "blue", "pink",
                         "red4", "darkorange4", "indianred4")

  ###SMS edit: assigning the "contrastColorList" object created by KA to the global environment for use in other functions. Namely, the plotTSNE.R function
  assign("contrastColorList", contrastColorList, envir = .GlobalEnv)

  ###SMS edit: Detecting what panel is present in the "ctInput" dataframe.
  # Detection<- any(colnames(ctInput) == "kmeans.cluster")
  # if (Detection == TRUE) {
  #   ColumnNames<- colnames(ctInput)[-c(1:8)]
  # } else {
  #   ColumnNames<- colnames(ctInput)[-c(1:7)]
  # }
  
  #Detection2<- any(colnames(ctInput) == "cohort")
  
  ###SMS edit: Detecting panels. Using the first gene in the Panels outlined above, we are going to scan the column names to see if we find a match.
  # if (Detection2 == FALSE) {
  #   FirstGene<- ColumnNames[1]
  #   LastGene<- ColumnNames[length(ColumnNames)]
  # } else {
  #   FirstGene<- ColumnNames[1]
  #   LastGene<- ColumnNames[(length(ColumnNames)-1)]
  # }
  
  #DetectFirstGene1<- any(Panel1 == FirstGene)
  DetectFirstGene1<- any(colnames(ctInput) == Panel1[1])
  if (Troubleshoot == TRUE) {print(DetectFirstGene1)}
  #DetectFirstGene2<- any(Panel3 == FirstGene)
  DetectFirstGene2<- any(colnames(ctInput) == Panel3[1])
  if (Troubleshoot == TRUE) {print(DetectFirstGene2)}
  
  #DetectLastGene1<- any(Panel1 == LastGene)
  DetectLastGene1<- any(colnames(ctInput) == Panel1[length(Panel1)])
  if (Troubleshoot == TRUE) {print(DetectLastGene1)}
  #DetectLastGene2<- any(Panel3 == LastGene)
  DetectLastGene2<- any(colnames(ctInput) == Panel3[length(Panel3)])
  if (Troubleshoot == TRUE) {print(DetectLastGene2)}
  
  if (DetectFirstGene1 == TRUE && DetectLastGene1 == TRUE && DetectFirstGene2 == FALSE && DetectLastGene2 == FALSE) {
    Panel <- "1"
  } else if (DetectFirstGene1 == TRUE && DetectLastGene2 == TRUE) {
    Panel <- "1 and 3"
  } else if (DetectFirstGene2 == TRUE && DetectLastGene1 == TRUE) {
    Panel<- "1 and 3"
  } else if (DetectFirstGene2 == TRUE && DetectLastGene2 == TRUE) {
    Panel <- "3"
  }
  
  print(paste0("The panel observed in the panel detection series in the 'clusterFilter.R' script is ", Panel))
  
  if (PanelNumberDetected == Panel) {
    print("The panel's found in ctInput are the same as the panel loaded by the User.")
  } else {
    warning(print("Warning! The panel detected and the panel number input by the user are not the same!"))
  }
  
  ###SMS edit, now that we know what panel is loaded into "ctInput," we can then jump into the HLADR test below but only if the "Species" setting is "HU".
  if (Species == "HU") {
    if (Panel == "1" | Panel == "1 and 3") {
      HLADRtest<- grepl("HLA.DRA", colnames(ctInput)[grep("HLA*", colnames(ctInput))])
      print(paste0("HLADR test is: ", HLADRtest))
    }


    if (Panel == "3") {
      print("Panel 3 is detected. No need to check for the 'HLA-DRA' spelling as it's not included in the panel.")
    } else {
      if (HLADRtest == TRUE) {
        names(ctInput)[grep("HLA.DRA", colnames(ctInput))] <- "HLA-DRA"
        print(paste0("HLADR test is ", HLADRtest, " meaning that the HLADRA gene is spelled as 'HLA.DRA' in this sheet. The code changed it back to 'HLA-DRA'."))
      } else {
        print(paste0("HLADR test is ", HLADRtest, " meaning that the HLADRA gene is spelled as 'HLA-DRA' in this sheet. No further action required of the code."))
      }
    }
  }
  if (Troubleshoot == TRUE) {print(colnames(ctInput))}
  
  assign("ctInputpostHLADR", ctInput, envir = .GlobalEnv)
  
  ###LG code below, K-means test begins here
  ###SMS edit, here we're testing to see if the "kmeans.cluster" column is in the 8th position because that means this function has been used the second time.
  kmeanstest<- str_detect(colnames(ctInput)[8], "kmeans.cluster")
  print(paste0("kmeanstest is "))
  print(kmeanstest)
  
  if (kmeanstest == TRUE) {
    Test1<- setequal(Panel1, colnames(ctInput)[9:(ncol(ctInput)-1)])
    print(paste0("Test 1 is ", Test1))
    Test2<- setequal(Panel3, colnames(ctInput)[9:(ncol(ctInput)-1)])
    print(paste0("Test 2 is ", Test2))
    Test1and2<- setequal(Panel1and3, colnames(ctInput)[9:(ncol(ctInput)-1)])
    print(paste0("Test1and2 is ", Test1and2))
    #print(setdiff(Panel1and3, colnames(ctInput)[9:(ncol(ctInput)-1)]))
    #print(colnames(ctInput[,9:(ncol(ctInput)-1)]))
    
  } else {
    Test1<- setequal(Panel1, colnames(ctInput)[8:ncol(ctInput)])
    print(paste0("Test 1 is ", Test1))
    Test2<- setequal(Panel3, colnames(ctInput)[8:ncol(ctInput)])
    print(paste0("Test 2 is ", Test2))
    Test1and2<- setequal(Panel1and3, colnames(ctInput)[8:ncol(ctInput)])
    print(paste0("Test1and2 is ", Test1and2))
    #print(setdiff(Panel1and3, colnames(ctInput)[8:ncol(ctInput)]))
    #print(colnames(ctInput[,8:(ncol(ctInput))]))
  }
  
  # Test1<- setequal(Panel1, colnames(ctInput)[9:ncol(ctInput)])
  # print(Test1)
  # Test2<- setequal(Panel2, colnames(ctInput)[9:ncol(ctInput)])
  # print(Test2)
  # Test1and2<- setequal(Panel1and3, colnames(ctInput)[9:ncol(ctInput)])
  # print(Test1and2)
  # print(setdiff(Panel1and3, colnames(ctInput)[9:ncol(ctInput)]))
  # 
  # print(colnames(ctInput[,8:(ncol(ctInput))]))
  # 
  
  ###SMS edit, trying to automate the subsetting of the ctInput columns
  if (kmeanstest == TRUE) {
    print(paste0("Kmeanstest is ", kmeanstest, ". We need to assesmble a list of column numbers for subsetting the dataframe while also keeping the prior established order."))
    ListofColumns<- c()
    if (Test1 == TRUE) {
     KmeanFirst<- (grep("kmeans.cluster", colnames(ctInput)) + 1)
     CohortLast<- (grep("cohort", colnames(ctInput))-1)
     print(colnames(ctInput))
     print(colnames(ctInput)[KmeanFirst])
     print(colnames(ctInput)[CohortLast])
    } else if (Test2 == TRUE) {
      KmeanFirst<- (grep("kmeans.cluster", colnames(ctInput)) + 1)
      CohortLast<- (grep("cohort", colnames(ctInput))-1)
      print(colnames(ctInput))
      print(colnames(ctInput)[KmeanFirst])
      print(colnames(ctInput)[CohortLast])
    } else if (Test1and2 == TRUE) {
      KmeanFirst<- (grep("kmeans.cluster", colnames(ctInput)) + 1)
      CohortLast<- (grep("cohort", colnames(ctInput))-1)
      print(colnames(ctInput))
      print(colnames(ctInput)[KmeanFirst])
      print(colnames(ctInput)[CohortLast])
    }
  } else {
    if (Test1 == TRUE) {
      TopNumber <- grep(Panel1[1], colnames(ctInput))
      BottomNumber <- grep(Panel1[length(Panel1)], colnames(ctInput))
      print(colnames(ctInput))
      print(TopNumber)
      print(BottomNumber)
      print(colnames(ctInput)[TopNumber:BottomNumber])
    } else if (Test2 == TRUE) {
      TopNumber <- grep(Panel3[1], colnames(ctInput))
      BottomNumber <- grep(Panel3[length(Panel3)], colnames(ctInput))
      print(colnames(ctInput))
      print(TopNumber)
      print(BottomNumber)
      print(colnames(ctInput)[TopNumber:BottomNumber])
    } else if (Test1and2 == TRUE) {
      TopNumber <- grep(Panel1and3[1], colnames(ctInput))
      BottomNumber <- grep(Panel1and3[length(Panel1and3)], colnames(ctInput))
      print(colnames(ctInput))
      print(TopNumber)
      print(BottomNumber)
      print(colnames(ctInput)[TopNumber:BottomNumber])
    }
  }

  ###Louis' original code
  #SMS edit, this includes the first gene depending on the gene panel detected
  if (kmeanstest == TRUE) {
    ctInput$ctSum <- rowSums(ctInput[,c(KmeanFirst:CohortLast)])
  } else {
    ctInput$ctSum <- rowSums(ctInput[,c(TopNumber:BottomNumber)])
  }
  
  assign("ctInput", ctInput, envir = .GlobalEnv)
  print(ctInput)
  #ctInput$ctSum <- rowSums(ctInput[,9:ncol(ctInput)]) #Louis' original, this splits the dataframe at column 9 (B2M) but does not include column 8 (AIM2)
  
  ctInput <- ctInput[order(ctInput$ctSum) ,]
  
  ctInput$ctSum <- NULL
  
  #print(colnames(ctInput)) #SMS edit, to figure out what the columns are before starting the Kmeans
  
  #ctTestK <- ctInput[,9:ncol(ctInput)] #Louis' orignal column numbers, it omitted AIM2
  
  ###SMS edit to detect the panel number and then adjust some of the column names that are going to be subsetted
  if (Species == "HU") {
    if (kmeanstest == TRUE) {
        ctTestK <- ctInput[,c(KmeanFirst:CohortLast)]
      } else {
        if (Test1 == TRUE) {
          FirstColumn<- grep("AIM2", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          LastColumn<- grep("ZEB2", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          ctTestK <- ctInput[,c(FirstColumn:LastColumn)]
        } else if (Test2 == TRUE) {
          FirstColumn<- grep("AKT1", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          LastColumn<- grep("ZBTB16", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          ctTestK <- ctInput[,c(FirstColumn:LastColumn)]
        } else if (Test1and2 == TRUE) {
          FirstColumn<- grep("AIM2", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          LastColumn<- grep("ZEB2", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          ctTestK <- ctInput[,c(FirstColumn:LastColumn)] 
        }
      }
    } else if (Species == "MS") {
      if (kmeanstest == TRUE) {
        ctTestK <- ctInput[,c(KmeanFirst:CohortLast)]
      } else {
        if (Test1 == TRUE) {
          FirstColumn<- grep("Aim2", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          LastColumn<- grep("Zeb2", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          ctTestK <- ctInput[,c(FirstColumn:LastColumn)]
          print(paste("Test1 is TRUE: columns to be sent for kmeans testing: "))
          print(colnames(ctTestK))
        } else if (Test2 == TRUE) {
          FirstColumn<- grep("ACTA2", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          LastColumn<- grep("ZAP70", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          ctTestK <- ctInput[,c(FirstColumn:LastColumn)]
          print(paste("Test2 is TRUE: columns to be sent for kmeans testing: "))
          print(colnames(ctTestK))
        } else if (Test1and2 == TRUE) {
          FirstColumn<- grep("ACTA2", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          LastColumn<- grep("Zeb2", colnames(ctInput)) #[8:ncol(ctInput)]) ###SMS - so using the subsetted columns leads to a wrong numbering for the colnames. So you're pulling the wrong column numbers.
          ctTestK <- ctInput[,c(FirstColumn:LastColumn)]
          print(paste("Test1and2 is TRUE: columns to be sent for kmeans testing: "))
          print(colnames(ctTestK))
        }
      }
    }
  #SMS edit, this includes AIM2
  #ctTestK <- ctInput[,8:(ncol(ctInput))]
  
  if (Troubleshoot == TRUE) {
    print("Column names after searching for the column pattern and after selecting the right columns. The following will be fed into the kmeans function.")
    print(colnames(ctTestK))
    #print(ctTestK)
    #print(paste("ctInput-line24; number of rows in ctInput;", nrow(ctTestK), sep=" "))
  }
  
  ## set seed
  set.seed(1)
  
  if(testK == T){
    ## run kmeans with 1-15 centers
    wss <- (nrow(ctTestK)-1)*sum(apply(ctTestK,2,var))
    for(i in 2:15){ 
      wss[i] <- sum(kmeans(ctTestK, centers=i)$withinss)
    }
    
    ## create within sum of squares dataframe and plot
    centers <- 1:15
    wssc <- as.data.frame(cbind(centers, wss))
    wssPlot <- ggplot(wssc, aes(centers, wss)) +
      geom_point(colour = "blue", fill = "cyan", size = 5, shape=21) +
      geom_path(size = 2.3) +
      ggtitle("kmeans scree plot") +
      ylab("total within-cluster sum of squares") +
      scale_x_continuous(breaks=1:15) +
      theme(text=element_text(size=25),
            axis.title.x=element_text(vjust=-0.5),
            plot.margin=unit(c(10,5,10,5),"cm"))
    print(wssPlot)
  }
  
  ## run kmeans with chosen number of centers
  set.seed(1)
  normFit <- kmeans(ctTestK, numCenters)
  
  ## attach results to dataframe, SMS edit: LG is adding the cluster ID's to the ctInput data frame
  ctClust <- data.frame(ctInput, normFit$cluster)
  
  if (Troubleshoot == TRUE) {
    print("Column Names for ctClust are: ")
    print(colnames(ctClust))
  }
  
  ## rename clustering results column
  names(ctClust)[which(names(ctClust)=="normFit.cluster")] <- "kmeans.cluster"
  
  ## move kmeans column to column 9 to keep with id variables, SMS Edit: changing from column 9 to column 8 for the kmeans column.
  #ctClust <- ctClust[, c(1:8, ncol(ctClust), 9:(ncol(ctClust)-1))]
  ctClust <- ctClust[, c(1:7, ncol(ctClust), 8:(ncol(ctClust)-1))]
  
  if (Troubleshoot == TRUE) {
    print("Column Numbers for ctClust after moving around the columns:")
    print(colnames(ctClust))
  }
  
  ##SMS edit, to prune out cohorts, tissues etc you are not interested in, PUT THE KMEANS TEST HERE!!!!!
  ctClust[grep("UCM", ctClust$probe, fixed=TRUE), "cohort"] <- "Adult"
  ctClust[grep("NBD", ctClust$probe, fixed=TRUE), "cohort"] <- "NBD"
  ctClust[grep("RUS", ctClust$probe, fixed=TRUE), "cohort"] <- "PED"
  ctClust[grep("OXF", ctClust$probe, fixed=TRUE), "cohort"] <- "OXFORD"
  ctClust[grep("GI4330", ctClust$probe, fixed=TRUE), "cohort"] <- "OXFORD-NBD"
  ctClust[grep("BigFoot", ctClust$probe, fixed=TRUE), "cohort"] <- "BigFoot"
  ctClust[grep("Astrios", ctClust$probe, fixed=TRUE), "cohort"] <- "Astrios"
  ctClust[grep("G", ctClust$probe, fixed=TRUE), "cohort"] <- "Glycopeptide"
  ctClust[grep("N", ctClust$probe, fixed=TRUE), "cohort"] <- "Naive"
  
  write.csv(ctClust, paste0(baseDir, "results/clusterFilter-002-ctInput.csv"))
  # ctClust<- subset(ctClust, ctClust$cohort != "OXFORD")
  # ctClust<- subset(ctClust, ctClust$cohort != "OXFORD-NBD")
  # ctClust<- subset(ctClust, ctClust$cellSource != "blood")
  # ctClust<- subset(ctClust, ctClust$age == "Spleen")
  write.csv(ctClust, paste0(baseDir, "results/clusterFilter-003.csv"))
  
  #### heatmap and cluster membership ####
  if(plotHeatmap == T){
    par(mar=c(1,1,1,1))
    
    ## account for kmeans.cluster column when in use
    if(heatmapFactor != "kmeans.cluster") {
      idCols <- as.numeric(grep("kmeans.cluster", colnames(ctClust))) #as.numeric(8) ###SMS edit, I believe Louis is trying to isolate the "kmeans.cluster" column
      ctClust <- subset(ctClust, select = -kmeans.cluster)
    } else{
      LenghtofKmeans<- grep("kmeans.cluster", colnames(ctClust))
      if (Troubleshoot == TRUE) {print(paste0("The values in lenghtofkmeans is: ", LenghtofKmeans))
        print(paste0("The length of lenghtofkmeans object is ", length(LenghtofKmeans)))}
      if (length(LenghtofKmeans) > 1) {
        idCols <- as.numeric(grep("kmeans.cluster", colnames(ctClust))[2]+1)
      } else {
        idCols <- as.numeric(grep("kmeans.cluster", colnames(ctClust))+1) #as.numeric(9) ###SMS edit, I believe Louis is trying to isolate the "kmeans.cluster" column
      }
    }
    
    ##SMS Edit, trying to figure out which column Louis is pulling for the Kruskal wallace rank sum test below##
    if (heatmapFactor != "kmeans.cluster") {
      if (Troubleshoot == TRUE) {print(paste0("When heatmapfactor is set to 'no kmeans.cluster', the first column being pulled is ", colnames(ctClust)[idCols]))}
      } else (
        if (Troubleshoot == TRUE) {print(paste0("When heatmapfactor is set to 'kmeans.cluster', the first column being pulled is ", colnames(ctClust)[idCols]))}
      )
    
    ## get p-values for differential expression for factors, SMS edit- added the "-1" to remove the last column from the kruskal wallace test ##
    print(paste0("Value laoded into idCols: ", idCols, " which corresponds to column ", colnames(ctClust)[idCols], ". The last column name is: ", colnames(ctClust)[ncol(ctClust)], ". The second to last column is ", colnames(ctClust)[ncol(ctClust)-1]))  ###SMS edit, I think LG is trying to exclude the kmeans column for the kruskal wallace test. So to confirm what idCols I'm printing the value.
    
    if (Troubleshoot == TRUE) {print(colnames(ctClust))}
    
    pvals <- NULL
    
    for(i in (idCols):(ncol(ctClust)-1)) {. ###SMS comment: HERE LG CODED A "+1", WHILE THIS MAKES SENSE INITIALLY, "IDCOLS" IS SET ABOVE WITH "+1" ADDED. SO YOU'RE DROPPING AIM2 HERE.
      pvals <- c(pvals, kruskal.test(ctClust[,i], factor(ctClust[, heatmapFactor]))$p.value)
    }
    
    if (Troubleshoot == TRUE) { 
      print(order(pvals))
      print(paste0("Length of pvals is ", length(pvals)))
    }
    
    ctClust <- ctClust[,c(1:(idCols-1), (order(pvals)+idCols), ncol(ctClust))] ##SMS edit, adding back the "cohort" column, adjusted the idCols function here to have a "-1"
    
    write.csv(ctClust, paste0(baseDir, "results/clusterFilter-004.csv"))
    
    AnnoColors <- NULL
    ## label cells by age
    # if("age" %in% heatmapColorBy){
    #   ctClust$ageColor <- NA
    #   ctClust[which(ctClust$age=="6"), "ageColor"] <- "deepskyblue2"
    #   ctClust[which(ctClust$age=="12"), "ageColor"] <- "peachpuff2"
    #   ctClust[which(ctClust$age=="infl"), "ageColor"] <- "deepskyblue2"
    #   ctClust[which(ctClust$age=="uninfl"), "ageColor"] <- "navy"
    #   ctClust[which(ctClust$age=="blood"), "ageColor"] <- "orangered"
    #   ctClust[which(ctClust$age=="Spleen"), "ageColor"] <- "cyan4"
    #   Ages <- ctClust$ageColor
    #   ctClust <- subset(ctClust, select=-c(ageColor))
    #   AnnoColors <- cbind(AnnoColors, Ages)
    # }
    
    ## KA edit:
    uniqueAges <- unique(ctClust[,"age"]) # make a list of the unique probe labels in the list
    colorsList <-contrastColorList[1:length(uniqueAges)]
    contrastColorList <- contrastColorList[-c(1:length(uniqueAges))] # removes the colors that are used for this matching
    ageColorKey <- cbind(uniqueAges, colorsList)
    if("age" %in% heatmapColorBy) {
      for (row in 1:nrow(ctClust)) {
        ctClust[row, "ageColor"] <- ageColorKey[which(ageColorKey[,"uniqueAges"]==ctClust[row, "age"]), "colorsList" ]
      }
      Ages <- ctClust$ageColor
      ctClust <- subset(ctClust, select=-c(ageColor))
      AnnoColors <- cbind(AnnoColors, Ages)
    }
    
    ## label cells by tissue
    # if("source" %in% heatmapColorBy){
    #   ctClust$sourceColor <- NA
    #   ctClust[which(ctClust$cellSource=="LN"), "sourceColor"] <- "blue3"
    #   ctClust[which(ctClust$cellSource=="I"), "sourceColor"] <- "orangered"
    #   ctClust[which(ctClust$cellSource=="S"), "sourceColor"] <- "turquoise3"
    #   ctClust[which(ctClust$cellSource=="PB"), "sourceColor"] <- "gold"
    #   ctClust[which(ctClust$cellSource=="UCM"), "sourceColor"] <- "turquoise3"
    #   ctClust[which(ctClust$cellSource=="NBD"), "sourceColor"] <- "gold"
    #   ctClust[which(ctClust$cellSource=="blood"), "sourceColor"] <- "orangered"
    #   ctClust[which(ctClust$cellSource=="tissue"), "sourceColor"] <- "blue"
    #   ctClust[which(ctClust$cellSource=="RUS"), "sourceColor"] <- "cyan4"
    #   Tissues <- ctClust$sourceColor
    #   ctClust <- subset(ctClust, select=-c(sourceColor))
    #   AnnoColors <- cbind(AnnoColors, Tissues)
    # }
    
    ## KA edit:
    uniqueSources <- unique(ctClust[,"cellSource"]) # make a list of the unique probe labels in the list
    colorsList <-contrastColorList[1:length(uniqueSources)]
    contrastColorList <- contrastColorList[-c(1:length(uniqueSources))] # removes the colors that are used for this matching
    sourceColorKey <- cbind(uniqueSources, colorsList)
    if("source" %in% heatmapColorBy) {
      for (row in 1:nrow(ctClust)) {
        ctClust[row, "sourceColor"] <- sourceColorKey[which(sourceColorKey[,"uniqueSources"]==ctClust[row, "cellSource"]), "colorsList" ]
      }
      Sources <- ctClust$sourceColor
      ctClust <- subset(ctClust, select=-c(sourceColor))
      AnnoColors <- cbind(AnnoColors, Sources)
    }
    
    ## label cells by probe
    # if("probe" %in% heatmapColorBy){
    #   ctClust$probeColor <- NA
    #   ctClust[which(ctClust$probe=="Obsc"), "probeColor"] <- "#FB9A99"
    #   ctClust[which(ctClust$probe=="12"), "probeColor"] <- "#A6CEE3"
    #   ctClust[which(ctClust$probe=="13"), "probeColor"] <- "#1F78B4"
    #   ctClust[which(ctClust$probe=="9D"), "probeColor"] <- "#B2DF8A"
    #   ctClust[which(ctClust$probe=="9Q"), "probeColor"] <- "#33A02C"
    #   ctClust[which(ctClust$probe=="IPDRp"), "probeColor"] <- "blue3"
    #   ctClust[which(ctClust$probe=="OXFGI4330"), "probeColor"] <- "orangered"
    #   ctClust[which(ctClust$probe=="UCM5"), "probeColor"] <- "#A6CEE3"
    #   ctClust[which(ctClust$probe=="UCM6"), "probeColor"] <- "#1F78B4"
    #   ctClust[which(ctClust$probe=="UCM9"), "probeColor"] <- "blue3"
    #   ctClust[which(ctClust$probe=="UCM10"), "probeColor"] <- "turquoise"
    #   ctClust[which(ctClust$probe=="UCM12"), "probeColor"] <- "navy"
    #   ctClust[which(ctClust$probe=="UCM13"), "probeColor"] <- "maroon"
    #   ctClust[which(ctClust$probe=="UCM14"), "probeColor"] <- "darkmagenta"
    #   ctClust[which(ctClust$probe=="UCM15"), "probeColor"] <- "thistle3"
    #   ctClust[which(ctClust$probe=="UCM16"), "probeColor"] <- "peru"
    #   ctClust[which(ctClust$probe=="UCM17"), "probeColor"] <- "purple4"
    #   ctClust[which(ctClust$probe=="NBD1"), "probeColor"] <- "#B2DF8A"
    #   ctClust[which(ctClust$probe=="NBD3"), "probeColor"] <- "#33A02C"
    #   ctClust[which(ctClust$probe=="NBD4"), "probeColor"] <- "yellow2"
    #   ctClust[which(ctClust$probe=="OXFGI4870"), "probeColor"] <- "plum4"
    #   ctClust[which(ctClust$probe=="RUS002"), "probeColor"] <- "#003333"
    #   ctClust[which(ctClust$probe=="RUS008"), "probeColor"] <- "#9933FF"
    #   ctClust[which(ctClust$probe=="OXFIBD1322"), "probeColor"] <- "violetred4"
    #   ctClust[which(ctClust$probe=="RUS007"), "probeColor"] <- "indianred1"
    #   ctClust[which(ctClust$probe=="RUS012"), "probeColor"] <- "cyan4"
    #   ctClust[which(ctClust$probe=="RUS013"), "probeColor"] <- "brown2"
    #   Probes <- ctClust$probeColor
    #   ctClust <- subset(ctClust, select=-c(probeColor))
    #   AnnoColors <- cbind(AnnoColors, Probes)
    # }
    
    ##KA edit: more efficient labeling cells by probe:
    uniqueProbes <- unique(ctClust[,"probe"]) # make a list of the unique probe labels in the list
    colorsList <-contrastColorList[1:length(uniqueProbes)]
    contrastColorList <- contrastColorList[-c(1:length(uniqueProbes))] # removes the colors that are used for this matching
    probeColorKey <- cbind(uniqueProbes, colorsList)
    if("probe" %in% heatmapColorBy) {
      for (row in 1:nrow(ctClust)) {
        ctClust[row, "probeColor"] <- probeColorKey[which(probeColorKey[,"uniqueProbes"]==ctClust[row, "probe"]), "colorsList" ]
      }

      Probes <- ctClust$probeColor
      ctClust <- subset(ctClust, select=-c(probeColor))
      AnnoColors <- cbind(AnnoColors, Probes)
    }
    
    ## label cells by cluster
    if(heatmapFactor == "kmeans.cluster"){
      
      ## for heatmap figure
      # clusterColorKey <- brewer.pal(8, "Dark2")
      # extraColorPalette <- clusterColorKey[c(2,3,4,6)]
      
      extraColorPalette <- brewer.pal(12, "Set3")
      clusterVals <- unique(ctClust$kmeans.cluster)
      clusterVals <- clusterVals[order(clusterVals)]
      
      ctClust$kmeansColor <- NULL
      for(i in 1:length(clusterVals)){
        ctClust[which(ctClust$kmeans.cluster == clusterVals[i]), 
                "kmeansColor"] <- extraColorPalette[i]
      }
      
      Kmeans.clusters <- ctClust$kmeansColor
      ctClust$kmeansColor <- NULL
      AnnoColors <- cbind(AnnoColors, Kmeans.clusters)
    }
    
    
    ## KA edit: commented out ageColorKey
    # ageColorKey <- rbind(c("6", "deepskyblue2"),
    #                      c("12", "peachpuff2"),
    #                      c("infl", "deepskyblue2"),
    #                      c("uninfl", "navy"),
    #                      c("blood", "orangered"),
    #                      c("Spleen", "#FF0000"))
    #ageColorKey <- ageColorKey[which(ageColorKey[,1] %in% ctClust$age),]
    
    ##SMS edit, if you only want "tissue" then the "which" function below turns the row into a vector, because it's a vector, Louis' tablefunction fails
    ###as it is expecting a two column matrix but instead is fed a vector. To ameliorate this, I've writen an "if" statement to maintain a matrix if the which function
    ####finds only one match
    
    if ("age" %in% heatmapColorBy) {
      LengthofDownSample <- which(ageColorKey[,1] %in% ctClust$age, arr.ind = TRUE)
      
      if (length(LengthofDownSample) == 1) {
        ageColorKey <- matrix(ageColorKey[LengthofDownSample,], nrow = 1, byrow = TRUE)
      } else {
        ageColorKey <- ageColorKey[LengthofDownSample, , drop = FALSE]
      }
      
      print(ageColorKey)
    }
    
    
    ## KA edit: commented out sourceColorKey
    # sourceColorKey <- rbind(c("I", "orangered"),
    #                         c("LN", "#2B65EC"),
    #                         c("S", "turquoise3"),
    #                         c("PB", "gold"),
    #                         c("UCM", "turquoise3"),
    #                         c("NBD", "gold"),
    #                         c("blood", "orangered"),
    #                         c("tissue", "blue"),
    #                         c("RUS", "cyan4"))
    # write.csv(sourceColorKey, paste0(baseDir, "results/clusterFilter-006.csv"))
    
    ##SMS edit, if you only want "tissue" then the "which" function below turns the row into a vector, because it's a vector, Louis' tablefunction fails
    ###as it is expecting a two column matrix but instead is fed a vector. To ameliorate this, I've writen an "if" statement to maintain a matrix if the which function
    ####finds only one match
    
    LengthofDownSample<- which(sourceColorKey[,1] %in% ctClust$cellSource, arr.ind = TRUE)
    write.csv(LengthofDownSample, paste0(baseDir, "results/clusterFilter-007.csv"))
    
    if(length(LengthofDownSample) == 1) {
      #print("ONE!!!!")
      print(sourceColorKey[LengthofDownSample,])
      sourceColorKey<- matrix(sourceColorKey[LengthofDownSample,], nrow = 1, byrow=TRUE)
    } else {
      #print("MORE THAN ONE!!!!")
      sourceColorKey <- sourceColorKey[which(sourceColorKey[,1] %in% ctClust$cellSource),]
    }
    
    #sourceColorKey <- sourceColorKey[which(sourceColorKey[,1] %in% ctClust$cellSource, arr.ind = TRUE),]
    
    write.csv(ctClust$cellSource, paste0(baseDir, "results/clusterFilter-008.csv"))
    #print(sourceColorKey) ##SMS Edit, trying to figure out what LG is doing.
    write.csv(sourceColorKey, paste0(baseDir, "results/clusterFilter-009.csv"))
    
    ## KA edit: automated above, so I commented this section out
    ###SMS Comment: Automate this. Pick 50 colors and have them automatically assigned to unique probes.
    # probeColorKey <- rbind(c("PI", "orange2"),
    #                        c("Obsc", "#FB9A99"),
    #                        c("12", "#A6CEE3"),
    #                        c("13", "#1F78B4"),
    #                        c("9D", "#B2DF8A"),
    #                        c("9Q", "#33A02C"),
    #                        c("IPDRp", "blue3"),
    #                        c("OXFGI4330", "orangered"),
    #                        c("UCM5", "#A6CEE3"),
    #                        c("UCM6", "#1F78B4"),
    #                        c("UCM9", "blue3"),
    #                        c("UCM10", "turquoise"),
    #                        c("UCM12", "navy"),
    #                        c("UCM13", "maroon"),
    #                        c("UCM14", "darkmagenta"),
    #                        c("UCM15", "thistle3"),
    #                        c("UCM16", "peru"),
    #                        c("UCM17", "purple4"),
    #                        c("NBD1", "#B2DF8A"),
    #                        c("NBD3", "#33A02C"),
    #                        c("NBD4", "yellow2"),
    #                        c("OXFGI4870", "plum4"),
    #                        c("RUS002", "#003333"),
    #                        c("RUS008", "#9933FF"),
    #                        c("OXFIBD1322", "violetred4"),
    #                        c("RUS007", "indianred1"),
    #                        c("RUS012", "cyan4"),
    #                        c("RUS013", "brown2"))
    # probeColorKey <- probeColorKey[which(probeColorKey[,1] %in% ctClust$probe),]
    # print(probeColorKey) ##SMS Edit
    
    insetX <- 0.256
    
      
#       clusterVals <- unique(ctClust$kmeans.cluster)
#       clusterVals <- clusterVals[order(clusterVals)]
#       
#       if(length(clusterVals) == 6){
#         extremeClusterVals <- c(-15, -9, -4, 4, 9, 15)
#       } else{
#         extremeClusterVals <- seq(-15, 15, (30/(length(clusterVals) - 1)))
#       }
#         
#       ctClust$extremeKmeans <- NULL
#       for(i in 1:length(clusterVals)){
#         ctClust[which(ctClust$kmeans.cluster == clusterVals[i]), 
#                 "extremeKmeans"] <- extremeClusterVals[i]
#       }
#       
#       kmeans.cluster <- ctClust$kmeans.cluster
#       ctClust$kmeans.cluster <- ctClust$extremeKmeans
#       ctClust$extremeKmeans <- NULL
#     }
    
    ##SMS Edit; trying to find which columns are not numeric##
    # for (i in 1:ncol(ctClust)) {
    #   x<- is.numeric(ctClust[,i])
    #   if (x == TRUE) {
    #     print(paste0("Column ", i, " is numeric"))
    #   } else {
    #     print(paste0("Column ", i,"-",colnames(ctClust)[i], " is NOT numeric"))
    #   }
    # }
    # 
    
    ## plot heatmap
    par(cex.main = 1.5)
#     heatmap.2(as.matrix(ctClust[(idCols + 1):ncol(ctClust)]), Colv=F, 
#               dendrogram="row", trace="none", 
#               col=colorpanel(15, "red", "white", "blue"),
#               RowSideColors = extraColors,
#               colRow=rowColors, adjRow = heatmapRowAdjust,
#               adjCol = c(1, NA), cexCol = 1.3, cexRow=heatmapRowSize, margins=c(10,10), 
#               density.info="density", denscol="black", key.title=NA, 
#               main=paste("heatmap for", heatmapTissueLabel, sep=" "), 
#               xlab="log2(gene expression / Gapdh expression)", 
#               ylab=heatmapYLabel)
    if(plotClustOnly == FALSE){
      print(colnames(ctClust))
      print(paste0("The value in idCols is ", idCols, " and the first column for the heatmap is ", colnames(ctClust[idCols]), " while the last column is ", colnames(ctClust)[ncol(ctClust)-1]))
      print(as.matrix(ctClust[(idCols):(ncol(ctClust)-1)]))
      heatmap.3(as.matrix(ctClust[(idCols):(ncol(ctClust)-1)]), Colv=F, 
                dendrogram="row", trace="none", srtCol = 75,
                col=colorpanel(15, "blue", "white", "red"),
                RowSideColors = AnnoColors,
                labRow=NA,
                adjCol = c(1, NA), cexCol = 1.00000000046, margins=c(10,10), ###SMS Edit: changed the cexCol value from 1.3 to what it is, to accomodate the text on the heatmap. '1.3' works for 96 genes but is too crowded for 192 genes.
                density.info="density", denscol="black", key.title=NA, 
                main=paste("heatmap for", heatmapTissueLabel, sep=" "), 
                xlab="log2(gene expression / Gapdh expression)", 
                ylab="Individual Cells\n\n\n")
      if("age" %in% heatmapColorBy){
        legend("topleft",
               legend=ageColorKey[,1],
               fill=ageColorKey[,2], 
               title="age",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.15))
        insetX <- insetX + 0.06
      }
      if("source" %in% heatmapColorBy){
        legend("topleft",
               legend=sourceColorKey[,1],
               fill=sourceColorKey[,2], 
               title="tissue",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.15))
        insetX <- insetX + 0.06
      }
      if("probe" %in% heatmapColorBy){
        legend("topleft",
               legend=probeColorKey[,1],
               fill=probeColorKey[,2], 
               title="probe",
             border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.15))
      insetX <- insetX + 0.075
        #         border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.1))
        # insetX <- insetX + 0.08
      }
      if(heatmapFactor == "kmeans.cluster"){
        legend("topleft",
               legend=clusterVals,
               fill=extraColorPalette[1:length(clusterVals)], 
               title="cluster",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(insetX,0.15))
      }
    }
      
      if(plotClustOnly == TRUE){
        heatmap.2(as.matrix(ctClust[(idCols + 1):(ncol(ctClust)-1)]), Colv=F, 
                  dendrogram="row", trace="none", 
                  col=colorpanel(15, "blue", "white", "red"),
                  RowSideColors = Kmeans.clusters,
                  labRow=NA,
                  adjCol = c(1, NA), cexCol = 1.3, margins=c(10,10), 
                  density.info="density", denscol="black", key.title=NA, 
                  main=paste("heatmap for", heatmapTissueLabel, sep=" "), 
                  xlab="log2(gene expression / Gapdh expression)", 
                  ylab="Individual Cells\n\n\n")
        legend("topleft",
               legend=clusterVals,
               fill=extraColorPalette[1:length(clusterVals)], 
               title="cluster",
               border=FALSE, bty="n", y.intersp = 0.7, cex=1.85, inset=c(0.26,0.15))
      }
    }
    par(cex.main = 1.2)
    
    
    ## remap clusters values to original
#     if(heatmapFactor == "kmeans.cluster"){
#       ctClust$kmeans.cluster <- kmeans.cluster
#     }
    
    
    #### fisher's exact test for all covariates ####
    if(heatmapFactor == "kmeans.cluster"){
      ## get cluster numbers
      clusterVals <- unique(ctClust$kmeans.cluster)
      clusterVals <- clusterVals[order(clusterVals)]
      
      ## probe
      if("probe" %in% fisherTests){
        probes <- unique(ctClust$probe)
        # probes <- probes[match(c("PI", "12", "13", "9D", "9Q", "Obsc", "IPDRp", "DRn", "UCM5", "UCM6", "UCM9", "UCM10", "UCM12", "UCM13", "UCM14", "UCM15", "UCM16", "UCM17", "NBD1", "NBD3", "NBD4", "RUS002", 
                                 # "RUS008", "OXFGI4330", "OXFGI4870", "OXFIBD1322", "RUS007", "RUS012", "RUS013"), probes, nomatch=F)]
        probeTable <- matrix(nrow=length(clusterVals), ncol=length(probes))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(probes)){
            probeTable[i, j] <- nrow(subset(ctClust, probe==probes[j] & 
                                                      kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(probeTable) <- paste("probe_", probes, sep="")
        rownames(probeTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("Probe vs. Cluster", quote=F)
        # print(probeTable)
        print(probeTable, quote = F)
        write.csv(probeTable, paste0(baseDir, "results/clusterFilter-010.csv")) ##SMS Edit
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(probeTable, simulate.p.value=T))
        
        ##SMS edit
        write.csv(probeTable, paste0(baseDir, "results/clusterFilter-011.csv"))
        write.csv(probeColorKey, paste0(baseDir, "results/clusterFilter-012.csv"))
        
        ## Bar Plots
        probeTablePlots <- tableBarPlots(plotFactor = "probe", plotFactorVals = probes, clusterVals = clusterVals, factorTable = probeTable, colorKey = probeColorKey)
      }
      
      ## tissue
      if("tissue" %in% fisherTests){
        cellSources <- unique(ctClust$cellSource)
        # cellSources <- cellSources[match(c("I", "LN", "S", "PB", "UCM", "NBD", "blood", "tissue"), cellSources, nomatch=F)]
        sourceTable <- matrix(nrow=length(clusterVals), ncol=length(cellSources))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(cellSources)){
            sourceTable[i, j] <- nrow(subset(ctClust, cellSource==cellSources[j] & 
                                               kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(sourceTable) <- paste("cellSource_", cellSources, sep="")
        rownames(sourceTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("Tissue vs. Cluster", quote=F)
        print(sourceTable)
        
        #SMS edit, printing sourceTable because the number of factors are not in alignment with the number of colors. Found the problem; the sourceColorKey
        ##does not have a column 2, it's all in one column....idk why that is....I'll print all the tables feeding into the tableBarPlots function to see how that's done
        write.csv(sourceTable, paste0(baseDir, "results/clusterFilter-013.csv"))
        write.csv(sourceColorKey, paste0(baseDir, "results/clusterFilter-014.csv"))
        
        # print(fisher.test(factor(ctClust$cellSource), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(sourceTable, simulate.p.value=T))
        
        ## Bar Plots
        tissueTablePlots <- tableBarPlots(plotFactor = "tissue", plotFactorVals = cellSources, clusterVals = clusterVals, factorTable = sourceTable, colorKey = sourceColorKey)
      }
      ## patient
      if("patient" %in% fisherTests){
        patients <- unique(ctClust$patient)
        patients <- patients[order(patients)]   # or keep as-is
        
        patientTable <- matrix(nrow = length(clusterVals),
                               ncol = length(patients))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(patients)){
            patientTable[i, j] <- nrow(subset(ctClust,
                                              patient == patients[j] &
                                                kmeans.cluster == clusterVals[i]))
          }
        }
        
        colnames(patientTable) <- paste("patient_", patients, sep = "")
        rownames(patientTable) <- paste("cluster_", clusterVals, sep = "")
        
        print("     ", quote = FALSE)
        print("     ", quote = FALSE)
        print("patient vs. Cluster", quote = FALSE)
        print(patientTable)
        print(chisq.test(patientTable, simulate.p.value = TRUE))
        
        ## simple color key for patients
        patientColors   <- brewer.pal(12, "Set3")[1:length(patients)]
        patientColorKey <- cbind(patients, patientColors)
        
        ## Bar plots/heatmap (same structure as probe/age)
        patientTablePlots <- tableBarPlots(
          plotFactor   = "patient",
          plotFactorVals = patients,
          clusterVals  = clusterVals,
          factorTable  = patientTable,
          colorKey     = patientColorKey
        )
        
        ## optionally draw them on their own
        grid.arrange(grobs = list(patientTablePlots[[3]],
                                  patientTablePlots[[4]]),
                     nrow = 1)
        grid.arrange(grobs = list(patientTablePlots[[1]],
                                  patientTablePlots[[2]]),
                     nrow = 1)
      }
      
      ## age
      if("age" %in% fisherTests){
        ages <- unique(ctClust$age)
        ages <- ages[order(as.numeric(ages))]
        ageTable <- matrix(nrow=length(clusterVals), ncol=length(ages))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(ages)){
            ageTable[i, j] <- nrow(subset(ctClust, age==ages[j] & 
                                                kmeans.cluster==clusterVals[i]))
          }
        }
        colnames(ageTable) <- paste("age_", ages, sep="")
        rownames(ageTable) <- paste("cluster_", clusterVals, sep="")
        print("     ", quote=F)
        print("     ", quote=F)
        print("age vs. Cluster", quote=F)
        print(ageTable)
        # print(fisher.test(factor(ctClust$age), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(ageTable, simulate.p.value=T))
        
        ## Bar Plots
        ageTablePlots <- tableBarPlots(plotFactor = "age", plotFactorVals = ages, clusterVals = clusterVals, factorTable = ageTable, colorKey = ageColorKey)
        
#         grid.arrange(grobs = list(probeTablePlots[[3]], probeTablePlots[[4]], tissueTablePlots[[3]], tissueTablePlots[[4]], ageTablePlots[[3]], ageTablePlots[[4]]), nrow=3)
#         
#         grid.arrange(grobs = list(probeTablePlots[[1]], probeTablePlots[[2]], tissueTablePlots[[1]], tissueTablePlots[[2]], ageTablePlots[[1]], ageTablePlots[[2]]), nrow=3)
      }
      ng <- nullGrob()
      if("probe" %in% fisherTests){}
#       if("age" %in% fisherTests){}
#       if("tissue" %in% fisherTests){}
      
      if("probe" %in% fisherTests & "tissue" %in% fisherTests & "age" %in% fisherTests){
        grid.arrange(grobs = list(probeTablePlots[[3]], probeTablePlots[[4]], tissueTablePlots[[3]], tissueTablePlots[[4]], ageTablePlots[[3]], ageTablePlots[[4]]), nrow=3)
        
        grid.arrange(grobs = list(probeTablePlots[[1]], probeTablePlots[[2]], tissueTablePlots[[1]], tissueTablePlots[[2]], ageTablePlots[[1]], ageTablePlots[[2]]), nrow=3)
      }
      
      if("probe" %in% fisherTests & "age" %in% fisherTests){
        grid.arrange(grobs = list(probeTablePlots[[3]], probeTablePlots[[4]], ageTablePlots[[3]], ageTablePlots[[4]], ng, ng), nrow=3)
        
        grid.arrange(grobs = list(probeTablePlots[[1]], probeTablePlots[[2]], ageTablePlots[[1]], ageTablePlots[[2]], ng, ng), nrow=3)
      }
      
      ## probe and age
      if("probe.age" %in% fisherTests){
        ctClust$probe.age <- paste(ctClust$probe, ctClust$age, sep = ".")
        probes <- unique(ctClust$probe)
        # probes <- probes[match(c("PI", "12", "13", "9D", "9Q", "Obsc", "UCM9", "UCM10", "UCM12", "UCM13", "UCM14", "UCM15", "UCM16", "UCM17", "NBD1", "NBD3", "NBD4", "RUS002", "RUS008","OXFGI4330", "OXFGI4870", "OXFIBD1322", "RUS007", "BigFoot", "Astrios"), 
                               # probes, nomatch=F)]
        ages <- unique(ctClust$age)
        ages <- ages[order(as.numeric(ages))]
        eg <- expand.grid(probes, ages)
        probes.ages <- sprintf('%s.%s', eg[,1], eg[,2])
        
        probe.ageTable <- matrix(nrow=length(clusterVals), ncol=length(probes.ages))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(probes.ages)){
            probe.ageTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                          probe.age==probes.ages[j]))
          }
        }
        colnames(probe.ageTable) <- probes.ages
        rownames(probe.ageTable) <- paste("cluster_", clusterVals, sep="")
        probe.ageTable <- probe.ageTable[, which(colSums(probe.ageTable) > 0)]
        probes.ages <- colnames(probe.ageTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Probe.Age vs. Cluster", quote=F)
        # print(probeTable)
        print(probe.ageTable, quote = F)
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(probe.ageTable, simulate.p.value=T))
        
        write.csv(probe.ageTable, paste0(baseDir, "results/clusterFilter-015.csv"))
        write.csv(probe.ageColorKey, paste0(baseDir, "results/clusterFilter-016.csv"))
        
        ## Bar Plots
        probe.ageColors <- brewer.pal(12, "Paired")
        probe.ageColorKey <- cbind(probes.ages, probe.ageColors[1:length(probes.ages)])
        
        probe.ageTablePlots <- tableBarPlots(plotFactor = "probe.age", plotFactorVals = probes.ages, clusterVals = clusterVals, factorTable = probe.ageTable, colorKey = probe.ageColorKey)
        
        grid.arrange(grobs = list(probe.ageTablePlots[[1]], probe.ageTablePlots[[2]],  ng, ng, ng, ng), nrow=3)
        
        ctClust <- subset(ctClust, select = -probe.age)
      }
      
      ## probe and tissue
      if("probe.tissue" %in% fisherTests){
        ctClust$probe.tissue <- paste(ctClust$probe, ctClust$cellSource, sep = ".")
        probes <- unique(ctClust$probe)
        # probes <- probes[match(c("PI", "12", "13", "9D", "9Q", "Obsc"), probes, nomatch=F)]
        tissues <- unique(ctClust$cellSource)
        # tissues <- tissues[match(c("I", "LN", "S", "PB"), tissues, nomatch=F)]
        eg <- expand.grid(probes, tissues)
        probes.tissues <- sprintf('%s.%s', eg[,1], eg[,2])
        
        probe.tissueTable <- matrix(nrow=length(clusterVals), ncol=length(probes.tissues))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(probes.tissues)){
            probe.tissueTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                     probe.tissue==probes.tissues[j]))
          }
        }
        colnames(probe.tissueTable) <- probes.tissues
        rownames(probe.tissueTable) <- paste("cluster_", clusterVals, sep="")
        probe.tissueTable <- probe.tissueTable[, which(colSums(probe.tissueTable) > 0)]
        probes.tissues <- colnames(probe.tissueTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Probe.Tissue vs. Cluster", quote=F)
        # print(probeTable)
        print(probe.tissueTable, quote = F)
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(probe.tissueTable, simulate.p.value=T))
        
        ## Bar Plots
        probe.tissueColors <- brewer.pal(12, "Paired")
        probe.tissueColorKey <- cbind(probes.tissues, probe.tissueColors[1:length(probes.tissues)])
        
        probe.tissueTablePlots <- tableBarPlots(plotFactor = "probe.tissue", plotFactorVals = probes.tissues, clusterVals = clusterVals, factorTable = probe.tissueTable, colorKey = probe.tissueColorKey)
        
        ctClust <- subset(ctClust, select = -probe.tissue)
      }
      
      ## age & tissue
      if("age.tissue" %in% fisherTests){
        ctClust$age.tissue <- paste(ctClust$age, ctClust$cellSource, sep = ".")
        
        cellSources <- unique(ctClust$cellSource)
        # cellSources <- cellSources[match(c("I", "LN", "S", "PB"), cellSources, nomatch=F)]
        ages <- unique(ctClust$age)
        ages <- ages[order(as.numeric(ages))]
        
        eg <- expand.grid(ages, cellSources)
        ages.tissues <- sprintf('%s.%s', eg[,1], eg[,2])
        
        age.tissueTable <- matrix(nrow=length(clusterVals), ncol=length(ages.tissues))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(ages.tissues)){
            age.tissueTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                         age.tissue==ages.tissues[j]))
          }
        }
        colnames(age.tissueTable) <- ages.tissues
        rownames(age.tissueTable) <- paste("cluster_", clusterVals, sep="")
        age.tissueTable <- age.tissueTable[, which(colSums(age.tissueTable) > 0)]
        ages.tissues <- colnames(age.tissueTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Age.Tissue vs. Cluster", quote=F)
        # print(probeTable)
        print(age.tissueTable, quote = F)
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(age.tissueTable, simulate.p.value=T))
        
        ## Bar Plots
        age.tissueColors <- brewer.pal(12, "Paired")
        age.tissueColorKey <- cbind(ages.tissues, age.tissueColors[1:length(ages.tissues)])
        
        age.tissueTablePlots <- tableBarPlots(plotFactor = "age.tissue", plotFactorVals = ages.tissues, clusterVals = clusterVals, factorTable = age.tissueTable, colorKey = age.tissueColorKey)
        
        grid.arrange(grobs = list(probe.ageTablePlots[[1]], probe.ageTablePlots[[2]], probe.tissueTablePlots[[1]], probe.tissueTablePlots[[2]], age.tissueTablePlots[[1]], age.tissueTablePlots[[2]]), nrow=3)
        
        ctClust <- subset(ctClust, select = -age.tissue)
      }
      
      
      ## probe & tissue & age
      if("probe.tissue.age" %in% fisherTests){
        ctClust$probe.tissue.age <- paste(ctClust$probe, ctClust$cellSource, ctClust$age, sep = ".")
      
        cellSources <- unique(ctClust$cellSource)
        # cellSources <- cellSources[match(c("I", "LN", "S", "PB"), cellSources, nomatch=F)]
        probes <- unique(ctClust$probe)
        # probes <- probes[match(c("PI", "12", "13", "9D", "9Q", "Obsc"), probes, nomatch=F)]
        ages <- unique(ctClust$age)
        ages <- ages[order(as.numeric(ages))]
        
        eg <- expand.grid(probes, cellSources, ages)
        probes.tissues.ages <- sprintf('%s.%s.%s', eg[,1], eg[,2], eg[,3])
        
        probe.tissue.ageTable <- matrix(nrow=length(clusterVals), ncol=length(probes.tissues.ages))
        for(i in 1:length(clusterVals)){
          for(j in 1:length(probes.tissues.ages)){
            probe.tissue.ageTable[i, j] <- nrow(subset(ctClust, kmeans.cluster==clusterVals[i] & 
                                                  probe.tissue.age==probes.tissues.ages[j]))
          }
        }
        colnames(probe.tissue.ageTable) <- probes.tissues.ages
        rownames(probe.tissue.ageTable) <- paste("cluster_", clusterVals, sep="")
        probe.tissue.ageTable <- probe.tissue.ageTable[, which(colSums(probe.tissue.ageTable) > 0)]
        probes.tissues.ages <- colnames(probe.tissue.ageTable)
        
        print("     ", quote=F)
        print("     ", quote=F)
        print("Probe.Tissue.Age vs. Cluster", quote=F)
        # print(probeTable)
        print(probe.tissue.ageTable, quote = F)
        # print(fisher.test(factor(ctClust$probe), factor(ctClust$kmeans.cluster), workspace = 1000000000))
        print(chisq.test(probe.tissue.ageTable, simulate.p.value=T))
        
        ## Bar Plots
        probe.tissue.ageColors <- brewer.pal(12, "Set3")
        probe.tissue.ageColorKey <- cbind(probes.tissues.ages, probe.tissue.ageColors[1:length(probes.tissues.ages)])
        
        probe.tissue.ageTablePlots <- tableBarPlots(plotFactor = "probe.tissue.age", plotFactorVals = probes.tissues.ages, clusterVals = clusterVals, factorTable = probe.tissue.ageTable, colorKey = probe.tissue.ageColorKey)
        
        ng <- nullGrob()
        grid.arrange(grobs = list(probe.tissue.ageTablePlots[[1]], probe.tissue.ageTablePlots[[2]], ng), nrow=3)
        
        ctClust <- subset(ctClust, select = -probe.tissue.age)
      }
      
      ## mouse
      # if("mouse" %in% fisherTests){
      #   mouses <- unique(ctClust$mouse)
      #   mouseTable <- matrix(nrow=length(clusterVals), ncol=length(mouses))
      #   for(i in 1:length(clusterVals)){
      #     for(j in 1:length(mouses)){
      #       mouseTable[i, j] <- nrow(subset(ctClust, mouse==mouses[j] & 
      #                                           kmeans.cluster==clusterVals[i]))
      #     }
      #   }
      #   colnames(mouseTable) <- paste("mouse_", mouses, sep="")
      #   rownames(mouseTable) <- paste("cluster_", clusterVals, sep="")
      #   print("     ", quote=F)
      #   print("     ", quote=F)
      #   print("Mouse vs. Cluster", quote=F)
      #   print(mouseTable)
      #   # print(fisher.test(factor(ctClust$mouse), factor(ctClust$kmeans.cluster), workspace = 1000000000))
      #   print(chisq.test(mouseTable, simulate.p.value=T))
      #   
      # } 
    }
#  }
  
  if(cumulativeExpHist == T){
    #### print histogram of cumulative expression for each cell, colored by cluster ####
    #ctClust$ctSum <- rowSums(ctClust[,10:ncol(ctClust)])
    ctExp <- (ctClust[,10:(ncol(ctClust)-1)])^2
    ctClust$ctSum <- log2(rowSums(ctExp))
    ctSumMin <- min(ctClust$ctSum) - (min(ctClust$ctSum) %% 50)
    ctSumMax <- 50 + max(ctClust$ctSum) - (max(ctClust$ctSum) %% 50)
    ctClust <- ctClust[order(ctClust$kmeans.cluster),]
    
    hist <- ggplot(ctClust, aes(ctSum, fill=factor(kmeans.cluster))) + 
      geom_histogram(bins = 150) +
      scale_x_continuous(breaks = seq(ctSumMin, ctSumMax, 50), 
                         minor_breaks = seq(ctSumMin, ctSumMax, 10)) +
      scale_fill_brewer(palette = "Set3") +
      guides(fill=guide_legend(title="clusters")) +
      ggtitle("histogram of cumulative expression values per cell") +
      xlab("sum of normalized expression values") +
      theme(text=element_text(size=25),
            axis.title.x=element_text(vjust=-0.5),
            plot.margin=unit(c(8,2,8,2),"cm"))
    print(hist)
    
    ## for insulin paper
    # date <- gsub("-", "_", Sys.Date())
    # ggsave(paste("cumulative_exp_hist_", date, ".png", sep=""), hist, device = "png", width = 14, height = 13, path = "~/Documents/abe/biomark/qpcr/figures/insulin/results")
    
    ctClust$ctSum <- NULL
  }
    
    ## filter out low expression clusters
  if(filterClusters == T){
    totalCells <- nrow(ctClust)
    
    ctClust <- subset(ctClust, !(kmeans.cluster %in% clustersToRemove))
    
    remCells <- totalCells - nrow(ctClust)
    remClust <- paste(clustersToRemove, collapse = ", ")
    
    if(removalText == T){
      cat(paste(remCells, " / ", totalCells, " cells were removed due to low gene expression (cluster ID of removed clusters: ", remClust, ")", sep=""))
    }
  }
  
  write.csv(ctClust, paste0(baseDir, "results/clusterFilter-017-output.csv"))
  return(ctClust)
  
}