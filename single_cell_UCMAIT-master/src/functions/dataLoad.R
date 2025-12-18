dataLoad <- function(Project,Species,FreshorFrozen,Panel,Troubleshoot,CombineCommonGenes,RemoveDonor,RemoveSample){
  
  library(plyr)
  library(dplyr)
  library(tidyverse)
  library(stringr)
  
  ## set data directory 
  dataDir <- paste(baseDir, "data/", sep="")
  #dataDir <- ("/Users/smsharma/Documents/OneDrive - Scripps Research/Laboratory Work/Projects MAIT- UC project/Log-BioMark_exported_csv_files/")
  
  assign("Troubleshoot", Troubleshoot, envir = .GlobalEnv)
  assign("CombineCommonGenes", CombineCommonGenes, envir = .GlobalEnv)
  assign("Species", Species, envir = .GlobalEnv)
  
  HumanPanels<- grepl("HU", Species)
  MousePanels<- grepl("MS", Species)
  if (HumanPanels == TRUE) {
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
   #Panel3<- c("ADAMTS2","ANG","APP","AREG","ATOH8","C3","C5","CAMP","CCL24","CCL3","CCL8","CDC42","CEACAM1","CFH","CMA1","CRISPLD2","CSF1",
   #           "CSF1R","CSF2","CTNNB1","CXCL10","CXCL2","DEFB1","ECM1","EGFL6","EGR1","EGR2","ENG","F13A1","F2R","F2RL2","FFAR2","FGF9","FGFR2",
   #           "FNDC4","FURIN","FZD1","GRN","HBEGF","HGF","HIF1A","HMGB1","HMOX1","IFNG","IFNGR1","IGF1","IL12RB1","IL17A","IL1A","IL1B","IL22","IL26",    
   #           "ITGA4","ITGAL","ITGB1","ITGB2","JAG1","LEF1","LGALS3","LRP5","LRP6","MMP2","MMP24","MMP25","MMP28","MMP9","NPNT","OSM","P3H4","PDCD1",
   #           "PDGFA","PDGFB","PI3","PTGES2","RAC1","RAC2","REG3G","RHO","SDC1","SYK","TCF7","TGFA","TGFB1","TGFB3","THBS1","TNF","TNFRSF1A","TNFRSF21",
   #           "TPSB2","TRAF3IP2","TSPAN2","VEGFA","VEGFB","VEGFC","VWF","WNT1")
  #CellID<- c("ACTA2", "ACVR1", "ADGRE1", "ANGPT1", "ANPEP", "BMP2", "CCL19", "CD14", "CD24", "CD36", "CD3E", "CD3G", "CD4", "CD44", "CD74", "CD80", "CD83", 
  #            "CD86", "CD8A", "CLEC7A", "COL11A1", "COL1A1", "COL2A1", "COL4A1", "CSF1R", "CSF3R", "CSPG4", "CXCL12", "CXCL13", "DES", "EGFR", "FCGR1A", "FGF2", 
  #            "FGFR1", "FGFR3", "FGR", "FLT1", "FLT4", "FN1", "GFAP", "HIF1A", "HLA-DMA", "HLA-DRA", "ICAM1", "ICAM2", "ICOSLG", "IFNG", "IGF1", "IGF2", "IL12A", 
  #            "IL1B", "IL6", "ITGAE", "ITGAM", "ITGAX", "ITGB1", "KDR", "LCK", "LY75", "LYVE1", "MMP13", "MMP2", "MMP3", "MMP9", "NFATC1", "NLRP3", "PDGFA", "PDGFRA", 
  #            "PDGFRB", "PDPN", "PECAM1", "PTGS2", "PTK2B", "RSPO1", "RSPO2", "SFRP1", "SMN1", "SPP1", "TEK", "TGFB1", "TIMP1", "TIMP2", "TLR3", "TLR4", "TLR5", "TLR7", 
  #            "TLR8", "TLR9", "TNC", "TNFSF11", "VCAM1", "VEGFA", "MR1-a1toa2", "VEGFC", "WNT2B", "WNT4")
     
   CellID<- c("ACTA2",	"CD19",	"CD74",	"COL2A1",	"EGFR",	"FLT1",	"ICAM1", "IL6",	"MMP2",	"PECAM1",	"TGFB1",	"TLR9",
              "ACVR1",	"CD1D",	"CD80",	"COL4A1",	"FCGR1A",	"FLT4",	"ICAM2",	"ITGAE",	"MMP3",	"PTGS2",	"TIMP1",	"TNC",
              "ADGRE1",	"CD24",	"CD83",	"CSF1R",	"FGF2",	"FN1",	"ICOSLG",	"ITGAM",	"MMP9",	"PTK2B",	"TIMP2",	"TNFSF11",
              "ANGPT1",	"CD36",	"CD86",	"CSF3R",	"FGFR1",	"GFAP",	"IFNG",	"ITGAX",	"NLRP3",	"RSPO1",	"TLR3",	"VCAM1",
              "ANPEP",	"CD3E",	"CD8A",	"CSPG4",	"FGFR3",	"HIF1A",	"IGF1",	"ITGB1",	"PDGFA", "SFRP1",	"TLR4",	"VEGFA",
              "BMP2",	"CD3G",	"CLEC7A",	"CXCL12",	"FGR",	"HLA-B",	"IGF2",	"KDR",	"PDGFRA",	"SMN1",	"TLR5",	"VEGFC",
              "CCL19",	"CD4",	"COL11A1",	"CXCL13",	"CD45",	"HLA-DMA",	"IL12A",	"LY75",	"PDGFRB",	"SPP1",	"TLR7",	"WNT2B",
              "CD14",	"CD44",	"COL1A1",	"DES",	"MR1",	"HLA-DRA",	"IL1B",	"LYVE1",	"PDPN",	"TEK",	"TLR8",	"WNT4")
   
    } else if (MousePanels == TRUE) {
   Panel1<- c("AIM2","BCL2","BCL6","CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CD28","CD3E","CD4","CD40","CD44","CD80","CD86","CD8A","CEACAM1","CTLA4",
              "CXCL10","CXCR3","CXCR4","FOXP3","FYN","GAPDH","GATA4","GSK3A","GSK3B","HPRT","ICAM1","ICOS","IFI44","IFI44L","IFIT1","IFIT3","IFNG","IFNGR1",
              "IL-21","IL10","IL12B","IL12RB","IL17A","IL18R1","IL1R2","IL2","IL25","IL27","IL27R","IL2RA","IL3","IL4","IL4RA","IL5","IL5RA","IL6","IL7","IL7R",
              "IRF1","IRF2","IRF4","IRF7","ISG15","JAK1","JAK2","LY6E","MAP2K6","MAPK8","MX1","NFKB1","NUR77","OAS1B","OAS2","OASL1","PD1","PDL-1","PPARA","PPARG",
              "PPARGC1A","PTEN","RSAD2","SOCS3","STAT1","STAT3","STAT4","STAT5","TBX21","TGFBR2","TNF","TNFAIP3","TNFRSF1A","TNFRSF1B","TRAF2","VAV1","ZAP70","ZEB2")

   CellID<- c("ACTA2","ACVR1","ADGRE1","ANGPT1","ANPEP","BMP5","BMP7","CD14","CD24A","CD36","CD3E","CD4","CD44","CD74","CD80","CD83","CD86","CD8A",   
              "CLEC7A","COL11A1","COL1A1","COL1A2","CSF1","CSF1R","CSF2RA","CSF2RB","CXCL13","DES","EGFR","FAP","FCGR1","FGFR1","FGFR3","FGR","FLT4","GCG",    
              "GFAP","GHRL","GM13889","H2-AA","H2-DMA","HIF1A","IAPP","ICAM1","ICAM2","ICOSL","IFNG","IGF1","IGF2","IL1A","IL1B","IL34","INS1","INS2",   
              "ITGAX","ITGB1","KDR","KLF5","LCK","LEPR","LY75","MMP1A","MMP2","MMP3","MMP9","NFATC1","NLRP3","PDGFA","PDGFB","PDGFRB","PDPN","PECAM1", 
              "PPY","PTGS2","PTK2","RSPO1","SELE","SFRP1","SPP1","SST","TEK","TGFB1","TIMP1","TIMP2","TLR3","TLR4","TLR7","TLR9","TNC","TNFSF11",
              "VCAM1","VEGFA","VEGFB","WNT2B","WNT4","ZAP70")
   Panel3<- c()
    #warning(print("You have yet to add the mouse panels into the 'dataLoad.R' function. Go to the function and add it in."))
  }
  
  if (Panel == "1 and 3") {
    CommonGenes <- intersect(Panel1, Panel3)
    assign("CommonGenes", CommonGenes, envir = .GlobalEnv)
    print(CommonGenes)
  } else if (Panel == "1 and CellID") {
    # For Mouse (and optionally Human) Panel1 + CellID overlap
    CommonGenes <- intersect(Panel1, CellID)
    assign("CommonGenes", CommonGenes, envir = .GlobalEnv)
    print(CommonGenes)
  }
  

  files_list <- list.files(dataDir, pattern = "F*csv")
  Pro <- files_list[grep(Project, files_list)]
  Pro_Species <- Pro[grep(Species, Pro)]
  FRvsFZ<- grepl("All", FreshorFrozen)
  if (FRvsFZ == TRUE) {
    Pro_Species_F<- Pro_Species
  } else {
    Pro_Species_F <- Pro_Species[grep(FreshorFrozen, Pro_Species)]
  }
  print(paste0("The data files listed below match the conditions in the 'dataload'function:"))
  print(Pro_Species_F)

  
  ###A cumbersome for-loop to do a myriad different things including: reading in the files we need. It also removes wells that are "EMPTY", 
  ###checks to see if the gene listed is Panel 1, 2 etc, and adds the plate name to the well name to give the full "cell name".
  ListofFiles<- list()
  ListofFilesNames <- c()
  for (i in 1:length(Pro_Species_F)) {
    #print(paste0(dataDir,Pro_Species_F[i]))
    all_content<- readLines(paste(dataDir, Pro_Species_F[i], sep=""))
    ListofFiles[[i]] <- read.csv(textConnection(all_content[-c(1:11)]), header=TRUE, stringsAsFactors=FALSE)
    ListofFiles[[i]]<- ListofFiles[[i]][grep("EMPTY", ListofFiles[[i]]$Name, invert = TRUE),]
    ListofFiles[[i]]<- ListofFiles[[i]][grep("Empty", ListofFiles[[i]]$Name, invert = TRUE),]
    
    PanelDetect1<- setequal(Panel1, unique(sort(ListofFiles[[i]]$Name.1)))
    PanelDetect2<- setequal(Panel3, unique(sort(ListofFiles[[i]]$Name.1)))
    PanelDetectCellID<- setequal(CellID, unique(sort(ListofFiles[[i]]$Name.1)))
    
    ListofFiles[[i]]$Comments <- ""
    
    ###Giving the elements in the list a name, the name consists of the Panel number and the element number. Below I accidentally wrote "PanelDetect2" with "Panel3". Excuse the error.
    
    if (PanelDetect1 == TRUE) {
      ListofFiles[[i]]$PanelNumber<- "Panel1"
      ListofFilesNames[i] <- paste0("Panel1","-",i)
    } else if (PanelDetect2 == TRUE) {
      ListofFiles[[i]]$PanelNumber<- "Panel3"
      ListofFilesNames[i] <- paste0("Panel3","-",i)
    } else if (PanelDetectCellID == TRUE) {
      ListofFiles[[i]]$PanelNumber<- "CellID"
      ListofFilesNames[i] <- paste0("CellID","-",i)
    }
    
    for (x in 1:nrow(ListofFiles[[i]])) {
      ListofFiles[[i]][x,"Name"]<- paste0(str_split(Pro_Species_F[i], "_[:digit:][:digit:][:digit:][:digit:][:digit:][:digit:][:digit:][:digit:][:digit:][:digit:].csv")[[1]][1],"_",ListofFiles[[i]][x,"Name"])
      ListofFiles[[i]][x,"PlateName"] <- paste0(str_split(Pro_Species_F[i], "_[:digit:][:digit:][:digit:][:digit:][:digit:][:digit:][:digit:][:digit:][:digit:][:digit:].csv")[[1]][1])
    }
  }
  names(ListofFiles) <- ListofFilesNames
  
  ###SMS edit: to put the panel number into an object in R
  assign("PanelNumberLoaded", Panel, envir = .GlobalEnv)
  
  ###Now we're going to separate out the .csv files based on the panel number (ie. it will only load the CSV files for specified panel number). If we want to do both panel 1 and 3,
  ###then we'll keep these separate and merge them later.
  
  PanelTest1 <- c()
  PanelTest2 <- c()
  PanelTest3 <- c()
  PanelTestCellID <- c()
  
  ListofFiles_Panel1<- list()
  ListofFiles_Panel2<- list()
  ListofFiles_Panel3<- list()
  ListofFiles_PanelCellID<- list()
  
  PanelTest1 <-  grep("Panel1", names(ListofFiles))
  PanelTest2 <-  grep("Panel2", names(ListofFiles))
  PanelTest3 <-  grep("Panel3", names(ListofFiles))
  PanelTestCellID <-  grep("CellID", names(ListofFiles))

  ListofFiles_Panel1 <- ListofFiles[PanelTest1]
  ListofFiles_Panel2 <- ListofFiles[PanelTest2]
  ListofFiles_Panel3 <- ListofFiles[PanelTest3]
  ListofFiles_PanelCellID <- ListofFiles[PanelTestCellID]
  
  # for (i in 1:length(ListofFiles)) {
  #   print(colnames(ListofFiles[[i]]))
  # }
  # 
  
  Panel1_SingleTest <- str_detect(Panel, "1")
  if (Troubleshoot == TRUE) {print(Panel1_SingleTest)}
  Panel2_SingleTest <- str_detect(Panel, "2")
  if (Troubleshoot == TRUE) {print(Panel2_SingleTest)}
  Panel3_SingleTest <- str_detect(Panel, "3")
  if (Troubleshoot == TRUE) {print(Panel3_SingleTest)}
  PanelCellID_SingleTest <- str_detect(Panel, "CellID")
  if (Troubleshoot == TRUE) {print(PanelCellID_SingleTest)}
  
  # if (Panel1_SingleTest == TRUE) {
  #   if (Panel2_SingleTest == FALSE) {
  #     print("User selected Panel 1 only")
  #     ctTableCombine<- do.call(rbind, ListofFiles_Panel1)
  #   } else if (Panel2_SingleTest == TRUE) {
  #     print("User selected Panel 1 and 2")
  #   } 
  # } else if (Panel1_SingleTest == FALSE) {
  #   if (Panel2_SingleTest == TRUE) {
  #     print("User selected Panel 2 only")
  #     ctTableCombine<- do.call(rbind, ListofFiles_Panel2)
  #   } else if (Panel1_SingleTest == TRUE) {
  #     print("User selected Panel 2 and 1")
  #   }
  # }
  
  PrintList <- c()
  
  # Panel 1 only
  if (Panel1_SingleTest == TRUE && Panel3_SingleTest == FALSE && PanelCellID_SingleTest == FALSE) {
    
    print("User selected Panel 1 only")
    ctTableCombine <- do.call(rbind, ListofFiles_Panel1)
    print("The following CSV files were loaded into the analysis:")
    for (i in 1:length(ListofFiles_Panel1)) {
      PrintList[i] <- ListofFiles_Panel1[[i]][1, "PlateName"]
    }
    print(PrintList)
    write.csv(ctTableCombine, paste0(baseDir, "results/dataLoad-001-ctTableCombine.csv"))
    
    # Panel 3 only
  } else if (Panel1_SingleTest == FALSE && Panel3_SingleTest == TRUE && PanelCellID_SingleTest == FALSE) {
    
    print("User selected Panel 3 only")
    ctTableCombine <- do.call(rbind, ListofFiles_Panel3)
    print("The following CSV files were loaded into the analysis:")
    for (i in 1:length(ListofFiles_Panel3)) {
      PrintList[i] <- ListofFiles_Panel3[[i]][1, "PlateName"]
    }
    print(PrintList)
    write.csv(ctTableCombine, paste0(baseDir, "results/dataLoad-001-ctTableCombine.csv"))
    
    # CellID only
  } else if (Panel1_SingleTest == FALSE && Panel3_SingleTest == FALSE && PanelCellID_SingleTest == TRUE) {
    
    print("User selected CellID only")
    ctTableCombine <- do.call(rbind, ListofFiles_PanelCellID)
    write.csv(ctTableCombine, paste0(baseDir, "results/dataLoad-001-ctTableCombine.csv"))
    print("The following CSV files were loaded into the analysis:")
    for (i in 1:length(ListofFiles_PanelCellID)) {
      PrintList[i] <- ListofFiles_PanelCellID[[i]][1, "PlateName"]
    }
    print(PrintList)
    
    # Panel 1 and 3 (existing behavior)
  } else if (Panel1_SingleTest == TRUE && Panel3_SingleTest == TRUE && PanelCellID_SingleTest == FALSE) {
    
    print("User selected Panel 1 and 3")
    PlateNames_Panel1 <- c()
    PlateNames_Panel3 <- c()
    for (i in 1:length(ListofFiles_Panel1)) {
      PlateNames_Panel1[i] <- ListofFiles_Panel1[[i]][1, "PlateName"]
    }
    names(ListofFiles_Panel1) <- PlateNames_Panel1
    
    for (i in 1:length(ListofFiles_Panel3)) {
      PlateNames_Panel3[i] <- ListofFiles_Panel3[[i]][1, "PlateName"]
    }
    names(ListofFiles_Panel3) <- PlateNames_Panel3
    
    CommonPlateNames <- intersect(PlateNames_Panel1, PlateNames_Panel3)
    print("These are the plates that were assessed by Panels 1 and 3. They will be loaded into R for analysis.")
    print(CommonPlateNames)
    
    ListofPanel1and2 <- list()
    for (i in CommonPlateNames) {
      ListofPanel1and2[[i]] <- bind_rows(ListofFiles_Panel1[[i]], ListofFiles_Panel3[[i]])
      ListofPanel1and2[[i]] <- ListofPanel1and2[[i]][order(ListofPanel1and2[[i]][, "Name"]), ]
    }
    
    ctTableCombine <- bind_rows(ListofPanel1and2)
    write.csv(ctTableCombine, paste0(baseDir, "results/dataLoad-001-ctTableCombine.csv"))
    
    # ✅ NEW: Panel 1 and CellID (for MS)
  } else if (Panel1_SingleTest == TRUE && Panel3_SingleTest == FALSE && PanelCellID_SingleTest == TRUE) {
    
    print("User selected Panel 1 and CellID")
    PlateNames_Panel1 <- c()
    PlateNames_CellID <- c()
    
    for (i in 1:length(ListofFiles_Panel1)) {
      PlateNames_Panel1[i] <- ListofFiles_Panel1[[i]][1, "PlateName"]
    }
    names(ListofFiles_Panel1) <- PlateNames_Panel1
    
    for (i in 1:length(ListofFiles_PanelCellID)) {
      PlateNames_CellID[i] <- ListofFiles_PanelCellID[[i]][1, "PlateName"]
    }
    names(ListofFiles_PanelCellID) <- PlateNames_CellID
    
    CommonPlateNames <- intersect(PlateNames_Panel1, PlateNames_CellID)
    print("These are the plates that were assessed by Panels 1 and CellID. They will be loaded into R for analysis.")
    print(CommonPlateNames)
    
    ListofPanel1andCellID <- list()
    for (i in CommonPlateNames) {
      ListofPanel1andCellID[[i]] <- bind_rows(ListofFiles_Panel1[[i]], ListofFiles_PanelCellID[[i]])
      ListofPanel1andCellID[[i]] <- ListofPanel1andCellID[[i]][order(ListofPanel1andCellID[[i]][, "Name"]), ]
    }
    
    ctTableCombine <- bind_rows(ListofPanel1andCellID)
    write.csv(ctTableCombine, paste0(baseDir, "results/dataLoad-001-ctTableCombine.csv"))
    
  } else {
    warning("Panel selection combination not explicitly handled in dataLoad(). Check the 'Panel' argument.")
  }
  
  
  ##KA edit
  #print("THIS IS WHATS GOING INTO THE TABLE WE EVENTUALLY USE FOR OTHER STUFF")
  #assign("ctTableCombine", ctTableCombine, envir = .GlobalEnv)
  
  #ctTableCombine<- do.call(rbind, ListofFiles)
 
   # --- sanity check / repair for Name column ---
  if (!"Name" %in% colnames(ctTableCombine)) {
    message("ctTableCombine is missing 'Name'. Columns are: ",
            paste(colnames(ctTableCombine), collapse = ", "))
    
    # Try to reconstruct it from likely columns
    if ("Sample" %in% colnames(ctTableCombine)) {
      ctTableCombine$Name <- ctTableCombine$Sample
    } else if ("Name.1" %in% colnames(ctTableCombine)) {
      ctTableCombine$Name <- ctTableCombine$Name.1
    } else {
      stop("No suitable column to use as Name in ctTableCombine.")
    }
  }
  
  ###We are breaking out most of the meta-data from the cell name or using the patterns in the cell name to put the proper metadata into the file.
  for(i in 1:nrow(ctTableCombine)) {
    ctTableCombine[i, "cell.type"] <- str_split(ctTableCombine[i, "Name"], "_")[[1]][12]
    ctTableCombine[i, "probe"] <- str_split(ctTableCombine[i, "Name"], "_")[[1]][13]
    ctTableCombine[i, "FRorFZ"] <-str_split(ctTableCombine[i, "Name"], "_")[[1]][1]
    ctTableCombine[i, "patient"] <-str_split(ctTableCombine[i, "Name"], "_")[[1]][10]
    ctTableCombine[i,"cellSource"]<- if (grepl("_BL_", ctTableCombine$Name[i]) == TRUE) {"blood"} else {"tissue"}
    ctTableCombine[i,"age"]<- if (grepl("_BL_", ctTableCombine$Name[i]) == TRUE) {
      "blood"
      } else if (grepl("-I-", ctTableCombine$Name[i]) == TRUE) {
        "infl"
      } else if (grepl("-UI-", ctTableCombine$Name[i]) == TRUE) {
          "uninfl"
      } else if (grepl("_IS_", ctTableCombine$Name[i]) == TRUE) {
          "islets"
      }
  }


  ###We are going to try to use the established columns for Kenna's metadata. We'll have to fix this later.
  if (Project == "TFH") {
    Kenna<- str_split(unique(ctTableCombine$Name), pattern = "P[:digit:]_")
    Kenna <- unlist(Kenna)
    KennaVector <- Kenna[seq(2,192,2)]
    ctTableCombine$cell.type<- "CD4_TFH"
    ctTableCombine$cellSource<- sapply(1:nrow(ctTableCombine), function (x) {as.factor(str_sub(str_split(ctTableCombine[x,"Name"], pattern = "_")[[1]][7],1,1))})
    ctTableCombine$probe <- sapply(1:nrow(ctTableCombine), function (x) {str_sub(str_split(ctTableCombine[x,"Name"], pattern = "_")[[1]][7],2,2)})
    ctTableCombine$age<- sapply(1:nrow(ctTableCombine), function (x) {str_sub(str_split(ctTableCombine[x,"Name"], pattern = "_")[[1]][7],3,4)})
  }
  
  ## change column names
  names(ctTableCombine)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  #print(colnames(ctTableCombine))
  print(paste0("The panel's pulled for analysis are: "))
  print(unique(ctTableCombine$PanelNumber))
  
  ###As there are a few genes that are common between our panels, you will need to change the names so that the recast function in the cleanCt function doesn't fail
  ###However, you may prefer to just average the common genes. If you'd like to average the common genes, even if there is a large difference between the two, use "TRUE" in "CombineCommonGenes".
  ###A "FALSE" setting will use put a "_#" next to the common genes.
  
  # Combine overlapping genes across panels (works for 1+3 or 1+CellID)
  if (CombineCommonGenes == TRUE && exists("CommonGenes")) {
    Test <- ctTableCombine
    print(paste0("Number of rows prior to duplicate genes, per cell, being combined: ", nrow(Test)))
    
    Test[Test$ct == 999, "ct"] <- NA
    Test <- Test %>%
      group_by(cellType, gene) %>%
      mutate(mean_ct = mean(ct, na.rm = TRUE)) %>%
      as.data.frame()
    
    for (i in CommonGenes) {
      if (Troubleshoot == TRUE) {
        print(paste0("Rows that contain data for ", i, ":"))
      }
      Rows <- which(Test$gene == i)
      even <- seq(2, by = 2, len = (length(Rows) / 2))
      Test <- Test[-c(Rows[even]), ]
      if (Troubleshoot == TRUE) {
        print(paste0("Number of rows AFTER duplicate genes have been removed: ", nrow(Test)))
      }
    }
    
    # Generalized prediction (no hard-coded 192)
    nCells <- length(unique(ctTableCombine$cellType))
    Predicted <- nrow(ctTableCombine) - length(CommonGenes) * nCells
    print(paste0("Predicted number of rows after removing duplicate genes, per cell: ", Predicted))
    
    if (Predicted != nrow(Test)) {
      print(warning("ERROR: The predicted number of rows does not match the number of rows, post duplicate gene removal"))
    } else {
      print("The predicted number of rows DOES match the number of rows, post duplicate gene removal")
    }
    
    ctTableCombine <- Test
    ctTableCombine[ctTableCombine$mean_ct == "NaN", "mean_ct"] <- 999
    ctTableCombine <- ctTableCombine %>%
      select(-"ct") %>%
      relocate("mean_ct", .after = "Type.1") %>%
      dplyr::rename("ct" = "mean_ct")
  }
  

  if (CombineCommonGenes == FALSE) {
    if (Panel1_SingleTest == TRUE && Panel3_SingleTest == TRUE) {
      CommonGenes
      for (i in CommonGenes) {
        RowNumber <- grep(paste0("^",i,"$"), ctTableCombine$gene)
        for (x in RowNumber) {
          Panelidentified <- ctTableCombine[x,"PanelNumber"]
          Panel1Test<- grepl("Panel1", Panelidentified)
          Panel3Test<- grepl("Panel3", Panelidentified)
          if (Panel1Test == TRUE) {
            ctTableCombine[x, "gene"] <- paste0(i,"_1")
          } else if (Panel3Test == TRUE) {
            ctTableCombine[x, "gene"] <- paste0(i,"_2")
          }
        }
      }
    }
  }
  
  ###if statement to remove donors. Intended to remove only data from donors that have technical issues that are appearing as separate clusters.
  '%notin%' <- Negate('%in%')
  if (length(RemoveDonor) > 0) {
    print(paste("Removed the following donor(s):", paste(RemoveDonor, collapse = ' & '), sep = " "))
    ctTableCombine<- ctTableCombine %>% filter(probe %notin% RemoveDonor)
  }
  
  if (length(RemoveSample) > 0) {
    print(paste("Removed the following samples(s):", paste(RemoveSample, collapse = ' & '), sep = " "))
    ctTableCombine<- ctTableCombine %>% filter(!grepl(RemoveSample, cellType))
  }
  
  write.csv(ctTableCombine, paste0(baseDir, "results/dataLoad-002-output.csv"))
  
  print(paste0("Are blood samples in this table? ", any(ctTableCombine$cellSource == "blood")))
  
  #print(unique(ctTableCombine$probe))
  #print(ctTableCombine)
  
  return(ctTableCombine)
}

###Ishan Taneja's code to help average CT values using dplyr functions
# tmp = ctTableCombine
# tmp[tmp$ct == 999,'ct'] = NA
# #tmp = tmp %>% group_by(cellType, gene) %>% summarise(count=n(), mean_ct = mean(ct, na.rm=TRUE), cv_ct = sd(ct, na.rm=TRUE)/mean(ct, na.rm=TRUE)) %>% as.data.frame()
# tmp = tmp %>% group_by(cellType, gene) %>% mutate(rep_num = row_number(), count=n(), mean_ct = mean(ct, na.rm=TRUE), cv_ct = sd(ct, na.rm=TRUE)/mean(ct, na.rm=TRUE)) %>% as.data.frame()
# #tmp = tmp %>% filter(rep_num == 1) %>% as.data.frame()

###Sidd building on Ishan's code above, adding a for-loop to remove rows that have the duplicate measurement of a gene,removing the original ct column and replacing it with the "mean-ct" column
# Test<- ctTableCombine
# Test[Test$ct == 999,'ct'] = NA
# Test = Test %>% group_by(cellType, gene) %>% mutate(mean_ct = mean(ct, na.rm=TRUE)) %>% as.data.frame()
# nrow(Test)
# 
# for (i in CommonGenes) {
#   print(nrow(Test))
#   print(paste0("Rows that contain data for ", i,":"))
#   print(which(Test$gene == i))
#   Rows<- which(Test$gene == i)
#   #odd <- seq(1,by=2, len=(length(Rows)/2))
#   even<- seq(2,by=2, len=(length(Rows)/2))
#   print(even)
#   print(Rows[even])
#   print(Test[-c(Rows[even]),])
#   Test<-Test[-c(Rows[even]),]
#   print(nrow(Test))
# }
# nrow(ctTableCombine)
# nrow(Test)
