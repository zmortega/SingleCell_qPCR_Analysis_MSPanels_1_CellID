#### load UC-MAIT data ####


dataLoadUCMAIT <- function(){
  ## set data directory
  dataDir <- paste(baseDir, "data/", sep="")
  
  ### read in data
  ## probes: NBD1 - UCM5
  ctTable1 <- read.csv(paste(dataDir, "1362292369-NBD1_UCM5-HeatMapResultsUpdated.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable1$Comments<- NA
  ctTable1$cellSource <- "blood"
  ctTable1$probe <- NA
  for(i in 1:nrow(ctTable1)) {
    if (grepl("UCM5", ctTable1[i,2], fixed=TRUE)){
      ctTable1[i, "probe"] <- "UCM5"
    }
    if (grepl("NBD1", ctTable1[i,2], fixed=TRUE)){
      ctTable1[i, "probe"] <- "NBD1"
    }
  }
  ctTable1$age <- "blood"
  ctTable1$cell.type<- "MAIT"
  ctTable1$FRorFZ<-"FR"

  ### read in data
  ## probes: NBD3 - UCM6
  ctTable2 <- read.csv(paste(dataDir, "1362292371-NBD3-UCM6-HeatMapResults.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable2$Comments<- NA
  ctTable2$cellSource <- "blood"
  ctTable2$probe <- NA
  for(i in 1:nrow(ctTable2)) {
    if (grepl("UCM6", ctTable2[i,2], fixed=TRUE)){
      ctTable2[i, "probe"] <- "UCM6"
    }
    if (grepl("NBD3", ctTable2[i,2], fixed=TRUE)){
      ctTable2[i, "probe"] <- "NBD3"
    }
  }
  ctTable2$age <- "blood"
  ctTable2$cell.type<- "MAIT"
  ctTable2$FRorFZ<-"FR"

  ### read in data
  ## probes: NBD4
  ctTable3 <- read.csv(paste(dataDir, "1362351425-NBD4.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable3$Comments<- NA
  ctTable3$cellSource <- NA
  ctTable3$probe <- "NBD4"
  ctTable3$age <- NA
  ctTable3$cell.type<- "MAIT"
  ctTable3$FRorFZ<-"FR"

  ### read in data
  ## probes: UCM14
  ctTable4 <- read.csv(paste(dataDir, "1362351426-UCM14.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable4$Comments<- NA
  ctTable4$cellSource <- NA
  ctTable4$probe <- "UCM14"
  ctTable4$age <- NA
  ctTable4$cell.type<- "MAIT"
  ctTable4$FRorFZ<-"FR"

  ### read in data
  ## probes: UCM13
  ctTable5 <- read.csv(paste(dataDir, "1362351427-UCM13.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable5$Comments<- NA
  ctTable5$cellSource <- NA
  ctTable5$probe <- "UCM13"
  ctTable5$age <- NA
  ctTable5$cell.type<- "MAIT"
  ctTable5$FRorFZ<-"FR"

  ### read in data
  ## probes: UCM12
  ctTable6 <- read.csv(paste(dataDir, "1362351428_UCM12.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable6$Comments<- NA
  ctTable6$cellSource <- NA
  ctTable6$probe <- "UCM12"
  ctTable6$age <- NA
  ctTable6$cell.type<- "MAIT"
  ctTable6$FRorFZ<-"FR"

  ### read in data
  ## probes: UCM10
  ctTable7 <- read.csv(paste(dataDir, "1362351431_UCM10.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable7$Comments<- NA
  ctTable7$cellSource <- NA
  ctTable7$probe <- "UCM10"
  ctTable7$age <- NA
  ctTable7$cell.type<- "MAIT"
  ctTable7$FRorFZ<-"FR"

  ### read in data
  ## probes: UCM15
  ctTable8 <- read.csv(paste(dataDir, "1362356328-UCM15.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable8$Comments<- NA
  ctTable8$cellSource <- NA
  ctTable8$probe <- "UCM15"
  ctTable8$age <- NA
  ctTable8$cell.type<- "MAIT"
  ctTable8$FRorFZ<-"FR"

  ### read in data
  ## probes: UCM16
  ctTable9 <- read.csv(paste(dataDir, "1362356329-UCM16.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable9$Comments<- NA
  ctTable9$cellSource <- NA
  ctTable9$probe <- "UCM16"
  ctTable9$age <- NA
  ctTable9$cell.type<- "MAIT"
  ctTable9$FRorFZ<-"FR"

  ### read in data
  ## probes: UCM17
  ctTable10 <- read.csv(paste(dataDir, "1362351435_Processed-TNET001-UCM17.csv", sep=""),
                       header=TRUE, stringsAsFactors=FALSE)
  ctTable10 <- subset(ctTable10, grepl("UCM17", ctTable10$Name, fixed = TRUE))
  ctTable10$Comments<- NA
  ctTable10$cellSource <- NA
  ctTable10$probe <- "UCM17"
  ctTable10$age <- NA
  ctTable10$cell.type<- "MAIT"
  ctTable10$FRorFZ<-"FR"

  ### read in data
  ## probes: RUS001
  all_content<- readLines(paste(dataDir, "FR_RNA_UCM_HU_20200610_P1_1362427403.csv", sep=""))
  ctTable11 <- read.csv(textConnection(all_content[-c(1:11)]),
                        header=TRUE, stringsAsFactors=FALSE)
  #ctTable11 <- ctTable11[-c(1:9),]
  ctTable11$Comments<- NA
  ctTable11$cellSource <- NA
  ctTable11$probe <- "RUS001"
  ctTable11$age <- NA
  ctTable11$cell.type<- "MAIT"
  ctTable11$FRorFZ<-"FR"

  ### read in data
  ## probes: RUS002
  all_content<- readLines(paste(dataDir, "FR_RNA_UCM_HU_20200903_P1_1362537186.csv", sep=""))
  ctTable12 <- read.csv(textConnection(all_content[-c(1:11)]),
                        header=TRUE, stringsAsFactors=FALSE)
  #ctTable12 <- ctTable12[-c(1:9),]
  ctTable12$Comments<- NA
  ctTable12$cellSource <- NA
  ctTable12$probe <- "RUS002"
  ctTable12$age <- NA
  ctTable12$cell.type<- "MAIT"
  ctTable12$FRorFZ<-"FR"

  ### read in data
  ## probes: RUS007_P1
  all_content<- readLines(paste(dataDir, "FR_RNA_UCM_HU_20201016_P1_1362537463.csv", sep=""))
  ctTable13 <- read.csv(textConnection(all_content[-c(1:11)]),
                        header=TRUE, stringsAsFactors=FALSE)
  #ctTable13 <- ctTable13[-c(1:9),]
  ctTable13$Comments<- NA
  ctTable13$cellSource <- NA
  ctTable13$probe <- "RUS007"
  ctTable13$age <- NA
  ctTable13$cell.type<- "MAIT"
  ctTable13$FRorFZ<-"FR"

  ### read in data
  ## probes: RUS008_Pl1
  all_content<- readLines(paste(dataDir, "FR_RNA_UCM_HU_20201113_P1_1362576147.csv", sep=""))
  ctTable14 <- read.csv(textConnection(all_content[-c(1:11)]),
                        header=TRUE, stringsAsFactors=FALSE)
  #ctTable14 <- ctTable14[-c(1:9),]
  ctTable14$Comments<- NA
  ctTable14$cellSource <- NA
  ctTable14$probe <- "RUS008"
  ctTable14$age <- NA
  for(i in 1:nrow(ctTable14)) {
    if (grepl("5OPRU", ctTable14[i,"Name"], fixed=TRUE)){
      ctTable14[i, "cell.type"] <- "MAIT"
    }
    if (grepl("PBS57", ctTable14[i,"Name"], fixed=TRUE)){
      ctTable14[i, "cell.type"] <- "NKT"
    }
  }
  ctTable14<- subset(ctTable14, cell.type != "NKT")
  ctTable14$FRorFZ<-"FR"

  ### read in data
  ## probes: RUS008_Plate 2
  all_content<- readLines(paste(dataDir, "FR_RNA_UCM_HU_20201113_P2_1362576022.csv", sep=""))
  ctTable15 <- read.csv(textConnection(all_content[-c(1:11)]),
                        header=TRUE, stringsAsFactors=FALSE)
  #ctTable15 <- ctTable15[-c(1:9),]
  ctTable15$Comments<- NA
  ctTable15$cellSource <- NA
  for(i in 1:nrow(ctTable15)) {
    if (grepl("RUS008", ctTable15[i,2], fixed=TRUE)){
      ctTable15[i, "probe"] <- "RUS008"
    }
    if (grepl("Empty", ctTable15[i,2], fixed=TRUE)){
      ctTable15[i, "probe"] <- "Empty"
    }
  }
  ctTable15<- subset(ctTable15, probe != "Empty")
  ctTable15$age <- NA
  for(i in 1:nrow(ctTable15)) {
    if (grepl("5OPRU", ctTable15[i,"Name"], fixed=TRUE)){
      ctTable15[i, "cell.type"] <- "MAIT"
    }
    if (grepl("PBS57", ctTable15[i,"Name"], fixed=TRUE)){
      ctTable15[i, "cell.type"] <- "NKT"
    }
  }
  ctTable15<- subset(ctTable15, cell.type != "NKT")
  ctTable15$FRorFZ<-"FR"

  ### read in data
  ## probes: RUS011_Plate 1
  all_content<- readLines(paste(dataDir, "FR_RNA_UCM_HU_20210312_P1_1362576007.csv", sep=""))
  ctTable16 <- read.csv(textConnection(all_content[-c(1:11)]),
                        header=TRUE, stringsAsFactors=FALSE)
  #ctTable16 <- ctTable16[-c(1:9),]
  ctTable16$Comments<- NA
  ctTable16$cellSource <- NA
  for(i in 1:nrow(ctTable16)) {
    if (grepl("RUS011", ctTable16[i,2], fixed=TRUE)){
      ctTable16[i, "probe"] <- "RUS011"
    }
    if (grepl("Empty", ctTable16[i,2], fixed=TRUE)){
      ctTable16[i, "probe"] <- "Empty"
    }
  }
  ctTable16<- subset(ctTable16, probe != "Empty")
  ctTable16$age <- NA
  for(i in 1:nrow(ctTable16)) {
    if (grepl("5OPRU", ctTable16[i,"Name"], fixed=TRUE)){
      ctTable16[i, "cell.type"] <- "MAIT"
    }
    if (grepl("PBS57", ctTable16[i,"Name"], fixed=TRUE)){
      ctTable16[i, "cell.type"] <- "NKT"
    }
  }
  ctTable16<- subset(ctTable16, cell.type != "NKT")
  ctTable16$FRorFZ<-"FR"

  ### read in data
  ## probes: OXFGI4330-OXFGFI4870
  all_content<- readLines(paste(dataDir, "FZ_RNA_UCM_HU_20200917_P1_1362537616.csv", sep=""))
  ctTable17 <- read.csv(textConnection(all_content[-c(1:11)]),
                        header=TRUE, stringsAsFactors=FALSE)
  #ctTable17 <- ctTable17[-c(1:9),]
  ctTable17$Comments<- NA
  ctTable17$cellSource <- NA
  for(i in 1:nrow(ctTable17)) {
    if (grepl("OXFGI4330", ctTable17[i,2], fixed=TRUE)){
      ctTable17[i, "probe"] <- "OXFGI4330"
    }
    if (grepl("OXFGI4870", ctTable17[i,2], fixed=TRUE)){
      ctTable17[i, "probe"] <- "OXFGI4870"
    }
    if (grepl("Empty", ctTable17[i,2], fixed=TRUE)){
      ctTable17[i, "probe"] <- "Empty"
    }
  }
  ctTable17<- subset(ctTable17, probe != "Empty")
  ctTable17$age <- NA
  for(i in 1:nrow(ctTable17)) {
    if (grepl("5OPRU", ctTable17[i,"Name"], fixed=TRUE)){
      ctTable17[i, "cell.type"] <- "MAIT"
    }
    if (grepl("PBS57", ctTable17[i,"Name"], fixed=TRUE)){
      ctTable17[i, "cell.type"] <- "NKT"
    }
  }
  ctTable17<- subset(ctTable17, cell.type != "NKT")
  ctTable17$FRorFZ<-"FZ"

  ### read in data
  ## probes: OXF-IBD1322
  all_content<- readLines(paste(dataDir, "FZ_RNA_UCM_HU_20201023_P1_1362537452.csv", sep=""))
  ctTable18 <- read.csv(textConnection(all_content[-c(1:11)]),
                        header=TRUE, stringsAsFactors=FALSE)
  #ctTable12 <- ctTable12[-c(1:9),]
  ctTable18$Comments<- NA
  ctTable18$cellSource <- NA
  for(i in 1:nrow(ctTable18)) {
    if (grepl("OXFIBD1322", ctTable18[i,2], fixed=TRUE)){
      ctTable18[i, "probe"] <- "OXFIBD1322"
    }
    if (grepl("Empty", ctTable18[i,2], fixed=TRUE)){
      ctTable18[i, "probe"] <- "Empty"
    }
  }
  ctTable18<- subset(ctTable18, probe != "Empty")
  ctTable18$age <- NA
  for(i in 1:nrow(ctTable18)) {
    if (grepl("5OPRU", ctTable18[i,"Name"], fixed=TRUE)){
      ctTable18[i, "cell.type"] <- "MAIT"
    }
    if (grepl("PBS57", ctTable18[i,"Name"], fixed=TRUE)){
      ctTable18[i, "cell.type"] <- "NKT"
    }
  }
  ctTable18<- subset(ctTable18, cell.type != "NKT")
  ctTable18$FRorFZ<-"FZ"

  ### read in data
  ## probes: RUS007
  all_content<- readLines(paste(dataDir, "FR_RNA_UCM_HU_20201016_P1_1362537463.csv", sep=""))
  ctTable19 <- read.csv(textConnection(all_content[-c(1:11)]),
                        header=TRUE, stringsAsFactors=FALSE)
  #ctTable19 <- ctTable19[-c(1:9),]
  ctTable19$Comments<- NA
  ctTable19$cellSource <- NA
  for(i in 1:nrow(ctTable19)) {
    if (grepl("RUS007", ctTable19[i,2], fixed=TRUE)){
      ctTable19[i, "probe"] <- "RUS007"
    }
    if (grepl("Empty", ctTable19[i,2], fixed=TRUE)){
      ctTable19[i, "probe"] <- "Empty"
    }
  }
  ctTable19<- subset(ctTable19, probe != "Empty")
  ctTable19$age <- NA
  for(i in 1:nrow(ctTable19)) {
    if (grepl("5OPRU", ctTable19[i,"Name"], fixed=TRUE)){
      ctTable19[i, "cell.type"] <- "MAIT"
    }
    if (grepl("PBS57", ctTable19[i,"Name"], fixed=TRUE)){
      ctTable19[i, "cell.type"] <- "NKT"
    }
  }
  ctTable19<- subset(ctTable19, cell.type != "NKT")
  ctTable19$FRorFZ<-"FR"

  # ### read in data
  # ## probes: FR_RNA_T1D_MS_20210719_P1
  # all_content<- readLines(paste(dataDir, "FR_RNA_T1D_MS_20210719_P1_1362626235.csv", sep=""))
  # ctTable20 <- read.csv(textConnection(all_content[-c(1:11)]),
  #                       header=TRUE, stringsAsFactors=FALSE)
  # #ctTable19 <- ctTable19[-c(1:9),]
  # ctTable20$Comments<- NA
  # ctTable20$cellSource <- NA
  # for(i in 1:nrow(ctTable20)) {
  #   if (grepl("SEEP05", ctTable20[i,2], fixed=TRUE)){
  #     ctTable20[i, "probe"] <- "BigFoot"
  #   }
  #   if (grepl("Empty", ctTable20[i,2], fixed=TRUE)){
  #     ctTable20[i, "probe"] <- "Empty"
  #   }
  # }
  # ctTable20<- subset(ctTable20, probe != "Empty")
  # ctTable20$age <- NA
  # for(i in 1:nrow(ctTable20)) {
  #   if (grepl("12-20", ctTable20[i,"Name"], fixed=TRUE)){
  #     ctTable20[i, "cell.type"] <- "CD4T"
  #   }
  #   if (grepl("13-21", ctTable20[i,"Name"], fixed=TRUE)){
  #     ctTable20[i, "cell.type"] <- "CD4T"
  #   }
  # }
  # #ctTable20<- subset(ctTable20, cell.type != "NKT")
  # ctTable20$FRorFZ<-"FR"
  # 
  # ### read in data
  # ## probes: FR_RNA_T1D_MS_20210719_P2
  # all_content<- readLines(paste(dataDir, "FR_RNA_T1D_MS_20210719_P2_1362626161.csv", sep=""))
  # ctTable21 <- read.csv(textConnection(all_content[-c(1:11)]),
  #                       header=TRUE, stringsAsFactors=FALSE)
  # #ctTable19 <- ctTable19[-c(1:9),]
  # ctTable21$Comments<- NA
  # ctTable21$cellSource <- NA
  # for(i in 1:nrow(ctTable21)) {
  #   if (grepl("SEEP05", ctTable21[i,2], fixed=TRUE)){
  #     ctTable21[i, "probe"] <- "Astrios"
  #   }
  #   if (grepl("Empty", ctTable21[i,2], fixed=TRUE)){
  #     ctTable21[i, "probe"] <- "Empty"
  #   }
  # }
  # ctTable21<- subset(ctTable21, probe != "Empty")
  # ctTable21$age <- NA
  # for(i in 1:nrow(ctTable21)) {
  #   if (grepl("12-20", ctTable21[i,"Name"], fixed=TRUE)){
  #     ctTable21[i, "cell.type"] <- "CD4T"
  #   }
  #   if (grepl("13-21", ctTable21[i,"Name"], fixed=TRUE)){
  #     ctTable21[i, "cell.type"] <- "CD4T"
  #   }
  # }
  # #ctTable21<- subset(ctTable21, cell.type != "NKT")
  # ctTable21$FRorFZ<-"FR"
  
  ## everything
  ctTableCombine <- rbind(ctTable1, ctTable2, ctTable3, ctTable4, ctTable5, ctTable6, ctTable7, ctTable8, ctTable9)
  ctTableCombine <- rbind(ctTableCombine, ctTable10, ctTable11, ctTable12, ctTable13, ctTable14, ctTable15, ctTable16)
  ctTableCombine <- rbind(ctTableCombine, ctTable17, ctTable18, ctTable19) #, ctTable20, ctTable21)
  
  ## change column names
  names(ctTableCombine)[c(2, 5, 7)] <- c("cellType", "gene", "ct")
  
  ## rename cellSource column value
  ctTableCombine[grepl("Infl", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "tissue"
  ctTableCombine[grepl("-I-", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "tissue"       #SMS edit
  ctTableCombine[grepl("Uninf", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "tissue"
  ctTableCombine[grepl("-UI-", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "tissue"      #SMS edit
  ctTableCombine[grepl("Blood", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "blood"
  ctTableCombine[grepl("BL", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "blood"       #SMS edit
  ctTableCombine[grepl("EMPTY", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "control"
  ctTableCombine[grepl("SP", ctTableCombine$cellType, fixed=TRUE), "cellSource"] <- "tissue"
  
  ## create age (inflamed/uninflamed) column
  ctTableCombine[grepl("Infl", ctTableCombine$cellType, fixed=TRUE), "age"] <- "infl"
  ctTableCombine[grepl("-I-", ctTableCombine$cellType, fixed=TRUE), "age"] <- "infl"                #SMS edit
  ctTableCombine[grepl("Uninf", ctTableCombine$cellType, fixed=TRUE), "age"] <- "uninfl"
  ctTableCombine[grepl("-UI-", ctTableCombine$cellType, fixed=TRUE), "age"] <- "uninfl"             #SMS edit
  ctTableCombine[grepl("Blood", ctTableCombine$cellType, fixed=TRUE), "age"] <- "blood"
  ctTableCombine[grepl("BL", ctTableCombine$cellType, fixed=TRUE), "age"] <- "blood"              #SMS edit
  ctTableCombine[grepl("EMPTY", ctTableCombine$cellType, fixed=TRUE), "age"] <- "controlIslet"
  ctTableCombine[grepl("SP", ctTableCombine$cellType, fixed=TRUE), "age"] <- "Spleen"
  #print(colnames(ctTableCombine))
  write.csv(ctTableCombine, "dataLoadUCMAIT_script_line400.csv")
  return(ctTableCombine)
}







