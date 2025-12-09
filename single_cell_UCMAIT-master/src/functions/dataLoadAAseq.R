dataLoadAAseq <- function(loadSheet) {
  value <- 3
  dataLoaded <- c()
  if (loadSheet == "Alpha") {
    dataLoaded <- list(dataLoadAAseq_subset("Alpha"))
  }else if (loadSheet == "Beta") {
    dataLoaded <- list(dataLoadAAseq_subset("Beta"))
  }else if (loadSheet == "Both") {
    dataLoaded <- list(dataLoadAAseq_subset("Alpha"), dataLoadAAseq_subset("Beta"))
  }

  ## combines the AA barcode with the sequence ID 
  allDataFullBC <- list() # create a list where the data frames will be stored once the combined barcode is made
  for (sheet in dataLoaded) {
    newBC <- c()
    currDF <- as.data.frame(sheet)
    for (i in 1:nrow(currDF)) {
      barcodeID <- currDF$AccessArray_Barcode[i]
      seqID <- sub("_.*", "", currDF$Sequence.ID[i]) # take only the "BC#" portion of the Sequence.ID
      fullBC <- paste(barcodeID, seqID, sep = "_") # combine the barcodeID with the BC# portion of the seqID
      newBC <- append(newBC, fullBC) #add it to a vector that will be used to replace the old column
    }
    
    ## replace the column that contained the barcode with the column of new barcodes
    currDF$AccessArray_Barcode <- newBC
    allDataFullBC <- append(allDataFullBC, currDF)
  }
  
  # combineCellName(dataLoaded)
  
  
  ## returns a list of data frames that contain the subsetted data for either the Alpha, Beta or Both
  return(allDataFullBC)
}

## a function that matches the access array barcode to the name from the plate and adds the full name to a new column in the loaded data
combineCellName <- function(dataLoaded, cellNamesPlate) { 
  matchVect <- match(dataLoaded[["AccessArray_Barcode"]], cellNamesPlate[["Barcode_ID"]])
  matchedNames <- c()
  for (i in 1:length(matchVect)) {
    if (is.na(matchVect[i])) { # carry over the NA if there is no matching name found in the list
      matchedNames <- append(matchedNames, NA) 
    }else { # if there was a match, use the index found to pair the name with the barcode
      matchedNames <- append(matchedNames, cellNamesPlate[matchVect[i], "Cell_Name"])
    }
  }
  dataLoaded$Cell_Name <- matchedNames
  return(dataLoaded)
}

dataLoadAAseq_subset <- function(loadHelp){
  ##libraries
  library("readxl")
  
  ## set up data directory
  dataDir <- paste(getwd(), "ALL_MAIT_TCRA_TCRB_Genes_Master_List-Hu.xlsx", sep="/")
  # print(dataDir)
  
  ## Read the alpha sheet
  all_data <-read_excel(dataDir, sheet=loadHelp )
  coln_names <- c("AccessArray_Barcode", "Sequence.ID", "V.domain.functionality", "J.GENE and allele", "AA.Junction")
  
  ## create a data set that only includes the relevant data
  refinedData <- data.frame(all_data[coln_names])
  # print(dim(refinedData)) #check the df
  ## subset only the productive TCRs
  prodTCRs <- data.frame(matrix(ncol=length(coln_names), nrow=0))
  colnames(prodTCRs) <- coln_names
  for (i in 1:nrow(refinedData)) {
    if (refinedData[i, "V.domain.functionality"] == "productive"){
      prodTCRs[nrow(prodTCRs) + 1,] = refinedData[i,]
    }
  }
  # print(dim(prodTCRs)) # check the df
  # print(is.data.frame(prodTCRs))
  return(prodTCRs)
}
# dataLoadAAseq("Alpha")
str(combineCellName(dataLoadAAseq("Alpha"), mergeAllSheets(c("1511-250-149")))) 
