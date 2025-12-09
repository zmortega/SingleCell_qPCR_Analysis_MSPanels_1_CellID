mergeAllSheets <- function(sheetList = c("1511-250-149", "1511-250-159", "1511-250-158", "1511-250-157", "1511-250-000")) {
  allSheetsDF <- data.frame()
  for (sheet in sheetList){
    merged <- mergeSheet(sheet)
    colnames(merged)[3] <- "Barcode_ID"
    allSheetsDF <- rbind(merged, allSheetsDF)
  }
  
  ## return a data frame with all of the barcodes from all of the sheets
  return(allSheetsDF)
}

mergeSheet <- function(BCID) {
  ## run readPlateXlsx to format the data from one sheet in a temporary csv and use plater library 
  nameBC <- as.data.frame(readPlateXlsx(BCID))
  ## separate the Barcode ID
  for (c in colnames(nameBC)) {
    if (grepl("Barcode_ID", c)) {
      fullcoln <- c
      barcodeID <- sub(".*Barcode_ID:", "", c)
    }
  }
  
  # add a column to the data frame that combines the barcode ID for the plate with each cell barcode
  newBC <- c()
  for (bc in nameBC[fullcoln]) {
    fullBC <- paste(barcodeID, bc, sep = "_")
    newBC <- append(newBC, fullBC)
    # print(bc)
  }
  
  ## replace the column that contained the barcode with the column of new barcodes
  nameBC[fullcoln] <- newBC
  
  ## returns a data frame with the cell names in one column and the full barcode in a second column
  return(nameBC)
}

readPlateXlsx <- function(sheetBC) {
  library("plater")
  library("readxl")
  
  ## loads one sheet of the excel spreadsheet
  dataDir <- paste(getwd(), "Logs-UCM-AccessArray_Layout.xlsx", sep="/")
  plateSheet <- as.data.frame(read_excel(dataDir, sheet=sheetBC), optional = TRUE)
  

  ## create csv file in correct format for the plater library
  csvDir <- paste(getwd(), "tempSheet.csv", sep="/")

  write.csv(plateSheet, csvDir, row.names = FALSE, na ="", quote = FALSE)
  
  cellNameDF <- as.data.frame(read_plate(file=csvDir, sep = ","))

  unlink(csvDir)
  return(cellNameDF)
}

mergeAllSheets()  