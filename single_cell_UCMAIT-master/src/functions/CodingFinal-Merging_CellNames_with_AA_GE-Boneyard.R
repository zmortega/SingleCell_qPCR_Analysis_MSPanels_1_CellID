####Testing out Plater R package on our dataset
library(plater)
library(stringr)
library(gtools)
library(gdata)
library(readxl)
library(dplyr)
library(plyr)
library(reshape)

#Load in AA plate layouts
setwd("/Users/sid/Documents/OneDrive - Lab files/OneDrive - The Scripps Research Institute/Laboratory Work/TCR_Sequencing_Results/")

AA_Directory<- list.files("./Log-Samples_Submitted_to_Genomics_Core/Log-T1D-All_AccessArray_PlateLayouts/", pattern = glob2rx("*.csv"))
#Test<- as.data.frame(read_plate(file = "/Users/sid/Documents/OneDrive - Lab files/OneDrive - The Scripps Research Institute/Laboratory Work/"))

AccessArray<- lapply(AA_Directory, function (x) {
  AA<- as.data.frame(read_plate(paste("./Log-Samples_Submitted_to_Genomics_Core/Log-T1D-All_AccessArray_PlateLayouts/",x,sep = "")))
  AA$AccessArray_ID<- str_split(x, ".csv")[[1]][1]
  columnnames<- colnames(AA)
  AA$Date<- str_split(columnnames[4], "_")[[1]][3]
  names(AA)<-c("Wells","AA_Updated_Cell_Name","AA-Original_Cell_Name","Cell_Barcode","AA_ChipBarcode","Date_of_AA")
  AA$AA_Combined_Barcode<- sapply(1:nrow(AA), function (y) {paste(AA[y,"AA_ChipBarcode"], AA[y,"Cell_Barcode"], sep = "_")})
  #print(AA)
  #names(AA)<-c("Wells","AA-Updated_Cell_name","AA-Original_Cell_Name","Cell_Barcode","AA_ChipBarcode","Date_of_AA", "Combined_Barcode")
  AA
})
AccessArrayLayout<- bind_rows(AccessArray)
AccessArrayLayout<- select(AccessArrayLayout, -"Wells")
#write.csv(AccessArrayLayout, paste0("Log-All_AccessArray-", Sys.Date(),"-Rearranged_into_a_DataFrame.csv"))

#Load in TCR list and merge the cell names and barcodes (using the AA ID and Barcodes) to put TCR segments with cell names
T1D_TCRa_data<- read.xls("/Users/sid/Documents/OneDrive - Lab files/OneDrive - The Scripps Research Institute/Laboratory Work/TCR_Sequencing_Results/T1D_project/Human_TCRs/ALL_T1D_TCRA_TCRB_Genes_Master_List_Hu.xlsx", sheet=1)
T1D_TCRb_data<- read.xls("/Users/sid/Documents/OneDrive - Lab files/OneDrive - The Scripps Research Institute/Laboratory Work/TCR_Sequencing_Results/T1D_project/Human_TCRs/ALL_T1D_TCRA_TCRB_Genes_Master_List_Hu.xlsx", sheet=2)

TCR<- list(T1D_TCRa_data, T1D_TCRb_data)
names(TCR)<- c("Alpha", "Beta")

TCR_updated<- lapply(1:length(TCR), function (x) {
  TCR[[x]]$Cell_Barcode<- sapply(1:nrow(TCR[[x]]), function (y) {str_split(TCR[[x]]$Sequence.ID[y], "_")[[1]][1]})
  TCR[[x]]$AA_Combined_Barcode<- sapply(1:nrow(TCR[[x]]), function (z) {paste(TCR[[x]][z,"AA_ChipBarcode"],TCR[[x]][z,"Cell_Barcode"], sep = "_")})
  TCR[[x]]
})
names(TCR_updated)<- c("Alpha", "Beta")

#Merge CellNames with TCR information
TCRa_merged<- merge(x=AccessArrayLayout, y=TCR_updated[["Alpha"]], by.x = "AA_Combined_Barcode", by.y = "AA_Combined_Barcode")
RemoveColumns<- grep("X", colnames(TCRa_merged))
TCRa_merged<- select(TCRa_merged, -c(RemoveColumns))

TCRb_merged<- merge(x=AccessArrayLayout, y=TCR_updated[["Beta"]], by.x = "AA_Combined_Barcode", by.y = "AA_Combined_Barcode")
RemoveColumns<- grep("X", colnames(TCRb_merged))
TCRb_merged<- select(TCRb_merged, -c(RemoveColumns))

###Adding in the cohort information into the appropriate column in the TCR list above. YOU SHOULD REWRITE THE ABOVE TO REDUCE THE AMOUNT OF HAND CURATING THAT IS CURRENTLY REQUIRED.
# Cohort_Info<- read.xls("/Users/sid/Documents/OneDrive - Lab files/OneDrive - The Scripps Research Institute/Laboratory Work/Projects T1D/Log-T1D-List_of_Donors_with_Cohort.xlsx")
# 
# for (i in 1:nrow(TCR_updated[["Alpha"]])) {
#   if (TCR_updated[["Alpha"]][i, ""])
# }

#Time to slim down the TCR data to just the basics needed for merging with transcriptomic data.
##Looking to only retain TRAV/TRBV, TRBD, TRAJ/TRBJ, CDR3A, CDR3B, with full cell name. There will be two TRAV columns (TRAV1 and TRAV2) to accoomodate cells with two alpha chains.
###Perhaps there will also be two TRBV columns too....this happens so rarely (and because suspicion is that two TRBVs is a doublet, this probably shouldn't be included)

####Starting with ALPHA chain, stripping out duplicates, rearranging, and then rbinding back to original DF
TCRa_merged_sub<- TCRa_merged[,c("AA_Combined_Barcode", "AA_Updated_Cell_Name", "V.domain.functionality", "Productive.after.fixing", "V.GENE.and.allele", "J.GENE.and.allele", "D.GENE.and.allele", "CDR3.IMGT.length", "AA.Junction")]
TCRa_merged_sub<- filter(TCRa_merged_sub, Productive.after.fixing == "Productive-after-fixing" | V.domain.functionality == "productive")
TCRa_merged_sub$TRAV_GENE<- sapply(1:nrow(TCRa_merged_sub), function (x) {
  str_split(str_split(TCRa_merged_sub[x,"V.GENE.and.allele"], " ")[[1]][2], "\\*")[[1]][1]
})

TCRa_merged_sub$TRAJ_GENE<- sapply(1:nrow(TCRa_merged_sub), function (x) {
  str_split(str_split(TCRa_merged_sub[x,"J.GENE.and.allele"], " ")[[1]][2], "\\*")[[1]][1]
})

TCRa_reduced<- TCRa_merged_sub[,c(1,2,8:11)]

####Below is the code for taking the duplicate TCRb that are in multiple rows and moving them into one row. The goal is to make the data from one cell
###all sit in one row...That's the goal anyways
TCRa_dupl<- data.frame(table(as.factor(TCRa_reduced$AA_Combined_Barcode)))
TCRa_dupl<- subset(TCRa_dupl, TCRa_dupl$Freq > 1)

ListDuplicates_A<- lapply(1:nrow(TCRa_dupl), function (x) {
  z<- str_which(TCRa_reduced$AA_Combined_Barcode, paste0(TCRa_dupl$Var1[x],"$"))
    Temp<- TCRa_reduced[c(z),]
    Temp
})

####put this in a for-loop
for (i in 1:length(ListDuplicates_A)) {
  x<- nrow(ListDuplicates_A[[i]])
  if (x<3) {
    ListDuplicates_A[[i]][1,c(7:12)]<- ListDuplicates_A[[i]][2,c(1:6)]
  } else if (x >= 3) {
    rows<- nrow(ListDuplicates_A[[i]])
    print(paste(rows, i, sep = "_"))
    for (y in 2:rows) {
      ColNum<- ncol(ListDuplicates_A[[1]])
      ListDuplicates_A[[i]][1,c(((6*y-5):(6*y)))]<- ListDuplicates_A[[i]][y,c(1:6)]
    }
  }
  if (x< 3) {
    ListDuplicates_A[[i]]<- ListDuplicates_A[[i]][-2,]
  } else {
    ListDuplicates_A[[i]]<- ListDuplicates_A[[i]][-c(2:x),]
  }
}
ListDuplicates_A<- do.call(rbind.fill, ListDuplicates_A)
TCRa_reduced<- rbind.fill(TCRa_reduced, ListDuplicates_A)
write.csv(TCRa_reduced, paste0("Log-TCRa_data_merged_with_Full_Cell_Names-", Sys.Date(),".csv"))

####TCRb-taking the duplicates out, rearranging the rows into columns and then rbinding it back to the original DF and then writing 
TCRb_merged_sub<- TCRb_merged[,c("AA_Combined_Barcode", "AA_Updated_Cell_Name", "V.domain.functionality", "Productive.after.fixing", "V.GENE.and.allele", "J.GENE.and.allele", "D.GENE.and.allele", "CDR3.IMGT.length", "AA.Junction")]
TCRb_merged_sub<- TCRb_merged_sub %>% filter(Productive.after.fixing == "Productive-after-fixing" | V.domain.functionality == "productive")
TCRb_merged_sub$TRBV_GENE<- sapply(1:nrow(TCRb_merged_sub), function (x) {
  str_split(str_split(TCRb_merged_sub[x,"V.GENE.and.allele"], " ")[[1]][2], "\\*")[[1]][1]
})
TCRb_merged_sub$TRBD_GENE<- sapply(1:nrow(TCRb_merged_sub), function (x) {
  str_split(str_split(TCRb_merged_sub[x,"D.GENE.and.allele"], " ")[[1]][2], "\\*")[[1]][1]
})
TCRb_merged_sub$TRBJ_GENE<- sapply(1:nrow(TCRb_merged_sub), function (x) {
  str_split(str_split(TCRb_merged_sub[x,"J.GENE.and.allele"], " ")[[1]][2], "\\*")[[1]][1]
})

TCRb_reduced<- TCRb_merged_sub[,c(1,2,8:12)]

####Below is the code for taking the duplicate TCRb that are in multiple rows and moving them into one row. The goal is to make the data from one cell
###all sit in one row...That's the goal anyways
TCRb_dupl<- data.frame(table(as.factor(TCRb_merged_sub$AA_Combined_Barcode)))
TCRb_dupl<- subset(TCRb_dupl, TCRb_dupl$Freq > 1)
 
ListDuplicates<- lapply(1:nrow(TCRb_dupl), function (x) {
   z<- str_which(TCRb_reduced$AA_Combined_Barcode, paste0(TCRb_dupl$Var1[x],"$"))
   Temp<- TCRb_reduced[c(z),]
   Temp
})
 
####put this in a for-loop
for (i in 1:length(ListDuplicates)) {
 x<- nrow(ListDuplicates[[i]])
 if (x<3) {
   ListDuplicates[[i]][1,c(8:14)]<- ListDuplicates[[i]][2,c(1:7)]
   } else if (x >= 3) {
   rows<- nrow(ListDuplicates[[i]])
   print(paste(rows, i, sep = "_"))
   for (y in 2:rows) {
     ColNum<- ncol(ListDuplicates[[1]])
     ListDuplicates[[i]][1,c(((7*y-6):(7*y)))]<- ListDuplicates[[i]][y,c(1:7)]
    }
   } 
 if (x< 3) {
   ListDuplicates[[i]]<- ListDuplicates[[i]][-2,]
   } else {
   ListDuplicates[[i]]<- ListDuplicates[[i]][-c(2:x),]
   }
}
ListDuplicates<- do.call(rbind.fill, ListDuplicates)
TCRb_reduced<- rbind.fill(TCRb_reduced, ListDuplicates)
write.csv(TCRb_reduced, paste0("Log-TCRb_data_merged_with_Full_Cell_Names-", Sys.Date(),".csv"))
