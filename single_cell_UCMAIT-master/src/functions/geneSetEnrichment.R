#### gene set enrichment ####
#### create gene sets and enrichment function ####
geneGroups <- list(chemokines=c("Ccr1", "Ccr7", "Ccr2", "Ccr3", "Ccr4", "Ccr5", "Ccr6", "Cxcr3", "Cxcr4", "Cxcl10"), 
                   cell_surface_receptors=c("Tnfrsf1a", "Tnfrsf1b", "Il.7r", "Il.5ra", "Il.4ra", "Ly6e", "Ifngr1", "Il.12rb", "Il.18r1", "Il.27r", "Icam1", "Cd44", "Ceacam1", "Il.1r2", "Cd28", "Cd3e", "Cd4", "Cd40", "Cd80", "Cd86", "Cd8a", "Ctla4", "Pdl.1"), 
                   signaling_molecules=c("Stat1", "Stat3", "Stat4", "Stat5", "Map2k6", "Mapk8", "Fyn", "Jak1", "Jak2", "Pten"), 
                   cytokines=c("Ifng", "Il.10", "Il.12b", "Il.17A", "Il.2", "IL.21", "Il.25", "Il.27", "Il.3", "Il.4", "Il.5", "Il.6", "Il.7", "Tnf"), 
                   transcription_factors=c("Bcl6", "Foxp3", "Gata4", "Ppara", "Pparg", "Ppargc1a", "Tbx21", "Bcl2", "Nfkb1", "Nur77"), 
                   metabolism=c("Gsk3a", "Gsk3b", "Hprt"), 
                   interferon_response=c("Ifi44", "Ifi44l", "Ifit1", "Ifit3", "Irf1", "Irf2", "Irf4", "Irf7", "Isg15", "Mx1"))

gse <- function(geneList, numSigGenes){
  geneSets <- as.vector(names(geneGroups))
  significantGenes <- geneList[1:numSigGenes]
  gsePvalues <- NULL
  overlaps <- NULL
  drawTest <- length(significantGenes)
  if(drawTest > 95){
    draw <- 95
  }else{
    draw <- drawTest
  }
  for(i in 1:length(names(geneGroups))){
    overlaps[i] <- paste(": ", sum(significantGenes %in% geneGroups[[i]]),
                         " / ", length(geneGroups[[i]]), " genes", sep="")
    gsePvalues[i] <- as.numeric(1 - phyper(sum(significantGenes %in% geneGroups[[i]]) - 1, 
                                           length(geneGroups[[i]]), 95-length(geneGroups[[i]]), draw))
  }
  gsePvalues <- round(gsePvalues, digits=3)
  gseResults <- cbind(geneSets, overlaps, gsePvalues)
  print("                 ", row.names=F, quote=F)
  print("                 ", row.names=F, quote=F)
  print("Gene set enrichment results", row.names=F, quote=F)
  print(gseResults, row.names=F, quote=F)
}