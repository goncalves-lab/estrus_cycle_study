#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This script uses a seurat object ({estrus_tissue_processing.R}) 


library(Seurat)
library(plyr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

setwd(args[1])

t1 <- loadRData(args[2])

metadata <- data.frame("cell_type"=t1$level_2, "cell_name"=colnames(t1), "condition"=t1$estrus_phase, "tissue"=substr(t1$orig.ident, nchar(t1$orig.ident)-1, nchar(t1$orig.ident)))
counts <- t1[["SCT"]]@counts



pro_cyto <- c("Il1b", "Il6",  "Il12a", "Tnf", "Ifng", "Csf2","Il17a") 
anti_cyto <- c("Il4", "Il10", "Il11", "Tgfb1", "Tgfb3")
fibrosis <- c("Hgf", "Egf", "Fgf2", "Tgfb2", "Pdgfa", "Pdgfb", "Pdgfc", "Pdgfd")
cyto <- c(pro_cyto, anti_cyto, fibrosis)




lr_pair_df <-data.frame("lr_pair"="", "tissue"= "", "expr"=1, "cell_type"="")

for (i in unique(metadata$cell_type)){
  for (j in seq(from=1, to=length(cyto))){
    print(j)
    ligand <- cyto[j]
    metadata_L <- metadata[metadata$cell_type==i,]
    counts_l <- counts[,colnames(counts)%in%metadata_L$cell_name]
    counts_l <- as.data.frame(counts_l[rownames(counts_l) == ligand,])
    if (nrow(counts_l)>0){
      counts_l <- as.data.frame(t(counts_l))
      rownames(counts_l) <- ligand
      metadata_L <- metadata_L[match(colnames(counts_l), metadata_L$cell_name),]
      cl <- i
      original <- rowSums(counts_l)
      res_df <- data.frame("lr_pair"=ligand, "tissue"=unique(metadata$tissue), "expr"=original, "cell_type"=i)
      lr_pair_df <- rbind(lr_pair_df, res_df)
    }}}





lr_pair_df <- lr_pair_df[2:nrow(lr_pair_df),]
write.csv(lr_pair_df, file=args[3], row.names = F, quote = F)