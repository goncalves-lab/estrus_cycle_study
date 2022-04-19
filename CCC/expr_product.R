#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This script uses LR databases of CellChat and CelltalkDB databases and a seurat object ({estrus_tissue_processing.R}) 


library(Seurat)
library(plyr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

setwd(args[1])
load(args[2])

t1 <- loadRData(args[3])

metadata <- data.frame("cell_type"=t1$level_2, "cell_name"=colnames(t1), "condition"=t1$estrus_phase, "tissue"=substr(t1$orig.ident, nchar(t1$orig.ident)-1, nchar(t1$orig.ident)))
counts <- t1[["SCT"]]@counts
rownames(counts) <- toupper(rownames(counts))

keeper <- c()
for (j in seq(nrow(mouse_db))){
  lr_pair <- mouse_db[j,]
  lr <- strsplit(rownames(lr_pair)[1], "_")
  receptor <- lr[[1]][2:length(lr[[1]])]
  ligand <- lr[[1]][1]
  a <- c(ligand%in%rownames(counts), receptor%in%rownames(counts))
  if (!F%in%a){
    keeper <- c(keeper, j)
  }
}

mouse_db <- mouse_db[keeper,]



lr_pair_df <-data.frame("lr_pair"="", "tissue"= "", "expr"=1)
for (j in seq(from=1, to=nrow(mouse_db))){
  print(j)
  lr_pair <- mouse_db[j,]
  lr <- strsplit(rownames(lr_pair)[1], "_")
  receptor <- lr[[1]][2:length(lr[[1]])]
  ligand <- lr[[1]][1]
  if (args[4] == "R"){
    metadata_LR <- metadata[metadata$cell_type=="F",]
    counts_r <- counts[,colnames(counts)%in%metadata_LR$cell_name]
    counts_r <- as.data.frame(counts_r[rownames(counts_r)%in%receptor,])
    if (length(lr[[1]]) == 2){
      counts_r <- as.data.frame(t(counts_r))
      rownames(counts_r) <- receptor
      
    }
    
    if (length(lr[[1]])>2){
      receptor_min <- names(rowMeans(counts_r))[which.min(rowMeans(counts_r))]
      counts_r <- as.data.frame(counts_r[rownames(counts_r) == receptor_min,])
    }
    counts_l <- as.data.frame(counts[rownames(counts) == ligand,])
    counts_l <- as.data.frame(t(counts_l))
    rownames(counts_l) <- ligand
    metadata_LR <- metadata_LR[match(colnames(counts_r), metadata_LR$cell_name),]
    metadata <- metadata[match(colnames(counts_l), metadata$cell_name),]
    cr <- "F"
    cl <- "other"
  }
  if (args[4] == "L"){
    counts_r <- as.data.frame(counts[rownames(counts)%in%receptor,])
    if (length(lr[[1]]) == 2){
      counts_r <- as.data.frame(t(counts_r))
      rownames(counts_r) <- receptor
    }
    
    if (length(lr[[1]])>2){
      receptor_min <- names(rowMeans(counts_r))[which.min(rowMeans(counts_r))]
      counts_r <- as.data.frame(counts_r[rownames(counts_r) == receptor_min,])
    }
    metadata_LR <- metadata[metadata$cell_type=="F",]
    counts_l <- counts[,colnames(counts)%in%metadata_LR$cell_name]
    counts_l <- as.data.frame(counts_l[rownames(counts_l) == ligand,])
    counts_l <- as.data.frame(t(counts_l))
    rownames(counts_l) <- ligand
    metadata_LR <- metadata_LR[match(colnames(counts_l), metadata_LR$cell_name),]
    metadata <- metadata[match(colnames(counts_r), metadata$cell_name),]
    cl <- "F"
    cr <- "other"
  }      
  
  original <- rowMeans(counts_r)*rowMeans(counts_l)
  res_df <- data.frame("lr_pair"=row.names(lr_pair)[1], "tissue"=unique(metadata$tissue), "expr"=original)
  lr_pair_df <- rbind(lr_pair_df, res_df)}  




lr_pair_df <- lr_pair_df[2:nrow(lr_pair_df),]
write.csv(lr_pair_df, file=args[5], row.names = F, quote = F)