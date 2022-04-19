#!/usr/bin/env Rscript

#this script uses two Seurat objects (young and old) (estrus_{tissue}_processing.R and 18m_{tissue}_processing.R)
args = commandArgs(trailingOnly=TRUE)
library(Seurat)
library(EntropyExplorer)
library(plyr)
library(data.table)
library(stringr)
library(dplyr)


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

setwd(args[1])
t_combined <- loadRData(args[2])
t_combined_old <- loadRData(args[3])
counts <- t_combined[["SCT"]]@counts
counts_old <- t_combined_old[["SCT"]]@counts

metadata <- data.frame("cell_types"=t_combined$level_2, "condition"=t_combined$estrus_phase)
metadata_old <-  data.frame("cell_types"=t_combined_old$level_2, "condition"="18m")

EntropeR <- function(df_1, df_2){
  ent <- EntropyExplorer(df_1, df_2, "dse", "bu", shift = c("auto","auto"))
  ent <- as.data.frame(ent)
  return(ent)
}


for (i in(unique(metadata$cell_types))){
  df_cells <- metadata[metadata$cell_types==i,]
  df_1_metadata <- df_cells[df_cells$condition=="diestrus",]
  df_cells_old <- metadata_old[metadata_old$cell_types==i,]
  df_2_metadata <- df_cells_old[df_cells_old$condition=="18m",]
  df_1 <- as.data.frame(counts[,colnames(counts)%in%rownames(df_1_metadata)])
  df_2 <-as.data.frame(counts_old[,colnames(counts_old)%in%rownames(df_2_metadata)])
  df_1 <- df_1[rownames(df_1)%in%rownames(df_2),]
  df_2 <- df_2[rownames(df_2)%in%rownames(df_1),]
  df_2 <- df_2[match(rownames(df_1), rownames(df_2)),]
  
  
  if (nrow(df_1)&nrow(df_2)>0){
    ent_dp <- EntropeR(df_1, df_2)
    write.csv(ent_dp, file=paste0(i, "_18m_ShE"), quote=F)
  }
}

