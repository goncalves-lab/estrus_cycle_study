#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#This script uses LR databases of CellChat and CelltalkDB databases and two seurat object ({estrus_tissue_processing.R}) for which LR
# significance is to be calculated


library(Seurat)
library(plyr)

permutation.test <- function(original, treatment_l, treatment_r, outcome_l, outcome_r, n, ct_1, ct_2){
  distribution=c()
  result=0
  for(i in 1:n){
    sample_l <- sample(treatment_l, length(treatment_l), FALSE)
    expression_l <- data.frame("phase" = sample_l, "expression" = outcome_l)
    expression_l_1 <- expression_l[expression_l$phase==ct_1,]
    expression_l_2 <- expression_l[expression_l$phase==ct_2,]
    sample_r <- sample(treatment_r, length(treatment_r), FALSE)
    expression_r <- data.frame("phase" = sample_r, "expression" = outcome_r)
    expression_r_1 <- expression_r[expression_r$phase==ct_1,]
    expression_r_2 <- expression_r[expression_r$phase==ct_2,]
    distribution[i]= log2(((mean(expression_r_1$expression)*mean(expression_l_1$expression))+0.01)/
                            ((mean(expression_r_2$expression)*mean(expression_l_2$expression))+0.01))
    
  }
  result=sum(abs(distribution) >= abs(original))/(n)
  
  return(list(result, distribution))
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

setwd(args[1])
load(args[2])

t1 <- loadRData(args[3])
t2 <- loadRData(args[4])
t2 <- RenameCells(t2, "t2")

print("before merging")
t_combined <- merge(t1, y=c(t2), project = "estrus", merge.data = T)
if (args[8] == "young"){
  metadata <- data.frame("cell_type"=t_combined$level_2, "cell_name"=colnames(t_combined), "condition"=t_combined$estrus_phase, "tissue"=substr(t_combined$orig.ident, nchar(t_combined$orig.ident)-1, nchar(t_combined$orig.ident)))
} else {metadata <- data.frame("cell_type"=t_combined$level_2, "cell_name"=colnames(t_combined), "condition"="18m", "tissue"=substr(t_combined$orig.ident, nchar(t_combined$orig.ident)-1, nchar(t_combined$orig.ident)))}



counts <- t_combined[["SCT"]]@counts
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

if (args[7]=="subset"){
  pro_cyto <- c("Il1b", "Il6", "Cxcl8", "Il12a", "Tnf", "Ifng", "Csf2") 
  anti_cyto <- c("Il4", "Il10", "Il11", "Tgfb1", "Tgfb3")
  
  cyto <- c(pro_cyto, anti_cyto)
  
  mouse_db <- mouse_db[mouse_db$ligand%in%cyto,]}


lr_pair_df <-data.frame("lr_pair"="", "tissue_1"= "", "tissue_2"="", 
                        "log2expr"=1, "p_value"=1)
for (j in seq(from=1, to=nrow(mouse_db))){
  print(j)
  lr_pair <- mouse_db[j,]
  lr <- strsplit(rownames(lr_pair)[1], "_")
  receptor <- lr[[1]][2:length(lr[[1]])]
  ligand <- lr[[1]][1]
  if (args[5] == "R"){
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
  if (args[5] == "L"){
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
  
  t <- unique(metadata$tissue)
  metadata_t1 <- metadata[metadata$tissue==t[1],]
  metadata_t2 <- metadata[metadata$tissue==t[2],]
  counts_r_t1 <- counts_r[,colnames(counts_r)%in%metadata_t1$cell_name]
  counts_r_t2 <- counts_r[,colnames(counts_r)%in%metadata_t2$cell_name]
  counts_l_t1 <- counts_l[,colnames(counts_l)%in%metadata_t1$cell_name]
  counts_l_t2 <- counts_l[,colnames(counts_l)%in%metadata_t2$cell_name]
  
  
  original <- log2(((rowMeans(counts_r_t1)*rowMeans(counts_l_t1))+0.01)/
                     ((rowMeans(counts_r_t2)*rowMeans(counts_l_t2))+0.01))
  if (original != 0){
    if (args[5] == "R"){
      res <- permutation.test(original, metadata$tissue, metadata_LR$tissue, as.matrix(counts_l)[1,], as.matrix(counts_r)[1,], 1000, t[1], t[2])}
    if (args[5] == "L"){
      res <- permutation.test(original, metadata_LR$tissue, metadata$tissue, as.matrix(counts_l)[1,], as.matrix(counts_r)[1,], 1000, t[1], t[2])}
    
    res_df <- data.frame("lr_pair"=row.names(lr_pair)[1], "tissue_1"=t[1], 
                         "tissue_2"=t[2], "log2expr"=original, "p_value"=res[[1]])
    lr_pair_df <- rbind(lr_pair_df, res_df)
  }
  else {res_df <- data.frame("lr_pair"=row.names(lr_pair)[1], "tissue_1"=t[1], 
                             "tissue_2"=t[2], "log2expr"=0, "p_value"=NA)
  lr_pair_df <- rbind(lr_pair_df, res_df)}  
} 



lr_pair_df <- lr_pair_df[2:nrow(lr_pair_df),]
lr_pair_df$p_value_adj <- p.adjust(lr_pair_df$p_value, method = "BH",n=length(lr_pair_df$p_value))
write.csv(lr_pair_df, file=args[6], row.names = F, quote = F)
