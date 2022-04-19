#this cript uses list of DEG genes generated in the script ml_glm_model.R and the gene sets retrieved from MsigDb
#https://www.gsea-msigdb.org/gsea/msigdb/

###############################################

library(dplyr)
library("RColorBrewer")
library(Seurat)
library(ComplexHeatmap)
library(biomaRt)
library(superheat)
library(data.table)
library(stringr)
library(biomaRt)
library(ggplot2)
library(AUCell)
library(GSEABase)
library(easyGgplot2)
library(RColorBrewer)
library(viridis)
library(circlize)
library(ggrepel)


#pathway analysis

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}


#functions used to load gmt files
GeneSetListFromGmt <- function(path2gmt) {
  gene_sets <- lapply(readLines(path2gmt), str_split, "\t")
  set_list=lapply(gene_sets, handle_gene_sets) 
  return(set_list)
}

handle_gene_sets <- function(gene_sets) {
  gene_set_vec <- gene_sets[[1]]
  if (length(gene_set_vec[2]) == 0) 
    description <- "None" 
  else
    description <- gene_set_vec[2]
  genes <- gene_set_vec[3:length(gene_set_vec)][str_length(gene_set_vec[3:length(gene_set_vec)]) > 0]
  gene_set_list <- list(id = gene_set_vec[1],
                        description = description,
                        genes = genes)
  return( gene_set_list )
}

#function used to sort genes as differentially regulated or not
thresholdeR <- function(x) {
  if (abs(x[1]) == Inf && !is.na(x[1])) {
    return(0)
  } 
  else {
    return(x[2])
  }
}

#function used to produce DF with ORA results
enricheR=function(gs_mmu_data, gene_df, bk){
  df_hyper <- data.frame("pathway"="","description"="","in_set_df"=1,"df"=1,"n_df"=1,"in_set"=1,"p_value"=1, check.names = F)
  for (j in seq(from=1, to =length(gs_mmu_data))){
    in_set_de_v=gene_df%in%gs_mmu_data[[j]][[3]]
    in_set_de_v=in_set_de_v[in_set_de_v==TRUE]
    num_vector_1=c(length(in_set_de_v), length(gene_df), bk-length(gene_df), length(gs_mmu_data[[j]][[3]]))
    #num_vector_2=c(length(in_set_de_v),length(df_1_de$entrezgene_id)-length(in_set_de_v),length(in_set_nde_v), length(df_1_nde$entrezgene_id)-length(in_set_nde_v))
    p_value_1=phyper(num_vector_1[1]-1, num_vector_1[2], num_vector_1[3], num_vector_1[4], lower.tail = FALSE, log.p = FALSE)
    num_vector_1[5]=p_value_1
    df_hyper=rbind(df_hyper,c(gs_mmu_data[[j]][[1]],gs_mmu_data[[j]][[2]],unname(as.list(num_vector_1))))
    df_hyper[,1]=as.character(df_hyper[,1])
    df_hyper[,2]=as.character(df_hyper[,2])
    
  }
  colnames(df_hyper)=c("pathway","description","in_set_df","df","n_df","in_set","p_value")
  df_hyper <- df_hyper[2:nrow(df_hyper),]
  df_hyper$p_value_adj=p.adjust(df_hyper$p_value, method = "BH", n = length(df_hyper$p_value))
  return(df_hyper)
  
}


gs_mmu_data <- GeneSetListFromGmt("c5.go.bp.v7.2.symbols.gmt")


file_names <- list.files(wd) 

genes_df <- data.frame()
for (i in file_names){
  genes <- read.csv(i)
  genes$genes <- rownames(genes)
  genes$tissue <- substr(i,1,2)
  genes_df <- rbind(genes, genes_df)
}

ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")

pathway_conserved_df <- data.frame()
for (i in unique(genes_df$tissue)){
  genes_tissue <- genes_df[genes_df$tissue==i,]
  genes_tissue <- genes_tissue[genes_tissue$p_value_adj_nb<0.05,]
  ensemnl_list <- getBM(attributes=c("hsapiens_homolog_associated_gene_name"), filters="external_gene_name", values = genes_tissue$genes, mart=ensembl)
  pathway_conserved <- enricheR(gs_mmu_data, ensemnl_list$hsapiens_homolog_associated_gene_name, 18500)
  pathway_conserved <- pathway_conserved[pathway_conserved$p_value_adj<0.05,]
  pathway_conserved$gene_ratio <- pathway_conserved$in_set_df/pathway_conserved$df
  pathway_conserved$tissue <- i
  pathway_conserved_df <- rbind(pathway_conserved, pathway_conserved_df)
}

