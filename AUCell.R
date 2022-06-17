#this script uses gene sets retrieved from MsigDb and seurat object generated in script {condition}_{tissue}_processing.R


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
library(gplots)


ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")

pathway_mmu <- read.table("/omics/groups/OE0433/internal/ivana/18m/DEG/pathways_mmu", sep=";", header=T)
pathway_mmu$genes <- toupper(pathway_mmu$genes)
geneSets_final <- list()
for (i in seq(nrow(pathway_mmu))){
  pathway=strsplit(pathway_mmu[i,2], ",")
  pathway=unique(pathway[[1]])
  ensemnl_list <- getBM(attributes=c("mmusculus_homolog_associated_gene_name", "external_gene_name"), filters="external_gene_name", values = pathway, mart=ensembl)
  for (j in seq(nrow(ensemnl_list))){
    if (ensemnl_list[j,1] == ""){
      ensemnl_list[j,1] = paste0(substr(ensemnl_list[j,2],1,1), substr(tolower(ensemnl_list[j,2]), 2, nchar(ensemnl_list[j,2])))
    }
  }
  geneSets <- GeneSet(unique(ensemnl_list$mmusculus_homolog_associated_gene_name), setName=pathway_mmu[i,1])
  geneSets_final <- c(geneSets, geneSets_final)
}


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


gs_mmu_data <- GeneSetListFromGmt("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/gene_sets/c5.go.bp.v7.2.symbols.gmt")
pathway_selected <- c("GO_NEGATIVE_REGULATION_OF_CELL_DEATH", "GO_EMBRYO_IMPLANTATION", "GO_POSITIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION",
                      "GO_RESPONSE_TO_CYTOKINE", "GO_WOUND_HEALING", "GO_EXTRACELLULAR_STRUCTURE_ORGANIZATION",
                      "GO_RESPONSE_TO_STEROID_HORMONE")
gs_mmu_data <- gs_mmu_data[which(as.character(lapply(gs_mmu_data, `[[`, 1))%in%pathway_selected)]

for (i in seq(length(gs_mmu_data))){
  pathway=gs_mmu_data[[i]][[3]]
  ensemnl_list <- getBM(attributes=c("mmusculus_homolog_associated_gene_name", "external_gene_name"), filters="external_gene_name", values = pathway, mart=ensembl)
  for (j in seq(nrow(ensemnl_list))){
    if (ensemnl_list[j,1] == ""){
      ensemnl_list[j,1] = paste0(substr(ensemnl_list[j,2],1,1), substr(tolower(ensemnl_list[j,2]), 2, nchar(ensemnl_list[j,2])))
    }
  }
  geneSets <- GeneSet(unique(ensemnl_list$mmusculus_homolog_associated_gene_name), setName=gs_mmu_data[[i]][[1]])
  geneSets_final <- c(geneSets, geneSets_final)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}


estrus <- loadRData("/omics/odcf/analysis/OE0538_projects/DO-0002/mmus/estrus/analysis_II/seurat_objects/ut_combined.Rdata")
decidua <- loadRData("/omics/groups/OE0433/internal/ivana/decidua/ut_combined_deci_merged_final.Rdata")
decidua$estrus_phase <- "decidua"
ut_list <- c(SplitObject(estrus, split.by = "orig.ident"), SplitObject(decidua, split.by = "orig.ident"))
ut_list <- lapply(X = ut_list, FUN = SCTransform, method = "glmGamPoi")
ut_combined <- merge(ut_list[[1]], y = unlist(ut_list[2:20], use.names=FALSE), project = "mouse_ut", merge.data=TRUE)
Idents(ut_combined) <- ut_combined$level_2
ut_combined_f <- subset(ut_combined, idents=c("F", "DeC"))

counts_m <- as.matrix(ut_combined_f[["SCT"]]@counts)
cells_rankings_m <- AUCell_buildRankings(counts_m, nCores=1, plotStats=TRUE)

null_rank_final_m <- c()

for (i in seq(from=1, to=ncol(counts_m))){
  rank <- as.data.frame(counts_m[,i])
  #print(i)
  colnames(rank) <- "rank"
  rank <- rank[order(rank$rank, decreasing = T),]
  null_rank <- min(which(rank==0))
  #print(null_rank)
  null_rank_final_m <- c(null_rank_final_m, null_rank)
}

auc_df_final_m <- data.frame()
for(i in geneSets_final){
  print(i)
  cells_AUC <- AUCell_calcAUC(i, cells_rankings_m, aucMaxRank = ceiling(0.045 * nrow(cells_rankings_m)))
  auc_df <- as.data.frame(getAUC(cells_AUC))
  #auc_df <- t(auc_df)
  auc_df_final_m <- rbind(auc_df, auc_df_final_m)
}

auc_df_final_m <- as.data.frame(t(auc_df_final_m))

metadata_m <- data.frame("cell_names" = colnames(ut_combined_f), "phase"= ut_combined_f$estrus_phase, "cell_type"=ut_combined_f$level_2, "batch"=ut_combined_f$orig.ident)
auc_df_m <- merge(auc_df_final_m, metadata_m, by=0)

auc_df_m_ecm <- data.frame(filler=1:5)

for (i in c("DeC", "F")){
  auc_df_ct <- auc_df_m[auc_df_m$cell_type==i,]
  auc_df_ct_ecm <- auc_df_ct %>%
    group_by(phase) %>%
    dplyr::summarize(Mean = mean(ECM_ORGANISATION, na.rm=TRUE))
  colnames(auc_df_ct_ecm)[2] <- i
  auc_df_m_ecm <- cbind(as.data.frame(auc_df_ct_ecm[,2]), auc_df_m_ecm)
  
  auc_df_ct_in <- auc_df_ct %>%
    group_by(phase) %>%
    dplyr::summarize(Mean = mean(INFLAMMATION, na.rm=TRUE))
  colnames(auc_df_ct_in)[2] <- i
  auc_df_m_in <- cbind(as.data.frame(auc_df_ct_in[,2]), auc_df_m_in)
  
}
