#!/usr/bin/env Rscript
#################################
# this script uses seurat objects generated in script estrus_ov_processing.R, estrus_od_processing.R, estrus_ut_processing.R, 
#estrus_ce_processing.R, estrus_va_processing.R, estrus_sp_processing.R
#################################

args = commandArgs(trailingOnly=TRUE)


if (length(args)<3) {
  stop("Three arguments must be supplied (directory containing seurat objects, spleen or not_spleen, output directory for UMAPS).n", call.=FALSE)}

library(Seurat)
library(dplyr)
library(ggplot2)
library("ggrepel")
ulimit::memory_limit(380000)

options(future.globals.maxSize = 80000 * 1024^2)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

integratoR <- function(rep_list, ref_anchor) {
  rep_list <- lapply(X = rep_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
  })
  
  # select features that are repeatedly variable across datasets for integration run PCA on each
  # dataset using these features
  features <- SelectIntegrationFeatures(object.list = rep_list, nfeatures = 5000)
  rep_list <- lapply(X = rep_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  t_anchors <- FindIntegrationAnchors(object.list = rep_list, anchor.features = features, reduction = "rpca", reference = ref_anchor)
  # this command creates an 'integrated' data assay
  all.combined <- IntegrateData(anchorset = t_anchors)
  return(all.combined)
}


data_dir <- args[1]
setwd(data_dir)


rep_list <-list()
max_cells <- c()
ref <- c()

file_names <- list.files(data_dir)

if (args[2]=="spleen"){
  file_names <- file_names
} else {file_names <- file_names[!file_names=="sp_combined_new.Rdata"]}

for (i in seq(from=1, to=length(file_names))){
  t_combined <- loadRData(file_names[i])
  DefaultAssay(t_combined) <- "RNA"
  #t_combined <- RenameCells(t_combined, add.cell.id= substr(file_names[1],1,2))
  t_combined[["SCT"]] <- NULL
  t_combined[["umap"]]<- NULL
  t_list <- SplitObject(t_combined, split.by = "orig.ident")
  for(j in seq(from=1, to=length(t_list))){
    max_cells[j] <- length(t_list[[j]]$orig.ident)
  }
  ref <- c(ref, names(t_list[which.max(max_cells)]))
  rep_list <- c(rep_list, t_list)
  
}  
rm (t_combined, t_list)

ref_anchor <- c()
ref_anchor_temp <- c()
for (i in seq(from=1, to=length(ref))){
  ref_anchor_temp <- which(names(rep_list) == ref[i])
  ref_anchor <- c(ref_anchor_temp, ref_anchor)
}

rep_combined_sct <- integratoR(rep_list, ref_anchor)
DefaultAssay(rep_combined_sct) <- "integrated"

# Run the standard workflow for visualization and clustering
rep_combined_sct <- FindVariableFeatures(rep_combined_sct, selection.method = "vst", nfeatures = 5000)
rep_combined_sct  <- ScaleData(rep_combined_sct , verbose = FALSE)
rep_combined_sct <- RunPCA(rep_combined_sct, verbose = FALSE)
rep_combined_sct <- RunUMAP(rep_combined_sct, reduction = "pca", dims = 1:30, metric ='euclidean',min.dist = 0.2)

metadata <- read.csv("~/scRNA-seq/Scripts/metadata_cell_types.csv")
color_vector_level2 <- metadata$Cell_type_color
names(color_vector_level2) <- metadata$Cell_type_abbreviation

umap_plot <- as.data.frame(rep_combined_sct[["umap"]]@cell.embeddings)
#umap_plot$level_1 <- unname(rep_combined_sct$level_1)
umap_plot$level_2 <- unname(rep_combined_sct$level_2)
label.cent_level_2 = umap_plot %>% group_by(level_2) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
#label.cent_level_1 = umap_plot %>% group_by(level_1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

#umap_plot <- umap_plot[order((umap_plot$level_2), decreasing = TRUE),]
ggplot(umap_plot) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill=level_2), shape=21, size = 1)+
  scale_fill_manual(values=color_vector_level2)+
  theme_classic()+ NoAxes()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+guides(fill = guide_legend(override.aes = list(size = 3), title=""))#+
#geom_label_repel(aes(x = UMAP_1, y = UMAP_2, label=level_2), size=5, force=10, data=label.cent_level_2)

label_df <- data.frame("cells"=colnames(rep_combined_sct), 
                       "organ" = substr(rep_combined_sct$orig.ident, nchar(rep_combined_sct$orig.ident)-1, nchar(rep_combined_sct$orig.ident)),
                       "cell_type" = rep_combined_sct$level_2)


organ_UMAPeR <- function(label_df, organ, umap_plot, level_order){
  color_label <- c()
  shape_label <- c()
  for (i in seq(from=1, to=nrow(label_df))){
    if (label_df[i,2]==organ){
      color_label[i] = label_df[i,3]
    } else {color_label[i] = "other"}
  }
  for (i in seq(from=1, to=nrow(label_df))){
    if (label_df[i,2]==organ){
      shape_label[i] = "shape_1"
    } else {shape_label[i] = "other"}
  }
  color_vector_young <- "grey85"
  names(color_vector_young) <- "other"
  color_vector_level2 <- c(color_vector_level2, color_vector_young)
  umap_plot$label <- color_label
  umap_plot$shape <- shape_label
  umap_plot$organ <- label_df$organ
  umap_plot <- umap_plot %>%
    arrange(factor(organ, levels=level_order))
  print(ggplot(umap_plot) + 
          geom_point(aes(x = UMAP_1, y = UMAP_2, color=label), shape=16, size=0.01)+
          scale_color_manual(breaks = c(levels(as.factor(umap_plot$label))[levels(as.factor(umap_plot$label)) != "other"]),values=color_vector_level2)+
          theme_classic()+ NoAxes()+
          theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
          theme(legend.text=element_text(size=20))+guides(color = guide_legend(override.aes = list(size = 3), title="")))
}

organ_UMAPeR(label_df, "ut", umap_plot, c("ov", "ce", "od","va", "ut"))
organ_UMAPeR(label_df, "ov", umap_plot, c("ce", "ut", "od","va", "ov"))
organ_UMAPeR(label_df, "od", umap_plot, c("ov", "ut", "ce","va", "od"))
organ_UMAPeR(label_df, "ce", umap_plot, c("ov", "ut", "od","va", "ce"))
organ_UMAPeR(label_df, "va", umap_plot,  c("ov", "ut", "od","ce", "va"))

umap_plot <- umap_plot %>%
  arrange(factor(organ, levels=c("ce", "ut", "od","va", "ov")))
print(ggplot(umap_plot) + 
        geom_point(aes(x = UMAP_1, y = UMAP_2, color=level_2), shape=16, size=0.01)+
        scale_color_manual(values=color_vector_level2)+
        theme_classic()+ NoAxes()+
        theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
        theme(legend.text=element_text(size=20))+guides(color = guide_legend(override.aes = list(size = 3), title="")))


setwd(args[3])
save(rep_combined_sct, file="estrus_integrated.Rdata")

