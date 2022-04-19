#################################
# this script uses raw filtered matrices deposited in Array Arrayexpress under accession number E-MTAB-11491, 
#marker list provided in Supplementary file 1 and seurat object generated in script estrus_sp_processing.R
#################################

library(ggrepel)

source("qc_normalization_cluster_annotation.R")



lotR_alt("matrix/spleen/")
t_combined_old_sp <- merge(`0002_18_mo_mm2_sp_filtered`, y = c(`0002_18_mo_mm3_sp_filtered`,`0002_18_mo_mm5_sp_filtered`), project = "mouse_sp", merge.data=TRUE)

t_combined_old_sp <- FindVariableFeatures(t_combined_old_sp, selection.method = "vst", nfeatures = 2000)
t_combined_old <- reductoR(t_combined_old_sp, 1)
save(t_combined_old,file="sp_combined_18m.Rdata")

t_lv_1 <- read.csv("sp_lv1_identity.csv")
t_lv_2 <- read.csv("sp_lv2_identity.csv")
t_combined <- loadRData("sp_combined_18m.Rdata")
seurat_clusters <- t_combined$seurat_clusters

t_combined$labels <- t_lv_1$identity_0.9
svm_cluster_lv1 <- garnett_extender(t_combined)

svm_annotation_ex_lv_1 <- data.frame("svm" = t_lv_1$identity_0.9, "seurat_clusters" = seurat_clusters)

t_combined$labels <- t_lv_2$identity_0.9
svm_cluster_lv2 <- garnett_extender(t_combined)

svm_annotation_ex_lv_2 <- data.frame("svm" = t_lv_2$identity_0.9, "seurat_clusters" = seurat_clusters)

markers <- list.load("spleen_markers.rdata")
markers_df <- data.frame(matrix(unlist(markers), nrow=length(markers), byrow=T),stringsAsFactors=FALSE)
markers_df <- as.data.frame(t(markers_df))
rownames(markers_df) <- NULL
colnames(markers_df) <- "markers"

list[t_combined_clusters, t_subset_clusters]=lotR_clusters(t_combined,markers_df)
print(DimPlot(t_combined, reduction = "umap", group.by = "seurat_clusters", label = T, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

svm_annotation_ex_lv_1 <- svm_annotation_ex_lv_1 %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 0  | seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 3 | 
    seurat_clusters == 4 | seurat_clusters == 5 | seurat_clusters == 6 | seurat_clusters == 7 |seurat_clusters == 8 | 
    seurat_clusters == 9  | seurat_clusters == 10 | seurat_clusters == 11 | seurat_clusters == 13 | seurat_clusters == 12 | 
    seurat_clusters == 14 | seurat_clusters == 15| seurat_clusters == 17 | seurat_clusters == 19  | seurat_clusters == 20 | 
    seurat_clusters == 21 | seurat_clusters == 22 | seurat_clusters == 24  | seurat_clusters == 25 | seurat_clusters == 26  ~ "IC",  
  seurat_clusters == 16 | seurat_clusters == 18 | seurat_clusters == 23 ~ "Er",
  
))

#26 out of 27 identified by SVM 
svm_annotation_ex_lv_2 <- svm_annotation_ex_lv_2 %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 0  | seurat_clusters == 3 | seurat_clusters == 4 | 
    seurat_clusters == 5 | seurat_clusters == 6 | seurat_clusters == 7 | seurat_clusters == 8 |
    seurat_clusters == 10 | seurat_clusters == 12  | seurat_clusters == 14 | seurat_clusters == 21 | seurat_clusters == 25~ "B-2",
  seurat_clusters == 1 | seurat_clusters == 13~ "iNKT",
  seurat_clusters == 24  | seurat_clusters == 15  ~ "MTC",
  seurat_clusters == 17  ~ "M1Mp",
  seurat_clusters == 9  | seurat_clusters == 22 | seurat_clusters == 26  ~ "DC",
  seurat_clusters == 11  ~ "NKC",
  seurat_clusters == 19 ~ "PC",
  seurat_clusters == 20 ~ "APC I",
  seurat_clusters == 2 ~ "CD4",
  seurat_clusters == 18  | seurat_clusters == 16 | seurat_clusters == 23~ "Er",
))


t_combined$level_1 <- svm_annotation_ex_lv_1$cell_labels_smv_ex
t_combined$level_2 <- svm_annotation_ex_lv_2$cell_labels_smv_ex


umap_plot <- as.data.frame(t_combined[["umap"]]@cell.embeddings)
umap_plot$level_1 <- unname(t_combined$level_1)
umap_plot$level_2 <- unname(t_combined$level_2)
label.cent_level_2 = umap_plot %>% group_by(level_2) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
label.cent_level_1 = umap_plot %>% group_by(level_1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

#metadata file is used only for color sheme (each cell types is associated with particular color)
metadata <- read.csv("metadata_cell_types.csv")

color_vector_level1 <- metadata$Cell_type_color_level_1
names(color_vector_level1) <- metadata$Cell_type_abbreviation_level_1

color_vector_level2 <- metadata$Cell_type_color
names(color_vector_level2) <- metadata$Cell_type_abbreviation

ggplot(umap_plot) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill=level_1), shape=21, color = 'black', size = 3)+
  scale_fill_manual(values=color_vector_level1)+
  theme_classic()+ NoAxes()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+guides(fill = guide_legend(override.aes = list(size = 3), title=""))

ggplot(umap_plot) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill=level_2), shape=21, color = 'black', size = 1.5)+
  scale_fill_manual(values=color_vector_level2)+
  theme_classic()+ NoAxes()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+guides(fill = guide_legend(override.aes = list(size = 3), title="", ncol=2))#+
#geom_label_repel(aes(x = UMAP_1, y = UMAP_2, label=level_2), size=5, force=10, data=label.cent_level_2)

t_list <- count_produceR(t_combined)
t_combined_identity <- t_list[[2]]
t_combined_identity <- t_combined_identity%>% mutate(condition = case_when(
  batch == "0002_18_mo_mm2_sp" | batch == "0002_18_mo_mm5_sp" | batch == "0002_18_mo_mm3_sp" ~ "18_m"))


save(t_combined,file="sp_combined_18m.Rdata")

write.table(t_combined_identity,"sp_combined_identity",sep=",",quote=F,row.names = F)