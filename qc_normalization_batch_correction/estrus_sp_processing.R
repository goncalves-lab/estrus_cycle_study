#################################
# this script uses raw filtered matrices deposited in Array Arrayexpress under accession number E-MTAB-11491 and 
#marker list provided in Supplementary file 1
#################################

library("ggrepel")


source("qc_normalization_cluster_annotation.R")


lotR("matrix/spleen")

t_combined <- merge(`0002_mouse1_sp_filtered`, y = c(`0002_mouse10_sp_filtered`,`0002_mouse11_sp_filtered`,`0002_mouse12_sp_filtered`,
                                                     `0002_mouse13_sp_filtered`,`0002_mouse14_sp_filtered`,`0002_mouse15_sp_filtered`,
                                                     `0002_mouse16_sp_filtered`,`0002_mouse2_sp_filtered`,`0002_mouse3_sp_filtered`,
                                                     `0002_mouse6_sp_filtered`,`0002_mouse8_sp_filtered`,`0002_mouse7_sp_filtered`,
                                                     `0002_mouse9_sp_filtered`), project = "mouse_sp", merge.data=TRUE)     

t_combined_cell_metadata <- data.frame("cell_names"=colnames(t_combined), "batch" = t_combined$orig.ident)
t_combined_cell_metadata <- t_combined_cell_metadata %>% mutate(condition = case_when(
  batch == "0002_mouse1_sp" | batch == "0002_mouse6_sp" | batch == "0002_mouse12_sp" | batch == "0002_mouse13_sp" ~ "proestrus",
  batch == "0002_mouse2_sp" | batch == "0002_mouse8_sp" | batch == "0002_mouse7_sp"  | batch == "0002_mouse16_sp"  ~ "diestrus", 
  batch == "0002_mouse3_sp" | batch == "0002_mouse9_sp" | batch == "0002_mouse14_sp" ~ "estrus",
  batch == "0002_mouse10_sp" | batch == "0002_mouse11_sp" | batch == "0002_mouse15_sp" ~ "metestrus"))


t_combined <- FindVariableFeatures(t_combined, selection.method = "vst", nfeatures = 2000)
t_combined <- reductoR(t_combined, 1.5)

print(DimPlot(t_combined, reduction = "umap", group.by = "seurat_clusters", label = T, pt.size = 0.7,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))


markers <- list.load("spleen_markers.rdata")
markers_df <- data.frame(matrix(unlist(markers), nrow=length(markers), byrow=T),stringsAsFactors=FALSE)
markers_df <- as.data.frame(t(markers_df))
rownames(markers_df) <- NULL
colnames(markers_df) <- "markers"

list[t_combined_clusters, t_subset_clusters]=lotR_clusters(t_combined,markers_df)
#writing out DF with information about cell identities and phase (condition) to be used with 
#cell_adundance script
write.table(t_combined_clusters,"sp_combined_clusters",sep=",",quote=F,row.names = F)
write.table(t_subset_clusters,"sp_combined_clusters_subset",sep=",",quote=F,row.names = F)

new_cluster_id_t_com=c("BC", "BC" , "APC", "APC", "BC", "TC", "TC", "BC", "APC","TC","APC", "TC", "TC", "?", "BC", "NKC", "TC", "Er", "APC", "Er", 
                       "TC", "APC", "APC","APC", "TC", "TC", "BC","BC", "TC", "BC", "APC", "APC", "APC", "APC", "Er")
t_combined <- lotR_ultimate(t_combined, new_cluster_id_t_com)

t <- t_combined$seurat_clusters
t <- data.frame("cells" = names(t), "type" = unname(t))

t <- t %>% mutate(cell_names = case_when(
  type == 13 ~ "?",
  type == 0 | type == 1 | type == 4 | type ==  7 | type == 14 | type == 26 | type == 27 | type ==  29 ~ "BC",
  type == 2 | type == 3 | type == 8 | type == 10 | type == 18 | type == 21 | type == 22 | type == 23 | type == 30 | type == 32 | type == 33 
  | type==31 ~ "APC",
  type == 5 | type == 6 | type == 9 | type == 11 | type == 12 | type == 16 | type == 20 | type == 24 | type == 25  | type == 28 ~ "TC",
  type == 15  ~ "NKC",
  type == 17 | type == 34 | type == 19 ~ "Er" 
))

t_label <- t$cell_names
names(t_label) <- t$cells
t_combined$top_level <- t_label


#garnett
rownames(t_combined_cell_metadata) <- t_combined_cell_metadata$cell_names
t_combined_cell_metadata$tissue <- "spleen"

t_combined_cds <- garnett_annotatoR_part1(t_combined, t_combined_cell_metadata, "spleen_markers_garnett")
t_combined <- garnett_annotatoR_part2(t_combined_cds, t_combined_cds, "spleen_markers_garnett", t_combined)

garnett_labels <- t_combined$labels
seurat_clusters <- t_combined$seurat_clusters

garnett_annotation_ex <- data.frame("garnett" = garnett_labels, "seurat_clusters" = seurat_clusters)

garnett_cluster <- garnett_extender(t_combined) 

garnett_annotation_ex <- garnett_annotation_ex %>% mutate(cell_labels_garnett_ex = case_when(
  seurat_clusters == 0  | seurat_clusters == 4 | seurat_clusters ==  7 | seurat_clusters == 14 | seurat_clusters == 26  
  | seurat_clusters ==  29  ~ "BC",
  seurat_clusters == 2 | seurat_clusters == 3 | seurat_clusters == 8 | seurat_clusters == 10 | seurat_clusters == 18 | seurat_clusters == 21 
  | seurat_clusters == 30 | seurat_clusters == 27 | seurat_clusters == 22~ "APC",
  seurat_clusters == 5 | seurat_clusters == 6 | seurat_clusters == 9 | seurat_clusters == 11 | seurat_clusters == 12 | seurat_clusters == 16 
  | seurat_clusters == 20 | seurat_clusters == 24 | seurat_clusters == 25  | seurat_clusters == 28 ~ "TC",
  seurat_clusters == 15  ~ "NKC",
  seurat_clusters == 17 | seurat_clusters == 34 | seurat_clusters == 19 | seurat_clusters == 1 | seurat_clusters == 13 | seurat_clusters == 23
  | seurat_clusters == 32 | seurat_clusters == 33 | seurat_clusters==31~ "Er" 
))

t_combined$garnett_cluster <- garnett_annotation_ex$cell_labels_garnett_ex
t_combined$labels <- NULL
print(DimPlot(t_combined, reduction = "umap", group.by = "garnett_cluster", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

garnett <- t_combined$garnett_cluster
manual <- t_combined$top_level
annotation_df <- data.frame("cells" = names(garnett), "garnett" = unname(garnett), "manual" = manual)
annotation_df$matched_cells = "?"

for (i in seq(from=1, to=nrow(annotation_df))){
  if (annotation_df$garnett[i] == annotation_df$manual[i]){
    annotation_df$matched_cells[i] = "match"
  } else (annotation_df$matched_cells[i] = "non-match")
}

matched_cells <- annotation_df$matched_cells
names(matched_cells) <- annotation_df$cells
t_combined$match <- matched_cells

list[t_non_match_ge, t_match_ge] <- svm_subseteR(t_combined)

pred.te_identity <- svm_returneR(t_match_ge, t_non_match_ge, 0.7, 0.2)
Idents(t_combined) <- t_combined$match
t_match <- subset(t_combined, idents = "match") 
t_non_match <- subset(t_combined, idents = "non-match") 

t_non_match$labels <- pred.te_identity
svm_cluster <- garnett_extender(t_non_match)
seurat_clusters <- t_non_match$seurat_clusters

svm_annotation_ex <- data.frame("svm" = pred.te_identity, "seurat_clusters" = seurat_clusters)
svm_annotation_ex <- svm_annotation_ex %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 1  | seurat_clusters == 13 | seurat_clusters == 27 ~ "BC" ,
  seurat_clusters == 23 | seurat_clusters == 31 | seurat_clusters == 32 | seurat_clusters ==  33 ~ "APC"
))


t_non_match$svm_cluster <- svm_annotation_ex$cell_labels_smv_ex

assigned.te <- as.data.frame(t_match$top_level) 
predicted.te <- as.data.frame(t_non_match$svm_cluster)
colnames(assigned.te) <- "identity"
colnames(predicted.te) <- "identity"

final_labels <- rbind (assigned.te, predicted.te)
final_labels$cells <- rownames(final_labels)

manual_labels <- as.data.frame(t_combined$top_level) 
manual_labels$cells <- rownames(manual_labels)
final_labels <- final_labels[manual_labels[["cells"]],]
final_labels_vector <- final_labels$identity
names(final_labels_vector) <- final_labels$cells
t_combined$final_labels <- final_labels_vector

print(DimPlot(t_combined, reduction = "umap", group.by = "final_labels", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

save(t_combined,file="sp_combined_new.Rdata")


########################
#TC
#######################

Idents(t_combined) <- t_combined$final_labels
list[t_tc, tc_cluster, tc_subcluster] <- subseteR_part_1(t_combined, "TC", 1)
new_cluster_id_tc=c("MTC", "MTC", "?","CD8", "CD4","iNKT", "CD8", "MTC", "?", "?", "?", "?", "MAIT", "iNKT", "iNKT", "MTC", "MTC")
tc <- subseteR_part_2(t_tc, new_cluster_id_tc)

tc <- tc %>% mutate(cell_names = case_when(
  type == 12 | type ==  15 ~ "MAIT",
  type == 5 | type == 14 | type == 13~ "iNKT",
  type == 1| type == 0 | type == 7 | type == 16    ~ "MTC",
  type ==  2 | type == 8 | type == 9 | type ==  10 | type == 11 ~ "?",
  type == 3 | type == 6~ "CD8",
  type == 4 ~ "CD4"))


t_tc$sublevel <- tc$cell_names

print(DimPlot(t_tc, reduction = "umap", label = T ,pt.size = 0.5,group.by="sublevel",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

t_tc_cell_metadata <- data.frame("cell_names"=colnames(t_tc), "batch" = t_tc$orig.ident)
t_tc_cell_metadata$tissue <- "spleen"
t_combined_cds_ep <- garnett_annotatoR_part1(t_tc, t_tc_cell_metadata, "sp_tc")
t_tc <- garnett_annotatoR_part2(t_combined_cds_ep,t_combined_cds_ep, "sp_tc", t_tc)


garnett_labels_tc <- t_tc$labels
seurat_clusters_tc <- t_tc$seurat_clusters

garnett_annotation_ex_tc <- data.frame("garnett" = garnett_labels_tc, "seurat_clusters" = seurat_clusters_tc)

garnett_cluster_tc <- garnett_extender(t_tc) 

garnett_annotation_ex_tc <- garnett_annotation_ex_tc %>% mutate(cell_labels_garnett_ex = case_when(
  seurat_clusters ==  0 | seurat_clusters == 4  | seurat_clusters == 9 | seurat_clusters == 10 ~ "CD4",
  seurat_clusters == 5 ~ "iNKT",
  seurat_clusters == 15 ~ "MAIT" ,
  seurat_clusters == 6 | seurat_clusters == 11 ~ "CD8",
  seurat_clusters == 3 | seurat_clusters == 8 | seurat_clusters == 12   ~ "Unknown",
  seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 7 | seurat_clusters == 13 | seurat_clusters == 14 | seurat_clusters == 16 ~ "MTC"
))

t_tc$garnett_cluster <- garnett_annotation_ex_tc$cell_labels_garnett_ex
t_tc$labels <- NULL

print(DimPlot(t_tc, reduction = "umap", label = T ,pt.size = 0.5,group.by="garnett_cluster",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

garnett <- t_tc$garnett_cluster
manual <- t_tc$sublevel
annotation_df_tc <- data.frame("cells" = names(garnett), "garnett" = unname(garnett), "manual" = manual)
annotation_df_tc$matched_cells = "?"

for (i in seq(from=1, to=nrow(annotation_df_tc))){
  if (annotation_df_tc$garnett[i] == annotation_df_tc$manual[i]){
    annotation_df_tc$matched_cells[i] = "match"
  } else (annotation_df_tc$matched_cells[i] = "non-match")
}

matched_cells_tc <- annotation_df_tc$matched_cells
names(matched_cells_tc) <- annotation_df_tc$cells
t_tc$match <- matched_cells_tc
print(DimPlot(t_tc, reduction = "umap", group.by = "match", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

list[t_non_match_ge_tc, t_match_ge_tc] <- svm_subseteR(t_tc)

pred.te_identity_tc <- svm_returneR(t_match_ge_tc, t_non_match_ge_tc, 0.7, 0.9)
Idents(t_tc) <- t_tc$match
t_match_tc <- subset(t_tc, idents = "match") 
t_non_match_tc <- subset(t_tc, idents = "non-match") 

t_non_match_tc$labels <- pred.te_identity_tc
svm_cluster_tc <- garnett_extender(t_non_match_tc)
seurat_clusters_tc <- t_non_match_tc$seurat_clusters

svm_annotation_ex_tc <- data.frame("svm" = pred.te_identity_tc, "seurat_clusters" = seurat_clusters_tc)

svm_annotation_ex_tc <- svm_annotation_ex_tc %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 0 | seurat_clusters == 8 | seurat_clusters == 9 | seurat_clusters == 10 | seurat_clusters == 11  ~ "MTC",
  seurat_clusters == 3 | seurat_clusters == 12 | seurat_clusters == 14 ~ "iNKT",
  seurat_clusters == 13 ~ "CD8",
  seurat_clusters == 2 ~ "TC I"
))

t_non_match_tc$svm_cluster <- svm_annotation_ex_tc$cell_labels_smv_ex

assigned.te_tc <- as.data.frame(t_match_tc$sublevel) 
predicted.te_tc <- as.data.frame(t_non_match_tc$svm_cluster)
colnames(assigned.te_tc) <- "identity"
colnames(predicted.te_tc) <- "identity"

final_labels_tc <- rbind (assigned.te_tc, predicted.te_tc)
final_labels_tc$cells <- rownames(final_labels_tc)

manual_labels_tc <- as.data.frame(t_tc$sublevel) 
manual_labels_tc$cells <- rownames(manual_labels_tc)
final_labels_tc <- final_labels_tc[manual_labels_tc[["cells"]],]
final_labels_tc_vector <- final_labels_tc$identity
names(final_labels_tc_vector) <- final_labels_tc$cells
t_tc$final_labels <- final_labels_tc_vector

print(DimPlot(t_tc, reduction = "umap", group.by = "final_labels", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

save(t_tc,file="sp_tc_final.Rdata")


########################
#APC
#######################
Idents(t_combined) <- t_combined$final_labels
list[t_apc, apc_cluster, apc_subcluster] <- subseteR_part_1(t_combined, "APC", 1)
new_cluster_id_apc=c("BC", "BC", "BC", "BC", "BC", "?", "MC", "?", "?", "Mp","M1Mp", "DC", "DC", "DC", "DC", "DC", "Mp", "DC", 
                     "?","M2Mp","M2Mp", "BC", "MTC", "EC")
apc <- subseteR_part_2(t_apc, new_cluster_id_apc)

apc <- apc %>% mutate(cell_names = case_when(
  type == 11 | type == 12 | type == 13 | type== 14 | type== 15 | type== 17 ~ "DC",
  type == 10 ~ "M1Mp",
  type == 0 | type== 1 | type== 2| type== 3 | type== 4 | type== 21~ "BC",
  type == 7 | type == 8 | type== 18 | type== 5 ~ "?",
  type== 19 | type== 20 ~ "M2Mp",
  type== 9 | type== 16 ~ "Mp",
  type== 6  ~ "MC",
  type== 22 ~ "MTC",
  type== 23 ~ "EC"
  
))

t_apc$sublevel <- apc$cell_names

print(DimPlot(t_apc, reduction = "umap", label = T ,pt.size = 0.5,group.by="sublevel",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

t_apc_cell_metadata <- data.frame("cell_names"=colnames(t_apc), "batch" = t_apc$orig.ident)
t_apc_cell_metadata$tissue <- "spleen"
t_combined_cds_ep <- garnett_annotatoR_part1(t_apc, t_apc_cell_metadata, "sp_markers_apc")
t_apc <- garnett_annotatoR_part2(t_combined_cds_ep,t_combined_cds_ep, "sp_markers_apc", t_apc)


garnett_labels_apc <- t_apc$labels
seurat_clusters_apc <- t_apc$seurat_clusters

garnett_annotation_ex_apc <- data.frame("garnett" = garnett_labels_apc, "seurat_clusters" = seurat_clusters_apc)

garnett_cluster_apc <- garnett_extender(t_apc) 

garnett_annotation_ex_apc <- garnett_annotation_ex_apc %>% mutate(cell_labels_garnett_ex = case_when(
  seurat_clusters == 0 | seurat_clusters == 5 | seurat_clusters == 13 | seurat_clusters== 14 | seurat_clusters== 15 
  | seurat_clusters== 17 | seurat_clusters== 18 | seurat_clusters== 2 ~ "DC",
  seurat_clusters == 10 | seurat_clusters==11 ~ "M1Mp",
  seurat_clusters== 1 | seurat_clusters== 4 | seurat_clusters== 21~ "BC",
  seurat_clusters == 9 | seurat_clusters == 8 | seurat_clusters== 12  ~ "Unknown",
  seurat_clusters== 3 | seurat_clusters== 20 | seurat_clusters== 7 | seurat_clusters== 19~ "M2Mp",
  seurat_clusters== 16 ~ "Mp",
  seurat_clusters== 6  ~ "MC",
  seurat_clusters== 22 ~ "MTC",
  seurat_clusters== 23 ~ "EC"
))

t_apc$garnett_cluster <- garnett_annotation_ex_apc$cell_labels_garnett_ex
t_apc$labels <- NULL

print(DimPlot(t_apc, reduction = "umap", label = T ,pt.size = 0.5,group.by="garnett_cluster",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

garnett <- t_apc$garnett_cluster
manual <- t_apc$sublevel
annotation_df_apc <- data.frame("cells" = names(garnett), "garnett" = unname(garnett), "manual" = manual)
annotation_df_apc$matched_cells = "?"

for (i in seq(from=1, to=nrow(annotation_df_apc))){
  if (annotation_df_apc$garnett[i] == annotation_df_apc$manual[i]){
    annotation_df_apc$matched_cells[i] = "match"
  } else (annotation_df_apc$matched_cells[i] = "non-match")
}

matched_cells_apc <- annotation_df_apc$matched_cells
names(matched_cells_apc) <- annotation_df_apc$cells
t_apc$match <- matched_cells_apc
print(DimPlot(t_apc, reduction = "umap", group.by = "match", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

list[t_non_match_ge_apc, t_match_ge_apc] <- svm_subseteR(t_apc)

pred.te_identity_apc <- svm_returneR(t_match_ge_apc, t_non_match_ge_apc, 0.7, 0.9)
Idents(t_apc) <- t_apc$match
t_match_apc <- subset(t_apc, idents = "match") 
t_non_match_apc <- subset(t_apc, idents = "non-match") 

t_non_match_apc$labels <- pred.te_identity_apc
svm_cluster_apc <- garnett_extender(t_non_match_apc)
seurat_clusters_apc <- t_non_match_apc$seurat_clusters

svm_annotation_ex_apc <- data.frame("svm" = pred.te_identity_apc, "seurat_clusters" = seurat_clusters_apc)

svm_annotation_ex_apc <- svm_annotation_ex_apc %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 0 | seurat_clusters== 2 | seurat_clusters== 3 | seurat_clusters== 5 | seurat_clusters== 7
  | seurat_clusters== 8~ "BC",
  seurat_clusters == 9 | seurat_clusters== 11 | seurat_clusters== 12  ~ "DC",
  seurat_clusters== 18 ~ "APC I"
  
))

t_non_match_apc$svm_cluster <- svm_annotation_ex_apc$cell_labels_smv_ex

assigned.te_apc <- as.data.frame(t_match_apc$sublevel) 
predicted.te_apc <- as.data.frame(t_non_match_apc$svm_cluster)
colnames(assigned.te_apc) <- "identity"
colnames(predicted.te_apc) <- "identity"

final_labels_apc <- rbind (assigned.te_apc, predicted.te_apc)
final_labels_apc$cells <- rownames(final_labels_apc)

manual_labels_apc <- as.data.frame(t_apc$sublevel) 
manual_labels_apc$cells <- rownames(manual_labels_apc)
final_labels_apc <- final_labels_apc[manual_labels_apc[["cells"]],]
final_labels_apc_vector <- final_labels_apc$identity
names(final_labels_apc_vector) <- final_labels_apc$cells
t_apc$final_labels <- final_labels_apc_vector

print(DimPlot(t_apc, reduction = "umap", group.by = "final_labels", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

save(t_apc,file="sp_apc_final.Rdata")


########################
#Level 2
#######################

Idents(t_combined) <- t_combined$final_labels
list[t_nkc, nkc_cluster, nkc_subcluster] <- subseteR_part_1(t_combined, "NKC",1)
final_labels_nkc <- data.frame("identity"="NKC", "cells"=names(t_nkc$seurat_clusters))
rownames(final_labels_nkc) <- final_labels_nkc$cells

list[t_bc, bc_cluster, bc_subcluster] <- subseteR_part_1(t_combined, "BC",1)
final_labels_bc <- data.frame("identity"="BC", "cells"=names(t_bc$seurat_clusters))
rownames(final_labels_bc) <- final_labels_bc$cells

list[t_er, er_cluster, er_subcluster] <- subseteR_part_1(t_combined, "Er",1)
final_labels_er <- data.frame("identity"="Er", "cells"=names(t_er$seurat_clusters))
rownames(final_labels_er) <- final_labels_er$cells

final_labels_level_2 <- rbind (final_labels_apc, final_labels_nkc, final_labels_bc, final_labels_tc, final_labels_er)

manual_labels <- as.data.frame(t_combined$final_labels) 
manual_labels$cells <- rownames(manual_labels)
final_labels_level_2 <- final_labels_level_2[manual_labels[["cells"]],]
final_labels_level_2_vector <- final_labels_level_2$identity
names(final_labels_level_2_vector) <- final_labels_level_2$cells
t_combined$level_2 <- final_labels_level_2_vector

level_1_label <- data.frame("label"=t_combined$level_2)

level_1_label <- level_1_label %>% mutate(level_1_final = case_when(
  label == "MTC" | label == "TC I"  | label == "CD8" | label == "CD4" | label == "iNKT" | label == "MAIT" 
  | label == "NKC" ~ "TC",
  label == "BC" ~ "BC",
  label == "MC" ~ "MC",
  label == "Er" ~ "Er",
  label == "EC" ~ "EC",
  label == "APC I" | label == "Mp" | label == "M2Mp" | label == "DC" | label == "M1Mp" ~ "APC"
))

t_combined$final_labels <- level_1_label$level_1_final
########################
#BC
#######################

Idents(t_combined) <- t_combined$level_2
list[t_bc, bc_cluster, bc_subcluster] <- subseteR_part_1(t_combined, "BC", 1)
new_cluster_id_bc=c("?", "?", "?", "?", "B-2", "B-2", "?", "?", "?", "?","?", "B-2","B-2", "?","B-2", "B-2", "MBC", "?","PC", "?")
bc <- subseteR_part_2(t_bc, new_cluster_id_bc)

bc <- bc %>% mutate(cell_names = case_when(
  type == 11 | type== 12 | type==14 | type==15 |type == 4 | type== 5 ~ "B-2",
  type== 18 ~ "PC",
  type== 16 ~ "MBC",
  type== 0 | type== 1 | type == 2 | type== 3 | type== 7 | type== 6 | type == 8 | type== 9 | type== 10 | 
    type == 17 | type== 13 | type==19 ~ "?"
  
))

t_bc$sublevel <- bc$cell_names

print(DimPlot(t_bc, reduction = "umap", label = T ,pt.size = 0.5,group.by="sublevel",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

t_bc_cell_metadata <- data.frame("cell_names"=colnames(t_bc), "batch" = t_bc$orig.ident)
t_bc_cell_metadata$tissue <- "spleen"
t_combined_cds_bc <- garnett_annotatoR_part1(t_bc, t_bc_cell_metadata, "sp_markers_bc")
t_bc <- garnett_annotatoR_part2(t_combined_cds_bc,t_combined_cds_bc, "sp_markers_bc", t_bc)


garnett_labels_bc <- t_bc$labels
seurat_clusters_bc <- t_bc$seurat_clusters

garnett_annotation_ex_bc <- data.frame("garnett" = garnett_labels_bc, "seurat_clusters" = seurat_clusters_bc)

garnett_cluster_bc <- garnett_extender(t_bc) 

garnett_annotation_ex_bc <- garnett_annotation_ex_bc %>% mutate(cell_labels_garnett_ex = case_when(
  seurat_clusters == 2 | seurat_clusters == 4 | seurat_clusters == 7 | seurat_clusters ==9
  | seurat_clusters ==12 | seurat_clusters ==14 | seurat_clusters ==15 | seurat_clusters ==19 ~ "B-2",
  seurat_clusters== 18 | seurat_clusters== 8 | seurat_clusters== 3 ~ "PC",
  seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters== 5 | seurat_clusters== 6 | seurat_clusters== 10
  | seurat_clusters== 11 | seurat_clusters==13 | seurat_clusters==17~ "Unknown",
  seurat_clusters== 16  ~ "MBC"
))

t_bc$garnett_cluster <- garnett_annotation_ex_bc$cell_labels_garnett_ex
t_bc$labels <- NULL

print(DimPlot(t_bc, reduction = "umap", label = T ,pt.size = 0.5,group.by="garnett_cluster",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

garnett <- t_bc$garnett_cluster
manual <- t_bc$sublevel
annotation_df_bc <- data.frame("cells" = names(garnett), "garnett" = unname(garnett), "manual" = manual)
annotation_df_bc$matched_cells = "?"

for (i in seq(from=1, to=nrow(annotation_df_bc))){
  if (annotation_df_bc$garnett[i] == annotation_df_bc$manual[i]){
    annotation_df_bc$matched_cells[i] = "match"
  } else (annotation_df_bc$matched_cells[i] = "non-match")
}

matched_cells_bc <- annotation_df_bc$matched_cells
names(matched_cells_bc) <- annotation_df_bc$cells
t_bc$match <- matched_cells_bc
print(DimPlot(t_bc, reduction = "umap", group.by = "match", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

list[t_non_match_ge_bc, t_match_ge_bc] <- svm_subseteR(t_bc)

pred.te_identity_bc <- svm_returneR(t_match_ge_bc, t_non_match_ge_bc, 0.7, 0.9)
Idents(t_bc) <- t_bc$match
t_match_bc <- subset(t_bc, idents = "match") 
t_non_match_bc <- subset(t_bc, idents = "non-match") 

t_non_match_bc$labels <- pred.te_identity_bc
svm_cluster_bc <- garnett_extender(t_non_match_bc)
seurat_clusters_bc <- t_non_match_bc$seurat_clusters

svm_annotation_ex_bc <- data.frame("svm" = pred.te_identity_bc, "seurat_clusters" = seurat_clusters_bc)

svm_annotation_ex_bc <- svm_annotation_ex_bc %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 0 | seurat_clusters== 2 | seurat_clusters== 3 | seurat_clusters== 1 | seurat_clusters== 6
  | seurat_clusters== 5  | seurat_clusters== 9  | seurat_clusters== 7  | seurat_clusters== 8  | seurat_clusters==10
  | seurat_clusters== 11|  seurat_clusters == 13 |  seurat_clusters == 17 |  seurat_clusters == 19  ~ "B-2"
  
))

t_non_match_bc$svm_cluster <- svm_annotation_ex_bc$cell_labels_smv_ex

assigned.te_bc <- as.data.frame(t_match_bc$sublevel) 
predicted.te_bc <- as.data.frame(t_non_match_bc$svm_cluster)
colnames(assigned.te_bc) <- "identity"
colnames(predicted.te_bc) <- "identity"

final_labels_bc <- rbind (assigned.te_bc, predicted.te_bc)
final_labels_bc$cells <- rownames(final_labels_bc)

manual_labels_bc <- as.data.frame(t_bc$sublevel) 
manual_labels_bc$cells <- rownames(manual_labels_bc)
final_labels_bc <- final_labels_bc[manual_labels_bc[["cells"]],]
final_labels_bc_vector <- final_labels_bc$identity
names(final_labels_bc_vector) <- final_labels_bc$cells
t_bc$final_labels <- final_labels_bc_vector

print(DimPlot(t_bc, reduction = "umap", group.by = "final_labels", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

save(t_bc,file="sp_bc_final.Rdata")

Idents(t_combined) <- t_combined$level_2
t_combined_no_bc <- subset(x = t_combined, idents = "BC", invert=T)
final_labels_no_bc <- data.frame("identity"=t_combined_no_bc$level_2, "cells"=names(t_combined_no_bc$seurat_clusters))
rownames(final_labels_no_bc) <- final_labels_no_bc$cells

final_labels_level_2 <- rbind (final_labels_no_bc, final_labels_bc)

manual_labels <- as.data.frame(t_combined$level_2) 
manual_labels$cells <- rownames(manual_labels)
final_labels_level_2 <- final_labels_level_2[manual_labels[["cells"]],]
final_labels_level_2_vector <- final_labels_level_2$identity
names(final_labels_level_2_vector) <- final_labels_level_2$cells
t_combined$level_2 <- final_labels_level_2_vector




#final plots
level_1 <- data.frame("cells"=colnames(t_combined), "old_labels"=t_combined$final_labels)

level_1 <- level_1 %>% mutate(new_labels = case_when(
  old_labels == "APC" | old_labels == "TC" | old_labels == "MC" | old_labels == "NKC" | old_labels == "BC" ~ "IC",
  old_labels == "EC" ~ "EC",
  old_labels == "EpC" ~ "EpC",
  old_labels == "SC" ~ "SC",
  old_labels == "Er" ~ "Er"))
t_combined$level_1 <- level_1$new_labels

t_combined$level_2 <- str_replace(t_combined$level_2, "MC", "MaC")
umap_plot <- as.data.frame(t_combined[["umap"]]@cell.embeddings)
umap_plot$level_1 <- unname(t_combined$level_1)
umap_plot$level_2 <- unname(t_combined$level_2)
label.cent_level_2 = umap_plot %>% group_by(level_2) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
label.cent_level_1 = umap_plot %>% group_by(level_1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

metadata <- read.csv("metadata_cell_types.csv")

color_vector_level1 <- metadata$Cell_type_color_level_1
names(color_vector_level1) <- metadata$Cell_type_abbreviation_level_1

color_vector_level2 <- metadata$Cell_type_color
names(color_vector_level2) <- metadata$Cell_type_abbreviation

ggplot(umap_plot) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill=level_1), shape=21, color = 'black', size = 1.5)+
  scale_fill_manual(values=color_vector_level1)+
  theme_classic()+ NoAxes()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+guides(fill = guide_legend(override.aes = list(size = 3), title=""))



ggplot(umap_plot) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill=level_2), shape=21, color = 'black', size = 1.5)+
  scale_fill_manual(values=color_vector_level2)+
  theme_classic()+ NoAxes()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+guides(fill = guide_legend(override.aes = list(size = 3), title=""))+
  geom_label_repel(aes(x = UMAP_1, y = UMAP_2, label=level_2), size=5, force=5, data=label.cent_level_2)

subplotteR=function(t_sub){
  umap_plot_sub <- as.data.frame(t_sub[["umap"]]@cell.embeddings)
  umap_plot_sub$level_2 <- unname(t_sub$final_labels)
  label.cent_sub_level_2 = umap_plot_sub %>% group_by(level_2) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
  
  print(ggplot(umap_plot_sub) + 
          geom_point(aes(x = UMAP_1, y = UMAP_2, fill=level_2), shape=21, color = 'black', size = 3)+
          scale_fill_manual(values=color_vector_level2)+
          theme_classic()+ NoAxes()+
          theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
          theme(legend.text=element_text(size=20))+guides(fill = guide_legend(override.aes = list(size = 3), title=""))+
          geom_label_repel(aes(x = UMAP_1, y = UMAP_2, label=level_2), size=5, data=label.cent_sub_level_2))
}

subplotteR(t_sc)
subplotteR(t_epc)
subplotteR(t_tc)
subplotteR(t_apc)

lotR_alt("matrix/spleen")

t_combined_alt <-  merge(`0002_mouse1_sp_filtered`, y = c(`0002_mouse10_sp_filtered`,`0002_mouse11_sp_filtered`,`0002_mouse12_sp_filtered`,
                                                          `0002_mouse13_sp_filtered`,`0002_mouse14_sp_filtered`,`0002_mouse15_sp_filtered`,
                                                          `0002_mouse16_sp_filtered`,`0002_mouse2_sp_filtered`,`0002_mouse3_sp_filtered`,
                                                          `0002_mouse6_sp_filtered`,`0002_mouse8_sp_filtered`,`0002_mouse7_sp_filtered`,
                                                          `0002_mouse9_sp_filtered`), project = "mouse_sp", merge.data=TRUE)     

cell_names <- colnames(t_combined)
t_combined_alt_subset <- subset(t_combined_alt, cells = cell_names)
t_combined[["SCT"]] <- t_combined_alt_subset[["SCT"]]


#adding to Seurat object phase (condition) information and producing DF with information about cell identities and phase (condition) to be used with 
#cell_adundance script
t_list <- count_produceR(t_combined)
t_combined_identity <- t_list[[2]]
t_combined_identity <- t_combined_identity%>% mutate(condition = case_when(
  batch == "0002_mouse1_sp" | batch == "0002_mouse6_sp" | batch == "0002_mouse12_sp" | batch == "0002_mouse13_sp" ~ "proestrus",
  batch == "0002_mouse2_sp" | batch == "0002_mouse8_sp" | batch == "0002_mouse7_sp"  | batch == "0002_mouse16_sp"  ~ "diestrus", 
  batch == "0002_mouse3_sp" | batch == "0002_mouse9_sp" | batch == "0002_mouse14_sp" ~ "estrus",
  batch == "0002_mouse10_sp" | batch == "0002_mouse11_sp" | batch == "0002_mouse15_sp" ~ "metestrus"))

save(t_combined,file="sp_combined_new.Rdata")

write.table(t_combined_identity,"sp_combined_identity",sep=",",quote=F,row.names = F)
