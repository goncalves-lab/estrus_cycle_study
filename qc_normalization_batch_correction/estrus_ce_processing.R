#################################
# this script uses raw filtered matrices deposited in Array Arrayexpress under accession number E-MTAB-11491 and 
#marker list provided in Supplementary file 1
#################################



library("ggrepel")

source("qc_normalization_cluster_annotation.R")
#manual anotation
lotR("matrix/cervix")

t_combined <- merge(`0002_mouse1_ce_filtered`, y = c(`0002_mouse10_ce_filtered`,`0002_mouse11_ce_filtered`,`0002_mouse12_ce_filtered`,
                                                     `0002_mouse13_ce_filtered`,`0002_mouse14_ce_filtered`,`0002_mouse15_ce_filtered`,
                                                     `0002_mouse16_ce_filtered`,`0002_mouse2_ce_filtered`,`0002_mouse3_ce_filtered`,
                                                     `0002_mouse6_ce_filtered`,`0002_mouse7_ce_filtered`,`0002_mouse8_ce_filtered`,
                                                     `0002_mouse9_ce_filtered`), project = "mouse_ce", merge.data=TRUE)       

t_combined_cell_metadata <- data.frame("cell_names"=colnames(t_combined), "batch" = t_combined$orig.ident)
t_combined_cell_metadata <- t_combined_cell_metadata %>% mutate(condition = case_when(
  batch == "0002_mouse1_ce" | batch == "0002_mouse6_ce" | batch == "0002_mouse12_ce" | batch == "0002_mouse13_ce" ~ "proestrus",
  batch == "0002_mouse2_ce" | batch == "0002_mouse8_ce" | batch == "0002_mouse7_ce" | batch == "0002_mouse16_ce"  ~ "diestrus", 
  batch == "0002_mouse3_ce" | batch == "0002_mouse9_ce" | batch == "0002_mouse14_ce" ~ "estrus",
  batch == "0002_mouse10_ce" | batch == "0002_mouse11_ce" | batch == "0002_mouse15_ce" ~ "metestrus"))


t_combined <- FindVariableFeatures(t_combined, selection.method = "vst", nfeatures = 2000)
t_combined <- reductoR(t_combined, 1.5)

print(DimPlot(t_combined, reduction = "umap", group.by = "seurat_clusters", label = T, pt.size = 0.7,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))


markers <- list.load("cervix_markers.rdata")
markers_df <- data.frame(matrix(unlist(markers), nrow=length(markers), byrow=T),stringsAsFactors=FALSE)
markers_df <- as.data.frame(t(markers_df))
rownames(markers_df) <- NULL
colnames(markers_df) <- "markers"

list[t_combined_clusters, t_subset_clusters]=lotR_clusters(t_combined,markers_df)
#writing out DF with information about cell identities and phase (condition) to be used with 
#cell_adundance script
write.table(t_combined_clusters,"ce_combined_clusters",sep=",",quote=F,row.names = F)
write.table(t_subset_clusters,"ce_combined_clusters_subset",sep=",",quote=F,row.names = F)

new_cluster_id_t_com=c("EpC", "EpC" , "EpC", "?", "EpC", "APC", "EpC", "APC", "EpC","EpC","?", "EpC", "EpC", "?", "?", "EpC", "SC", "EpC", "EpC", "APC", 
                       "?", "EpC", "EpC","TC", "EpC", "EpC", "?","SC", "?", "EpC", "?", "?", "EpC", "EpC", "SC")
t_combined <- lotR_ultimate(t_combined, new_cluster_id_t_com)

t <- t_combined$seurat_clusters
t <- data.frame("cells" = names(t), "type" = unname(t))

t <- t %>% mutate(cell_names = case_when(
  type == 5 | type == 7 | type == 19  ~ "APC",
  type == 16 | type == 27 | type == 34 ~ "SC",
  type == 0 | type == 1 | type == 2 | type ==  4 | type == 6 | type == 8 | type == 9 | type ==  11 | type == 12 | type == 15
  | type == 17 | type == 18 | type == 21 | type == 22 | type == 24 | type == 25 | type == 29 | type == 32 | type == 33 ~ "EpC",
  type == 3 | type == 10 | type == 13 | type == 14 | type == 20 | type == 28 | type == 30 | type == 31 | type == 26 ~ "?",
  type == 23  ~ "TC"
))

t_label <- t$cell_names
names(t_label) <- t$cells
t_combined$top_level <- t_label


#garnett
rownames(t_combined_cell_metadata) <- t_combined_cell_metadata$cell_names
t_combined_cell_metadata$tissue <- "cervix"


t_combined_cds <- garnett_annotatoR_part1(t_combined, t_combined_cell_metadata, "cervix_markers_garnett")
t_combined <- garnett_annotatoR_part2(t_combined_cds, t_combined_cds, "cervix_markers_garnett", t_combined)

garnett_labels <- t_combined$labels
seurat_clusters <- t_combined$seurat_clusters

garnett_annotation_ex <- data.frame("garnett" = garnett_labels, "seurat_clusters" = seurat_clusters)

garnett_cluster <- garnett_extender(t_combined) 

garnett_annotation_ex <- garnett_annotation_ex %>% mutate(cell_labels_garnett_ex = case_when(
  seurat_clusters == 5 | seurat_clusters == 7 | seurat_clusters == 19  ~ "APC",
  seurat_clusters == 16  | seurat_clusters == 34 | seurat_clusters == 27 ~ "SC",
  seurat_clusters ==  0 | seurat_clusters == 2 | seurat_clusters == 3 | seurat_clusters ==  6 | seurat_clusters == 8 |
    seurat_clusters == 13 | seurat_clusters == 15 | seurat_clusters == 18 | seurat_clusters == 21 | seurat_clusters == 22 | seurat_clusters == 24  ~ "EpC",
  seurat_clusters == 2  | seurat_clusters == 1 | seurat_clusters == 4 | seurat_clusters == 9 | seurat_clusters == 10 | seurat_clusters == 11 | seurat_clusters == 12 
  | seurat_clusters == 14 | seurat_clusters == 17 | seurat_clusters == 20 | seurat_clusters == 21 | seurat_clusters == 25 | seurat_clusters == 26
  | seurat_clusters == 27 | seurat_clusters == 28 | seurat_clusters == 30 | seurat_clusters == 31 | seurat_clusters == 32 | seurat_clusters == 33
  | seurat_clusters == 29 ~ "Unknown",
  seurat_clusters == 23 ~ "TC"
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


#writing out Seurat object
save(t_combined,file="ov_combined_new.Rdata")

list[t_non_match_ge, t_match_ge] <- svm_subseteR(t_combined)

pred.te_identity <- svm_returneR(t_match_ge, t_non_match_ge, 0.9, 0.5)
Idents(t_combined) <- t_combined$match
t_match <- subset(t_combined, idents = "match") 
t_non_match <- subset(t_combined, idents = "non-match") 

t_non_match$labels <- pred.te_identity
svm_cluster <- garnett_extender(t_non_match)
seurat_clusters <- t_non_match$seurat_clusters

svm_annotation_ex <- data.frame("svm" = pred.te_identity, "seurat_clusters" = seurat_clusters)
svm_annotation_ex <- svm_annotation_ex %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 1  | seurat_clusters == 3 | seurat_clusters == 4 | seurat_clusters == 9 | seurat_clusters == 10 | seurat_clusters == 11 
  | seurat_clusters ==  12 | seurat_clusters == 13 | seurat_clusters == 14 | seurat_clusters == 17 | seurat_clusters == 20| seurat_clusters == 25 
  | seurat_clusters ==  26 | seurat_clusters == 28 | seurat_clusters == 29 | seurat_clusters == 30 | seurat_clusters == 31| seurat_clusters == 32| seurat_clusters == 33 ~ "EpC"
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

save(t_combined,file="ce_combined_new.Rdata")

########################
#Ep
#######################

Idents(t_combined) <- t_combined$final_labels
list[t_epc, epc_cluster, epc_subcluster] <- subseteR_part_1(t_combined, "EpC", 1)
new_cluster_id_epc=c("PB/BEpC", "PB/BEpC", "PB/BEpC", "?", "PB/BEpC", "PB/BEpC", "S/IEpC", "?", "PB/BEpC", "?","?", "S/IEpC", "?", "PB/BEpC", "CEpC", "?",
                     "PB/BEpC", "PB/BEpC", "?", "PB/BEpC", "S/IEpC", "?", "?", "GC")
epc <- subseteR_part_2(t_epc, new_cluster_id_epc)

epc <- epc %>% mutate(cell_names = case_when(
  type == 6 | type == 11  | type == 20 | type == 3 ~ "S/IEpC",
  type == 0 | type == 1 | type == 2 | type ==4 | type==5  | type ==8 | type==13| type ==16 | type==17 | type ==19   ~ "P/BEpC",
  type == 14 ~ "CEpC",
  type == 23 ~ "GC",
  type== 15 | type == 7 | type == 9 | type == 10 | type == 12 | type ==18 | type==21| type ==22 ~ "?"
  
))

t_epc$sublevel <- epc$cell_names

print(DimPlot(t_epc, reduction = "umap", label = T ,pt.size = 0.5,group.by="sublevel",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

t_epc_cell_metadata <- data.frame("cell_names"=colnames(t_epc), "batch" = t_epc$orig.ident)
t_epc_cell_metadata$tissue <- "cervix"
t_combined_cds_ep <- garnett_annotatoR_part1(t_epc, t_epc_cell_metadata, "ce_markers_ep")
t_epc <- garnett_annotatoR_part2(t_combined_cds_ep,t_combined_cds_ep, "ce_markers_ep", t_epc)

garnett_labels_epc <- t_epc$labels
seurat_clusters_epc <- t_epc$seurat_clusters

garnett_annotation_ex_epc <- data.frame("garnett" = garnett_labels_epc, "seurat_clusters" = seurat_clusters_epc)

garnett_cluster_epc <- garnett_extender(t_epc) 

garnett_annotation_ex_epc <- garnett_annotation_ex_epc %>% mutate(cell_labels_garnett_ex = case_when(
  seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 5 | seurat_clusters == 6 | seurat_clusters == 11 | seurat_clusters == 13 
  | seurat_clusters == 15 | seurat_clusters == 16 | seurat_clusters == 17 | seurat_clusters == 19 | seurat_clusters == 22 | seurat_clusters == 21~ "P/BEpC",
  seurat_clusters == 23 ~ "GC",
  seurat_clusters ==  3 | seurat_clusters == 4 | seurat_clusters == 8 | seurat_clusters == 20  ~ "S/IEpC",
  seurat_clusters == 7 | seurat_clusters == 9 | seurat_clusters ==  10 | seurat_clusters == 12 | seurat_clusters == 14 | seurat_clusters == 18   ~ "CEpC", 
))

t_epc$garnett_cluster <- garnett_annotation_ex_epc$cell_labels_garnett_ex
t_epc$labels <- NULL

print(DimPlot(t_epc, reduction = "umap", label = T ,pt.size = 0.5,group.by="garnett_cluster",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

garnett <- t_epc$garnett_cluster
manual <- t_epc$sublevel
annotation_df_epc <- data.frame("cells" = names(garnett), "garnett" = unname(garnett), "manual" = manual)
annotation_df_epc$matched_cells = "?"

for (i in seq(from=1, to=nrow(annotation_df_epc))){
  if (annotation_df_epc$garnett[i] == annotation_df_epc$manual[i]){
    annotation_df_epc$matched_cells[i] = "match"
  } else (annotation_df_epc$matched_cells[i] = "non-match")
}

matched_cells_epc <- annotation_df_epc$matched_cells
names(matched_cells_epc) <- annotation_df_epc$cells
t_epc$match <- matched_cells_epc
print(DimPlot(t_epc, reduction = "umap", group.by = "match", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

list[t_non_match_ge_epc, t_match_ge_epc] <- svm_subseteR(t_epc)

pred.te_identity_epc <- svm_returneR(t_match_ge_epc, t_non_match_ge_epc, 0.7, 0.2)
Idents(t_epc) <- t_epc$match
t_match_epc <- subset(t_epc, idents = "match") 
t_non_match_epc <- subset(t_epc, idents = "non-match") 

t_non_match_epc$labels <- pred.te_identity_epc
svm_cluster_epc <- garnett_extender(t_non_match_epc)
seurat_clusters_epc <- t_non_match_epc$seurat_clusters

svm_annotation_ex_epc <- data.frame("svm" = pred.te_identity_epc, "seurat_clusters" = seurat_clusters_epc)

svm_annotation_ex_epc <- svm_annotation_ex_epc %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 4 | seurat_clusters == 6 |  seurat_clusters == 7 | seurat_clusters == 15 | seurat_clusters == 8 | seurat_clusters == 9 
  | seurat_clusters == 10 | seurat_clusters == 11 | seurat_clusters == 12 | seurat_clusters == 18 | seurat_clusters == 21 | seurat_clusters == 22 ~ "P/BEpC"
))

t_non_match_epc$svm_cluster <- svm_annotation_ex_epc$cell_labels_smv_ex

assigned.te_epc <- as.data.frame(t_match_epc$sublevel) 
predicted.te_epc <- as.data.frame(t_non_match_epc$svm_cluster)
colnames(assigned.te_epc) <- "identity"
colnames(predicted.te_epc) <- "identity"

final_labels_epc <- rbind (assigned.te_epc, predicted.te_epc)
final_labels_epc$cells <- rownames(final_labels_epc)

manual_labels_epc <- as.data.frame(t_epc$sublevel) 
manual_labels_epc$cells <- rownames(manual_labels_epc)
final_labels_epc <- final_labels_epc[manual_labels_epc[["cells"]],]
final_labels_epc_vector <- final_labels_epc$identity
names(final_labels_epc_vector) <- final_labels_epc$cells
t_epc$final_labels <- final_labels_epc_vector

print(DimPlot(t_epc, reduction = "umap", group.by = "final_labels", label = T, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

save(t_epc,file="ce_epc_final.Rdata")

########################
#SC
#######################

Idents(t_combined) <- t_combined$final_labels
list[t_sc, sc_cluster, sc_subcluster] <- subseteR_part_1(t_combined, "SC", 1)
new_cluster_id_sc=c("F", "?", "?", "F", "F", "MC","MC", "F", "EC",  "F", "?", "F", "?")
sc <- subseteR_part_2(t_sc, new_cluster_id_sc)

sc <- sc %>% mutate(cell_names = case_when(
  type == 0 | type == 3  | type == 4 | type == 7 | type ==  9 | type ==  11 ~ "F",
  type == 5  | type == 6  ~ "MC",
  type == 1 | type == 2 | type == 10 | type == 12~ "?",
  type == 8 ~ "EC"
  
))

t_sc$sublevel <- sc$cell_names

print(DimPlot(t_sc, reduction = "umap", label = T ,pt.size = 0.5,group.by="sublevel",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

t_sc_cell_metadata <- data.frame("cell_names"=colnames(t_sc), "batch" = t_sc$orig.ident)
t_sc_cell_metadata$tissue <- "cervix"
t_combined_cds_ep <- garnett_annotatoR_part1(t_sc, t_sc_cell_metadata, "ce_markers_sc")
t_sc <- garnett_annotatoR_part2(t_combined_cds_ep,t_combined_cds_ep, "ce_markers_sc", t_sc)

garnett_labels_sc <- t_sc$labels
seurat_clusters_sc <- t_sc$seurat_clusters

garnett_annotation_ex_sc <- data.frame("garnett" = garnett_labels_sc, "seurat_clusters" = seurat_clusters_sc)

garnett_cluster_sc <- garnett_extender(t_sc) 

garnett_annotation_ex_sc <- garnett_annotation_ex_sc %>% mutate(cell_labels_garnett_ex = case_when(
  seurat_clusters ==  0 | seurat_clusters == 1 | seurat_clusters == 3 | seurat_clusters ==  4 | seurat_clusters == 9 | 
    seurat_clusters == 7 | seurat_clusters == 9 | seurat_clusters == 11  ~ "F",
  seurat_clusters == 5 | seurat_clusters == 6 ~ "MC",
  seurat_clusters == 2 | seurat_clusters == 10 | seurat_clusters == 12 ~ "Unknown",
  seurat_clusters == 8 ~ "EC"
  
))

t_sc$garnett_cluster <- garnett_annotation_ex_sc$cell_labels_garnett_ex
t_sc$labels <- NULL

print(DimPlot(t_sc, reduction = "umap", label = T ,pt.size = 0.5,group.by="garnett_cluster",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

garnett <- t_sc$garnett_cluster
manual <- t_sc$sublevel
annotation_df_sc <- data.frame("cells" = names(garnett), "garnett" = unname(garnett), "manual" = manual)
annotation_df_sc$matched_cells = "?"

for (i in seq(from=1, to=nrow(annotation_df_sc))){
  if (annotation_df_sc$garnett[i] == annotation_df_sc$manual[i]){
    annotation_df_sc$matched_cells[i] = "match"
  } else (annotation_df_sc$matched_cells[i] = "non-match")
}

matched_cells_sc <- annotation_df_sc$matched_cells
names(matched_cells_sc) <- annotation_df_sc$cells
t_sc$match <- matched_cells_sc
print(DimPlot(t_sc, reduction = "umap", group.by = "match", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

list[t_non_match_ge_sc, t_match_ge_sc] <- svm_subseteR(t_sc)

pred.te_identity_sc <- svm_returneR(t_match_ge_sc, t_non_match_ge_sc, 0.7, 0.9)
Idents(t_sc) <- t_sc$match
t_match_sc <- subset(t_sc, idents = "match") 
t_non_match_sc <- subset(t_sc, idents = "non-match") 

t_non_match_sc$labels <- pred.te_identity_sc
svm_cluster_sc <- garnett_extender(t_non_match_sc)
seurat_clusters_sc <- t_non_match_sc$seurat_clusters

svm_annotation_ex_sc <- data.frame("svm" = pred.te_identity_sc, "seurat_clusters" = seurat_clusters_sc)

svm_annotation_ex_sc <- svm_annotation_ex_sc %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 12 ~ "F",
  seurat_clusters == 10 ~ "SC I"
))

t_non_match_sc$svm_cluster <- svm_annotation_ex_sc$cell_labels_smv_ex

assigned.te_sc <- as.data.frame(t_match_sc$sublevel) 
predicted.te_sc <- as.data.frame(t_non_match_sc$svm_cluster)
colnames(assigned.te_sc) <- "identity"
colnames(predicted.te_sc) <- "identity"

final_labels_sc <- rbind (assigned.te_sc, predicted.te_sc)
final_labels_sc$cells <- rownames(final_labels_sc)

manual_labels_sc <- as.data.frame(t_sc$sublevel) 
manual_labels_sc$cells <- rownames(manual_labels_sc)
final_labels_sc <- final_labels_sc[manual_labels_sc[["cells"]],]
final_labels_sc_vector <- final_labels_sc$identity
names(final_labels_sc_vector) <- final_labels_sc$cells
t_sc$final_labels <- final_labels_sc_vector

print(DimPlot(t_sc, reduction = "umap", group.by = "final_labels", label = F, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

save(t_sc,file="ce_sc_final.Rdata")

########################
#TC
#######################

Idents(t_combined) <- t_combined$final_labels
list[t_tc, tc_cluster, tc_subcluster] <- subseteR_part_1(t_combined, "TC", 1)
new_cluster_id_tc=c("MTC", "MTC", "MTC","MTC", "?","MTC", "iNKT", "?")
tc <- subseteR_part_2(t_tc, new_cluster_id_tc)

tc <- tc %>% mutate(cell_names = case_when(
  type == 0 ~ "MAIT",
  type == 6  ~ "iNKT",
  type == 1  | type == 2 | type == 3 | type == 5 ~ "MTC",
  type == 4 | type == 7 ~ "?"))


t_tc$sublevel <- tc$cell_names

print(DimPlot(t_tc, reduction = "umap", label = T ,pt.size = 0.5,group.by="sublevel",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

t_tc_cell_metadata <- data.frame("cell_names"=colnames(t_tc), "batch" = t_tc$orig.ident)
t_tc_cell_metadata$tissue <- "cervix"
t_combined_cds_ep <- garnett_annotatoR_part1(t_tc, t_tc_cell_metadata, "ce_markers_tc")
t_tc <- garnett_annotatoR_part2(t_combined_cds_ep,t_combined_cds_ep, "ce_markers_tc", t_tc)


garnett_labels_tc <- t_tc$labels
seurat_clusters_tc <- t_tc$seurat_clusters

garnett_annotation_ex_tc <- data.frame("garnett" = garnett_labels_tc, "seurat_clusters" = seurat_clusters_tc)

garnett_cluster_tc <- garnett_extender(t_tc) 

garnett_annotation_ex_tc <- garnett_annotation_ex_tc %>% mutate(cell_labels_garnett_ex = case_when(
  seurat_clusters ==  0 | seurat_clusters == 2  | seurat_clusters == 5 ~ "MAIT",
  seurat_clusters == 3 | seurat_clusters == 6 ~ "iNKT",
  seurat_clusters == 4 | seurat_clusters == 7   ~ "Unknown",
  seurat_clusters == 1 ~ "MTC"
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
  seurat_clusters == 2 | seurat_clusters == 3 | seurat_clusters == 4 | seurat_clusters == 7  ~ "MTC",
  seurat_clusters == 5 ~ "MAIT"
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

save(t_tc,file="ce_tc_final.Rdata")


########################
#APC
#######################

Idents(t_combined) <- t_combined$final_labels
list[t_apc, apc_cluster, apc_subcluster] <- subseteR_part_1(t_combined, "APC", 1)
new_cluster_id_apc=c("N", "M1Mp", "M1Mp", "?", "M1Mp", "DC", "M1Mp", "DC", "?","DC","DC")
apc <- subseteR_part_2(t_apc, new_cluster_id_apc)

apc <- apc %>% mutate(cell_names = case_when(
  type == 5 | type == 7 | type == 9 | type== 10 ~ "DC",
  type == 1 | type == 2  | type == 4 | type == 6 ~ "M1Mp",
  type == 0 ~ "N",
  type == 3 | type == 8 ~ "?"
  
))

t_apc$sublevel <- apc$cell_names

print(DimPlot(t_apc, reduction = "umap", label = T ,pt.size = 0.5,group.by="sublevel",cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

t_apc_cell_metadata <- data.frame("cell_names"=colnames(t_apc), "batch" = t_apc$orig.ident)
t_apc_cell_metadata$tissue <- "cervix"
t_combined_cds_ep <- garnett_annotatoR_part1(t_apc, t_apc_cell_metadata, "ce_markers_apc")
t_apc <- garnett_annotatoR_part2(t_combined_cds_ep,t_combined_cds_ep, "ce_markers_apc", t_apc)


garnett_labels_apc <- t_apc$labels
seurat_clusters_apc <- t_apc$seurat_clusters

garnett_annotation_ex_apc <- data.frame("garnett" = garnett_labels_apc, "seurat_clusters" = seurat_clusters_apc)

garnett_cluster_apc <- garnett_extender(t_apc) 

garnett_annotation_ex_apc <- garnett_annotation_ex_apc %>% mutate(cell_labels_garnett_ex = case_when(
  seurat_clusters == 2 | seurat_clusters == 4 ~ "M1Mp",
  seurat_clusters == 0 | seurat_clusters == 6 ~ "N",
  seurat_clusters == 5 | seurat_clusters == 7 | seurat_clusters == 10 | seurat_clusters == 9 | seurat_clusters == 8  ~ "DC",
  seurat_clusters == 1 | seurat_clusters ==  3 ~ "Unknown" 
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

pred.te_identity_apc <- svm_returneR(t_match_ge_apc, t_non_match_ge_apc, 0.8, 0.9)
Idents(t_apc) <- t_apc$match
t_match_apc <- subset(t_apc, idents = "match") 
t_non_match_apc <- subset(t_apc, idents = "non-match") 

t_non_match_apc$labels <- pred.te_identity_apc
svm_cluster_apc <- garnett_extender(t_non_match_apc)
seurat_clusters_apc <- t_non_match_apc$seurat_clusters

svm_annotation_ex_apc <- data.frame("svm" = pred.te_identity_apc, "seurat_clusters" = seurat_clusters_apc)

svm_annotation_ex_apc <- svm_annotation_ex_apc %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 3 ~ "N",
  seurat_clusters == 1 | seurat_clusters == 6~ "M1Mp",
  seurat_clusters == 8 ~ "DC"
  
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

save(t_apc,file="ce_apc_final.Rdata")


################
#Unknown
##############

########################
#Level 2
#######################

final_labels_level_2 <- rbind (final_labels_apc, final_labels_epc, final_labels_sc, final_labels_tc)

manual_labels <- as.data.frame(t_combined$final_labels) 
manual_labels$cells <- rownames(manual_labels)
final_labels_level_2 <- final_labels_level_2[manual_labels[["cells"]],]
final_labels_level_2_vector <- final_labels_level_2$identity
names(final_labels_level_2_vector) <- final_labels_level_2$cells
t_combined$level_2 <- final_labels_level_2_vector

save(t_combined,file="ce_combined_new.Rdata")


#final plots
level_1 <- data.frame("cells"=colnames(t_combined), "old_labels"=t_combined$final_labels)

level_1 <- level_1 %>% mutate(new_labels = case_when(
  old_labels == "APC" | old_labels == "TC" | old_labels == "MaC" | old_labels == "NKC" ~ "IC",
  old_labels == "EC" ~ "EC",
  old_labels == "EpC" ~ "EpC",
  old_labels == "SC" ~ "SC"))
t_combined$level_1 <- level_1$new_labels


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
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill=level_2), shape=21, color = 'black', size = 5)+
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

lotR_alt("matrix/cervix")

t_combined_alt <-  merge(`0002_mouse1_ce_filtered`, y = c(`0002_mouse10_ce_filtered`,`0002_mouse11_ce_filtered`,`0002_mouse12_ce_filtered`,
                                                          `0002_mouse13_ce_filtered`,`0002_mouse14_ce_filtered`,`0002_mouse15_ce_filtered`,
                                                          `0002_mouse16_ce_filtered`,`0002_mouse2_ce_filtered`,`0002_mouse3_ce_filtered`,
                                                          `0002_mouse6_ce_filtered`,`0002_mouse7_ce_filtered`,`0002_mouse8_ce_filtered`,
                                                          `0002_mouse9_ce_filtered`), project = "mouse_ce", merge.data=TRUE)       

cell_names <- colnames(t_combined)
t_combined_alt_subset <- subset(t_combined_alt, cells = cell_names)
t_combined[["SCT"]] <- t_combined_alt_subset[["SCT"]]

#adding to Seurat object phase (condition) information and producing DF with information about cell identities and phase (condition) to be used with 
#cell_adundance script
t_list <- count_produceR(t_combined)
t_combined_identity <- t_list[[2]]
t_combined_identity <- t_combined_identity%>% mutate(condition = case_when(
  batch == "0002_mouse1_ce" | batch == "0002_mouse6_ce" | batch == "0002_mouse12_ce" | batch == "0002_mouse13_ce" ~ "proestrus",
  batch == "0002_mouse2_ce" | batch == "0002_mouse8_ce" | batch == "0002_mouse7_ce" | batch == "0002_mouse16_ce"  ~ "diestrus", 
  batch == "0002_mouse3_ce" | batch == "0002_mouse9_ce" | batch == "0002_mouse14_ce" ~ "estrus",
  batch == "0002_mouse10_c" | batch == "0002_mouse11_c" | batch == "0002_mouse15_ce" ~ "metestrus"))

save(t_combined,file="ce_combined_new.Rdata")

write.table(t_combined_identity,"ce_combined_identity",sep=",",quote=F,row.names = F)
