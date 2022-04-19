#################################
# this script uses raw filtered matrices deposited in Array Arrayexpress under accession number E-MTAB-11491, 
#marker list provided in Supplementary file 1 and seurat object generated in script estrus_ut_processing.R
#################################


library(ggrepel)

source("qc_normalization_cluster_annotation.R")

lotR_alt("/decidua/matrix/")
t_combined_deci = merge(`0002_DC_mmus01_ut_run_filtered`, y = c(`0002_DC_mmus03_ut_run_filtered`,
                                                                `OE0538_DO-0002_mmus_D06__uterus01_run_filtered`,
                                                                `OE0538_DO-0002_mmus_D07__uterus01_run_filtered`,
                                                                `OE0538_DO-0002_mmus_D08__uterus01_run_filtered`,
                                                                `OE0538_DO-0002_mmus_D09__uterus01_run_filtered`,
                                                                `OE0538_DO-0002_mmus_D10__uterus01_run_filtered`), 
                        project = "mouse_ov", merge.data=TRUE)
VariableFeatures(t_combined_deci) <- Reduce(intersect, list(`0002_DC_mmus03_ut_run_filtered`$SCT@var.features,
                                                            `OE0538_DO-0002_mmus_D06__uterus01_run_filtered`$SCT@var.features,
                                                            `OE0538_DO-0002_mmus_D07__uterus01_run_filtered`$SCT@var.features,
                                                            `OE0538_DO-0002_mmus_D08__uterus01_run_filtered`$SCT@var.features,
                                                            `OE0538_DO-0002_mmus_D09__uterus01_run_filtered`$SCT@var.features,
                                                            `OE0538_DO-0002_mmus_D10__uterus01_run_filtered`$SCT@var.features))
t_combined_deci <- reductoR(t_combined_deci, 1)

DimPlot(t_combined_deci, reduction = "umap",group.by = "orig.ident", split.by = "orig.ident")
save(t_combined_deci,file="ut_combined_deci_merged.Rdata")

svm_subseteR <- function(t_combined, y_var){
  t_ge <- t_combined[["SCT"]]@counts 
  #t_combined <- FindVariableFeatures(t_combined, selection.method = "vst", nfeatures = 2000)
  var_genes <- t_combined@assays$SCT@var.features
  t_ge <- t_ge[rownames(t_ge)%in%var_genes,]
  t_ge <- t(t_ge)  
  y <- as.factor(y_var)
  t_ge <- data.frame(x=t_ge, y=y)
  list(t_ge)}

svm_returneR <- function(t_ge, t_deci, subset_value){
  prob <- c(0.7, 0.8, 0.9)
  pred.te_identity_final <- data.frame(row.names = row.names(t_deci))
  set.seed(1)
  t_train <- data.frame()
  t_ge$y <- as.factor(t_ge$y)
  for (i in levels(t_ge$y)){
    subset_size <- round(table(t_ge$y)[i]*subset_value)
    t_ge_subset <- t_ge[grep(i, t_ge$y),]
    t_ge_train_subset <- t_ge_subset[sample.int(subset_size),]
    t_train <- rbind(t_train, t_ge_train_subset)
  }
  
  
  out <- svm(y~., data = t_train, kernel ="linear", cost=10, scale = F, probability = T)
  pred.te <- predict (out, newdata = t_deci, probability = T)
  pred.te_probabilites <- attr(pred.te,"probabilities")
  
  
  for (j in prob){
    pred.te_identity <- c()
    for (i in seq(from=1, to=nrow(pred.te_probabilites))){
      if (T %in% (pred.te_probabilites[i,] > j)){
        pred.te_identity[i] <- names(pred.te_probabilites[i,][pred.te_probabilites[i,]==max(pred.te_probabilites[i,])])
      } else{pred.te_identity[i] <- "unknown"}
    }
    pred.te_identity <- as.data.frame(pred.te_identity)
    colnames(pred.te_identity) <- paste("identity", j, sep="_") 
    row.names(pred.te_identity) <- row.names(pred.te_probabilites)
    row.names(pred.te_identity_final) <- row.names(pred.te_identity)
    pred.te_identity_final <- cbind(pred.te_identity, pred.te_identity_final)
  }
  return(pred.te_identity_final)}


t_combined <- loadRData("ut_combined.Rdata")


Idents(t_combined_deci) <-t_combined_deci$seurat_clusters
list[t_ge_lv_1] <- svm_subseteR(t_combined, t_combined$level_1)
list[t_ge_lv_2] <- svm_subseteR(t_combined, t_combined$level_2)
list[t_deci] <- svm_subseteR(t_combined_deci, t_combined_deci$seurat_clusters)
t_ge_lv_1 <- t_ge_lv_1[,colnames(t_ge_lv_1)%in%colnames(t_deci)]
t_ge_lv_2 <- t_ge_lv_2[,colnames(t_ge_lv_2)%in%colnames(t_deci)]


markers <- list.load("uterus_markers.rdata")
markers_df <- data.frame(matrix(unlist(markers), nrow=length(markers), byrow=T),stringsAsFactors=FALSE)
markers_df <- as.data.frame(t(markers_df))
rownames(markers_df) <- NULL
colnames(markers_df) <- "markers"
list[t_combined_clusters, t_subset_clusters]=lotR_clusters(t_combined_deci,markers_df)

t_combined <- RunUMAP(t_combined, dim=1:30)

print(DimPlot(t_combined_deci, reduction = "umap", group.by = "seurat_clusters", label= T,pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))


pred.te_identity_lv_1 <- svm_returneR(t_ge_lv_1, t_deci, 0.5)
t_combined_deci$labels <- pred.te_identity_lv_1$identity_0.9
svm_cluster_lv1 <- garnett_extender(t_combined_deci)

svm_annotation_ex_lv_1 <- data.frame("svm" = pred.te_identity_lv_1$identity_0.9, "seurat_clusters" = t_combined_deci$seurat_clusters)

svm_annotation_ex_lv_1 <- svm_annotation_ex_lv_1 %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 0  | seurat_clusters == 1 | seurat_clusters == 2 | seurat_clusters == 4 | 
    seurat_clusters == 5 | seurat_clusters == 7 | seurat_clusters == 6 |
    seurat_clusters == 8  | seurat_clusters == 9  | seurat_clusters == 11 | seurat_clusters == 12 |
    seurat_clusters == 13  | seurat_clusters == 14  | seurat_clusters == 16 | seurat_clusters == 17
  | seurat_clusters == 20 | seurat_clusters == 24 | seurat_clusters == 27~ "SC",
  seurat_clusters == 3 | seurat_clusters == 10  | seurat_clusters == 15 | seurat_clusters == 18 |
    seurat_clusters == 21  | seurat_clusters == 22 | seurat_clusters == 25 | seurat_clusters == 26  ~ "IC",  
  seurat_clusters == 19  ~ "EpC",
  seurat_clusters == 23  ~ "EC"
))

t_combined_deci$level_1 <- svm_annotation_ex_lv_1$cell_labels_smv_ex
#############
##EPC

Idents(t_combined) <- t_combined$level_1
#t_combined <- FindVariableFeatures(t_combined)
t_epc <- subset(x = t_combined, idents = "EpC")
t_epc <- normalizeR_alt(t_epc)
t_epc <- reductoR(t_epc, 1)

Idents(t_combined_deci) <- t_combined_deci$level_1
list[t_epc_deci, epc_cluster_deci, epc_subcluster_deci] <- subseteR_part_1(t_combined_deci, "EpC", 1)

list[t_ge_epc] <- svm_subseteR(t_epc, t_epc$level_2)
list[t_ge_epc_deci] <- svm_subseteR(t_epc_deci, t_epc_deci$seurat_clusters)
t_ge_epc <- t_ge_epc[,colnames(t_ge_epc)%in%colnames(t_ge_epc_deci)]

pred.te_identity_lv_1_epc <- svm_returneR(t_ge_epc, t_ge_epc_deci, 0.3)

t_epc_deci$labels <- pred.te_identity_lv_1_epc$identity_0.9
svm_cluster_epc <- garnett_extender(t_epc_deci)

svm_annotation_ex_epc <- data.frame("svm" = pred.te_identity_lv_1_epc$identity_0.9, "seurat_clusters" = t_epc_deci$seurat_clusters)

svm_annotation_ex_epc <- svm_annotation_ex_epc %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 0  | seurat_clusters == 1 ~ "CEpC",
  seurat_clusters ==2   ~ "GlC",
  seurat_clusters == 3 ~ "CC"
  
))

t_epc_deci$level_2 <- svm_annotation_ex_epc$cell_labels_smv_ex

print(DimPlot(t_epc_deci, reduction = "umap", group.by = "level_2", label = T, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))


#############
##SC

Idents(t_combined) <- t_combined$level_1
t_sc <- subset(x = t_combined, idents = "SC")
t_sc <- normalizeR_alt(t_sc)
t_sc <- reductoR(t_sc, 1)

Idents(t_combined_deci) <- t_combined_deci$level_1
list[t_sc_deci, sc_cluster_deci, sc_subcluster_deci] <- subseteR_part_1(t_combined_deci, "SC", 1)

list[t_ge_sc] <- svm_subseteR(t_sc, t_sc$level_2)
list[t_ge_sc_deci] <- svm_subseteR(t_sc_deci, t_sc_deci$seurat_clusters)
t_ge_sc <- t_ge_sc[,colnames(t_ge_sc)%in%colnames(t_ge_sc_deci)]

print(DimPlot(t_sc_deci, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident",label = T, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

pred.te_identity_lv_1_sc <- svm_returneR(t_ge_sc, t_ge_sc_deci, 0.3)

#CCA correction
#############################################################
sc.list <- SplitObject(t_sc_deci, split.by = "orig.ident")
sc.list <- lapply(X = sc.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 3000)
sc.list <- PrepSCTIntegration(object.list = sc.list, anchor.features = features)
sc.anchors <- FindIntegrationAnchors(object.list = sc.list, normalization.method = "SCT", 
                                     anchor.features = features)
sc.combined.sct <- IntegrateData(anchorset = sc.anchors, normalization.method = "SCT")
sc.combined.sct <- RunPCA(sc.combined.sct, verbose = FALSE)
sc.combined.sct <- RunUMAP(sc.combined.sct, reduction = "pca", dims = 1:30)
sc.combined.sct <- FindNeighbors(sc.combined.sct, reduction = "pca", dims = 1:30)
sc.combined.sct <- FindClusters(sc.combined.sct, resolution = 0.5)
print(DimPlot(sc.combined.sct, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident",label = T, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

#############################################################################

t_sc_deci$labels <- pred.te_identity_lv_1_sc$identity_0.9
svm_cluster_sc <- garnett_extender(t_sc_deci)

svm_annotation_ex_sc <- data.frame("svm" = pred.te_identity_lv_1_sc$identity_0.9, "seurat_clusters" = t_sc_deci$seurat_clusters)

svm_annotation_ex_sc <- svm_annotation_ex_sc %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 0 | seurat_clusters == 1  | seurat_clusters == 2 | seurat_clusters == 3 |
    seurat_clusters == 4  | seurat_clusters == 5 | seurat_clusters == 6 | seurat_clusters == 7 | seurat_clusters == 8 | seurat_clusters == 10 | seurat_clusters == 11
  | seurat_clusters == 13 | seurat_clusters == 15 | seurat_clusters == 16 | seurat_clusters == 17 ~ "F",  
  seurat_clusters == 18 ~ "MC",
  seurat_clusters == 19 ~ "SC I - DC",
  seurat_clusters == 9 | seurat_clusters == 12 | seurat_clusters == 14  ~ "DeC"
))

t_sc_deci$level_2 <- svm_annotation_ex_sc$cell_labels_smv_ex

print(DimPlot(t_sc_deci, reduction = "umap", group.by = "level_2", label = T, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

#############
#estrus and dc combined
t_combined_de <- merge(t_sc_deci, y=t_sc)
t_combined_de <- FindVariableFeatures(t_combined_de)
t_combined_de <- RunPCA(t_combined_de, verbose = FALSE)
t_combined_de <- RunUMAP(t_combined_de, reduction = "pca", dims = 1:30)
t_combined_de <- FindNeighbors(t_combined_de, reduction = "pca", dims = 1:30)
t_combined_de <- FindClusters(t_combined_de, resolution = 0.5)

print(DimPlot(t_combined_de, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident" , label = T, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))

##IC

Idents(t_combined) <- t_combined$level_1
t_ic <- subset(x = t_combined, idents = "IC")
t_ic <- normalizeR_alt(t_ic)
t_ic <- reductoR(t_ic, 1)

Idents(t_combined_deci) <- t_combined_deci$level_1
list[t_ic_deci, ic_cluster_deci, ic_subcluster_deci] <- subseteR_part_1(t_combined_deci, "IC", 1)

list[t_ge_ic] <- svm_subseteR(t_ic, t_ic$level_2)
list[t_ge_ic_deci] <- svm_subseteR(t_ic_deci, t_ic_deci$seurat_clusters)
t_ge_ic <- t_ge_ic[,colnames(t_ge_ic)%in%colnames(t_ge_ic_deci)]

pred.te_identity_lv_1_ic <- svm_returneR(t_ge_ic, t_ge_ic_deci, 0.9)

t_ic_deci$labels <- pred.te_identity_lv_1_ic$identity_0.7
svm_cluster_ic <- garnett_extender(t_ic_deci)

svm_annotation_ex_ic <- data.frame("svm" = pred.te_identity_lv_1_ic$identity_0.9, "seurat_clusters" = t_ic_deci$seurat_clusters)

svm_annotation_ex_ic <- svm_annotation_ex_ic %>% mutate(cell_labels_smv_ex = case_when(
  seurat_clusters == 2 | seurat_clusters == 11 | seurat_clusters == 15 ~ "DC",
  seurat_clusters == 0 | seurat_clusters == 1 | seurat_clusters == 5 | seurat_clusters == 6 | seurat_clusters == 9 ~ "NKC",  
  seurat_clusters ==12 ~ "MTC",
  seurat_clusters == 10 ~ "MAIT",
  seurat_clusters == 3 | seurat_clusters == 4 | seurat_clusters == 7 | seurat_clusters==8 ~ "M2Mp",
  seurat_clusters == 13 ~ "M1Mp",
  seurat_clusters == 14 ~ "MaC"
  
))


t_ic_deci$level_2 <- svm_annotation_ex_ic$cell_labels_smv_ex
print(DimPlot(t_ic_deci, reduction = "umap", group.by = "level_2", label = T, pt.size = 0.5, cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))



list[t_ec, ec_cluster, ec_subcluster] <- subseteR_part_1(t_combined_deci, "EC",1)
final_labels_ec <- data.frame("identity"="EC", "cells"=names(t_ec$seurat_clusters))
rownames(final_labels_ec) <- final_labels_ec$cells



final_labels_sc <- as.data.frame(svm_annotation_ex_sc$cell_labels_smv_ex)
rownames(final_labels_sc) <- colnames(t_sc_deci)
final_labels_sc$cells <- rownames(final_labels_sc)
colnames(final_labels_sc) <- c("identity", "cells")

final_labels_epc <- as.data.frame(svm_annotation_ex_epc$cell_labels_smv_ex)
rownames(final_labels_epc) <- colnames(t_epc_deci)
final_labels_epc$cells <- rownames(final_labels_epc)
colnames(final_labels_epc) <- c("identity", "cells")

final_labels_ic <- as.data.frame(svm_annotation_ex_ic$cell_labels_smv_ex)
rownames(final_labels_ic) <- colnames(t_ic_deci)
final_labels_ic$cells <- rownames(final_labels_ic)
colnames(final_labels_ic) <- c("identity", "cells")


final_labels_level_2 <- rbind (final_labels_ec, final_labels_sc, final_labels_ic, final_labels_epc )

manual_labels <- as.data.frame(t_combined_deci$seurat_clusters) 
manual_labels$cells <- rownames(manual_labels)
final_labels_level_2 <- final_labels_level_2[manual_labels[["cells"]],]
final_labels_level_2_vector <- final_labels_level_2$identity
names(final_labels_level_2_vector) <- final_labels_level_2$cells
t_combined_deci$level_2 <- final_labels_level_2_vector

metadata <- read.csv("/omics/groups/OE0433/internal/ivana/scripts/EstrusProject/metadata_cell_types.csv")

color_vector_level1 <- metadata$Cell_type_color_level_1
names(color_vector_level1) <- metadata$Cell_type_abbreviation_level_1

color_vector_level2 <- metadata$Cell_type_color
names(color_vector_level2) <- metadata$Cell_type_abbreviation

#final plots
t_combined_deci <- RunUMAP(t_combined_deci, dims=1:30, metric = "euclidean")

umap_plot <- as.data.frame(t_combined_deci[["umap"]]@cell.embeddings)
umap_plot$level_1 <- unname(t_combined_deci$level_1)
umap_plot$level_2 <- unname(t_combined_deci$level_2)
label.cent_level_2 = umap_plot %>% group_by(level_2) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)
label.cent_level_1 = umap_plot %>% group_by(level_1) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

color_vector_level2 <- color_vector_level2[names(color_vector_level2)%in%unique(t_combined_deci$level_2)]

ggplot(umap_plot) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, color=level_1), shape=16, size = 0.7, alpha=0.3)+
  scale_color_manual(values=color_vector_level1)+
  theme_classic()+ NoAxes()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+guides(color = guide_legend(override.aes = list(size = 3), title=""))
#+geom_label_repel(aes(x = UMAP_1, y = UMAP_2, label=level_1), size=5, data=label.cent_level_1)

ggplot(umap_plot) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, color=level_2), shape=16, alpha=0.3, size = 0.7)+
  scale_color_manual(values=color_vector_level2)+
  theme_classic()+ NoAxes()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+guides(color = guide_legend(override.aes = list(size = 3), title=""))+
  geom_label_repel(aes(x = UMAP_1, y = UMAP_2, label=level_2), size=5, force=10, data=label.cent_level_2)

ggplot(umap_plot) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, color=level_2), shape=16, alpha=0.3, size = 0.7)+
  scale_color_manual(values=color_vector_level2)+
  theme_classic()+ NoAxes()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+guides(color = guide_legend(override.aes = list(size = 3), title=""))

#adding to Seurat object phase (condition) information and producing DF with information about cell identities and phase (condition) to be used with 
#cell_adundance script
t_list <- count_produceR(t_combined_deci)
t_combined_identity <- t_list[[2]]
t_combined_identity$condition <- "dc"
save(t_combined_deci,file="ut_combined_deci_merged_final.Rdata")

write.table(t_combined_identity,"ut_combined_identity",sep=",",quote=F,row.names = F)
