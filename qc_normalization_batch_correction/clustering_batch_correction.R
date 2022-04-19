#################################
# this script uses seurat objects generated in scripts estrus_ov_processing.R, 18m_ov_processing.R, estrus_od_processing.R, 18m_od_processing.R, estrus_ut_processing.R, 18m_ut_processing.R,
#estrus_ce_processing.R, 18m_ce_processing.R, estrus_va_processing.R, 18m_va_processing.R, estrus_sp_processing.R, 18m_sp_processing.R
#################################


library(RColorBrewer)
library(GISTools)
library(SeuratWrappers)
library("harmony")
library(Seurat)

options(future.globals.maxSize= 2621440000)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

#color scheme for UMAPs
color_vector=c(brewer.pal(9,"Reds"),brewer.pal(9,"Blues"),brewer.pal(9,"Greens"),brewer.pal(9,"Purples"))
color_vector=color_vector[c(4,13,22,7,16,25,2,11,20,3,12,21,5,14,23,6,15,24,8,17,26,9,18,27,29,30,31,31,33,34,35)]
color_vector_tr=add.alpha(color_vector,0.7)


#plot UMAPs of Seurat combined object (no batch correction)
umap_plotteR=function(ds_combined){
  ds_combined <- FindVariableFeatures(ds_combined, selection.method = "vst", nfeatures = 2000)
  ds_combined <- RunPCA(ds_combined, features = VariableFeatures(object = ds_combined))
  ds_combined <- FindNeighbors(ds_combined, dims = 1:30)
  ds_combined <- FindClusters(ds_combined, resolution = 0.5)
  ds_combined <- RunUMAP(ds_combined, dims = 1:30)
  print(DimPlot(ds_combined, reduction = "umap", label = TRUE, pt.size = 0.5,cols=color_vector))
  print(DimPlot(ds_combined, reduction = "umap",group.by ="orig.ident",cols = color_vector_tr, label=T))
  return(ds_combined)
}


selectoR <- function(data_dir, young){
  setwd(data_dir)
  rep_list <-list()
  file_names <- list.files(data_dir)
  file_names <- file_names[!file_names=="sp_combined_new.Rdata"]
  file_names <- file_names[!file_names=="sp_combined_18m.Rdata"]
  file_names <- file_names[grep("combined", file_names)]
  for (i in seq(from=1, to=length(file_names))){
    t_combined <- loadRData(file_names[i])
    DefaultAssay(t_combined) <- "RNA"
    t_combined[["SCT"]] <- NULL
    t_combined[["umap"]]<- NULL
    Idents(t_combined) <- t_combined$level_2
    t_combined <- subset(t_combined, idents="F")
    if (young == "young"){
      Idents(t_combined) <- t_combined$estrus_phase
      t_combined <- subset(t_combined, idents="diestrus")
    }
    t_list <- SplitObject(t_combined, split.by = "orig.ident")
    rep_list <- c(rep_list, t_list)
    
  }  
  return(rep_list)
}


mergeR <- function(data_dir, young){
  setwd(data_dir)
  rep_list <-list()
  file_names <- list.files(data_dir)
  file_names <- file_names[!file_names=="sp_combined_new.Rdata"]
  file_names <- file_names[!file_names=="sp_combined_18m.Rdata"]
  file_names <- file_names[grep("combined", file_names)]
  for (i in seq(from=1, to=length(file_names))){
    t_combined <- loadRData(file_names[i])
    Idents(t_combined) <- t_combined$level_2
    t_combined <- subset(t_combined, idents="F")
    if (young == "young"){
      Idents(t_combined) <- t_combined$estrus_phase
      t_combined <- subset(t_combined, idents="diestrus")
    }
    rep_list <- c(t_combined, rep_list)
  }  
  return(rep_list)
}

#Harmony integration-function; plot UMAPs of Seurat combined object (Harmony corrected)
harmonizeR=function(ds_combined){
  ds_combined <- SCTransform(ds_combined)
  ds_combined <- FindVariableFeatures(ds_combined, selection.method = "vst", nfeatures = 2000)
  ds_combined <- RunPCA(ds_combined, verbose = FALSE)
  ds_combined_har <- RunHarmony(ds_combined, group.by.vars = "orig.ident", assay="SCT")
  ds_combined_har <- RunUMAP(ds_combined_har, reduction = "harmony", dims = 1:30)
  ds_combined_har <- FindNeighbors(ds_combined_har, reduction = "harmony", dims = 1:30)
  return(ds_combined_har)}


rep_list_young <- selectoR("seurat_objects/", "young")
rep_list_old <- selectoR("seurat_objects/", "old")

rep_list_young <-mergeR("seurat_objects/", "young")
rep_list_old <-mergeR("seurat_objects/", "old")

for (i in seq(length(rep_list_old))){
  rep_list_old[[i]] <- RenameCells(rep_list_old[[i]], i)
}

t_combined <- merge(rep_list_young[[1]], y=c(rep_list_young[[2]], rep_list_young[[3]], rep_list_young[[4]], rep_list_young[[5]],
                                             rep_list_old[[1]], rep_list_old[[2]], rep_list_old[[3]], rep_list_old[[4]], rep_list_old[[5]]))

t_combined <- ScaleData(t_combined)
t_combined <- umap_plotteR(t_combined)

t_combined$estrus_phase[is.na(t_combined$estrus_phase)] <- "18_m"

t_combined <- FindVariableFeatures(t_combined, selection.method = "vst", nfeatures = 2000)
t_combined <- RunPCA(t_combined, features = VariableFeatures(object = t_combined))
t_combined <- FindNeighbors(t_combined, dims = 1:30)
t_combined <- FindClusters(t_combined, resolution = 0.5)
t_combined <- RunUMAP(t_combined, dims = 1:30)

t_combined$tissue <- substr(t_combined$orig.ident, nchar(t_combined$orig.ident)-1, nchar(t_combined$orig.ident))
print(DimPlot(t_combined, reduction = "umap",group.by ="estrus_phase",cols=color_vector_tr, pt.size=0.1))
print(DimPlot(t_combined, reduction = "umap",group.by ="tissue", pt.size = 0.1, cols=color_vector_tr))


t_har <- harmonizeR(t_combined)

print(DimPlot(t_har, reduction = "umap",group.by ="estrus_phase",cols=color_vector_tr, pt.size=0.1))
t_har$tissue <- substr(t_har$orig.ident, nchar(t_har$orig.ident)-1, nchar(t_har$orig.ident))
print(DimPlot(t_har, reduction = "umap",group.by ="tissue", pt.size = 0.1, cols=color_vector_tr))
print(DimPlot(t_har, reduction = "umap",group.by ="tissue", pt.size = 0.1, cols=color_vector_tr, split.by="tissue"))


metadata_F <- merge(metadata_F, auc_df_final, by=0)
metadata_F$label <- paste0(as.character(lapply(str_split(metadata_F$cell_names,"-"), function(l) l[[1]])), "_", metadata_F$sample)
t_har$label <- paste0(as.character(lapply(str_split(colnames(t_har),"-"), function(l) l[[1]])), "_", t_har$orig.ident)

t_har <- FindClusters(t_har, graph.name = "SCT_nn")

umap_df <- data.frame("UMAP_1"=t_har[["umap"]]@cell.embeddings[,1], "UMAP_2"=t_har[["umap"]]@cell.embeddings[,2], "phase"=t_har$estrus_phase, "label"=t_har$label)
umap_df_18m <- umap_df[umap_df$phase=="18_m",]
umap_df_18m$label <- substr(umap_df_18m$label, 3, nchar(umap_df_18m$label)) 
umap_df <- umap_df[umap_df$phase != "18_m",]
umap_df <- rbind(umap_df, umap_df_18m)

umap_df <- merge(umap_df, metadata_F, by="label")

print(DimPlot(t_har, reduction = "umap",group.by ="estrus_phase",cols=color_vector_tr, pt.size=0.1))
print(DimPlot(t_har, reduction = "umap",group.by ="tissue", pt.size = 0.1, cols=color_vector_tr))
ggplot(umap_df) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, color=INFLAMMATION), shape=20, size = 0.5)+
  theme_classic()+ NoAxes()+ scale_color_viridis(option = "D")+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+facet_wrap(~phase.x)

ggplot(umap_df) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, color=ECM_ORGANISATION), shape=20, size = 0.5)+
  theme_classic()+ NoAxes()+ scale_color_viridis(option = "D")+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+facet_wrap(~phase.x)
print(DimPlot(t_har, reduction = "umap",group.by ="seurat_clusters", pt.size = 0.1, cols=color_vector_tr, label=T))

cluster2.markers <- FindMarkers(t_har, ident.1 = 2, min.pct = 0.25)
gene_inflammation <- cluster2.markers[cluster2.markers$genes%in%geneSets_final[[6]]@geneIds,]
p_value_inflammation=phyper(23-1, 322, 20000-322, 197, lower.tail = FALSE, log.p = FALSE)
##########################################

umap_df <- data.frame("UMAP_1"=t_har[["umap"]]@cell.embeddings[,1], "UMAP_2"=t_har[["umap"]]@cell.embeddings[,2], 
                      auc=metadata_har$auc, "phase"=t_har$estrus_phase, "sample"=t_har$orig.ident, "cluster"=t_har$seurat_clusters)

umap_df[,7] <- 1:nrow(umap_df)
colnames(umap_df)[7] <- "cluster_all"
for (i in seq(nrow(umap_df))){
  if (umap_df[i,6]==2){umap_df[i,7]="cluster 2"} else {umap_df[i,7]="remaining clusters"}
}

umap_df$tissue <- substr(umap_df$sample, nchar(umap_df$sample)-1, nchar(umap_df$sample))

library(vioplot)

pdf("violin_plot_tissue.pdf", height=6, width=7)
for(i in unique(umap_df$tissue)){
  umap_df_i <- umap_df[umap_df$tissue==i,]
  umap_df_1 <- umap_df_i[umap_df_i$cluster_all=="remaining clusters",]
  umap_df_2 <- umap_df_i[umap_df_i$cluster_all=="cluster 2",]
  vioplot(auc~phase, data=umap_df_1, col = "grey", plotCentre = "line", side = "right", xlab = "", ylab = "AUCell score",  axes=F, main=i)
  vioplot(auc~phase, data=umap_df_2, col = "white", plotCentre = "line", side = "left", add = T, bty="n", frame.plot=FALSE, axes=F)
  legend("topright", inset=c(0,-0.1) ,fill = c("grey", "white"), legend = c("other", "cluster_2"), title = "", bty = "n")
}

dev.off()
umap_df_1 <- umap_df[umap_df$cluster_all=="remaining clusters",]
umap_df_2 <- umap_df[umap_df$cluster_all=="cluster 2",]
vioplot(auc~phase, data=umap_df_1, col = "grey", plotCentre = "line", side = "right", xlab = "", ylab = "AUCell score",  axes=F)
vioplot(auc~phase, data=umap_df_2, col = "white", plotCentre = "line", side = "left", add = T, bty="n", frame.plot=FALSE, axes=F)
legend("topright", inset=c(0,-0.1) ,fill = c("grey", "white"), legend = c("other", "cluster 2"), title = "", bty = "n")

umap_df_1_count <- umap_df_1 %>% group_by(sample) %>% dplyr::count(sample)
umap_df_2_count <- umap_df_2 %>% group_by(sample) %>% dplyr::count(sample)

umap_df_count <- merge(umap_df_1_count, umap_df_2_count, by="sample")
umap_df_count$perc <- umap_df_count$n.y/(umap_df_count$n.x+umap_df_count$n.y)
umap_df_count$phase <- c(rep("18m",16), rep("diestrus",19))
umap_df_count$tissue <- substr(umap_df_count$sample, nchar(umap_df_count$sample)-1, nchar(umap_df_count$sample))

ggplot(umap_df_count, aes(x=tissue, y=perc, fill=phase)) + 
  geom_boxplot()+theme_classic()+  scale_fill_manual(values = c("white", "grey"))+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=20))+ylab("% of cells in cluster 2")

