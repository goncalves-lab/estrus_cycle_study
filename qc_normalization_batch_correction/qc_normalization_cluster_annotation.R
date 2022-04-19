library(garnett)
library(dplyr)
library(Seurat)
library(ggplot2)
library(rlist)
library(SingleCellExperiment)
library(Matrix)
library(scater)
library(RColorBrewer)
library(GISTools)
library(gsubfn)
library(scran)
library(e1071)

###################
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}



############
#Script containing functions to run Seurat analysis on filtered matrices produced by Cellranger
############

#color scheme for UMAPs
color_vector=c(brewer.pal(9,"Reds"),brewer.pal(9,"Blues"),brewer.pal(9,"Greens"),brewer.pal(9,"Purples"))
color_vector=color_vector[c(4,13,22,7,16,25,2,11,20,3,12,21,5,14,23,6,15,24,8,17,26,9,18,27,29,30,31,31,33,34,35)]

###############
#Functions
##############

#detection of doublets (by clustering-identify clusters that have unusually low N using an outlier-based approach (Bach et al 2017) and by simulation)
doubledetectoR=function(mouse_seurat,file_name){
  mouse_seurat=umap_plotteR(mouse_seurat)
  mouse_sce=as.SingleCellExperiment(mouse_seurat)
  # Like 'findMarkers', this function will automatically
  # retrieve cluster assignments from 'colLabels'.
  dbl.out <- doubletCluster(mouse_sce,clusters=mouse_sce$seurat_clusters)
  chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$N, 
                                                type="lower", log=TRUE)]
  set.seed(100)
  # Setting up the parameters for consistency with denoisePCA();
  # this can be changed depending on your feature selection scheme.
  dbl.dens <- doubletCells(mouse_sce, subset.row=mouse_seurat@assays$SCT@var.features, 
                           d=ncol(reducedDim(mouse_sce)))
  mouse_sce$DoubletScore <- log10(dbl.dens+1)
  pdf(file=file_name, width=8.80, height=6.80)
  print(plotUMAP(mouse_sce, colour_by="DoubletScore"))
  print(plotUMAP(mouse_sce, colour_by="seurat_clusters"))
  print(plotColData(mouse_sce, x="seurat_clusters", y="DoubletScore", colour_by="seurat_clusters"))
  dev.off()
  return(chosen.doublet)
}

#normalizes reads using SCTranform method
normalizeR=function(read_matrix_analysis){
  read_matrix_analysis <- SCTransform(read_matrix_analysis, assay = 'RNA', new.assay.name = 'SCT')
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  read_matrix_analysis <- CellCycleScoring(read_matrix_analysis, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay = 'SCT')
  read_matrix_analysis <- SCTransform(read_matrix_analysis, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE,return.only.var.genes = FALSE,
                                      assay = 'RNA', new.assay.name = 'SCT')
  return(read_matrix_analysis)
}

normalizeR_alt=function(read_matrix_analysis){
  read_matrix_analysis <- SCTransform(read_matrix_analysis, assay = 'RNA', new.assay.name = 'SCT')
  return(read_matrix_analysis)
}


#performes PCA and UMAP reduction and clustering step  
reductoR=function(read_matrix_analysis,res_par){
  read_matrix_analysis <- RunPCA(read_matrix_analysis)
  read_matrix_analysis <- RunUMAP(read_matrix_analysis, dims = 1:30)
  #VizDimLoadings(read m_matrix_analysis, dims = 1:2, reduction = "pca")
  read_matrix_analysis <- FindNeighbors(read_matrix_analysis, dims = 1:30)
  read_matrix_analysis <- FindClusters(read_matrix_analysis, resolution = res_par)
  return(read_matrix_analysis) 
}

#filtering based on the median absolute deviation (MAD) from the median value of high mitochondrial RNA content, 
#extreme numbers of counts (count depth) and extreme numbers of genes per barcode.

filtereR <- function(read_matrix_analysis){
  sce <- SingleCellExperiment(assays=list(counts=as.matrix(read_matrix_analysis[["RNA"]]@data)))
  is.mito <- grep("^mt-",rownames(counts(sce)))
  names(is.mito) <- rownames(counts(sce))[grep("^mt-",rownames(counts(sce)))]
  df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
  reasons <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent"))
  filtered <- sce[,!reasons$discard]
  filtered <- as(assay(filtered,"counts"),"dgCMatrix")
  return(filtered)
}

#loading matrices
file_loader <- function(data_dir){
  dirnames <- list.dirs(data_dir)
  for (i in seq(from=2, to=length(dirnames))){
    filenames <- list.files(dirnames[i], full.names=TRUE)
    read_matrix <- Read10X_h5(filename = filenames[1])
    sample_name <- strsplit(filenames[1],"/")
    sample_name <- sample_name[[1]][[(length(sample_name[[1]])-1)]]
    read_matrix_analysis<- CreateSeuratObject(counts = read_matrix, project = sample_name, min.cells = 3, min.features = 200)
    filtered <- filtereR(read_matrix_analysis)
    read_matrix_analysis_filtered <- CreateSeuratObject(counts = filtered, project = sample_name, min.cells = 3, min.features = 200)
    read_matrix_analysis_filtered[["percent.mt"]] <- PercentageFeatureSet(read_matrix_analysis_filtered, pattern = "^mt-")
    read_matrix_analysis_filtered <- normalizeR(read_matrix_analysis_filtered)
    assign(paste0(sample_name,"_filtered"),read_matrix_analysis_filtered,envir = .GlobalEnv)
  }
}

#identifying positive markers of cells in each single cluster
marker_defineR <- function(data_matrix_analysis){
  cluster_markers <- FindAllMarkers(data_matrix_analysis, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  cluster_markers %>% group_by(cluster) %>% top_n(n = 5)
  return(cluster_markers)
}


lotR_clusters <- function(read_matrix_analysis, markers_df){
  clusters <- marker_defineR(read_matrix_analysis)
  clusters_markers <- clusters[clusters$gene%in%markers_df$markers,] 
  list(clusters,clusters_markers)
}

lotR <- function(data_dir,res_par){
  filenames <- list.dirs(data_dir, full.names=TRUE)
  for (i in seq(from=2, to=length(filenames))){
    read_matrix <- Read10X(data.dir =filenames[i])
    sample_name <- strsplit(filenames[i],"/")
    read_matrix_analysis <- CreateSeuratObject(counts = read_matrix, project = sample_name[[1]][lengths(sample_name)], min.cells = 3, min.features = 200)
    read_matrix_analysis[["percent.mt"]] <- PercentageFeatureSet(read_matrix_analysis, pattern = "^mt-")
    filtered=filtereR(read_matrix_analysis)
    read_matrix_analysis_filtered <- CreateSeuratObject(counts = filtered, project = sample_name[[1]][lengths(sample_name)], min.cells = 3, min.features = 200)
    read_matrix_analysis_filtered[["percent.mt"]] <- PercentageFeatureSet(read_matrix_analysis_filtered, pattern = "^mt-")
    read_matrix_analysis_filtered <- normalizeR(read_matrix_analysis_filtered)
    assign(paste0(sample_name[[1]][lengths(sample_name)],"_filtered"),read_matrix_analysis_filtered,envir = .GlobalEnv)
    chosen.doublet_mouse=doubledetectoR(read_matrix_analysis_filtered, paste0(sample_name[[1]][lengths(sample_name)],"_doublets"))
    assign(paste0(sample_name[[1]][lengths(sample_name)],"_doublet"),chosen.doublet_mouse,envir = .GlobalEnv)
  }
  
}

lotR_alt <- function(data_dir,res_par){
  filenames <- list.dirs(data_dir, full.names=TRUE)
  for (i in seq(from=2, to=length(filenames))){
    read_matrix <- Read10X(data.dir =filenames[i])
    sample_name <- strsplit(filenames[i],"/")
    read_matrix_analysis <- CreateSeuratObject(counts = read_matrix, project = sample_name[[1]][lengths(sample_name)], min.cells = 3, min.features = 200)
    read_matrix_analysis[["percent.mt"]] <- PercentageFeatureSet(read_matrix_analysis, pattern = "^mt-")
    filtered=filtereR(read_matrix_analysis)
    read_matrix_analysis_filtered <- CreateSeuratObject(counts = filtered, project = sample_name[[1]][lengths(sample_name)], min.cells = 3, min.features = 200)
    read_matrix_analysis_filtered[["percent.mt"]] <- PercentageFeatureSet(read_matrix_analysis_filtered, pattern = "^mt-")
    read_matrix_analysis_filtered <- normalizeR_alt(read_matrix_analysis_filtered)
    assign(paste0(sample_name[[1]][lengths(sample_name)],"_filtered"),read_matrix_analysis_filtered,envir = .GlobalEnv)
    
  }
}
lotR_ultimate <- function(data_matrix_analysis, new_cluster_id){
  names(new_cluster_id) <- levels(data_matrix_analysis)
  data_matrix_analysis <- RenameIdents(data_matrix_analysis, new_cluster_id)
  print(DimPlot(data_matrix_analysis, reduction = "umap", label = F, pt.size = 0.1,cols = color_vector,repel=T) + NoAxes()+
          theme(legend.text = element_text(size = 16)))
  return(data_matrix_analysis)
}


##################
#annotating the cell clusters
subseteR_part_1 <- function(seurat_obj, identities, res){
  tissue_combined_sub <- subset(x = seurat_obj, idents = identities)
  #tissue_combined_sub <- FindVariableFeatures(tissue_combined_sub, selection.method = "vst", nfeatures = 2000)
  tissue_combined_sub <- reductoR(tissue_combined_sub, res)
  print(DimPlot(tissue_combined_sub, reduction = "umap", label = T, pt.size = 0.7,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))
  list[subset_clusters, subset_subclusters] <- lotR_clusters(tissue_combined_sub, markers_df)
  list(tissue_combined_sub, subset_clusters, subset_subclusters)
}



subseteR_part_2 <- function(tissue_combined_sub, new_cluster_id){
  tissue_combined_sub <- lotR_ultimate(tissue_combined_sub, new_cluster_id)
  cell_sub <- tissue_combined_sub$seurat_clusters
  cell_sub <- data.frame("cells" = names(cell_sub), "type" = unname(cell_sub))
  return(cell_sub)
}


#Garnett
garnett_annotatoR_part1 <- function(t_sub, t_combined_cell_metadata, marker_file){
  t_combined_mat <- t_sub[["RNA"]]@counts
  t_combined_genes <- data.frame("gene_short_name" = rownames(t_combined_mat))
  rownames(t_combined_genes) <- rownames(t_combined_mat)
  t_combined_cell_metadata_sub <- t_combined_cell_metadata[t_combined_cell_metadata$cell_names%in%colnames(t_sub),]
  t_combined_cds_sub <- new_cell_data_set(as(t_combined_mat, "dgCMatrix"),
                                          cell_metadata = t_combined_cell_metadata_sub,
                                          gene_metadata = t_combined_genes)
  
  
  marker_check <- check_markers(t_combined_cds_sub, marker_file,
                                db="none",
                                cds_gene_id_type = "SYMBOL",
                                marker_file_gene_id_type = "SYMBOL")
  
  print(plot_markers(marker_check))
  return(t_combined_cds_sub)
}


garnett_annotatoR_part2 <- function(t_combined_cds_sub_tr, t_combined_cds_sub, marker_file, t_sub){
  t_classifier_sub <- train_cell_classifier(cds = t_combined_cds_sub_tr,
                                            marker_file = marker_file,
                                            db="none",
                                            cds_gene_id_type = "SYMBOL",
                                            marker_file_gene_id_type = "SYMBOL")
  t_combined_cds_sub <- classify_cells(t_combined_cds_sub, t_classifier_sub,
                                       db = "none",
                                       cluster_extend = TRUE,
                                       cds_gene_id_type = "SYMBOL")
  cells_clas_sub <- as.data.frame(pData(t_combined_cds_sub))
  t_sub$labels <- cells_clas_sub$cell_type
  print(DimPlot(t_sub, reduction = "umap", label = F, pt.size = 0.5, group.by="labels",repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))
  return(t_sub)
}

garnett_extender <- function(t_combined){
  garnett_labels <- t_combined$labels
  seurat_clusters <- t_combined$seurat_clusters
  
  garnett_annotation_ex <- data.frame("garnett" = garnett_labels, "seurat_clusters" = seurat_clusters)
  garnett_annotation_ex$seurat_clusters <- as.character(garnett_annotation_ex$seurat_clusters)
  garnett_annotation_ex$seurat_clusters <- as.factor(garnett_annotation_ex$seurat_clusters) 
  garnett_cluster <- data.frame()
  for(i in seq(from = 1, to = length(levels(garnett_annotation_ex$seurat_clusters)))){
    
    garnett_cluster_temp <- garnett_annotation_ex[grep(paste0("\\b",levels(garnett_annotation_ex$seurat_clusters),"\\b")[i], garnett_annotation_ex$seurat_clusters),]
    ex <- data.frame("cell_type" = names(table(garnett_cluster_temp$garnett)/sum(table(garnett_cluster_temp$garnett))*100), 
                     "percentage" = unname(table(garnett_cluster_temp$garnett)/sum(table(garnett_cluster_temp$garnett))*100))
    ex <- ex[ex$cell_type != "Unknown",]
    if (T %in% (ex$percentage.Freq>20)){
      garnett_temp <- ex[ex$percentage.Freq==max(ex$percentage.Freq),][1,1]
    } else {garnett_temp = "Unknown"}
    garnett_temp_df <- data.frame ("cluster" = garnett_cluster_temp[1,2], "cell"=garnett_temp)
    garnett_cluster <- rbind(garnett_cluster, garnett_temp_df)
  }
  
  garnett_cluster$cluster <- as.numeric(levels(garnett_cluster$cluster))
  garnett_cluster <- garnett_cluster [order(garnett_cluster$cluster),]
  return(garnett_cluster)
}

#SVM
svm_subseteR <- function(t_combined){
  Idents(t_combined) <- t_combined$match
  t_match <- subset(t_combined, idents = "match") 
  t_match_ge <- t_match[["SCT"]]@scale.data 
  t_match_ge <- t(t_match_ge)  
  y <- as.factor(t_match$garnett_cluster)
  t_match_ge <- data.frame(x=t_match_ge, y=y)
  
  t_non_match <- subset(t_combined, idents = "non-match") 
  t_non_match_ge <- t_non_match[["SCT"]]@scale.data 
  t_non_match_ge  <- t(t_non_match_ge)  
  t_non_match_ge <- data.frame(x=t_non_match_ge)
  
  t_match_ge <- t_match_ge[, colnames(t_match_ge)%in%colnames(t_non_match_ge)]
  t_match_ge$y <- y
  list(t_non_match_ge, t_match_ge)}

svm_returneR <- function(t_match_ge, t_non_match_ge, prob, subset_value){
  set.seed(1)
  t_match_train <- data.frame()
  t_match_ge$y <- as.factor(t_match_ge$y)
  for (i in levels(t_match_ge$y)){
    subset_size <- round(table(t_match_ge$y)[i]*subset_value)
    t_match_ge_subset <- t_match_ge[grep(i, t_match_ge$y),]
    t_match_ge_train_subset <- t_match_ge_subset[sample.int(subset_size),c(sample.int(10000),ncol(t_match_ge))]
    t_match_train <- rbind(t_match_train, t_match_ge_train_subset)
  }
  
  t_match_test <- t_match_ge[!rownames(t_match_ge)%in%rownames(t_match_train), colnames(t_match_ge)%in%colnames(t_match_train)]
  t_match_test_final <- data.frame()
  for (i in levels(t_match_test$y)){
    subset_size <- round(table(t_match_test$y)[i]*0.1)
    t_match_test_subset <- t_match_test[grep(i, t_match_test$y),]
    t_match_test_subset <- t_match_test_subset[sample.int(subset_size),c(sample.int(10000),10001)]
    t_match_test_final <- rbind(t_match_test_final, t_match_test_subset)
  }
  
  #tune.out <- tune(svm, y~., data=t_match_train, kernel="linear", ranges = list(cost=c(0.1,10,1000)))
  out <- svm(y~., data = t_match_train, kernel ="linear", cost=10, scale = F, probability = T)
  pred.te <- predict (out, newdata = t_non_match_ge, probability = T)
  pred.te_probabilites <- attr(pred.te,"probabilities")
  
  pred.te_identity <- c()
  for (i in seq(from=1, to=nrow(pred.te_probabilites))){
    if (T %in% (pred.te_probabilites[i,] > prob)){
      pred.te_identity[i] <- names(pred.te_probabilites[i,][pred.te_probabilites[i,]==max(pred.te_probabilites[i,])])
    } else{pred.te_identity[i] <- "unknown"}
  }
  names(pred.te_identity) <- row.names(pred.te_probabilites)
  return(pred.te_identity)}

#UMAP generation
count_produceR=function(tissue_combined){
  tissue_combined_counts=as.matrix(tissue_combined[["SCT"]]@counts)
  tissue_combined_counts=tissue_combined_counts[!!rowSums(abs(tissue_combined_counts)),]
  tissue_combined_identity=data.frame(names(tissue_combined$level_2),unname(tissue_combined$level_2))
  colnames(tissue_combined_identity)=c("cells_names","identity")
  tissue_combined_identity$batch=tissue_combined$orig.ident
  rownames(tissue_combined_identity)=tissue_combined_identity$cells_names
  tissue_list=list("counts"=tissue_combined_counts,"identity"=tissue_combined_identity)
  return(tissue_list)
}


UMAP_plotteR=function(tissue_combined,file_name){
  pdf(file=paste0(file_name ,"_uncorrected_combined"), width=8.80, height=6.80)
  print(DimPlot(tissue_combined, reduction = "umap", label = F, pt.size = 0.7,cols = color_vector,repel=T) + NoAxes()+theme(legend.text = element_text(size = 16)))
  dev.off()
  pdf(file=paste0(file_name ,"_uncorrected_sample"), width=8.80, height=6.80)
  print(DimPlot(tissue_combined, reduction = "umap",group.by ="orig.ident",cols = color_vector_tr,label = F, pt.size = 0.1) + NoAxes()+theme(legend.text = element_text(size = 16)))
  dev.off()
  pdf(file=paste0(file_name ,"_uncorrected_sample_separated"), width=30, height=6.80)
  print(DimPlot(tissue_combined, reduction = "umap",split.by="orig.ident",group.by ="orig.ident",cols=color_vector)+ NoAxes()+theme(legend.text = element_text(size = 16)))
  dev.off()
} 
