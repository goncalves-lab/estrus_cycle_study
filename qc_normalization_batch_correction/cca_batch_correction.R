#################################
# this script uses seurat objects generated in scripts estrus_ut_processing.R and pregnant_ut_processing.R 
#################################

source("qc_normalization_cluster_annotation.R")


# takes a Seurat object and outputs a sparse matrix 
# of smoothed expression values
get_smooth_expression <- function( seu ) {
  cat("setting default assay to 'RNA'...\n")
  DefaultAssay(seu) <- "RNA"
  cat("running SCTransform...\n")
  seu = SCTransform(seu, vars.to.regress = "percent.mt", 
                    method = "glmGamPoi", verbose = TRUE)  
  
  seu = RunPCA(seu, assay="SCT", npcs = 50)
  
  if( packageVersion("Seurat") != "4.0.3" ) {
    seu = FindNeighbors(seu, k.param = 50 )
  } else {
    seu = FindNeighbors(seu, k.param = 50, reduction = "pca", 
                        assay="SCT", dims = 1:50 )
  }
  on = names(seu@graphs)[grepl("_nn",names(seu@graphs))]
  
  cat("using", on, "...\n")
  nb = as(seu@graphs[[on]],"Matrix")
  raw = seu@assays$RNA@counts
  
  # the dot product gives the smoothed expression
  cat("dot product...\n")
  nb = Matrix::t(nb) 
  dp = raw %*% nb
  
  cat("dividing by total reads...\n")
  ts = as.numeric(colSums(dp))
  dp = dp %*% Matrix::Diagonal( x=1/ts )
  dp = dp*1000
  cat("collecting garbage...\n")
  gc()
  
  return( dp )
}


t_combined <- loadRData("ut_combined.Rdata")
decidua <- loadRData("ut_combined_deci_merged_final.Rdata")

decidua$estrus_phase <- "decidua"

Idents(t_combined) <- t_combined$estrus_phase
t_combined <- subset(t_combined, idents=c("metestrus"))
decidua <- decidua[, sample(colnames(decidua), size = ncol(t_combined), replace=F)]

# select only estrus
t_list <- SplitObject(t_combined, split.by = "orig.ident")
decidua_list <- SplitObject(decidua, split.by = "orig.ident") 
seul <- c(t_list, decidua_list)

seul <- lapply(X = seul, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x, method = "glmGamPoi")
})


features <- SelectIntegrationFeatures(object.list = seul, nfeatures = 3000)
seul <- PrepSCTIntegration(object.list = seul, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seul, normalization.method = "SCT", anchor.features = features)

combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

stab = get_smooth_expression(combined)
save(combined, stab, file="ut_metestrus_phase_decidua_integrated.RData")

###########################

load("ut_metestrus_phase_decidua_integrated.RData")
genes = c("Cebpb","Alpl","Bmp2")
genes = genes[genes %in% rownames(stab)]
exp = c()
for( gene in genes ) {
  exp = cbind(exp,stab[gene,])
}
colnames(exp) = paste( "gene_", genes, sep="" )
tab = cbind( combined@meta.data, exp, combined@reductions$umap@cell.embeddings )

pdf( "decidua_mestrus_phase_integrated_UMAP.pdf", width=5, height=5 )
ggplot(tab,aes(UMAP_1,UMAP_2,color=level_1)) + geom_point(size=0.1) + theme_void()
ggplot(tab,aes(UMAP_1,UMAP_2,color=level_2)) + geom_point(size=0.1) + theme_void()
ggplot(tab,aes(UMAP_1,UMAP_2,color=estrus_phase)) + geom_point(size=0.1) + theme_void()+scale_color_manual(values = c(brewer.pal(9, "Purples")[3], brewer.pal(9, "Greys")[5]))
ggplot(tab,aes(UMAP_1,UMAP_2,color=orig.ident)) + geom_point(size=0.1) + theme_void()
dev.off()

# separate color scale for estrus and decidua
a = tab
pdf("decidua_estrus_phase_integrated_smoothed_expression_pcpc.pdf", width=5, height=5 )
for( gene in genes ) {
  gene = paste("gene_",gene,sep="")
  a = a[order(a[,gene]),]
  for( phase in c("estrus","decidua") ) {
    p = ggplot( a[a$estrus_phase == phase & a$level_1 == "SC",], aes(UMAP_1,UMAP_2,color=get(gene)) ) + geom_point(size=0.1) + theme_void() +  scale_color_viridis() + ggtitle(paste(gene,phase)) 
    print(p)
  }
}
dev.off()

# same color scale for estrus and decidua
a = tab
pdf("decidua_estrus_phase_integrated_smoothed_expression_pcpc_wrap.pdf", width=7, height=5)
for( gene in genes ) {
  gene = paste("gene_",gene,sep="")
  a = a[order(a[,gene]),]
  p = ggplot( a[a$level_1 == "SC",], aes(UMAP_1,UMAP_2,color=get(gene)) ) + geom_point(size=0.1) + theme_void() +  scale_color_viridis(option = "G") + ggtitle(paste(gene,phase)) + facet_wrap(~ estrus_phase)
  print(p)
}
dev.off()