#this script requires seurat object (estrus_ut_processing.R), and list of cycle-associated genes (MI.R)

library(tidyverse)
library(caret)
library(MASS)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(uwot)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(biomaRt)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
library(dplyr)
library(mda)
#library(dyno)


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}

subseteR<- function(seurat_obj, gene_subset){
  cell_counts <-  as.matrix(seurat_obj[["SCT"]]@counts)
  cell_counts <- cell_counts[rownames(cell_counts)%in%gene_subset,]
  #sce <- SingleCellExperiment(list(logcounts = cell_counts))
  return(cell_counts)
}


ut_combined <- loadRData("ut_combined.Rdata")


ut_combined_f <- subset(ut_combined, idents="F")

dge_f=read.csv("ut_MI_genes_mouse")
dge_f=dge_f[order(dge_f$p_val_adj),]
dge_f=dge_f[1:2000,]

sce_f_pt_dge <- subseteR(ut_combined_f, dge_f$genes)
sce_f_pt_dge <- t(sce_f_pt_dge)
sce_f_pt_dge <- as.data.frame(sce_f_pt_dge)
metadata <- data.frame(group = ut_combined_f$group)

sce_f_pt_dge <- merge(sce_f_pt_dge, metadata, by=0)
set.seed(123)
training.samples <- sce_f_pt_dge$group %>%
  createDataPartition(p =0.1, list = FALSE)

train.data <- sce_f_pt_dge[training.samples, ]
train.data <- train.data[,2:2000]
sce_f_pt_dge <- sce_f_pt_dge[,2:2000]
sce_f_pt_dge$group <- as.factor(sce_f_pt_dge$group)
train.data$group <- as.factor(train.data$group)
# Fit the model
estrus_model <- fda(group~., data = sce_f_pt_dge)

plot(estrus_model, data=sce_f_pt_dge)

object <- estrus_model                        # generic/method
ff <- terms(object)
attr(ff, "intercept") <- 0
m <- model.frame(ff, sce_f_pt_dge)
x <- model.matrix(ff, m)
vars <- predict(object, x, type = "var")
g <- model.extract(m, "response")
means <- object$means
g <- as.factor(g)
cc <- as.numeric(g)
np <- seq(levels(g))
coords = c(1, 2)
if (!is.matrix(coords)) {
  coords <- matrix(coords, length(coords), length(coords))
  tt <- lower.tri(coords)
  coords <- cbind(t(coords)[tt], coords[tt])
}

coord.pair <- coords[1, ]
fda_matrix <- as.data.frame(vars[, coord.pair])
fda_matrix$group <- g

color_vec <- c(brewer.pal(9, "Greys")[c(3, 5, 6)], "black")
print(ggplot(fda_matrix, aes(x=V1, y=V2, color= group)) + geom_point(size=1)+
        theme_bw(base_size = 14) + 
        theme(legend.position = "right")+
        theme(legend.text=element_text(size=28),legend.title = element_blank())+labs(x="LD1",y="LD2")+
        theme(axis.text=element_text(size=34), axis.title=element_text(size=36))+
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())+
        guides(fill = guide_legend(override.aes = list(size = 3), title=""))+
        scale_color_manual(values=c(color_vec)))

rd1 <- as.matrix(fda_matrix[,1:2])
colnames(rd1) <- c('LD1', 'LD2')
sce <- SingleCellExperiment(list(logcounts = as.matrix(t(sce_f_pt_dge))))
reducedDims(sce) <- list(LDA=rd1)
pca_pt <- as.data.frame(reducedDims(sce)$LDA)
colnames(pca_pt) <- c("LD1", "LD2")
cl2 <- kmeans(rd1, centers = 2)$cluster
colData(sce)$kmeans <- cl2
print(plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1))
sce <- slingshot(sce, clusterLabels = 'kmeans', reducedDim = 'LDA',start.clus=1)


