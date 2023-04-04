#######################
#this script calculates cluster coefficient. It uses files produced in deconvolution script
########################

library("RColorBrewer")
library(Seurat)
library(data.table)
library(stringr)
library(ggplot2)
library(AUCell)
library(GSEABase)
library(rjson)
library(ggforce)
library(png)
library(dplyr)
library(viridis)
library(lme4)
library(igraph)
library(ggdist)
library(ggridges)
library(MuMIn)

load("/omics/odcf/analysis/OE0538_projects/DO-0002/mmus/aucell/geneSets_AUCell.Rdata")
geneSets_final <- geneSets_final[12]
AUCElleR <- function(counts, cell){
  cells_rankings_m <- AUCell_buildRankings(t(counts), nCores=1, plotStats=TRUE)
  auc_df_final_m <- data.frame()
  for(i in geneSets_final){
    print(i)
    cells_AUC <- AUCell_calcAUC(i, cells_rankings_m, aucMaxRank = ceiling(1 * nrow(cells_rankings_m)))
    auc_df <- as.data.frame(getAUC(cells_AUC))
    #auc_df <- t(auc_df)
    auc_df_final_m <- rbind(auc_df, auc_df_final_m)
  }
  auc_df_final_m <- as.data.frame(t(auc_df_final_m))
  auc_df_final_m$keys <- rownames(auc_df_final_m)
  #colnames(auc_df_final_m) <- c(paste0("INFLAMMATION_", cell), "keys")
  return(auc_df_final_m)
  
}


cluster_index <- function(spatial_loc, wd1, prop_loc, sc, wd2, sample_name, threshold_fixed){
  spatial <- Load10X_Spatial(spatial_loc)
  setwd(wd1)
  print("a")
  prop <- read.csv(prop_loc, row.names=1,sep="\t")
  coord <- GetTissueCoordinates(spatial)
  print("b")
  distances <- as.data.frame(as.matrix(dist(coord)))
  
  distances[distances == 0] <- 100
  distances[distances < 10] <- 1
  distances[distances > 10] <- 0
  counts_sc <- read.csv(sc, row.names=1,sep="\t")
  auc_sc <- AUCElleR(counts_sc, "SC")
  dev.off()
  g1c <- graph_from_adjacency_matrix( as.matrix(distances), mode="undirected")
  components <- igraph::clusters(g1c, mode="weak")
  cluster_id <- components[[1]]
  cluster_id_sc <- cluster_id[names(cluster_id)%in%rownames(auc_sc)]
  cluster_name <- as.numeric(names(table(cluster_id_sc))[unname(table(cluster_id_sc))>2])
  cluster_id <- cluster_id[cluster_id%in%cluster_name]
  # Define two datasets and store them to disk
  setwd(wd2)
  coords = read.csv("outs/spatial/tissue_positions_list.csv", header=F, row.names = 1)
  colnames(coords) = c("tissue", "row", "col", "imagerow", "imagecol")
  coords = coords[coords$tissue==1,]
  coords$imagerow = -coords$imagerow
  coords$keys = rownames(coords)
  
  scalef = fromJSON(file = "outs/spatial/scalefactors_json.json")$tissue_lowres_scalef
  spot = fromJSON(file = "outs/spatial/scalefactors_json.json")$spot_diameter_fullres / 2
  coords$scaled_row = coords$imagerow * scalef
  coords$scaled_col = coords$imagecol * scalef
  
  img <- readPNG("outs/spatial/tissue_lowres_image.png")
  selected_points <- coords[0, ]
  auc <- data.frame("keys"=rownames(auc_sc), "auc"=auc_sc$INFLAMMATION)
  coords_ab <- left_join(coords, auc, by="keys")
  setwd("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/ut/3m/")
  pdf(paste0(sample_name, "_INFLAMMATION.pdf"), height=6, width=6)
  print(ggplot(coords_ab, aes(x = scaled_col, y = scaled_row, r = spot * scalef, key=keys)) +
          geom_point(aes(color=auc),size = 2)+
          scale_color_gradient(low = brewer.pal(9, "YlOrRd")[1], high = brewer.pal(9, "YlOrRd")[9], na.value = "light grey")+
          coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
          theme_void())
  cl_c_all <- data.frame()
  for (j in unique(cluster_id)){
    cluster_names <- cluster_id[cluster_id==j]
    auc_ov <- auc[auc$keys%in%names(cluster_names),]
    distances_c <- distances[,colnames(distances)%in%names(cluster_names)]
    distances_c <- distances_c[rownames(distances_c)%in%names(cluster_names),]
    distances_sc <- distances_c[,colnames(distances_c)%in%auc_ov$keys]
    distances_sc <- distances_sc[rownames(distances_sc)%in%auc_ov$keys,]
    
    g1_sc <- graph_from_adjacency_matrix( as.matrix(distances_sc), mode="undirected")
    components_sc <- igraph::clusters(g1_sc, mode="weak")
    cluster_id_sc <- components_sc[[1]]
    big_cluster <- as.numeric(names(sort(table(cluster_id_sc),decreasing=TRUE)[1]))
    cluster_id_sc <- cluster_id_sc[cluster_id_sc%in%big_cluster]
    if (length(cluster_id_sc)>50){ #50 for ut, 20 for ov
      auc_ov <- auc_ov[auc_ov$keys%in%names(cluster_id_sc),]
      auc_ov <- auc_ov[order(auc_ov$auc),]
      threshold <- auc_ov$auc[round(nrow(auc_ov)*threshold_fixed)]
      
      auc_ov$auc[auc_ov$auc >= threshold] <- 1
      auc_ov$auc[auc_ov$auc < threshold] <- 0
      
      coords_ov <- left_join(coords, auc_ov, by="keys")
      coords_bin <- coords_ov
      coords_bin$auc <- as.character(coords_bin$auc)
      print(ggplot(coords_bin, aes(x = scaled_col, y = scaled_row, r = spot * scalef, key=keys)) +
              geom_point(aes(color=auc),size = 2)+scale_color_manual(values=c("black", brewer.pal(9, "YlOrRd")[7]))+
              coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
              theme_void())
      
      auc_in <- auc_ov[auc_ov$auc==1,]
      distances_auc <- distances_c[,colnames(distances_c)%in%auc_in$keys]
      distances_auc <- distances_auc[rownames(distances_auc)%in%auc_in$keys,]
      
      g1 <- graph_from_adjacency_matrix( as.matrix(distances_auc), mode="undirected")
      print(plot(g1, vertex.label=NA, vertex.size=3))
      a <- transitivity(g1, type="average",isolates="zero")
      cl_c <- data.frame(cluster_size=a, sample=sample_name,ovary=j)
      cl_c_all <- rbind(cl_c, cl_c_all) 
      
    }
    
     
  } 
  dev.off()
  return(cl_c_all)
  
}
##############
#ov
#young
sample_list <- list.dirs("./")
sample_list <- gsub(".//", "", sample_list)
sample_list <- sample_list[4:length(sample_list)]

cl_c_final_y <- data.frame()
for (i in sample_list){
  cl_c <- cluster_index(paste0("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/ha/3M/",i,"_manual_repmerged/outs/"),
                             "/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/ha/deconvoluted/3M/",
                             paste0(i,"/", i, "_cell_proportions.csv"),
                             paste0(i,"/", i, "_SC_gene_expression_thresholded.csv"),
                             paste0("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/spatial_tr/spaceranger_out/", i,"_manual_repmerged/"),
                             i, 0.75)
  cl_c_final_y <- rbind(cl_c, cl_c_final_y)
}

cl_c_final_y$group <- "young"

#old
sample_list <- list.dirs("./")
sample_list <- gsub(".//", "", sample_list)
sample_list <- sample_list[4:length(sample_list)]

cl_c_final <- data.frame()
for (i in sample_list){
  cl_c <- cluster_index(paste0("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/spatial_tr/spaceranger_out/", i,"_manual_repmerged/outs"),
                        "/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/ha/deconvoluted/18M/",
                        paste0(i,"/", i, "_cell_proportions.csv"),
                        paste0(i,"/", i, "_SC_gene_expression_thresholded.csv"),
                        paste0("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/follicle_project/spatial_tr/spaceranger_out/", i,"_manual_repmerged/"),
                        i, 0.75)
  cl_c_final <- rbind(cl_c, cl_c_final)
}
cl_c_final$group="old"
##############
#ut
setwd("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/ut/3m")
sample_list <- list.dirs("./")
sample_list <- gsub(".//", "", sample_list)
sample_list <- sample_list[4:length(sample_list)]

cl_c_final_y <- data.frame()
for (i in sample_list){
  cl_c <- cluster_index(paste0("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/",i,"_manual/outs/"),
                        "/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/ut/3m/",
                        paste0(i,"/", i, "_cell_proportions.csv"),
                        paste0(i,"/", i, "_SC_gene_expression_thresholded.csv"),
                        paste0("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/",i,"_manual/"),
                        i, 0.75)
  cl_c_final_y <- rbind(cl_c, cl_c_final_y)
}
cl_c_final_y$group="young"

setwd("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/ut/18m")
sample_list <- list.dirs("./")
sample_list <- gsub(".//", "", sample_list)
sample_list <- sample_list[3:length(sample_list)]

auc_all_o <- data.frame()
for (i in sample_list){
  counts_sc <- read.csv(paste0(i,"/", i, "_SC_gene_expression_thresholded.csv"), row.names=1,sep="\t")
  auc_sc <- AUCElleR(counts_sc, "SC")
  auc_sc$slide <- i
  auc_all_o <- rbind(auc_sc, auc_all_o)}


auc_all_o_in <- auc_all_o %>% group_by(slide) %>% summarize(Mean=mean(INFLAMMATION))
auc_all_o_in <- auc_all_o_in %>% group_by(sample) %>%summarize(Mean=mean(Mean))
cl_c_final <- data.frame()
for (i in sample_list){
  cl_c <- cluster_index(paste0("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/",i,"_manual/outs/"),
                        "/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/ut/18m/",
                        paste0(i,"/", i, "_cell_proportions.csv"),
                        paste0(i,"/", i, "_SC_gene_expression_thresholded.csv"),
                        paste0("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/",i,"_manual/"),
                        i, 0.75)
  cl_c_final <- rbind(cl_c, cl_c_final)
}
cl_c_final$group="old"


################
cl_y <- cl_c_final_y %>% group_by(sample) %>% summarise(Mean=mean(cluster_size))
cl_y$group <- "y"
cl_o <- cl_c_final %>% group_by(sample) %>% summarise(Mean=mean(cluster_size))
cl_o$group <- "o"
cl=rbind(cl_c_final_y,cl_c_final)
cl$group <- factor(cl$group,levels = c("young", "old"))
ggplot(cl, aes(x=group, y=cluster_size))+ geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', stackratio = .7)+
  theme_classic()

res <- wilcox.test(cluster_size ~ group, data = cl,
                   exact = F, correct=F)





