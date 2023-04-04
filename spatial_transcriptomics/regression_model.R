###############################
#this script calculates best neighborhood predictor of stromal cell inflammation, it takes as imoput files
#produced in deconvolution script"
##############################


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
library(MuMIn)
load("/omics/odcf/analysis/OE0538_projects/DO-0002/mmus/aucell/geneSets_AUCell.Rdata")
geneSets_final <- geneSets_final[10]
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
  colnames(auc_df_final_m) <- c(paste0("INFLAMMATION_", cell), "keys")
  return(auc_df_final_m)
  
}

spot_neigbouR <- function(coord,prop, auc_sc, auc_other, other_label, inflammation, cell, infl,lab){
  prop$keys <- rownames(prop)
  prop <- na.omit(prop)
  coord <- coord[rownames(coord)%in%rownames(prop),]
  distances <- as.data.frame(as.matrix(dist(coord)))
  prop <- left_join(prop, auc_sc, by="keys")
  if (other_label==T){
    prop <- left_join(prop, auc_other, by="keys")
  }
  prop[is.na(prop)] <- 0
  cor_df_final <- data.frame()
  for (i in colnames(distances)){
    spot <- data.frame(spots=distances[,i], spot_names=rownames(distances))
    spot <- spot[spot$spots<10,]
    spot <- spot[order(spot$spots),]
    fib <- prop[prop$keys%in%spot[1,2],]
    other <- prop[prop$keys%in%spot[2:nrow(spot),2],]
    
    if (inflammation==T){
      if (other_label==T){
        a=other[other[, infl]>0,]
        #b=c(nrow(a),b)
        if(nrow(a)==0){cor_df <- data.frame(sc=fib$INFLAMMATION_SC, other=0)}
        else{cor_df <- data.frame(sc=fib$INFLAMMATION_SC, other=sum(other[, infl])/nrow(a))}
      }
      else {
        a=other[other$INFLAMMATION_SC>0,]
        cor_df <- data.frame(sc=fib$INFLAMMATION_SC, other=sum(other$INFLAMMATION_SC)/nrow(a))}
    } else {cor_df <- data.frame(sc=fib$INFLAMMATION_SC, other=mean(other[,cell]))}
    cor_df_final <- rbind(cor_df, cor_df_final)
  }
  cor_df_final=cor_df_final[cor_df_final$sc>0,]
  cor_df_final=cor_df_final[cor_df_final$other>0,]
    print(ggplot(cor_df_final, aes(x = sc, y = other)) + 
          geom_point()+theme_classic()+xlab("INFLAMATION-center")+ylab(lab)
  )
  return(cor_df_final)
}

organizeR <- function(spatial_loc, prop_loc, sc, apc, tc, bc, slide){
  spatial <- Load10X_Spatial(spatial_loc)
  setwd("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/ut/3m")
  prop <- read.csv(prop_loc, row.names=1,sep=",")
  coord <- GetTissueCoordinates(spatial)
  
  counts_sc <- read.csv(sc, row.names=1,sep="\t")
  auc_sc <- AUCElleR(counts_sc, "SC")
  
  counts_apc <- read.csv(apc, row.names=1,sep="\t")
  auc_apc <- AUCElleR(counts_apc, "APC")
  
  counts_tc <- read.csv(tc, row.names=1,sep="\t")
  auc_tc <- AUCElleR(counts_tc, "TC")
  
  counts_bc <- read.csv(bc, row.names=1,sep="\t")
  auc_bc <- AUCElleR(counts_bc, "NKC")
  
  pdf(paste0(slide, "_sc_inflammation.pdf"), height=5, width=5)
  cor_sc_inflammation <- spot_neigbouR(coord, prop, auc_sc, auc_sc,F, T, "SC", "INFLAMMATION_SC", "INFLAMMATION-Neighbour")
  cor_sc_prop <- spot_neigbouR(coord, prop, auc_sc, auc_sc,F, F, "SC", "INFLAMMATION_SC","Proportion-Neighbour")
  cor_sc_inflammation$cell_type <- "sc"
  cor_sc_prop$cell_type <- "sc"
  dev.off()
  
  pdf(paste0(slide, "_apc_inflammation.pdf"), height=5, width=5)
  cor_apc_inflammation <- spot_neigbouR(coord, prop, auc_sc, auc_apc,T, T, "APC", "INFLAMMATION_APC","INFLAMMATION-Neighbour")
  cor_apc_prop <- spot_neigbouR(coord, prop, auc_sc, auc_apc ,T, F, "APC","INFLAMMATION_APC" ,"Proportion-Neighbour")
  cor_apc_inflammation$cell_type <- "apc"
  cor_apc_prop$cell_type <- "apc"
  dev.off()
  
  pdf(paste0(slide, "_tc_inflammation.pdf"), height=5, width=5)
  cor_tc_inflammation <- spot_neigbouR(coord, prop, auc_sc, auc_tc,T, T, "TC", "INFLAMMATION_TC","INFLAMMATION-Neighbour")
  cor_tc_prop <- spot_neigbouR(coord, prop, auc_sc, auc_tc ,T, F, "TC","INFLAMMATION_TC" ,"Proportion-Neighbour")
  cor_tc_inflammation$cell_type <- "tc"
  cor_tc_prop$cell_type <- "tc"
  dev.off()
  
  pdf(paste0(slide, "_nkc_inflammation.pdf"), height=5, width=5)
  cor_bc_inflammation <- spot_neigbouR(coord, prop, auc_sc, auc_bc,T, T, "NKC", "INFLAMMATION_NKC","INFLAMMATION-Neighbour")
  cor_bc_prop <- spot_neigbouR(coord, prop, auc_sc, auc_bc ,T, F, "BC","INFLAMMATION_NKC" ,"Proportion-Neighbour")
  cor_bc_inflammation$cell_type <- "nkc"
  cor_bc_prop$cell_type <- "nkc"
  dev.off()
  
  pdf(paste0(slide, "_il6_inflammation.pdf"), height=5, width=5)
  cor_sc_il6 <- spot_neigbouR_gene(coord, prop, auc_sc, counts_sc, "Il6", "Il6 expression")
  cor_sc_il6$cell_type <- "Il6"
  cor_sc_il6$feature <- "gene"
  dev.off()
  
  cor_inflammation <- rbind(cor_sc_inflammation, cor_apc_inflammation, cor_tc_inflammation, cor_bc_inflammation)
  cor_inflammation$feature <- "auc"
  cor_prop <- rbind(cor_sc_prop, cor_apc_prop, cor_tc_prop, cor_bc_prop)
  cor_prop$feature <- "prop"
  cor_df <- rbind(cor_sc_il6, cor_prop, cor_inflammation)
  return(cor_df)
}



V11M25_311_B1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/ha/3M/V11M25-311_B1_manual_repmerged/outs/",
                           "V11M25-311_B1_corrected_cell_proportions.csv", 
                           "V11M25-311_B1/V11M25-311_B1_SC_gene_expression_thresholded.csv",
                           "V11M25-311_B1/V11M25-311_B1_APC_gene_expression_thresholded.csv",
                           "V11M25-311_B1/V11M25-311_B1_TC_gene_expression_thresholded.csv",
                           "V11M25-311_B1/V11M25-311_B1_BC_gene_expression_thresholded.csv",
                           "V11M25-311_B1")

V11M25_312_A1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/ha/3M/V11M25-312_A1_manual_repmerged/outs/",
                           "V11M25-312_A1_corrected_cell_proportions.csv", 
                           "V11M25-312_A1/V11M25-312_A1_SC_gene_expression_thresholded.csv",
                           "V11M25-312_A1/V11M25-312_A1_APC_gene_expression_thresholded.csv",
                           "V11M25-312_A1/V11M25-312_A1_TC_gene_expression_thresholded.csv",
                           "V11M25-312_A1/V11M25-312_A1_BC_gene_expression_thresholded.csv",
                           "V11M25-312_A1")

V11M25_312_B1<- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/ha/3M/V11M25-312_B1_manual_repmerged/outs/",
                           "V11M25-312_B1_corrected_cell_proportions.csv", 
                           "V11M25-312_B1/V11M25-312_B1_SC_gene_expression_thresholded.csv",
                           "V11M25-312_B1/V11M25-312_B1_APC_gene_expression_thresholded.csv",
                           "V11M25-312_B1/V11M25-312_B1_TC_gene_expression_thresholded.csv",
                           "V11M25-312_B1/V11M25-312_B1_BC_gene_expression_thresholded.csv",
                           "V11M25-312_B1")

V11M25_312_C1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/ha/3M/V11M25-312_C1_manual_repmerged/outs/",
                           "V11M25-312_C1_corrected_cell_proportions.csv", 
                           "V11M25-312_C1/V11M25-312_C1_SC_gene_expression_thresholded.csv",
                           "V11M25-312_C1/V11M25-312_C1_APC_gene_expression_thresholded.csv",
                           "V11M25-312_C1/V11M25-312_C1_TC_gene_expression_thresholded.csv",
                           "V11M25-312_C1/V11M25-312_C1_BC_gene_expression_thresholded.csv",
                           "V11M25-312_C1")

V11M25_312_D1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/ha/3M/V11M25-312_D1_manual_repmerged/outs/",
                           "V11M25-312_D1_corrected_cell_proportions.csv", 
                           "V11M25-312_D1/V11M25-312_D1_SC_gene_expression_thresholded.csv",
                           "V11M25-312_D1/V11M25-312_D1_APC_gene_expression_thresholded.csv",
                           "V11M25-312_D1/V11M25-312_D1_TC_gene_expression_thresholded.csv",
                           "V11M25-312_D1/V11M25-312_D1_BC_gene_expression_thresholded.csv",
                           "V11M25-312_D1")
V11M25_311_B1$slide <- "V11M25_311_B1"
V11M25_312_A1$slide <- "V11M25_312_A1"
V11M25_312_B1$slide <- "V11M25_312_B1"
V11M25_312_C1$slide <- "V11M25_312_C1"
V11M25_312_D1$slide <- "V11M25_312_D1"

cor_df <- rbind(V11M25_311_B1, V11M25_312_A1, V11M25_312_B1, V11M25_312_C1, V11M25_312_D1)

model_fit_df <- data.frame()
for(i in unique(cor_df$feature)){
  print(i)
  cor_df_feature <- cor_df[cor_df$feature==i,] 
  for(j in unique(cor_df_feature$cell_type)){
    print(j)
    cor_df_cell <- cor_df_feature[cor_df_feature$cell_type==j,]
    model_pd=lmer(sc~log(other)+(1 | slide), data=cor_df_cell)
    r2 <- r.squaredGLMM(model_pd)
    model_fit <- data.frame(feature=i, cell_type=j, fit=r2)
    model_fit_df <- rbind(model_fit, model_fit_df)
    print(ggplot(cor_df_cell, aes(x = log(other), y = sc, color=slide))+ 
            scale_color_manual(values=c(brewer.pal(9,"Reds")[c(2,4,6,8,9)]))+
            geom_point()+theme_classic()+xlab(paste0(i,"_",j))+ylab("INFLAMATION-center")
            #stat_smooth(method = "lm", col = "blue")
            
    )
  }
}

model_fit_df$label <- paste0(model_fit_df$feature, "_", model_fit_df$cell_type)
model_fit_df <- model_fit_df[c(1:8),]
model_fit_df$label <- factor(model_fit_df$label, levels=c("auc_sc", "auc_apc",
                                                          "auc_nkc", "auc_tc",
                                                          "prop_sc", "prop_apc",
                                                          "prop_nkc", "prop_tc"))
ggplot(model_fit_df, aes(x=label, y=fit.R2m)) +
  geom_bar(stat="identity", fill="black")+theme_classic()


V12N14_363_A1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/V12N14-363_A1_manual/outs/",
                           "V12N14-363_A1_corrected_cell_proportions.csv", 
                           "V12N14-363_A1/V12N14-363_A1_SC_gene_expression_thresholded.csv",
                           "V12N14-363_A1/V12N14-363_A1_APC_gene_expression_thresholded.csv",
                           "V12N14-363_A1/V12N14-363_A1_TC_gene_expression_thresholded.csv",
                           "V12N14-363_A1/V12N14-363_A1_NKC_gene_expression_thresholded.csv",
                           "V12N14-363_A1")

V12N14_364_A1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/V12N14-364_A1_manual/outs/",
                           "V12N14-364_A1_corrected_cell_proportions.csv", 
                           "V12N14-364_A1/V12N14-364_A1_SC_gene_expression_thresholded.csv",
                           "V12N14-364_A1/V12N14-364_A1_APC_gene_expression_thresholded.csv",
                           "V12N14-364_A1/V12N14-364_A1_TC_gene_expression_thresholded.csv",
                           "V12N14-364_A1/V12N14-364_A1_NKC_gene_expression_thresholded.csv",
                           "V12N14-364_A1")

V12N14_364_B1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/V12N14-364_B1_manual/outs/",
                           "V12N14-364_B1_corrected_cell_proportions.csv", 
                           "V12N14-364_B1/V12N14-364_B1_SC_gene_expression_thresholded.csv",
                           "V12N14-364_B1/V12N14-364_B1_APC_gene_expression_thresholded.csv",
                           "V12N14-364_B1/V12N14-364_B1_TC_gene_expression_thresholded.csv",
                           "V12N14-364_B1/V12N14-364_B1_NKC_gene_expression_thresholded.csv",
                           "V12N14-364_B1")

V12N14_364_C1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/V12N14-364_C1_manual/outs/",
                           "V12N14-364_C1_corrected_cell_proportions.csv", 
                           "V12N14-364_C1/V12N14-364_C1_SC_gene_expression_thresholded.csv",
                           "V12N14-364_C1/V12N14-364_C1_APC_gene_expression_thresholded.csv",
                           "V12N14-364_C1/V12N14-364_C1_TC_gene_expression_thresholded.csv",
                           "V12N14-364_C1/V12N14-364_C1_NKC_gene_expression_thresholded.csv",
                           "V12N14-364_C1")

V12N14_364_D1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/V12N14-364_D1_manual/outs/",
                           "V12N14-364_D1_corrected_cell_proportions.csv", 
                           "V12N14-364_D1/V12N14-364_D1_SC_gene_expression_thresholded.csv",
                           "V12N14-364_D1/V12N14-364_D1_APC_gene_expression_thresholded.csv",
                           "V12N14-364_D1/V12N14-364_D1_TC_gene_expression_thresholded.csv",
                           "V12N14-364_D1/V12N14-364_D1_NKC_gene_expression_thresholded.csv",
                           "V12N14-364_D1")
V12N14_363_B1 <- organizeR("/omics/odcf/analysis/OE0538_projects/DO-0007/mmus/rso_ha/spatial_tr/data/spaceranger_out/V12N14-363_B1_manual/outs/",
                           "V12N14-363_B1_corrected_cell_proportions.csv", 
                           "V12N14-363_B1/V12N14-363_B1_SC_gene_expression_thresholded.csv",
                           "V12N14-363_B1/V12N14-363_B1_APC_gene_expression_thresholded.csv",
                           "V12N14-363_B1/V12N14-363_B1_TC_gene_expression_thresholded.csv",
                           "V12N14-363_B1/V12N14-363_B1_NKC_gene_expression_thresholded.csv",
                           "V12N14-363_B1")

V12N14_363_A1$slide <- "V12N14_363_A1"
V12N14_364_A1$slide <- "V12N14_364_A1"
V12N14_364_B1$slide <- "V12N14_364_B1"
V12N14_364_C1$slide <- "V12N14_364_C1"
V12N14_364_D1$slide <- "V12N14_364_D1"
V12N14_363_B1$slide <- "V12N14_363_B1"

cor_df <- rbind(V12N14_363_A1, V12N14_364_A1, V12N14_364_B1, V12N14_364_C1, V12N14_364_D1, V12N14_363_B1)

model_fit_df <- data.frame()
for(i in unique(cor_df$feature)){
  print(i)
  cor_df_feature <- cor_df[cor_df$feature==i,] 
  for(j in unique(cor_df_feature$cell_type)){
    print(j)
    cor_df_cell <- cor_df_feature[cor_df_feature$cell_type==j,]
    model_pd=lmer(sc~other+(1 | slide), data=cor_df_cell)
    r2 <- r.squaredGLMM(model_pd)
    model_fit <- data.frame(feature=i, cell_type=j, fit=r2)
    model_fit_df <- rbind(model_fit, model_fit_df)
    print(ggplot(cor_df_cell, aes(x = other, y = sc, color=slide))+ 
            scale_color_manual(values=c(brewer.pal(9,"Reds")[c(2,4,6,7,8,9)]))+
            geom_point()+theme_classic()+xlab(paste0(i,"_",j))+ylab("INFLAMATION-center")
          #stat_smooth(method = "lm", col = "blue")
          
    )
  }
}

model_fit_df$label <- paste0(model_fit_df$feature, "_", model_fit_df$cell_type)
model_fit_df <- model_fit_df[c(1:8),]
model_fit_df$label <- factor(model_fit_df$label, levels=c("auc_sc", "auc_apc",
                                                          "auc_nkc", "auc_tc",
                                                          "prop_sc", "prop_apc",
                                                          "prop_nkc", "prop_tc"))
ggplot(model_fit_df, aes(x=label, y=fit.R2m)) +
  geom_bar(stat="identity", fill="black")+theme_classic()

