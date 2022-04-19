library(dplyr)
library("RColorBrewer")
library(Seurat)
library(ComplexHeatmap)
library(biomaRt)
library(superheat)
library(data.table)
library(stringr)
library(biomaRt)
library(ggplot2)
library(AUCell)
library(GSEABase)
library(easyGgplot2)
library(RColorBrewer)
library(viridis)
library(circlize)
library(gplots)

#this script uses seurat objects generated in scripts estrus_ut_processing.R, pregnant_ut_processing.R, count matrix 
# retrieved from NCBI’s Gene Expression Omnibus (accession code GSE111976) and MI gene list calculated in script MI.R

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


##################
#loading Seurat object, generating counts and preparing metadata
#################

#######
#mouse
estrus <- loadRData("ut_combined.Rdata")
decidua <- loadRData("ut_combined_deci_merged_final.Rdata")
decidua$estrus_phase <- "decidua"
ut_list <- c(SplitObject(estrus, split.by = "orig.ident"), SplitObject(decidua, split.by = "orig.ident"))
ut_list <- lapply(X = ut_list, FUN = SCTransform, method = "glmGamPoi")
ut_combined <- merge(ut_list[[1]], y = unlist(ut_list[2:20], use.names=FALSE), project = "mouse_ut", merge.data=TRUE)
Idents(ut_combined) <- ut_combined$level_2
ut_combined_f <- subset(ut_combined, idents=c("F", "DeC"))


#####
#human
counts <- read.csv("GSE111976_ct.csv", header=T, row.names = 1, check.names = F)
endometrium_fibroblast <- CreateSeuratObject(counts = counts, project = "endometrium")
endometrium_fibroblast <- PercentageFeatureSet(endometrium_fibroblast, pattern = "^MT-", col.name = "percent.mt")
endometrium_fibroblast <- SCTransform(endometrium_fibroblast, vars.to.regress = "percent.mt", verbose = FALSE)

metadata <- read.csv("GSE111976_summary_C1_day_donor_ctype.csv", row.names = 1)

metadata <- metadata %>% mutate(group = case_when(
  donor == 7 | donor == 4 | donor == 5 | donor == 11 | donor ==6~ "menstrual",
  donor == 13 | donor == 40 | donor == 15 ~ "proliferative_early",
  donor == 33 | donor == 30 | donor == 38 ~ "proliferative_late",
  donor == 12 | donor == 14 | donor == 26 ~ "secretory_early",
  donor == 20 | donor == 56 | donor == 8 ~ "secretory_mid",
  donor == 62 | donor == 59  ~ "secretory_late"
))

metadata <- metadata %>% mutate(level_1 = case_when(
  cell_type == "Unciliated epithelia" | cell_type == "Ciliated" ~ "epithelial",
  cell_type == "Lymphocytes" | cell_type == "Macrophages" ~ "immune",
  cell_type == "Stromal fibroblasts" ~ "F",
  cell_type == "Endothelia" ~ "Endothelia"
))

endometrium_fibroblast$level_2 <- metadata$cell_type
endometrium_fibroblast$level_1 <- metadata$level_1
endometrium_fibroblast$estrus_phase <- metadata$group
Idents(endometrium_fibroblast) <- endometrium_fibroblast$level_2
endometrium_fibroblast <- subset(endometrium_fibroblast, idents="Endothelia", invert=T)
metadata <- metadata[!metadata$cell_type=="Endothelia",]




######################################
#scatter plot
mi_genes_mouse <- read.csv("ut_MI_genes_mouse")
colnames(mi_genes_mouse)[2] <- "p_adj"
mi_gene_human <- read.csv("ut_MI_genes_human", row.names = 1)
mi_gene_human <- mi_gene_human[,c(1,3)]
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="mmusculus_gene_ensembl")
ensemnl_list <- getBM(attributes=c("hsapiens_homolog_associated_gene_name", "external_gene_name"), filters="external_gene_name", values = mi_genes_mouse$genes, mart=ensembl)
ensemnl_list <- ensemnl_list[ensemnl_list$hsapiens_homolog_associated_gene_name != "",]
ensemnl_list_dp <- ensemnl_list$hsapiens_homolog_associated_gene_name[duplicated(ensemnl_list$hsapiens_homolog_associated_gene_name)]

ensemnl_list <- ensemnl_list[!ensemnl_list$hsapiens_homolog_associated_gene_name%in%ensemnl_list_dp,]
mi_genes_mouse <- merge(ensemnl_list,mi_genes_mouse, by.x="external_gene_name", by.y="genes")
mi_gene_human <- mi_gene_human[mi_gene_human$genes%in%mi_genes_mouse$hsapiens_homolog_associated_gene_name,]
mi_genes_mouse <- mi_genes_mouse[mi_genes_mouse$hsapiens_homolog_associated_gene_name%in%mi_gene_human$genes,]

mouse_sig <- mi_genes_mouse[mi_genes_mouse$p_adj<0.05,]
human_sig <- mi_gene_human[mi_gene_human$p_adj<0.05,]

hm_subset_genes <- unique(union(mouse_sig$hsapiens_homolog_associated_gene_name, human_sig$genes))

hm_subset <- merge(mi_genes_mouse, mi_gene_human, by.x="hsapiens_homolog_associated_gene_name", by.y="genes")
hm_subset <- hm_subset[hm_subset$hsapiens_homolog_associated_gene_name%in%hm_subset_genes,]

color_v <- c()
for (i in seq(from=1, to=nrow(hm_subset))){
  if (hm_subset[i,3]<0.05&hm_subset[i,4]<0.05)
  {color_v[i] <- "both"}
  if (hm_subset[i,3]<0.05&hm_subset[i,4]>0.05)
  {color_v[i] <- "mouse"}
  if (hm_subset[i,3]>0.05&hm_subset[i,4]<0.05)
  {color_v[i] <- "human"}
  
}

hm_subset$group <- color_v
hm_subset$p_adj.x <- hm_subset$p_adj.x+10^-5
hm_subset$p_adj.y <- hm_subset$p_adj.y+10^-5
hm_subset$p_adj.x <- -log(hm_subset$p_adj.x)
hm_subset$p_adj.y <- -log(hm_subset$p_adj.y)



print(ggplot(hm_subset, aes(x=p_adj.x, y=p_adj.y, color=group)) + 
        geom_point()+
        theme_classic()+
        scale_color_manual(values=c(brewer.pal(11, "RdGy")[2], brewer.pal(11, "RdGy")[8], brewer.pal(11, "RdGy")[10]))+
        theme(axis.text=element_text(size=20), axis.title=element_text(size=25))+
        labs(x=paste("padj mouse"), y=paste("padj human"))+
        theme(legend.position="right", 
              legend.title = element_text(colour="Black", size=18),
              legend.text = element_text(colour="black", size = 20))+
        guides(color=guide_legend(override.aes=list(fill=NA),title=""))+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "white")))

colnames(hm_subset) <- c("hs_gene_name", "mm_gene_name", "p_adj_mm", "p_adj_hs", "label")
write.csv(hm_subset, "figure_4_p_values", row.names = F, quote=F)


####################################
# overlap log_fc
###################################
h_counts <- as.data.frame(endometrium_fibroblast[["SCT"]]@data)
ensembl=useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
h_counts <- h_counts[rownames(h_counts)%in%mi_gene_human$genes,]
ensemnl_list=getBM(attributes=c("mmusculus_homolog_associated_gene_name", "external_gene_name"), filters="external_gene_name", values = rownames(h_counts), mart=ensembl)

ensemnl_list <- ensemnl_list[ensemnl_list$mmusculus_homolog_associated_gene_name != "",]
ensemnl_list_dp <- ensemnl_list$mmusculus_homolog_associated_gene_name[duplicated(ensemnl_list$mmusculus_homolog_associated_gene_name)]
ensemnl_list <- ensemnl_list[!ensemnl_list$mmusculus_homolog_associated_gene_name%in%ensemnl_list_dp,]

h_counts <- merge(ensemnl_list, h_counts, by.x="external_gene_name", by.y=0)
h_counts <- h_counts[h_counts$external_gene_name%in%mi_gene_human$genes,]
rownames(h_counts) <- make.unique(h_counts$mmusculus_homolog_associated_gene_name)
h_counts <- h_counts[,3:ncol(h_counts)]

m_counts <- as.data.frame(ut_combined_f[["SCT"]]@data)
save(m_counts, file="m_counts.Rdata")
m_counts <- m_counts[rownames(m_counts)%in%rownames(h_counts),]
h_counts <- h_counts[rownames(h_counts)%in%rownames(m_counts),] 

genes_phase_final_av <- data.frame(holder=1:nrow(m_counts))
for (i in unique(mouse_metadata)){
  mouse_metadata_phase <- mouse_metadata[mouse_metadata == i]
  mi_genes_mouse_phase_average <- m_counts[,colnames(m_counts)%in%names(mouse_metadata_phase),]
  genes_phase_av <- data.frame(phase=rowMeans(mi_genes_mouse_phase_average))
  rownames(genes_phase_av) <- rownames(mi_genes_mouse_phase_average)
  colnames(genes_phase_av) <- i
  genes_phase_final_av <- cbind(genes_phase_final_av, genes_phase_av)
}

genes_phase_final_av <- genes_phase_final_av[,2:7]

genes_phase_final_av <- genes_phase_final_av+0.01
genes_phase_final_av_c <- genes_phase_final_av
genes_phase_final_av$proestrus <- genes_phase_final_av$proestrus/genes_phase_final_av_c$diestrus
genes_phase_final_av$estrus <- genes_phase_final_av$estrus/genes_phase_final_av_c$proestrus
genes_phase_final_av$metestrus <- genes_phase_final_av$metestrus/genes_phase_final_av_c$estrus
genes_phase_final_av$decidua <- genes_phase_final_av$decidua/genes_phase_final_av_c$metestrus
genes_phase_final_av$`decidua-DeC` <- genes_phase_final_av$`decidua-DeC`/genes_phase_final_av_c$metestrus
genes_phase_final_av$diestrus <- genes_phase_final_av$diestrus/genes_phase_final_av_c$`decidua-DeC`

genes_phase_final_av <- log2(genes_phase_final_av)

genes_phase_final_av_h <- data.frame(holder=1:nrow(h_counts))
for (i in unique(metadata$group)){
  metadata_phase <- metadata[metadata$group == i,]
  h_subset_phase_average <- h_counts[,colnames(h_counts)%in%metadata_phase$cell_name]
  genes_phase_av <- data.frame(phase=rowMeans(h_subset_phase_average))
  rownames(genes_phase_av) <- rownames(h_subset_phase_average)
  colnames(genes_phase_av) <- i
  genes_phase_final_av_h <- cbind(genes_phase_final_av_h, genes_phase_av)
}


genes_phase_final_av_h <- genes_phase_final_av_h[,2:7]

genes_phase_final_av_h <- genes_phase_final_av_h+0.01
genes_phase_final_av_h_c <- genes_phase_final_av_h

genes_phase_final_av_h$proliferative_early <- genes_phase_final_av_h$proliferative_early/genes_phase_final_av_h_c$secretory_late
genes_phase_final_av_h$proliferative_late <- genes_phase_final_av_h$proliferative_late/genes_phase_final_av_h_c$proliferative_early
genes_phase_final_av_h$secretory_early <- genes_phase_final_av_h$secretory_early/genes_phase_final_av_h_c$proliferative_late
genes_phase_final_av_h$secretory_mid <- genes_phase_final_av_h$secretory_mid/genes_phase_final_av_h_c$secretory_early
genes_phase_final_av_h$secretory_late <- genes_phase_final_av_h$secretory_late/genes_phase_final_av_h_c$secretory_mid

genes_phase_final_av_h <- log2(genes_phase_final_av_h)

genes_phase_final_av$genes <- rownames(genes_phase_final_av)
genes_phase_final_av_h$genes <- rownames(genes_phase_final_av_h)

genes_final <- merge(genes_phase_final_av, genes_phase_final_av_h, by="genes")


hm_intersect <- intersect(human_sig$genes, mouse_sig$hsapiens_homolog_associated_gene_name)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
hm_intersect <- getBM(attributes=c("mmusculus_homolog_associated_gene_name"), filters="external_gene_name", values = hm_intersect, mart=ensembl)
hm_union <- getBM(attributes=c("mmusculus_homolog_associated_gene_name"), filters="external_gene_name", values = hm_subset_genes, mart=ensembl)

genes_final_union <- genes_final[genes_final$genes%in%hm_union$mmusculus_homolog_associated_gene_name,]
genes_final_intersect <- genes_final[genes_final$genes%in%hm_intersect$mmusculus_homolog_associated_gene_name,]


compareR_summary <- function(genes_final, column_choice){
  phase_final <- genes_final[,c(column_choice)]
  res <- cor(phase_final)
  label_v <- c()
  for(i in seq(nrow(phase_final))){
    if (phase_final[i,1] < 0 & phase_final[i,2]< 0){
      label_v[i] = "conserved_down"
    }
    if (phase_final[i,1] > 0 & phase_final[i,2]>0){
      label_v[i] = "conserved_up"
    }
    if (phase_final[i,1] > 0 & phase_final[i,2]<0){
      label_v[i] = "divergent"
    }
    if (phase_final[i,1] < 0 & phase_final[i,2]>0){
      label_v[i] = "divergent"
    }
    
  }
  
  phase_final$label <- label_v
  a <- as.vector(table(phase_final$label))
  names(a)=names(table(phase_final$label))
  a <- a[order(names(a))]
  return(list(res[1,2], a))
  
}


p_pe <- compareR_summary(genes_final_intersect, c(2,8))
cor_p_pe <- p_pe[[1]]
prop_p_pe <- p_pe[[2]]

e_pl <- compareR_summary(genes_final_intersect, c(4,12))
cor_e_pl <- e_pl[[1]]
prop_e_pl <- e_pl[[2]]

m_se <- compareR_summary(genes_final_intersect, c(3,10))
cor_m_se <- m_se[[1]]
prop_m_se <- m_se[[2]]


de_sm <- compareR_summary(genes_final_intersect, c(6,11))
cor_de_sm <- de_sm[[1]]
prop_de_sm <- de_sm[[2]]

dec_sm <- compareR_summary(genes_final_intersect, c(7,11))
cor_dec_sm <- dec_sm[[1]]
prop_dec_sm <- dec_sm[[2]]

d_sl <- compareR_summary(genes_final_intersect, c(5,13))
cor_d_sl <- d_sl[[1]]
prop_d_sl <- d_sl[[2]]

cor_i <- data.frame(p_pe=cor_p_pe, e_pl=cor_e_pl, m_se=cor_m_se, de_sm=cor_de_sm, dec_sm=cor_dec_sm ,d_sl=cor_d_sl)
rownames(cor_i) <- "cor"
cor_i <- as.data.frame(t(cor_i))
cor_i$group <- rownames(cor_i)
cor_i$label=""

cor_i$group <- factor(cor_i$group, levels=c("p_pe", "e_pl", "m_se", "de_sm", "dec_sm","d_sl"))

col_1 = colorRampPalette(rev(brewer.pal(9, "Reds")[9:2]))
print(ggplot(cor_i, aes(group, label))+
        geom_point(aes(colour=cor), size=10)+
        scale_colour_viridis(option = "G")+
        theme_void()+
        theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
        labs(x="", y="",size="")+
        theme(legend.position="right", 
              legend.title = element_text(colour="Black", size=18),
              legend.text = element_text(colour="black", size = 20))+
        scale_size_continuous(limits = c(1, 30), range = c(1,20), breaks = c(1,5,10, 15,20)))



write.csv(genes_final_intersect, "figure_4_log2fc_intersect", row.names = F, quote=F)
write.csv(genes_final_union, "figure_4_log2fc_union", row.names = F, quote=F)

prop_p_pe <- prop_p_pe/sum(prop_p_pe)
prop_e_pl <- prop_e_pl/sum(prop_e_pl)
prop_m_se <- prop_m_se/sum(prop_m_se)
prop_de_sm <- prop_de_sm/sum(prop_de_sm)
prop_dec_sm <- prop_dec_sm/sum(prop_dec_sm)
prop_d_sl <- prop_d_sl/sum(prop_d_sl)

prop <- as.data.frame(rbind(prop_p_pe, prop_e_pl, prop_m_se, prop_de_sm, prop_dec_sm, prop_d_sl))
prop <- stack(as.data.frame(prop))
prop$condition <- rep(c("PD", "EP", "ME", "DeM", "DecM","DDec"), 3)

prop$ind <- factor(prop$ind, levels=c("divergent", "conserved_down", "conserved_up"))
prop$condition <- factor(prop$condition, levels=c("PD", "EP", "ME", "DeM", "DecM","DDec"))


write.csv(prop, "figure_4_proportions", row.names = F, quote=F)
write.csv(cor, "figure_4_cor", row.names = F, quote=F)



####
#simulation

mouse_sig_total <- genes_phase_final_av[rownames(genes_phase_final_av)%in%mouse_sig$external_gene_name,]



ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
human_sig_m <- getBM(attributes=c("mmusculus_homolog_associated_gene_name"), filters="external_gene_name", values = human_sig$genes, mart=ensembl)

human_sig_total <- genes_phase_final_av_h[rownames(genes_phase_final_av_h)%in%human_sig_m$mmusculus_homolog_associated_gene_name,]




simulatoR <- function(prob_vector_m, prob_vector_h, n){
  distribution_up <- c()
  distribution_dw <- c()
  for (j in seq(n)){
    print(j)
    mouse <- sample(c(-1,0,1), 12589, replace=T, prob=prob_vector_m)
    human <- sample(c(-1,0,1), 12589, replace=T, prob=prob_vector_h)
    if (length(mouse) == 12589 & length(human ==12589)){
      phase_final <- data.frame(mouse=mouse, human=human)
      label_v <- c()
      for(i in seq(nrow(phase_final))){
        if (phase_final[i,1] < 0 & phase_final[i,2]< 0){
          label_v[i] = "conserved_down"
        }
        if (phase_final[i,1] > 0 & phase_final[i,2]>0){
          label_v[i] = "conserved_up"
        }
        if (phase_final[i,1] > 0 & phase_final[i,2]<0){
          label_v[i] = "divergent"
        }
        if (phase_final[i,1] < 0 & phase_final[i,2]>0){
          label_v[i] = "divergent"
        }
        if (phase_final[i,1] == 0 & phase_final[i,2]==0){
          label_v[i] = "ns"
        }
        if (phase_final[i,1] != 0 & phase_final[i,2]==0){
          label_v[i] = "divergent"
        }
        if (phase_final[i,1] == 0 & phase_final[i,2]!=0){
          label_v[i] = "divergent"
        }
      }
      
      phase_final$label <- label_v
      a <- as.vector(table(phase_final$label))
      names(a)=names(table(phase_final$label))
      up <- a["conserved_up"]
      down <- a["conserved_down"]
      distribution_up <- c(distribution_up, up)
      distribution_dw <- c(distribution_dw, down)
    }
  }
  na.zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
  }
  distribution_up <- na.zero(distribution_up)
  distribution_dw <- na.zero(distribution_dw)
  return(list("up"=distribution_up, "dw"=distribution_dw))
}


mouse_p_up <- mouse_sig_total$proestrus[mouse_sig_total$proestrus>0]
mouse_p_dw <- mouse_sig_total$proestrus[mouse_sig_total$proestrus<0]

human_pe_up <- human_sig_total$proliferative_early[human_sig_total$proliferative_early>0]
human_pe_dw <- human_sig_total$proliferative_early[human_sig_total$proliferative_early<0]

mouse_e_up <- mouse_sig_total$estrus[mouse_sig_total$estrus>0]
mouse_e_dw <- mouse_sig_total$estrus[mouse_sig_total$estrus<0]

human_pl_up <- human_sig_total$proliferative_late[human_sig_total$proliferative_late>0]
human_pl_dw <- human_sig_total$proliferative_late[human_sig_total$proliferative_late<0]

mouse_m_up <- mouse_sig_total$metestrus[mouse_sig_total$metestrus>0]
mouse_m_dw <- mouse_sig_total$metestrus[mouse_sig_total$metestrus<0]

human_se_up <- human_sig_total$secretory_early[human_sig_total$secretory_early>0]
human_se_dw <- human_sig_total$secretory_early[human_sig_total$secretory_early<0]

mouse_de_up <- mouse_sig_total$decidua[mouse_sig_total$decidua>0]
mouse_de_dw <- mouse_sig_total$decidua[mouse_sig_total$decidua<0]

mouse_dec_up <- mouse_sig_total$`decidua-DeC`[mouse_sig_total$`decidua-DeC`>0]
mouse_dec_dw <- mouse_sig_total$`decidua-DeC`[mouse_sig_total$`decidua-DeC`<0]

human_sm_up <- human_sig_total$secretory_mid[human_sig_total$secretory_mid>0]
human_sm_dw <- human_sig_total$secretory_mid[human_sig_total$secretory_mid<0]

mouse_d_up <- mouse_sig_total$diestrus[mouse_sig_total$diestrus>0]
mouse_d_dw <- mouse_sig_total$diestrus[mouse_sig_total$diestrus<0]

human_sl_up <- human_sig_total$secretory_late[human_sig_total$secretory_late>0]
human_sl_dw <- human_sig_total$secretory_late[human_sig_total$secretory_late<0]


total=12589
p_pe_p_value_s <- simulatoR(c(length(mouse_p_dw)/total, (total-length(mouse_p_dw)-length(mouse_p_up))/total,length(mouse_p_up)/total), 
                            c(length(human_pe_dw)/total, (total-length(human_pe_dw)-length(human_pe_up))/total, length(human_pe_up)/total),
                            100)

de_sm_p_value_s <- simulatoR(c(length(mouse_de_dw)/total, (total-length(mouse_de_dw)-length(mouse_de_up))/total,length(mouse_de_up)/total), 
                             c(length(human_sm_dw)/total, (total-length(human_sm_dw)-length(human_sm_up))/total, length(human_sm_up)/total),
                             100)


e_pl_p_value_s <- simulatoR(c(length(mouse_e_dw)/total, (total-length(mouse_e_dw)-length(mouse_e_up))/total,length(mouse_e_up)/total), 
                            c(length(human_pl_dw)/total, (total-length(human_pl_dw)-length(human_pl_up))/total, length(human_pl_up)/total),
                            100)

m_se_p_value_s <- simulatoR(c(length(mouse_m_dw)/total, (total-length(mouse_m_dw)-length(mouse_m_up))/total,length(mouse_m_up)/total), 
                            c(length(human_se_dw)/total, (total-length(human_se_dw)-length(human_se_up))/total, length(human_se_up)/total),
                            100)

dec_sm_p_value_s <- simulatoR(c(length(mouse_dec_dw)/total, (total-length(mouse_dec_dw)-length(mouse_dec_up))/total,length(mouse_dec_up)/total), 
                              c(length(human_sm_dw)/total, (total-length(human_sm_dw)-length(human_sm_up))/total, length(human_sm_up)/total),
                              100)
d_sl_p_value_s <- simulatoR(c(length(mouse_d_dw)/total, (total-length(mouse_d_dw)-length(mouse_d_up))/total,length(mouse_d_up)/total), 
                            c(length(human_sl_dw)/total, (total-length(human_sl_dw)-length(human_sl_up))/total, length(human_sl_up)/total),
                            100)


exp_value <- c((mean(p_pe_p_value_s[[1]])+mean(p_pe_p_value_s[[2]]))/941,
               (mean(e_pl_p_value_s[[1]])+mean(e_pl_p_value_s[[2]]))/941,
               (mean(m_se_p_value_s[[1]])+mean(m_se_p_value_s[[2]]))/941,
               (mean(de_sm_p_value_s[[1]])+mean(de_sm_p_value_s[[2]]))/941,
               (mean(dec_sm_p_value_s[[1]])+mean(dec_sm_p_value_s[[2]]))/941,
               (mean(d_sl_p_value_s[[1]])+mean(d_sl_p_value_s[[2]]))/941)

ggplot(prop, aes(fill=ind, y=values, x=condition)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(labels = c("divergent", "conserved-dw", "conserved-up"), values = c("grey",brewer.pal(9, "Reds")[c(7,8)]))+
  theme_void()+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  labs(x="", y="",size="")+
  theme(legend.position="right", 
        legend.title = element_text(colour="Black", size=18),
        legend.text = element_text(colour="black", size = 20))+
  geom_segment(aes(x=0.55,xend=1.45,y=exp_value[1],yend=exp_value[1]), size=1)+
  geom_segment(aes(x=1.55,xend=2.45,y=exp_value[2],yend=exp_value[2]), size=1)+
  geom_segment(aes(x=2.55,xend=3.45,y=exp_value[3],yend=exp_value[3]), size=1)+
  geom_segment(aes(x=3.55,xend=4.45,y=exp_value[4],yend=exp_value[4]), size=1)+
  geom_segment(aes(x=4.55,xend=5.45,y=exp_value[5],yend=exp_value[5]), size=1)+
  geom_segment(aes(x=5.55,xend=6.45,y=exp_value[6],yend=exp_value[6]), size=1)
