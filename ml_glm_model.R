#################################
# this script uses {tissue}_combined_identity files and seurat objects generated in scripts estrus_ov_processing.R, estrus_od_processing.R, estrus_ut_processing.R, 
#estrus_ce_processing.R, estrus_va_processing.R, estrus_sp_processing.R
#################################


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(dplyr)
library("MASS")
library(lme4)
library(parallel)

if (length(args)<6) {
  stop("six arguments must be supplied (path to seurat object, path to identity csv, output directory, condition 1, condition 2,output string).n", call.=FALSE)}

source("qc_normalization_cluster_annotation.R")


t_combined <- loadRData(args[1])

t_list= count_produceR(t_combined)

t_combined_counts=t_list[[1]]
t_combined_identity=read.csv(args[2])

overdisp_fun <- function(model) {
  rdf <- unname(model$AICtab[5])
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


DEGeR_nb=function(x,df_cell,condition_1,condition_2){
  gene_test=x
  tissue_combined_identity_cell_test=df_cell
  tissue_combined_identity_cell_test$gene=unname(gene_test)
  tissue_combined_identity_cell_test$cells_gene=names(gene_test)
  tissue_combined_identity_cell_test_pd_1=data.frame()
  tissue_combined_identity_cell_test_pd_2=data.frame()
  tissue_combined_identity_cell_test_pd_1=tissue_combined_identity_cell_test[tissue_combined_identity_cell_test$condition==condition_1,]
  tissue_combined_identity_cell_test_pd_2=tissue_combined_identity_cell_test[tissue_combined_identity_cell_test$condition==condition_2,]
  tissue_combined_identity_cell_test_pd=rbind(tissue_combined_identity_cell_test_pd_1,tissue_combined_identity_cell_test_pd_2)
  tissue_combined_identity_cell_test_pd$sample=as.numeric(as.factor(tissue_combined_identity_cell_test_pd$batch))
    model_pd=try(glmer.nb(gene~condition+(1 | sample), data=tissue_combined_identity_cell_test_pd))
    if (class(model_pd)!="try-error") {
      summary_nb=summary(model_pd)
      if(length(summary_nb$optinfo[4][1]$conv$lme4)!=0){
        fit_info <- "bad"
      } else {fit_info <- "good"}
    } else {fit_info <- "bad"}
    
    if ((class(model_pd)=="try-error")|(fit_info=="bad")){
      model_pd <- try(glmer.nb(
        gene~condition+(1 | sample),
        data = tissue_combined_identity_cell_test_pd,
        control = glmerControl(
          optimizer = "optimx",
          calc.derivs = FALSE,
          optCtrl = list(
            method = "nlminb",
            starttests = FALSE,
            kkt = FALSE
          )
        )  
      ))
      }
      if (class(model_pd)=="try-error"){
        p_value=NA
        overdispesion_nb=NA
        AIC=NA
        z_value=NA
      }
    else {summary_nb=summary(model_pd)
      overdispesion_nb=overdisp_fun(summary_nb)
      AIC=unname(summary_nb$AICtab[1])
        if(length(summary_nb$optinfo[4][1]$conv$lme4)==0){
          p_value=summary_nb$coefficients[8]
          z_value=summary_nb$coefficients[6]
        } 
        else {p_value=NA
        z_value=NA}}
  fold_change=log2((mean(tissue_combined_identity_cell_test_pd[grep(condition_1,tissue_combined_identity_cell_test_pd$condition),]$gene)+0.001)/
                     (mean(tissue_combined_identity_cell_test_pd[grep(condition_2,tissue_combined_identity_cell_test_pd $condition),]$gene)+0.001))
  my_list <- list("overdispersion"=unname(overdispesion_nb[4]),"AIC" = AIC, "fold_change" = fold_change,"z_value"= z_value,"p_value" = p_value)
  return(my_list)
}

setwd(args[3])#output directory (e.g."~/scRNA-seq/17270/dge")

t_combined_identity$identity <- as.factor(t_combined_identity$identity)

tissue_DGER=function(tissue_combined_identity,tissue_combined_counts, con_1, con_2,out_string,folder){
  gene_df_nb_final =data.frame()
  for (cell_pop in levels(tissue_combined_identity$identity)){
    file_names = list.files(folder)
    if (!paste0(out_string,"_",cell_pop,"_",con_1, "_",con_2)%in%file_names){
    tissue_combined_identity_cell = tissue_combined_identity[tissue_combined_identity$identity == cell_pop,]
    tissue_combined_counts_cell = tissue_combined_counts[,colnames(tissue_combined_counts) %in% tissue_combined_identity_cell$cells_names]
    tissue_combined_counts_cell = tissue_combined_counts_cell[!!rowSums(abs(tissue_combined_counts_cell)),]
    non_zero=c()
    for (i in seq(from=1, to=nrow(tissue_combined_counts_cell))){
      if ((names(table(tissue_combined_counts_cell[i,])[1])!="0")|(unname(table(tissue_combined_counts_cell[i,])[1])<ncol(tissue_combined_counts_cell)-round(ncol(tissue_combined_counts_cell)*0.05)))
      {non_zero[i]=i}
    }
    non_zero=non_zero[!is.na(non_zero)]
    tissue_combined_counts_cell=tissue_combined_counts_cell[c(non_zero),]
    gene_list_nb=list()
    gene_df_nb=data.frame()
    tissue_combined_identity_cell <- tissue_combined_identity_cell[tissue_combined_identity_cell$cells_names%in%colnames(tissue_combined_counts_cell),]
    tissue_combined_counts_cell_list <- setNames(split(tissue_combined_counts_cell, seq(nrow(tissue_combined_counts_cell))), rownames(tissue_combined_counts_cell))
    cores=detectCores()
    
    
    gene_list_nb=mclapply(tissue_combined_counts_cell_list[1:100], FUN=DEGeR_nb, df_cell=tissue_combined_identity_cell, condition_1= con_1, condition_2= con_2, mc.cores=5)
    gene_df_nb <- data.frame(matrix(unlist(gene_list_nb), nrow=length(gene_list_nb), byrow=T))
    rownames(gene_df_nb)=names(gene_list_nb)
    colnames(gene_df_nb)=c("Overdispersion","AIC","log2foldchange", "z_value","p_value")
    p_values_adj_nb=p.adjust(gene_df_nb$p_value, method = "BH", n = length(gene_df_nb$p_value))
    gene_df_nb$p_value_adj_nb=p_values_adj_nb
    write.table(gene_df_nb,file=paste0(out_string,"_",cell_pop,"_",con_1, "_",con_2),quote = F,sep=",",row.names = T)
    gene_df_nb$cell_type=cell_pop
    gene_df_nb_final =rbind(gene_df_nb_final, gene_df_nb)
  }}
  return(gene_df_nb_final)
}
tissue_DGER(t_combined_identity, t_combined_counts, args[5], args[6], args[4], args[3])
