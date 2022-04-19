#this script uses Seurat objects generated in scripts pregnant_ut_processing and estrus_ut_processing 

library(Seurat)
require(gtools)
library(dplyr)
library(aggregation)

### prepare data for MI calculation
mi_formatxy<-function(ct,wd,prefix,n_permute,time){
  print("permuting day...")
  day.pmt<-permute_day(pca.out_rownames=colnames(ct), day_master=time, n=n_permute)
  sum(names(day.pmt$ori)==colnames(day.pmt$pmt))
  
  print("generating data x: permutated data and y: original data...")
  x<-structure(rbind(day.pmt$ori,day.pmt$pmt), dimnames=list(c("X",paste0("NULL",1:n_permute)),names(day.pmt$ori)))
  y<-as.matrix(ct)
  print(dim(x))
  print(dim(y))
  
  print("output x, y for MI calculation...")
  setwd(wd)
  print(paste0("x_",prefix,"_miTest.dat"))
  print(paste0("y_",prefix,"_miTest.dat"))
  format.for.mi(x, paste0("x_",prefix,"_miTest.dat"))
  format.for.mi(y, paste0("y_",prefix,"_miTest.dat"))
  out<-list(x,y);names(out)<-c("x","y")
  return(out)
}

# permute data, such that cell labels are permutated agst day for n times
permute_day<-function(pca.out_rownames, day_master, n){
  day.mi<-day_master[,pca.out_rownames]
  day.tmp<-day.mi;names(day.tmp)<-NULL
  day.pmt<-replicate(n,permute(day.tmp))
  print(dim(day.pmt))
  rownames(day.pmt)<-names(day.mi)
  structure(list(day.mi,t(day.pmt)),names=c("ori","pmt"))
}

# format data for mi
format.for.mi<-function(dt, filename){
  tmp <- cbind(genes=rownames(dt), round(dt, 5))
  tmp <- rbind(colnames(tmp), tmp)
  cat(unlist(apply(tmp, 1, paste, collapse="\t"), use.names=FALSE), sep="\n", file=paste(filename, sep=""))
} 

### load MI results
mi_load<-function(miOut.path,mi_xy.out){
  print(paste0("loading miOut:",miOut.path))
  
  miOut<-as.matrix(read.csv(miOut.path, row.names=1, sep="",check.names = F))
  print(dim(miOut))
  print(length(rownames(mi_xy.out$x)))
  print(length(rownames(mi_xy.out$y)))
  miOut<-miOut[rownames(mi_xy.out$x), rownames(mi_xy.out$y)]
  print(dim(miOut))
  
  return(miOut)
}

### ecdf
ecdf.cal.sig<-function(miOut){
  print("calculating ecdf")
  ecdf_out<-apply(miOut[-1,],2,ecdf)
  p_gene<-sapply(1:ncol(miOut),function(a) ecdf_out[[a]](miOut[1,a]));names(p_gene)<-names(ecdf_out)
  
  print("selecting temporal genes")
  p.sig_genes<-p_gene[p_gene==1]
  p.sig_genes.sortedMI<-sort(miOut[1,names(p.sig_genes)],decreasing = T)
  
  genes_sorted.by.MI<-names(sort(miOut[1,],decreasing=T))
  
  ecdf.out<-list(ecdf_out,p_gene,p.sig_genes,p.sig_genes.sortedMI,genes_sorted.by.MI)
  names(ecdf.out)<-c("ecdf_out","p.val_gene","p.sig_genes","p.sig_genes.sortedMI","genes_sorted.by.MI")
  return(ecdf.out)
}

### correct p value based on FDR correction
fdr.adjustment<-function(miOut,fdr){
  library(sgof)
  print(fdr)
  ecdf.out<-ecdf(c(miOut[-1,]))
  p.val<-1-ecdf.out(miOut[1,]);names(p.val)<-colnames(miOut)
  fdr.correction_out<-BH(p.val,fdr)
  out<-list(fdr.correction_out,sort(p.val));names(out)<-c("FDR.correction","p.val")
  out
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}


### Load function and data
estrus <- loadRData("ut_combined.Rdata")
decidua <- loadRData("ut_combined_deci_merged_final.Rdata")
ut_list <- c(SplitObject(estrus, split.by = "orig.ident"), SplitObject(decidua, split.by = "orig.ident"))
ut_list <- lapply(X = ut_list, FUN = SCTransform, method = "glmGamPoi")
ut_list[[14]]$estrus_phase <- "decidua"
ut_list[[15]]$estrus_phase <- "decidua"
ut_list[[16]]$estrus_phase <- "decidua"
ut_list[[17]]$estrus_phase <- "decidua"
ut_list[[18]]$estrus_phase <- "decidua"
ut_list[[19]]$estrus_phase <- "decidua"
ut_list[[20]]$estrus_phase <- "decidua"

ut_combined <- merge(ut_list[[1]], y = unlist(ut_list[2:20], use.names=FALSE), project = "mouse_ut", merge.data=TRUE)

ut_metadata <- data.frame(condition=ut_combined$estrus_phase, cell_type=ut_combined$level_2, cell_names=colnames(ut_combined))

ut_metadata <- ut_metadata[ut_metadata$cell_type%in%c("F", "DeC"),]


for (j in seq(100)){
  metadata_phase <- data.frame()
  for (i in unique(ut_metadata$condition)){
    metadata_phase_temp <- ut_metadata[ut_metadata$condition==i,]
    subset_size <- round(2000*nrow(metadata_phase_temp)/nrow(ut_metadata))
    subset_sample <- sample(1:nrow(metadata_phase_temp), subset_size)
    metadata_phase_temp <- metadata_phase_temp[subset_sample,]
    metadata_phase <- rbind(metadata_phase, metadata_phase_temp)
  }
  
  save(metadata_phase, file=paste0("metadata_subset_",j,".Rdata"))
  j=j+1
}

ut_combined_cell_metadata <- data.frame("cell_names"=colnames(ut_combined), "batch" = ut_combined$orig.ident, "cell_type"=ut_combined$level_2)
ut_combined_cell_metadata <- ut_combined_cell_metadata %>% mutate(condition = case_when(
  batch == "0002_mouse1_ut" | batch == "0002_mouse6_ut" | batch == "0002_mouse12_ut" | batch == "0002_mouse13_ut" ~ 1,
  batch == "0002_mouse2_ut" | batch == "0002_mouse8_ut" | batch == "0002_mouse16_ut"  ~ 4, 
  batch == "0002_mouse3_ut" | batch == "0002_mouse9_ut" | batch == "0002_mouse14_ut" ~ 2,
  batch == "0002_mouse10_ut" | batch == "0002_mouse11_ut" | batch == "0002_mouse15_ut" ~ 3,
  batch == "0002_DC_mmus01_ut_run" | batch == "0002_DC_mmus03_ut_run" | batch == "OE0538_DO-0002_mmus_D06__uterus01_run" | batch == "OE0538_DO-0002_mmus_D07__uterus01_run" |
    batch == "OE0538_DO-0002_mmus_D08__uterus01_run" | batch == "OE0538_DO-0002_mmus_D09__uterus01_run" | batch == "OE0538_DO-0002_mmus_D10__uterus01_run"~ 5 ))


for (i in seq(nrow(ut_combined_cell_metadata))){
  if(ut_combined_cell_metadata[i,3]=="DeC"){ut_combined_cell_metadata[i,4]=6}
}


ut_combined$group <- ut_combined_cell_metadata$condition


ct <-  as.data.frame(ut_combined[["SCT"]]@counts)
day <- data.frame("day"=ut_combined$group)

for (j in seq(from=1, to=100)){
  print(j)
  load(paste0("metadata_subset_",j,".Rdata"))
  ct_t <- ct[,colnames(ct)%in%metadata_phase$cell_names]
  day_t <- as.data.frame(day[rownames(day)%in%metadata_phase$cell_names,])
  wd="/omics/groups/OE0433/internal/ivana/decidua/DEG-MI/re-run/"
  ### Prepare data for mutual information (MI) calculation
  # For each gene, permutate expression against day for n=1000 times
  # Format and output both permutated and un-permutated data to wd
  day_t <- t(day_t)
  colnames(day_t) <- colnames(ct_t)
  rownames(day_t) <- "day"
  mi_xy.out<-mi_formatxy(ct=ct_t,wd=wd,prefix=paste0("ut_mouse_", j),n_permute = 1000,time=day_t)
  ### Calculate MI using java script with the following command line (prefix=mi):
  # java -Xmx10000M -jar interCellMI.jar -e x_mi_miTest.dat -d y_mi_miTest.dat -o mi_miOut.txt
  ### import MI output 
  miOut<-mi_load(paste0("mi_ut_", j, ".txt"), mi_xy.out)
  ### for each gene, calculate significance level of MI results
  # and fdr correct
  ecdfOut<-ecdf.cal.sig(miOut)
  fdr.correction<-fdr.adjustment(miOut,0.05)
  
  ### select genes with significance level that rejects null
  sorted.by.MI<-names(sort(miOut[1,],decreasing = T))
  selected<-sorted.by.MI[sorted.by.MI%in%names(sort(fdr.correction$p.val))[1:fdr.correction$FDR.correction$Rejections]]
  
  
  
  sig_genes <- data.frame(genes=names(fdr.correction$p.val), p_value=fdr.correction$p.val)
  #sig_genes <- sig_genes[sig_genes$p_value<0.05,]
  sig_genes$p_adj <- p.adjust(sig_genes$p_value, method="BH")
  write.csv(sig_genes, paste0("ut_MI_genes_mouse_", j))
}


files <- list.files("/re-run") 
files <- files[grep("ut_MI_genes_mouse", files)]

ut_final <- data.frame()
for (i in files){
  ut <- read.csv(i, row.names = 1)
  ut_final <- rbind(ut, ut_final)
}

ut_final_count <- ut_final %>%
  group_by(genes) %>% count(genes)


ut_final_mi <- ut_final %>%
  group_by(genes) %>% summarize(p_val_adj = fisher(p_adj))

ut_final_mi_sig <- ut_final_mi[ut_final_mi$p_val_adj<0.05,]

write.csv(ut_final_mi, "ut_MI_genes_mouse", quote=F, row.names=F)
