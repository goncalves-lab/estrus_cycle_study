#################################
# this script uses {tissue}_combined_identity files generated in scripts estrus_ov_processing.R, estrus_od_processing.R, estrus_ut_processing.R, 
#estrus_ce_processing.R, estrus_va_processing.R, estrus_sp_processing.R
#################################


source("cell_abundance.R")

params <- list()
params$organ.label <- commandArgs(trailingOnly=TRUE)
params$data.dir <- paste0("cell_abundance/")
params$out.dir <- paste0("cell_abundance/results/", params$organ.label)
params$organ.label.short <- case_when(
  params$organ.label == "oviduct" ~ "od",
  params$organ.label == "vagina" ~ "va",
  params$organ.label == "ovary" ~ "ov",
  params$organ.label == "spleen" ~ "sp",
  params$organ.label == "cervix" ~ "ce",
  params$organ.label == "uterus" ~ "ut")

# create output dir in case it doesn't exists
dir.create(params$out.dir, recursive = T, showWarning=F)

# read data
tissue_combined_identity= read.table(file.path(params$data.dir, paste0(params$organ.label.short , "_combined_identity")),header=T,sep=",") # csv containing cell_names, cell identities, batch and phase information (specific for each organ)

# producing DF (which contains mean, sd, se and ci of cell populations)-named tgc_sd
tissue_combined_identity$condition=as.factor(tissue_combined_identity$condition)
tissue_combined_identity$identity=as.character(tissue_combined_identity$identity)
prepareR(tissue_combined_identity)

#adding 2.1 estus cycles to the DF for plotting
tgc_sd$phase=gsub("\\bestrus\\b","estrus_1",tgc_sd$phase)
tgc_sd$phase=factor(tgc_sd$phase,levels=c("proestrus","estrus_1","metestrus","diestrus"))
tgc_sd$phase_labels=tgc_sd$phase
tgc_sd$cell_types=as.factor(tgc_sd$cell_types)
tgc_sd=tgc_sd[c(1:length(levels(tgc_sd$cell_types)),(length(levels(tgc_sd$cell_types))*2):((length(levels(tgc_sd$cell_types))*3)-1),
                (length(levels(tgc_sd$cell_types))+1):((length(levels(tgc_sd$cell_types))*2)-1),
                (length(levels(tgc_sd$cell_types))*3):nrow(tgc_sd)),]
tgc_sd_2=tgc_sd
tgc_sd_2$phase_labels=paste0(tgc_sd_2$phase,"_2")
tgc_sd_3=rbind(tgc_sd,tgc_sd_2)
tgc_sd_3a=tgc_sd_3[grep("proestrus_2",tgc_sd_3$phase_labels),]
tgc_sd_3a$phase_labels=paste0("proestrus","_3")
tgc_sd_3=rbind(tgc_sd_3,tgc_sd_3a)


write.table(tgc_sd[,1:7],
  file.path(params$out.dir,paste0(params$organ.label.short, "_mean_sd_abundance.csv")),
  row.names = F,
  quote = F,
  sep=",")

cell_abundance = dcast(plot_df_final_model_com, sample~cell_types,value.var = "counts.Freq")

metadata=paste(tissue_combined_identity$batch, tissue_combined_identity$condition, sep="-")
metadata=unique(metadata)
sample_v=c()
phase_v=c()
for(i in seq(from=1, to=length(metadata))){
  sample_v[i]=str_split(metadata[i],"-")[[1]][1]
  phase_v[i]=str_split(metadata[i],"-")[[1]][2]
}
metadata=data.frame("sample"=sample_v, "condition"=phase_v)
cell_abundance=merge(cell_abundance,metadata,by="sample")

# performing compositional regression
cell_abundance$condition=factor(cell_abundance$condition)
cell_abundance[is.na(cell_abundance)] <- 0
phases_combinations=combinations(4, 2, v=levels(cell_abundance$condition), set=TRUE, repeats.allowed=FALSE)
cell_df=data.frame()
cell_vector=c()
for (i in seq(from=1, to=nrow(phases_combinations))){
  cell_subset=cell_abundance[cell_abundance$condition==phases_combinations[i,1]|cell_abundance$condition==phases_combinations[i,2],]
  
  for (j in seq(from=2, to=ncol(cell_subset)-1)){
    cell_subset$total=1-cell_subset[,j]
    Y=acomp(cell_subset[,c(j,ncol(cell_subset))])
    X=factor(cell_subset$condition)
    cell_model = lm(ilr(Y)~X)
    model=try(summary(cell_model))
    if (class(model)=="try-error"){
      cell_vector[1]=paste0(phases_combinations[i,1],"_",phases_combinations[i,2])
      names(cell_vector)[1]="contrast"
      cell_vector[j]=NA
      names(cell_vector)[j]=colnames(cell_subset)[j]} 
    else{
    cell_vector[1]=paste0(phases_combinations[i,1],"_",phases_combinations[i,2])
    names(cell_vector)[1]="contrast"
    cell_vector[j]=summary(cell_model)[[4]][8]
    names(cell_vector)[j]=colnames(cell_subset)[j]}
  }
  cell_df_temp=t(as.data.frame(unname(cell_vector)))
  cell_df=rbind(cell_df,cell_df_temp)
}

colnames(cell_df)=names(cell_vector)
rownames(cell_df)=NULL

write.table(cell_df,
            file.path(params$out.dir, paste0(params$organ.label.short, "_cell_abundance.csv")),
            row.names = F,quote = F,dec=".",sep=",") # writing out results of compositional regression (specific for each organ)

plotteR(tgc_sd, tgc_sd_3,cell_abundance,params$organ.label.short,output_dir = params$out.dir, show_value_labels = F)

