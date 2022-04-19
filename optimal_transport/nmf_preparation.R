#################################
# this script uses seurat objects generated in scripts estrus_ov_processing.R, estrus_od_processing.R, estrus_ut_processing.R, 
#estrus_ce_processing.R, estrus_va_processing.R, estrus_sp_processing.R
#################################


library(dplyr)
library(Seurat)
library(gsubfn)
library(dplyr)
library(mgsub)


reads_returneR=function(seurat_object){
  reads=seurat_object[["SCT"]]@data
  #reads=t(as.matrix(reads))
  return(reads)
}

reads_hvg_returneR=function(reads,seurat_object){
  seurat_object = FindVariableFeatures(seurat_object)
  hvg=VariableFeatures(seurat_object)
  reads_hvg=reads[rownames(reads) %in% hvg,]
  reads_hvg=t(as.matrix(reads_hvg))
  return(reads_hvg)
}

metadata_returneR=function(seurat_object){
  metadata=data.frame(id=seurat_object$estrus)
  metadata=metadata%>% mutate(day = case_when(
    id == "0002_mouse1_"| id == "0002_mouse12_"| id =="0002_mouse13_"| id == "0002_mouse6_" ~ "proestrus",
    id == "0002_mouse14_" | id == "0002_mouse3_" | id == "0002_mouse9_"   ~ "estrus", 
    id == "0002_mouse10_" | id == "0002_mouse15_" | id == "0002_mouse11_" ~ "metestrus",
    id == "0002_mouse16_"  | id == "0002_mouse2_" | id == "0002_mouse7_" | id == "0002_mouse8_" ~ "diestrus"))
  row.names(metadata)=colnames(seurat_object)
  return(metadata)
}

cell_ident_returneR=function(seurat_object){
  cell_ident=data.frame(cells=seurat_object$cell_ident)
  row.names(cell_ident)=colnames(seurat_object)
  return(cell_ident)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])}


prepareR <- function(t_combined){
  reads <- reads_returneR(t_combined)
  reads_hvg <- reads_hvg_returneR(reads, t_combined)
  metadata <- metadata_returneR(t_combined)
  cell_ident <- cell_ident_returneR(t_combined)
  metadata <- cbind(metadata, cell_ident)
  list(metadata, reads_hvg)
}

metadata_subseteR <- function(metadata){
    cell_counts <- metadata %>% group_by(day)%>%dplyr::count(cells)
    cell_counts <- cell_counts[cell_counts$n>100,]
    cell_counts <- cell_counts[cell_counts$cells%in%names(which(table(cell_counts$cells)==4)),]
    metadata <- metadata[metadata$cells%in%cell_counts$cells,]
}

file_names <- list.files("seurat_objects/")
file_names <- file_names[grep("combined_new", file_names)]
for (i in file_names){
  name <- substr(i, 1, nchar(i)-10)
  assign(name, loadRData(i))
}

DefaultAssay(ce_combined) <- "SCT"
ce_combined[["SCT_nr"]] <- NULL

DefaultAssay(sp_combined) <- "SCT"
sp_combined[["SCT_nr"]] <- NULL

t_combined <- merge(ov_combined, y=c(od_combined, ut_combined, ce_combined, va_combined, sp_combined), merge.data=T)
t_combined$estrus <- substr(t_combined$orig.ident, 1 , nchar(t_combined$orig.ident)-2)
t_combined$cell_ident <- paste0(t_combined$level_2, substr(t_combined$orig.ident, nchar(t_combined$orig.ident)-2 , nchar(t_combined$orig.ident)))
list[metadata, reads_hvg] <- prepareR(t_combined)

metadata$cell_type <- substr(metadata$cells,1, nchar(metadata$cells)-3)
metadata$cell_type <- mgsub::mgsub(metadata$cell_type, c("\\bB-2\\b", "\\bPC\\b", "\\bMBC\\b", "\\bCD4\\b","\\bCD8\\b","\\biNKT\\b","\\bMAIT\\b","\\bMTC\\b","\\bTC I\\b","\\bTC I-ov\\b","\\bTC I-sp\\b", "\\bM1Mp\\b","\\bM2Mp\\b","\\bMp\\b"), 
                                   c("BC", "BC", "BC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "TC", "Mp", "Mp", "Mp"))
metadata$organ <- substr(metadata$cells, nchar(metadata$cells)-1, nchar(metadata$cells))
metadata$cells <- paste(metadata$cell_type, metadata$organ, sep="_")
metadata <- metadata[,1:3]

metadata <- metadata_subseteR(metadata)
reads_hvg <- reads_hvg[rownames(reads_hvg)%in%rownames(metadata),]

write.table(reads_hvg, "reads_hvg_agg", sep=",", quote = F, row.names = F, col.names = F)
write.table(metadata, "metadata_hvg_agg", sep=",", quote = F)


