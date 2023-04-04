#this script calcutes iqr shown in figure 2

# meta contains all the cells with cell type assignments
load("metadata_matrix_13_3_2023.Rdata")
library("reshape")
library("ggrepel")
library("dplyr")
library(tidyr)

meta=matrix_final
meta=meta[grep("estrus", meta$condition),]
a_final = data.frame()


for( tissue in unique(meta$tissue) ) {
  cat( "\tat", tissue, "\n" )
  a=meta[meta$tissue==tissue,]
  a = a%>%group_by(condition,individual, identity) %>%dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  a=as.data.frame(a)
  sample_max <- names(table(a$individual)[table(a$individual)==max(table(a$individual))])[1]
  cell_type <- a[a$individual==sample_max,]$identity
  a_temp_final <- data.frame()
  for (i in names(table(a$individual)[table(a$individual)<max(table(a$individual))])){
    a_temp <- a[a$individual==i,]
    a_add <- data.frame(condition=a_temp$condition[1], individual=a_temp$individual[1],
                        identity=cell_type[!cell_type%in%a_temp$identity], n=0, freq=0)
    a_temp <- rbind(a_temp, a_add)
    a_temp_final <- rbind(a_temp, a_temp_final)
  }
  a <- a[!a$individual%in%names(table(a$individual)[table(a$individual)<max(table(a$individual))]),]
  a <- rbind(a, a_temp_final)
  range= a%>%group_by(condition,identity) %>% dplyr::summarize(Range = IQR(freq, na.rm=TRUE))
  a = a%>%group_by(condition,identity) %>% dplyr::summarize(Mean = mean(freq, na.rm=TRUE))
  a=as.data.frame(a)
  a = reshape(a, idvar = "identity", timevar = "condition", direction = "wide")
  a[is.na(a)] <- 0
  rownames(a) <- a$identity
  a=a[,2:5]
  
  a = data.frame(iqr=apply( a, 1, IQR ), label = "average", tisssue = tissue, cell_type = rownames(a))
  range = data.frame(iqr=range$Range, label = "phase", tisssue = tissue, cell_type = range$identity)
  a =rbind(a, range)
  a_final = rbind(a, a_final)
}

a_final <- as.data.frame(a_final)
a_final$tisssue <- factor(a_final$tisssue, levels=c("va", "ce", "ut", "od", "ov", "sp"))


##########################
#alternative range permuation on sample level, constrained in phase

iqr_range_final=data.frame()
for( tissue in unique(meta$tissue) ) {
  cat( "\tat", tissue, "\n" )
  a_p=meta[meta$tissue==tissue,]
  a_p = a_p%>%group_by(condition,individual, identity) %>%dplyr::summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  a_p=as.data.frame(a_p)
  sample_max <- names(table(a_p$individual)[table(a_p$individual)==max(table(a_p$individual))])[1]
  cell_type <- a_p[a_p$individual==sample_max,]$identity
  a_temp_final <- data.frame()
  for (i in names(table(a_p$individual)[table(a_p$individual)<max(table(a_p$individual))])){
    a_temp <- a_p[a_p$individual==i,]
    a_add <- data.frame(condition=a_temp$condition[1], individual=a_temp$individual[1],
                        identity=cell_type[!cell_type%in%a_temp$identity], n=0, freq=0)
    a_temp <- rbind(a_temp, a_add)
    a_temp_final <- rbind(a_temp, a_temp_final)
  }
  a_p <- a_p[!a_p$individual%in%names(table(a_p$individual)[table(a_p$individual)<max(table(a_p$individual))]),]
  a_p <- rbind(a_p, a_temp_final)
  for (i in unique(a_p$identity)){
    print(i)
    a_ct <- a_p[a_p$identity==i,]
    a_com <- crossing(var1 = as.vector(a_ct[a_ct$condition=="proestrus",]$freq), 
                    var2 = as.vector(a_ct[a_ct$condition=="estrus",]$freq), 
                    var3 = as.vector(a_ct[a_ct$condition=="metestrus",]$freq),
                    var4 = as.vector(a_ct[a_ct$condition=="diestrus",]$freq))
    all_com_iqr = apply( a_com, 1, IQR )
    iqr_range = data.frame(iqr=all_com_iqr, label="sim", tisssue=tissue, cell_type=i)
    iqr_range_final=rbind(iqr_range, iqr_range_final)
    }
}

