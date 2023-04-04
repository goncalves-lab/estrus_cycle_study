#this script fit glm to calculate relation between estrous cycle phase and AUCell activity scores

library(dplyr)
library(fitdistrplus)
library("MASS")
library(glmmTMB)
library("ggplot2")

load("metadata_matrix_2_12_2022.Rdata")

ab_matrix <- matrix_final[,c(1:5,26)]

for (i in seq(nrow(ab_matrix))){
  if (ab_matrix[i,1] == "SC I") {ab_matrix[i,1] = paste0("SC I-", ab_matrix[i,3])}
  if (ab_matrix[i,1] == "EpC I") {ab_matrix[i,1] = paste0("EpC I-", ab_matrix[i,3])}
  if (ab_matrix[i,1] == "TC I") {ab_matrix[i,1] = paste0("TC I-", ab_matrix[i,3])}
  if (ab_matrix[i,1] == "APC I") {ab_matrix[i,1] = paste0("APC I-", ab_matrix[i,3])}
}

ab_matrix <- ab_matrix[grep("estrus", ab_matrix$condition),]

ecm <- ab_matrix$INFLAMMATION[!is.na(ab_matrix$INFLAMMATION)]

descdist(ecm, discrete = FALSE)
fit.beta <- fitdist(ecm[ecm>0], "beta")
#plot(fit.beta)

ab_matrix <- ab_matrix[!is.na(ab_matrix$INFLAMMATION),]
ab_matrix <- ab_matrix[ab_matrix$INFLAMMATION>0,]

phase <- c("proestrus", "estrus", "metestrus", "diestrus")
phase_comb <- combn(phase, 2)
ecm_df_final <- data.frame()
ab_matrix <- ab_matrix[ab_matrix$identity=="F",]
ab_matrix <- ab_matrix[!is.na(ab_matrix$tissue),]


for (i in unique(ab_matrix$tissue)){
  ab_matrix_t <- ab_matrix[ab_matrix$tissue==i,]
  for (j in seq(ncol(phase_comb))){
    ab_matrix_p <- ab_matrix_t[ab_matrix_t$condition%in%phase_comb[,j],]
    ecm_model <- glmmTMB(INFLAMMATION~condition+(1 | batch), data = ab_matrix_p, family=beta_family())
    ecm_df <- data.frame("tissue"=i, "phase_1"=phase_comb[1,j], "phase_2"=phase_comb[2,j],
                         "p_value"=summary(ecm_model)[6][[1]][[1]][8])
    ecm_df_final <- rbind(ecm_df, ecm_df_final)
  }
    
}
ecm_df_final$p_adj <- p.adjust(ecm_df_final$p_value, method = "fdr")
write.csv(ecm_df_final, "INFLAMMATION_ic_aucell_gmm.csv" ,row.names = F, quote = F)


ab_matrix <- ab_matrix[ab_matrix$identity=="F",]
ab_matrix <- ab_matrix[!is.na(ab_matrix$tissue),]
ab_matrix <- ab_matrix[ab_matrix$condition%in%c("diestrus"),]

tissue <- c("ov", "od", "ut", "ce", "va")
phase_comb <- combn(tissue, 2)

ecm_df_final <- data.frame()


for (j in seq(ncol(phase_comb))){
  ab_matrix_p <- ab_matrix[ab_matrix$tissue%in%phase_comb[,j],]
  ecm_model <- glmmTMB(INFLAMMATION~tissue+(1 | batch), data = ab_matrix_p, family=beta_family())
  ecm_df <- data.frame("tissue_1"=phase_comb[1,j], "tissue_2"=phase_comb[2,j],
                         "p_value"=summary(ecm_model)[6][[1]][[1]][8])
  ecm_df_final <- rbind(ecm_df, ecm_df_final)
}
  

ecm_df_final$p_adj <- p.adjust(ecm_df_final$p_value, method = "fdr")
write.csv(ecm_df_final, "INFLAMMATION_diestrus_aucell_gmm.csv" ,row.names = F, quote = F)

