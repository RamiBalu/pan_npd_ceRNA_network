library(broom)
library(GSVA)
setwd("D:/CeRNA_Modules/PanData/Post_COMBAT/")
#setwd("/home/rami/PanData")
load("ASD_SCZ_BD_CTL_normalized_CR_Cleaned.RData")
load("LMSM_ASD_SCZ_BPD_noctl.RData")
datMeta=datMeta[datMeta$Diagnosis!=0,]
datExpr=datExpr[,row.names(datMeta)]
#DataPreparation
mrna_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "protein_coding"],]
lncRNA_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "lncRNA"],]
miRNAProbes=datProbes[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!="",]

mir_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!=""],]

miRNAProbes$mirbase_id=gsub("mir","miR",miRNAProbes$mirbase_id)
row.names(mir_expr)=miRNAProbes$mirbase_id[miRNAProbes$external_gene_id %in% rownames(mir_expr)]

datExpr=t(rbind(mrna_expr,lncRNA_expr))

expr <- t(datExpr)
gsva_es <- gsva(expr, LMSM_WGCNA_Modulegenes, mx.diff = FALSE)

datMeta$sample_id=rownames(datMeta)
subtype=datMeta[c("sample_id","Diagnosis")]


up_subtype_specific_module <- c()
down_subtype_specific_module <- c()
for (i in seq_len(dim(gsva_es)[1])){
  Response <- gsva_es[i, ]
  Treatment <- subtype$Diagnosis 
  RT <- tidy(pairwise.t.test(Response, Treatment, p.adjust.method = "BH"))
  mu <- unlist(lapply(seq(unique(subtype$Diagnosis)), function(i) mean(Response[Treatment==unique(subtype$Diagnosis)[i]])))
  up_subtype_specific <- unique(subtype$Diagnosis)[which(mu == max(mu))]
  down_subtype_specific <- unique(subtype$Diagnosis)[which(mu == min(mu))]
  down_subtype_specific_p.value <- RT$p.value[c(which(RT$group1 == down_subtype_specific), 
                                                which(RT$group2 == down_subtype_specific))]
  
  up_subtype_specific_p.value <- RT$p.value[c(which(RT$group1 == up_subtype_specific), 
                                              which(RT$group2 == up_subtype_specific))]
  up_subtype_specific_module[i] <- up_subtype_specific
  down_subtype_specific_module[i] <- down_subtype_specific
  if(all(up_subtype_specific_p.value<0.5)){
    up_subtype_specific_module[i] <- up_subtype_specific}
  else{
    up_subtype_specific_module[i] <- "uncertain"}
  if(all(down_subtype_specific_p.value<0.5)){
    down_subtype_specific_module[i] <- down_subtype_specific}
  else{
    down_subtype_specific_module[i] <- "uncertain"}
  print(paste(up_subtype_specific,down_subtype_specific,sep = "-"))
}
  
LMSM_subtype_specific_module_down <- cbind(rownames(LMSM_WGCNA_modules_p.values), down_subtype_specific_module)
LMSM_subtype_specific_module_up <- cbind(rownames(LMSM_WGCNA_modules_p.values), up_subtype_specific_module)

colnames(LMSM_subtype_specific_module_up) <- colnames(LMSM_subtype_specific_module_down) <- c("LMSM_modules", "subtype specific")
save.image(file="LMSM_GSVA_enrichment.RData")
