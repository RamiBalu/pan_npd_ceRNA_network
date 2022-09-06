setwd("D:/CeRNA_Modules/PanData/Post_COMBAT/")
load("ASD_SCZ_BD_normalized_CR_Cleaned.RData")
load("LMSM_ASD_SCZ_BPD_noctl.RData")
load("Pan_normCRcleaned_moduleData.RData")

datMeta=datMeta[datMeta$Diagnosis!=0,]
datExpr=datExpr[,row.names(datMeta)]
# expression_matrix

mrna_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "protein_coding"],]
lncRNA_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "lncRNA"],]
miRNAProbes=datProbes[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!="",]

geneProbes=datProbes[datProbes$gene_biotype == "protein_coding" | datProbes$gene_biotype == "lncRNA",]

mir_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!=""],]

miRNAProbes$mirbase_id=gsub("mir","miR",miRNAProbes$mirbase_id)
row.names(mir_expr)=miRNAProbes$mirbase_id[miRNAProbes$external_gene_id %in% rownames(mir_expr)]

datExpr=t(rbind(mrna_expr,lncRNA_expr))

sign_mod_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj<=0.05,]))
sponge_interactions=c()

for (s in sign_mod_index){
  int_mod=data.frame(LMSM_WGCNA_miRSponge[[s]])
  sponge_interactions=rbind.data.frame(sponge_interactions,int_mod)
}
colnames(sponge_interactions)=c("LncRNA","mRNA")
lncRNAs=unique(sponge_interactions$LncRNA)
mRNAs=unique(sponge_interactions$mRNA)

lncRNA_expr=data.frame((datExpr[,lncRNAs]))
LMSM_Modulegenes_expr=data.frame((datExpr[,unique(unlist(LMSM_WGCNA_Modulegenes[sign_mod_index]))]))

names(dynamicColors)=colnames(datExpr)
colorLabels=dynamicColors
save(lncRNA_expr,LMSM_Modulegenes_expr,colorLabels,file="D:/CeRNA_Modules/Validation/mod_presev_datExpr.RData")