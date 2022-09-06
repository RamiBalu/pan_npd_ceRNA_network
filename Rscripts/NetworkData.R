library(WGCNA)

library(nlme)
library(ggplot2)
library(reshape2)
setwd("D:/CeRNA_Modules/PanData/")
#setwd("/home/rami/PanData")

load("Post_COMBAT/ASD_SCZ_BD_CTL_normalized_CR_Cleaned.RData")

datMeta=datMeta[datMeta$Diagnosis!=0,]
datExpr=datExpr[,row.names(datMeta)]

mrna_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "protein_coding"],]
lncRNA_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "lncRNA"],]
miRNAProbes=datProbes[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!="",]

geneProbes=datProbes[datProbes$gene_biotype == "protein_coding" | datProbes$gene_biotype == "lncRNA",]

mir_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!=""],]

miRNAProbes$mirbase_id=gsub("mir","miR",miRNAProbes$mirbase_id)
row.names(mir_expr)=miRNAProbes$mirbase_id[miRNAProbes$external_gene_id %in% rownames(mir_expr)]

datExpr=t(rbind(mrna_expr,lncRNA_expr))


#----#####

load("Post_COMBAT/LMSM_ASD_SCZ_BPD_noctl.RData")
load("Post_COMBAT/Pan_normCRcleaned_moduleData.RData")
sign_module_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj<=0.05,]))
names(dynamicColors)=colnames(datExpr)
dynamicColors_list=c()
i=1
for (c in unique(dynamicColors)){
  dynamicColors_list[[i]]=list(colnames(datExpr)[dynamicColors==c])[[1]]
  
  i=i+1
}
names(dynamicColors_list)=unique(dynamicColors)

LMSM_modules_color=c()
for (c in sign_module_index){ #names(LMSM_WGCNA_Modulegenes)){ #
  for (m in names(dynamicColors_list)){
    if (LMSM_WGCNA_Modulegenes[[c]]==dynamicColors_list[[m]]){ #paste0(c("LMSM ",c),collapse = "")
      LMSM_modules_color[[m]]=dynamicColors_list[[m]]
    }
  }
}
MEs_row=paste0("ME",names(LMSM_modules_color))
MEs_LMSM=MEs[,MEs_row]

LMSM_gene_Expr=datExpr[,unname(unlist(LMSM_modules_color))]
LMSM_colors=dynamicColors[unlist(LMSM_modules_color)]

kMEtable = signedKME((LMSM_gene_Expr),MEs_LMSM)
tableS1 = data.frame(kMEtable[,paste0("kME", names(LMSM_modules_color))])
colnames(tableS1) = paste0("kME.Sponge", 1:10, ".", names(LMSM_modules_color))

LMSM_datProbes=datProbes[datProbes$external_gene_id %in% unlist(LMSM_modules_color),]
LMSM_datProbes=LMSM_datProbes[colnames(LMSM_gene_Expr),]
tableS1 = cbind(LMSM_datProbes, data.frame(Module.Color=LMSM_colors, Module.name = paste0("Sponge ",LMSM_colors)), tableS1)

write.csv(file="Post_COMBAT/kME_table.csv", tableS1)
save(file="Post_COMBAT/FinalizedNetwork.RData", LMSM_gene_Expr, datMeta, LMSM_datProbes, geneTree, LMSM_colors, MEs_LMSM, kMEtable)
