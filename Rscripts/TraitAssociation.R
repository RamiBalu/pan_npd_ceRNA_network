library(WGCNA)

library(nlme)
library(ggplot2)
library(reshape2)
setwd("D:/1_CeRNA_Modules/PanData/")
#setwd("/home/rami/PanData")

load("Post_COMBAT/ASD_SCZ_BD_normalized_CR_Cleaned.RData")
datExpr=datExpr[,rownames(datMeta)]
# expression_matrix
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



# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

datMeta=datMeta[order(datMeta$Diagnosis),]
MEs_LMSM=MEs_LMSM[rownames(datMeta),]
datMeta$SampleID=rownames(datMeta)
MEs_toPlot=MEs_LMSM[datMeta$SampleID,]
MEs_toPlot$Diagnosis=datMeta$Diagnosis[rownames(MEs_toPlot) %in% rownames(datMeta)]

pdf("Eigengene_Expression_LMSM_Modules_noctl.pdf")

barplot(MEs_toPlot$MEgreenyellow,col = MEs_toPlot$Diagnosis, main = "M1", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

barplot(MEs_toPlot$MEturquoise,col = MEs_toPlot$Diagnosis, main = "M2", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

barplot(MEs_toPlot$MEpurple,col = MEs_toPlot$Diagnosis, main = "M3", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

barplot(MEs_toPlot$MEmagenta,col = MEs_toPlot$Diagnosis, main = "M4", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

barplot(MEs_toPlot$MEpink,col = MEs_toPlot$Diagnosis, main = "M5", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

barplot(MEs_toPlot$MEbrown,col = MEs_toPlot$Diagnosis, main = "M6", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

barplot(MEs_toPlot$MEcyan,col = MEs_toPlot$Diagnosis, main = "M7", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

barplot(MEs_toPlot$MEred,col = MEs_toPlot$Diagnosis, main = "M8", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

barplot(MEs_toPlot$MEorange,col = MEs_toPlot$Diagnosis, main = "M9", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

barplot(MEs_toPlot$MEgreen,col = MEs_toPlot$Diagnosis, main = "M10", xlab = "Samples", ylab = "ME Expression", family = "serif")
legend("topright",legend = c("ASD","SCZ","BPD"),fill=c(1,2,3))

dev.off()

datMeta=datMeta[c("Diagnosis","ASD")]
colnames(datMeta)=c("Disorders","ASD vs SCZ BD")

moduleTraitCor = cor(MEs_LMSM, datMeta, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

tiff("Post_COMBAT/Figure_2.tiff",width = 8,height = 8,units = "in",res =200)

par(margin(4,4,4,4))
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datMeta),
               yLabels = c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10"),#names(MEs_LMSM)
               ySymbols = names(MEs_LMSM),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"), family="serif")

dev.off()

