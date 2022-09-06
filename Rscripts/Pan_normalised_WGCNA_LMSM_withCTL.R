library(WGCNA)
library(igraph)
library(SPONGE)
library(parallel)
library(doParallel)

library(SummarizedExperiment)
library(GSEABase)
library(flashClust)


library(PMA)
#library(genefu) #
library(varhandle)
library(broom)
#library(GSVA) #
library(pheatmap)
library(ggplot2)
library(survival)
library(SPONGE)
library(mldr) #
library(utiml) #
library(e1071)
library(corpcor)
library(miRspongeR)
library(WGCNA)
setwd("D:/CeRNA_Modules/PanData/")
#setwd("/home/rami/CeRNA_Modules/PanData")

load("ASD_SCZ_BD_normalized_CR_Cleaned.RData")

load("Data/ASD_SCZ_BD_CTL_normalized_CR_Cleaned.RData")
datMeta=datMeta[datMeta$Diagnosis!=0,]
datExpr=datExpr[,row.names(datMeta)]
save(datExpr,datProbes,datMeta,file = "Data/ASD_SCZ_BD_normalized_CR_Cleaned.RData")

# expression_matrix

mrna_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "protein_coding"],]
lncRNA_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "lncRNA"],]
miRNAProbes=datProbes[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!="",]

geneProbes=datProbes[datProbes$gene_biotype == "protein_coding" | datProbes$gene_biotype == "lncRNA",]

mir_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!=""],]

miRNAProbes$mirbase_id=gsub("mir","miR",miRNAProbes$mirbase_id)
row.names(mir_expr)=miRNAProbes$mirbase_id[miRNAProbes$external_gene_id %in% rownames(mir_expr)]

datExpr=t(rbind(mrna_expr,lncRNA_expr))

## NETWORK ANALYSIS
## ----------------
bsize = 5000
nSets = 1
powers = c(seq(1,9,by=1),seq(10,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()

pdf("WGCNA/WGCNA-softthresh.pdf")
par(mfrow=c(1,2))
n = 1
softThresh = pickSoftThreshold(data= datExpr, networkType = "signed", corFnc="bicor",verbose=5,powerVector=powers,blockSize = bsize)

sft = softThresh
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2")
abline(h=0.9, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
abline(h=100, col="red")

dev.off()


adjacency = adjacency(datExpr,type="signed", power = 18)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

save(dissTOM,file="WGCNA/Pan_Normalized_CRCleaned_TOM.Rdata")
load("WGCNA/Pan_Normalized_CRCleaned_TOM.Rdata")
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram) 
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)


# We like large modules, so we set the minimum module size relatively high:
minModule=50
# Module identification using dynamic tree cut: 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModule)
table(dynamicMods)

# Convert numeric lables into colors 
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath 
sizeGrWindow(8,6) 
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

# Calculate eigengenes 
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs); # Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result 
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")


# Save module colors and labels for use in subsequent parts

save(MEs, dynamicColors, geneTree, file ="Data_Preparation/Pan_normCRcleaned_moduleData.RData")
########

load("Data_Preparation/Pan_normCRcleaned_moduleData.RData")

mrna_expr=t(mrna_expr)
lncrna_expr=t(lncRNA_expr)
mir_expr=t(mir_expr)

############# loading Modules
source("LMSM/LMSM.R")

colorlevels <- unique(dynamicColors)
colorlevels <- colorlevels[-which(colorlevels=="grey")]

Modulegenes <- lapply(seq_len(length(colorlevels)), function(i)
  colnames(datExpr)[ which(dynamicColors==colorlevels[i]) ])

names(Modulegenes)=colorlevels


ceR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in% colnames(lncrna_expr))) )
mR_Num <- lapply(seq_along(Modulegenes), function(i) length(which(Modulegenes[[i]] %in% colnames(mrna_expr))) )

index <- which(ceR_Num >= 3 & mR_Num >=3)

CandidateModulegenes <- lapply(index, function(i) Modulegenes[[i]])
names(CandidateModulegenes) = names(Modulegenes)[index]

########## Targets

##### lncRNA interactions

lncbase=read.csv("D:/DataSource/miRNA_interactions/LncBasev2_download.csv",sep="\t")
lncbase_hsa=lncbase[(lncbase$species=="Homo sapiens" & lncbase$geneName %in% geneProbes$external_gene_id[geneProbes$gene_biotype == "lncRNA"] & lncbase$mirna %in% miRNAProbes$mirbase_id),][,c('mirna','geneName')]
lncbase_hsa=lncbase_hsa[!is.na(lncbase_hsa$geneName),]

npinter=read.csv("D:/DataSource/miRNA_interactions/NPInter/lncRNA_interaction.txt",sep="\t")
npinter_hsa=npinter[npinter$organism=="Homo sapiens" & npinter$ncName %in% geneProbes$external_gene_id[geneProbes$gene_biotype == "lncRNA"] & npinter$tarName %in% miRNAProbes$mirbase_id,][,c('tarName','ncName')]
colnames(npinter_hsa)=c('mirna','geneName')

mircode_lnc=read.csv("D:/DataSource/miRNA_interactions/miRcode/miRcode_lncRNA_miRNA_intrn.csv")
mircode_lnc_hsa=mircode_lnc[mircode_lnc$mirna %in% miRNAProbes$mirbase_id & mircode_lnc$geneName %in% geneProbes$external_gene_id[geneProbes$gene_biotype == "lncRNA"],]

lncInter_1=unique(rbind(lncbase_hsa,npinter_hsa))
lncInter=unique(rbind(lncInter_1,mircode_lnc_hsa))

###### mRNA interaction

mirtarbase=read.csv("D:/DataSource/miRNA_interactions/hsa_MTI.csv")
mirtarbase_hsa=mirtarbase[mirtarbase$miRNA %in% unique(lncInter$mirna) & mirtarbase$Target.Gene %in% geneProbes$external_gene_id[geneProbes$gene_biotype == "protein_coding"],][,c('miRNA','Target.Gene')]
colnames(mirtarbase_hsa)=c('mirna','geneName')

tarbase=read.csv("D:/DataSource/miRNA_interactions/TarBase_v8_download.txt",sep="\t")
tarbase_hsa=tarbase[tarbase$species == "Homo sapiens" & tarbase$mirna %in% unique(lncInter$mirna) & tarbase$geneName %in% geneProbes$external_gene_id[geneProbes$gene_biotype == "protein_coding"] ,][,c('mirna','geneName')]

tgtscn=read.csv("D:/DataSource/miRNA_interactions/TargetScan/Predicted_Targets_TargetScan.txt",sep='\t')
tgtscn_hsa=tgtscn[tgtscn$Gene.Tax.ID==9606 & tgtscn$miRNA %in% unique(lncInter$mirna) & tgtscn$Gene.Symbol %in% geneProbes$external_gene_id[geneProbes$gene_biotype == "protein_coding"],][,c('miRNA','Gene.Symbol')]
colnames(tgtscn_hsa)=c('mirna','geneName')

mircode_cdg=read.csv("D:/DataSource/miRNA_interactions/miRcode/miRcode_mRNA_miRNA_intrn.csv")
mircode_cdg_hsa=mircode_cdg[mircode_cdg$mirna %in% unique(lncInter$mirna) & mircode_cdg$geneName %in% geneProbes$external_gene_id[geneProbes$gene_biotype == "protein_coding"],]

mrnaInter_1=rbind(mirtarbase_hsa,tarbase_hsa)
mrnaInter_2=rbind(tgtscn_hsa,mircode_cdg_hsa)
mrnaInter=unique(rbind(mrnaInter_1,mrnaInter_2))


###### miR target
miRTarget=rbind(lncInter,mrnaInter)


num.cores=6
cl <- makeCluster(num.cores)
registerDoParallel(cl)
# LMSM method for miRNA sponge modules
CandidateModulegenes_WGCNA <- CandidateModulegenes
CommonmiRs_WGCNA <- share_miRs(mir_expr, lncrna_expr, mrna_expr, miRTarget, CandidateModulegenes_WGCNA)

LMSM_WGCNA <- LMSM(mir_expr, lncrna_expr, mrna_expr, miRTarget, CandidateModulegenes_WGCNA)

LMSM_WGCNA_Filter_modules <- LMSM_WGCNA[which(LMSM_WGCNA[, 3] >=3 & LMSM_WGCNA[, 5] < 0.05 & 
                                                LMSM_WGCNA[, 6] > 0.6), ]
LMSM_WGCNA_Modulegenes <- lapply(which(LMSM_WGCNA[, 3] >=3 & LMSM_WGCNA[, 5] < 0.05 & 
                                         LMSM_WGCNA[, 6] > 0.6), 
                                 function(i) CandidateModulegenes_WGCNA[[i]])
LMSM_WGCNA_CommonmiRs <- lapply(which(LMSM_WGCNA[, 3] >=3 & LMSM_WGCNA[, 5] < 0.05 & 
                                        LMSM_WGCNA[, 6] > 0.6), #
                                function(i) CommonmiRs_WGCNA[[i]])
rownames(LMSM_WGCNA_Filter_modules) <- names(LMSM_WGCNA_Modulegenes) <- names(LMSM_WGCNA_CommonmiRs) <- paste("LMSM", seq_along(LMSM_WGCNA_Modulegenes), sep=" ")


## Evaluate the significance of each LMSM module by using null model


num.cores=6
cl <- makeCluster(num.cores)
registerDoParallel(cl)
LMSM_WGCNA_precomputed_cov_matrices <- precomputed_cov_matrices
LMSM_WGCNA_null_model <- sponge_build_null_model(number_of_datasets = 1e+06, number_of_samples = 500, 
                                                 cov_matrices = LMSM_WGCNA_precomputed_cov_matrices,  
                                                 ks = seq(0.8, 0.9, 0.1), m_max = 1)
LMSM_WGCNA_modules <- data.frame(geneA = paste("ceR_module", seq(nrow(LMSM_WGCNA_Filter_modules))), 
                                 geneB = paste("mR_module", seq(nrow(LMSM_WGCNA_Filter_modules))), 
                                 df = replicate(nrow(LMSM_WGCNA_Filter_modules), 1), 
                                 cor = LMSM_WGCNA_Filter_modules[, 6], 
                                 pcor = LMSM_WGCNA_Filter_modules[, 9], 
                                 mscor = LMSM_WGCNA_Filter_modules[, 10])
LMSM_WGCNA_modules_p.values <- sponge_compute_p_values(sponge_result = LMSM_WGCNA_modules, null_model = LMSM_WGCNA_null_model)


LMSM_WGCNA_miRSponge <- lapply(seq(LMSM_WGCNA_Modulegenes), function(i) Extract.miRSponge(lncrna_expr, mrna_expr, LMSM_WGCNA_Modulegenes[[i]]))

## Understand predicted lncRNA-related miRNA sponge interactions, predicted and putative miRNA-target interactions of each LMSM module
LMSM_WGCNA_understand_miRSpongeTarget <- lapply(seq(LMSM_WGCNA_CommonmiRs), function(i) 
  Understand.miRSpongeTarget(lncrna_expr, mrna_expr, miRTarget, LMSM_WGCNA_CommonmiRs[[i]], 
                             LMSM_WGCNA_Modulegenes[[i]]))
LMSM_WGCNA_miRSpongeTarget <- do.call("rbind", LMSM_WGCNA_understand_miRSpongeTarget)
rownames(LMSM_WGCNA_miRSpongeTarget) <- names(LMSM_WGCNA_CommonmiRs)

LMSM_WGCNA_miRTarget <- lapply(seq(LMSM_WGCNA_CommonmiRs), function(i) Extract.miRTarget(LMSM_WGCNA_CommonmiRs[[i]], LMSM_WGCNA_Modulegenes[[i]]))

## miRNAs distribution in LMSM modules
LMSM_WGCNA_miR_distribution <- miR.distribution(LMSM_WGCNA_CommonmiRs)

save(LMSM_WGCNA_miRSponge,LMSM_WGCNA_miR_distribution,LMSM_WGCNA_miRSpongeTarget,LMSM_WGCNA_Filter_modules,LMSM_WGCNA_modules_p.values,LMSM_WGCNA_CommonmiRs,LMSM_WGCNA_Modulegenes,file="WGCNA/LMSM_ASD_SCZ_BPD_noctl.RData")

write.csv(LMSM_WGCNA_miRSpongeTarget,"LMSM_Pan_norm_noctl_miRSpongeTarget.csv",quote = F)

