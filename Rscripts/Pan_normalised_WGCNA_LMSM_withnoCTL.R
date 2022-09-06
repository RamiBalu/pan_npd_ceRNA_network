library(dplyr)
library(ggplot2)
library(WGCNA)

setwd("D:/1_CeRNA_Modules/PanData/Post_COMBAT/")
load("Data/ASD_SCZ_BD_normalized_CR_Cleaned.RData")

# expression_matrix

mrna_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "protein_coding"],]
lncRNA_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "lncRNA"],]
miRNAProbes=datProbes[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!="",]

geneProbes=datProbes[datProbes$gene_biotype == "protein_coding" | datProbes$gene_biotype == "lncRNA",]

mir_expr=datExpr[datProbes$external_gene_id[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!=""],]

miRNAProbes$mirbase_id=gsub("mir","miR",miRNAProbes$mirbase_id)
row.names(mir_expr)=miRNAProbes$mirbase_id[miRNAProbes$external_gene_id %in% rownames(mir_expr)]

datExpr=t(rbind(mrna_expr,lncRNA_expr))

save(datExpr,datProbes,datMeta,mir_expr,file="Data/ASD_SCZ_BD_normalized.RData")

## NETWORK ANALYSIS
## ----------------
bsize = 5000
nSets = 1
powers = c(seq(1,9,by=1),seq(10,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()

pdf("WGCNA-softthresh_noCTL.pdf")
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

save(dissTOM,file="Pan_Normalized_CRCleaned_TOM.Rdata")
load("Pan_Normalized_CRCleaned_TOM.Rdata")
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

load("Pan_normCRcleaned_moduleData.RData")

mrna_expr=t(mrna_expr)
lncrna_expr=t(lncRNA_expr)
mir_expr=t(mir_expr)

############# loading Modules
source("../LMSM/LMSM.R")

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
miRTarget=rbind(lncInter,mrnaInter)

colnames(lncInter)=c("miRNA","lncRNA")
colnames(mrnaInter)=c("miRNA","mRNA")
putative_interactions=merge.data.frame(lncInter,mrnaInter,by="miRNA")
write.csv(putative_interactions,"Putative_Interactions.csv",row.names = F,quote = F)
putative_interactions=read.csv("Putative_Interactions.csv")
###### miR target


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

save(LMSM_WGCNA_miRSponge,LMSM_WGCNA_miR_distribution,LMSM_WGCNA_miRSpongeTarget,LMSM_WGCNA_Filter_modules,LMSM_WGCNA_modules_p.values,LMSM_WGCNA_CommonmiRs,LMSM_WGCNA_Modulegenes,file="LMSM_ASD_SCZ_BPD_noctl.RData")

write.csv(LMSM_WGCNA_miRSpongeTarget,"LMSM_Pan_norm_noctl_miRSpongeTarget.csv",quote = F)

LMSM_sign_modulegenes=c()
i=1
for(s in sign_module_index){
  m=paste(c("M",i),collapse = "-")
  write.table(LMSM_WGCNA_Modulegenes[[s]],paste(c("Modules/",m,".txt"),collapse = ""),row.names = F,col.names = F,quote = F)
  i=i+1
  LMSM_sign_modulegenes[[i]]=LMSM_WGCNA_Modulegenes[[s]]
}
sign_Modulegenes=unlist(LMSM_sign_modulegenes)

sign_module_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj
                                                                  < 0.05,]))

write.csv(LMSM_WGCNA_miRSpongeTarget,"LMSM_Pan_norm_noctl_miRSpongeTarget.csv",quote = F)
write.csv(LMSM_WGCNA_modules_p.values,"LMSM_Pan_norm_miRSpongeP_value.csv",row.names = F,quote = F)


#######----------------- Putative_Interactions CeRNA network
###
load("LMSM_ASD_SCZ_BPD_noctl.RData")
sign_module_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj < 0.05,]))

sponge_int=read.csv("SPONGE/Pan_norm_batch_crcted_CeRNA_interactions.csv")
sponge_int=sponge_int[sponge_int$p.val<0.05,][c("geneA","geneB")]
#sponge_int=sponge_int[sponge_int$p.val<0.05,][c("geneA","geneB")]
colnames(sponge_int)=c("lncRNA","mRNA")
library(dplyr)
cmn_sponge_int=c()
all_int=c()
put_int=putative_interactions[c("lncRNA","mRNA")]
lncRNAs=c()
all_put_int_mod=c()
put_int_mod=c()
i=1
for (s in sign_module_index){
  int_mod=data.frame(LMSM_WGCNA_miRSponge[[s]])
  colnames(int_mod)=colnames(put_int)
  all_int=rbind(all_int,int_mod)
  lncRNAs[[i]]=int_mod$lncRNA
  put_int_mod[[i]]=unique(inner_join(put_int,int_mod))
  cmn_sponge_int=rbind(cmn_sponge_int,inner_join(sponge_int,int_mod))
  all_put_int_mod=rbind(all_put_int_mod,put_int_mod[[i]])
  i=i+1
}
all_put_int_mod=dplyr::bind_rows(put_int_mod, .id="Module")
lncRNAs=unique(unlist(lncRNAs))
write.csv((all_put_int_mod),"Putative_Interactions_miRNASpongeModules.csv",row.names = F,quote = F)

write.table(lncRNAs,"LncRNAs.txt",col.names = F,row.names = F,quote = F)

### DEGs

load("../Results/DEGs_from_SNM_TWAS_Gandal_2018.RData")

asd_de_lncrna=asd_degs[asd_degs$geneBiotype=="lncRNA",][c("geneName","Direction")]
scz_de_lncrna=scz_degs[scz_degs$geneBiotype=="lncRNA",][c("geneName","Direction")]
bpd_de_lncrna=bpd_degs[bpd_degs$geneBiotype=="lncRNA",][c("geneName","Direction")]

asd_de_mrna=asd_degs[asd_degs$geneBiotype=="protein_coding",][c("geneName","Direction")]
scz_de_mrna=scz_degs[scz_degs$geneBiotype=="protein_coding",][c("geneName","Direction")]
bpd_de_mrna=bpd_degs[bpd_degs$geneBiotype=="protein_coding",][c("geneName","Direction")]

colnames(bpd_de_lncrna) <- colnames(scz_de_lncrna) <- colnames(asd_de_lncrna) <- c("lncRNA","Direction")
colnames(bpd_de_mrna) <- colnames(scz_de_mrna) <- colnames(asd_de_mrna) <- c("mRNA","Direction")

### write_to_csv
load("LMSM_ASD_SCZ_BPD_noctl.RData")
sign_module_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj
                                                                  < 0.05,]))
encori_lncRNA=read.csv("Known_CeRNAs.csv")
colnames(encori_lncRNA)=c("lncRNA","mRNA")

i=1
for (s in sign_module_index){
  int_mod=data.frame(LMSM_WGCNA_miRSponge[[s]])
  colnames(int_mod)=c("lncRNA","mRNA")
  
  int_mod=merge(int_mod,asd_de_lncrna,by="lncRNA",all.x = TRUE)
  int_mod=merge(int_mod,scz_de_lncrna,by="lncRNA",all.x = TRUE)
  int_mod=merge(int_mod,bpd_de_lncrna,by="lncRNA",all.x = TRUE)
  
  colnames(int_mod) <- c("lncRNA","mRNA","ASD_DElncRNA","SCZ_DElncRNA","BPD_DElncRNA")
  
  int_mod=merge(int_mod,asd_de_mrna,by="mRNA",all.x = TRUE)
  int_mod=merge(int_mod,scz_de_mrna,by="mRNA",all.x = TRUE)
  int_mod=merge(int_mod,bpd_de_mrna,by="mRNA",all.x = TRUE)
  
  colnames(int_mod) <- c("mRNA","lncRNA","ASD_DElncRNA","SCZ_DElncRNA","BPD_DElncRNA","ASD_DEmRNA","SCZ_DEmRNA","BPD_DEmRNA")
  
  cmn_int=inner_join(put_int,int_mod)
  cmn_int$putative=rep(1,dim(cmn_int)[1])
  int_mod=merge(int_mod,cmn_int[c("lncRNA","mRNA","putative")],by=c("lncRNA","mRNA"),all.x = TRUE)
  
  cmn_int=inner_join(encori_lncRNA,int_mod)
  cmn_int$encori=rep(1,dim(cmn_int)[1])
  int_mod=merge(int_mod,cmn_int[c("lncRNA","mRNA","encori")],by=c("lncRNA","mRNA"),all.x = TRUE)
  int_mod[is.na(int_mod)] <- "-"
  write.csv(int_mod, file=paste("ModularInteractions/M",i,"_Interactions.csv",sep=""),row.names = F)
  i=i+1
}

### RNAs KME_Table

kme=read.csv("kME_table.csv",row.names = 1)[,-c(4:7,9:11,13:23)]
kme_lncRNA=kme[kme$gene_biotype=="lncRNA",]
kme_mRNA=kme[kme$gene_biotype=="protein_coding",]
colnames(bpd_de_lncrna) <- colnames(scz_de_lncrna) <- colnames(asd_de_lncrna) <- 
  colnames(bpd_de_mrna) <- colnames(scz_de_mrna) <- colnames(asd_de_mrna) <- c("external_gene_id","Direction")

kme_lncRNA=merge(kme_lncRNA,asd_de_lncrna,by="external_gene_id",all.x = TRUE)
kme_lncRNA=merge(kme_lncRNA,scz_de_lncrna,by="external_gene_id",all.x = TRUE)
kme_lncRNA=merge(kme_lncRNA,bpd_de_lncrna,by="external_gene_id",all.x = TRUE)

kme_mRNA=merge(kme_mRNA,asd_de_mrna,by="external_gene_id",all.x = TRUE)
kme_mRNA=merge(kme_mRNA,scz_de_mrna,by="external_gene_id",all.x = TRUE)
kme_mRNA=merge(kme_mRNA,bpd_de_mrna,by="external_gene_id",all.x = TRUE)

kme_summary=rbind.data.frame(kme_mRNA, kme_lncRNA)
colnames(kme_summary) <- colnames(kme_mRNA) <- colnames(kme_lncRNA) <- c(colnames(kme_summary)[1:5],"ASD_DEG","SCZ_DEG","BD_DEG")


### ROCs

asd_rocs=read.csv("D:/CeRNA_Modules/Validation/ROC_validation/ASD/AUCs_of_ASD_Data.csv")
scz_rocs=read.csv("D:/CeRNA_Modules/Validation/ROC_validation/SCZ/AUCs_of_SCZ_data.csv")
bpd_rocs=read.csv("D:/CeRNA_Modules/Validation/ROC_validation/BPD/AUCs_of_BPD_data.csv")
colnames(asd_rocs)[1] <- colnames(scz_rocs)[1] <- colnames(bpd_rocs)[1] <- "dot"


lncrna_dot=read.csv("lncRNA_dot.csv",header = T)
row.names(lncrna_dot)=lncrna_dot$dot

asd_rocs=merge(asd_rocs,lncrna_dot,by="dot",all = T)
scz_rocs=merge(scz_rocs,lncrna_dot,by="dot",all = T)
bpd_rocs=merge(bpd_rocs,lncrna_dot,by="dot",all = T)

asd_rocs=asd_rocs[-c(1)]
scz_rocs=scz_rocs[-c(1)]
bpd_rocs=bpd_rocs[-c(1)]



kme_lncRNA_summary=merge(kme_lncRNA,asd_rocs,by="external_gene_id",all.x = T)
kme_lncRNA_summary=merge(kme_lncRNA_summary,scz_rocs,by="external_gene_id",all.x = T)
kme_lncRNA_summary=merge(kme_lncRNA_summary,bpd_rocs,by="external_gene_id",all.x = T)

kme_summary[is.na(kme_summary)]<- 
  kme_lncRNA_summary[is.na(kme_lncRNA_summary)] <- "-"
write.csv(kme_summary,"Summary_ceRNAs.csv",row.names = F)
write.csv(kme_lncRNA_summary,"Summary_lncRNAs.csv",row.names = F)