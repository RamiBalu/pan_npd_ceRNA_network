library(dplyr)
library(ggplot2);library(sva); library(WGCNA)

setwd("D:/CeRNA_Modules/PanData/Data_Preparation/")
load("D:/NPD_Data/ASD/Data/ASD_4region_normalized_CR_cleaned.Rdata")
asd_datMeta=datMeta[datMeta$Diagnosis_=="ASD",]
asd_datExpr=datExpr[rownames(asd_datMeta)]
asd_datProbes=datProbes

load("D:/NPD_Data/GVEX/Data/GVEX_SCZ_BPD_Normalized_CR_cleaned.RData")
gvex_datMeta=datMeta[datMeta$Diagnosis!="Control",]
gvex_datExpr=datExpr[rownames(gvex_datMeta)]
gvex_datProbes=datProbes

length(intersect(rownames(asd_datExpr),rownames(gvex_datExpr)))


######MetaData

colnames(gvex_datMeta)
colnames(asd_datMeta)

gvex_datMeta$Study=rep("GVEX",length(gvex_datMeta$Individual_ID..RNAseq.library.BID.))
asd_datMeta$Study=rep("ASD",length(asd_datMeta$X))

gvex_datMeta$RNAseq_ID=rownames(gvex_datMeta)
gvex_datTrait=gvex_datMeta[,c("Study","BrainWeight","PMI","pH","Diagnosis","AgeDeath","Sex")]

asd_datTrait=asd_datMeta[,c("Study","Brain_Weight","PMI","pH","Diagnosis_","Age","Sex")]
colnames(asd_datTrait)=c("Study","BrainWeight","PMI","pH","Diagnosis","AgeDeath","Sex")

datTrait=rbind.data.frame(gvex_datTrait,asd_datTrait)

datTrait$Diagnosis[datTrait$Diagnosis=="ASD"]=1
datTrait$Diagnosis[datTrait$Diagnosis=="SCZ"]=2
datTrait$Diagnosis[datTrait$Diagnosis=="BP"]=3
datTrait$Diagnosis[datTrait$Diagnosis=="Control" | datTrait$Diagnosis=="CTL"]=0


datTrait$BrainWeight=as.numeric(datTrait$BrainWeight)
datTrait$PMI=as.numeric(datTrait$PMI)
datTrait$pH=as.numeric(datTrait$pH)
datTrait$AgeDeath=as.numeric(datTrait$AgeDeath)
datTrait$Diagnosis=as.numeric(datTrait$Diagnosis)
datTrait$Sex=as.numeric(datTrait$Sex)
datTrait$Study=as.factor(datTrait$Study)

datTrait$ASD[datTrait$Diagnosis==1]=1
datTrait$SCZ[datTrait$Diagnosis==2]=1
datTrait$BPD[datTrait$Diagnosis==3]=1

datTrait$ASD[datTrait$Diagnosis!=1]=0
datTrait$SCZ[datTrait$Diagnosis!=2]=0
datTrait$BPD[datTrait$Diagnosis!=3]=0

datMeta=datTrait


#####################################################

#####################################################

library(ggplot2);library(sva); library(WGCNA)

multiExpr = vector(mode="list",length = 2)
multiExpr[[1]]$datExpr= asd_datExpr
multiExpr[[1]]$datMeta= asd_datMeta
multiExpr[[2]]$datExpr= gvex_datExpr
multiExpr[[2]]$datMeta= gvex_datMeta

genes = rownames(multiExpr[[1]]$datExpr)
for(i in 2:length(multiExpr)) genes = intersect(genes, rownames(multiExpr[[i]]$datExpr))

all_datExpr = data.frame(row.names = genes)

for(i in 1:length(multiExpr)) {
all_datExpr = cbind(all_datExpr, multiExpr[[i]]$datExpr[match(genes, rownames(multiExpr[[i]]$datExpr)),])}

####### 

##QC Pre-Combat
sex_col = numbers2colors(datMeta$Sex, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$Sex, na.rm=T),max(datMeta$Sex, na.rm=T)))
age_col = numbers2colors(datMeta$AgeDeath, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$AgeDeath, na.rm=T),max(datMeta$AgeDeath, na.rm=T)))
pmi_col = numbers2colors(datMeta$PMI, blueWhiteRed(100), signed=F, centered=T, lim=c(min(datMeta$PMI, na.rm=T),max(datMeta$PMI, na.rm=T)))

plot(density(all_datExpr[,1]), xlim=c(-5,20), ylim=c(0, 0.5), col = as.numeric(datMeta$Study[1]), xlab="Intensity (log2)", ylab="Density", main="Mega-Analysis: Pre-Combat")
for(i in 2:dim(all_datExpr)[[2]])
  lines(density(all_datExpr[,i]), xlim=c(0,20), col = as.numeric(datMeta$Study[i]))  
legend("topleft", (levels(datMeta$Study)), col=c(1:8), pch=16, cex=0.5)

par(mfrow=c(1,1))
mds = cmdscale(dist(t(all_datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Study)), pch=20, main="Multidimensional Scaling Plot\nPre-ComBat", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
legend("topleft", (levels(datMeta$Study)), col=c(1:8), pch=16, cex=0.5)

tree = hclust(dist(t(all_datExpr)), method = "average")
par(mfrow=c(1,1))
plotDendroAndColors(tree, cbind(as.numeric(datMeta$Diagnosis), as.numeric(datMeta$Study), sex_col, age_col, pmi_col), 
                    groupLabels = c("Group", "Study", "Sex", "Age", "pH", "PMI", "RIN", "RNA"), cex.colorLabels=0.8, cex.dendroLabels=0.15,
                    main="Dendrogram\nPre-Combat")



#Normalize by Study
mod = model.matrix(~Diagnosis+AgeDeath+Sex, data=datMeta)
batch = as.factor(datMeta$Study)
datExpr = ComBat(all_datExpr, batch=batch, mod=mod, prior.plots = F)

pdf("QC_postcombat_withno_CTL.pdf")

##QC - PostCombat
par(mfrow=c(2,2))
plot(density(datExpr[,1]),xlim=c(0,17), ylim=c(0, 0.4), col = as.numeric(datMeta$Study[1]), xlab="", ylab="", main="")
for(i in 2:dim(datExpr)[[2]])
  lines(density(datExpr[,i]), xlim=c(0,16), ylim=c(0,0.3), col = as.numeric(datMeta$Study[i]))  
legend("topright", levels(datMeta$Study), col=c(1:8), pch=16,cex=0.7)

#MDS Plot
mds = cmdscale(dist(t(datExpr)), eig = T);   pc1 = mds$eig[1]^2 / sum(mds$eig^2);   pc2 = mds$eig[2]^2 / sum(mds$eig^2)
plot(mds$points, col=as.numeric(as.factor(datMeta$Study)), pch=20, main="MDS: Study", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC2 (", signif(100*pc2,3),"%)",sep="")); 
plot(mds$points, col=as.numeric(as.factor(datMeta$Diagnosis)), pch=16, main="MDS: Group",  xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
plot(mds$points, col=sex_col, pch=16, main="MDS - Sex", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 
plot(mds$points, col=age_col, pch=16, main="MDS - Age", xlab = paste("PC1 (", signif(100*pc1,3), "%)", sep=""), ylab = paste("PC1 (", signif(100*pc2,3),"%)",sep="")); 

#Dendrogram
par(mfrow=c(1,1))
tree = hclust(dist(t(datExpr)), method = "average")
plotDendroAndColors(tree, cbind(as.numeric(datMeta$Diagnosis), as.numeric(datMeta$Study), sex_col, age_col), 
                    groupLabels = c("Group", "Study", "Sex", "Age"), cex.colorLabels=0.6, cex.dendroLabels=0.2,
                    main="Dendrogram\nPost-Combat")

dev.off()

datProbes=intersect(asd_datProbes,gvex_datProbes)
lncRNA=c("processed_transcript","lincRNA","antisense","sense_overlapping","sense_intronic","3prime_overlapping_ncrna")
datProbes$gene_biotype[datProbes$gene_biotype %in% lncRNA] = "lncRNA"

gene_freq=as.data.frame(table(datProbes$external_gene_id))
datProbes=datProbes[datProbes$external_gene_id %in% gene_freq$Var1[gene_freq$Freq==1],]

datExpr=datExpr[rownames(datProbes),]
rownames(datProbes)=datProbes$external_gene_id
rownames(datExpr)=rownames(datProbes)
save(file="ASD_SCZ_BD_normalized_CR_Cleaned.RData",datExpr,datProbes,datMeta)
