library(dplyr)
library(ggplot2)
library(WGCNA)

setwd("D:/1_CeRNA_Modules/PanData/Post_COMBAT/")
## expression_matrix
load("Data/Pan_NPD/ASD_SCZ_BD_normalized.RData")

## NETWORK ANALYSIS parameters
bsize = 5000
nSets = 1
powers = c(seq(1,9,by=1),seq(10,30,by=2))
enableWGCNAThreads()
allowWGCNAThreads()

if(FALSE) { 
  # compute soft threshold for the network analysis
  pdf("WGCNA-softthresh.pdf")
  par(mfrow=c(1,2))
  n = 1
  softThresh = pickSoftThreshold(data= datExpr, networkType = "signed", corFnc="bicor",
                                 verbose=5,powerVector=powers,blockSize = bsize)
  
  sft = softThresh
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Thresh Power", ylab="Scale free R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       labels = powers, cex = 0.7, col="red",  xlab="Soft Thresh Power", ylab="Scale free R^2")
  abline(h=0.9, col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold power", ylab = "Mean connectivity", type = "n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.7, col="black")
  abline(h=100, col="red")
  
  dev.off()
}

##WGCNA

adjacency = adjacency(datExpr,type="signed", power = 18)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

#hierarchical clustering 
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram) 
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)


# We like large modules, so we set the minimum module size relatively high:
minModule=50
# Module identification using dynamic tree cut: 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = minModule)
table(dynamicMods)

# Convert numeric lables into colors 
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath 
sizeGrWindow(8,6) 
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

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

save(MEs, dynamicColors, geneTree, file ="Data/Pan_NPD/Pan_npd_moduleData.RData")
