library(WGCNA);library(ggplot2); library(reshape); library(igraph); library(RColorBrewer); library(WGCNA)
setwd("D:/CeRNA_Modules/PanData/Post_COMBAT/")
load("FinalizedNetwork.RData")


## Make module eigengene-MDS plot
eigmat = MEs_LMSM
eigmat = eigmat[,paste0("ME", unique(LMSM_colors))]
adj = bicor(eigmat)
mds = cmdscale(dist(t(eigmat)), eig = T);   
mds$points[,1] = -1* mds$points[,1]
g1 <- graph.adjacency(as.matrix(adj),mode="undirected",weighted=T,diag=FALSE)
layoutFR <- mds$points
c = paste0("M", 1:10)

pdf("ME_MDSPlot.pdf",height=4,width=4,useDingbats=FALSE)
edgecolors = numbers2colors(E(g1)$weight, colors = redWhiteGreen(9, gamma=3), signed=T, centered=T, lim=c(-1,1))
plot.igraph(g1, vertex.label = c,
            vertex.label.dist=0.3, 
            vertex.size=6,
            vertex.label.color="black",
            vertex.label.family = "sans",
            vertex.color = unique(LMSM_colors),
            vertex.label.cex=0.6,
            layout=layoutFR,
            edge.color=edgecolors,
            edge.width=2,asp=1, main="Module Eigengene MDS")
labs = seq(1,-1,by=-.25)
str = paste(labs,sep="\n")
#text(-1.25,-1, labels = paste(labs,collapse='\n'),pos = 4,cex = .5)
p = matrix(NA,nrow=9,ncol=4)
p[,1]=-1.35; p[,2]=-1.3
p[,3]=p[,4] = -.6-.6*seq(0,1,by=.12)
for(i in 1:9) {
  lines(x=p[i,1:2],y=p[i,3:4], lwd = 2, col=numbers2colors(as.numeric(labs[i]), colors = redWhiteGreen(9, gamma=3), signed=T, centered=T, lim=c(-1,1)))
  text(x=-1.3,y=p[i,3],labs[i],cex=.4,pos=4)
}
dev.off()