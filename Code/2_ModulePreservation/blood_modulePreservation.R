setwd("D:/1_CeRNA_Modules/PanData/Post_COMBAT/")
library(WGCNA)
source("Code/2_ModulePreservation/preservation_plot.R")


###### ASD input
load("Data/Validation/Bloodbased/ASD_Blood_GSE18123_normalized.RData")
asd_datExpr=datExpr

###### Module information
load("Data/Pan_NPD/mod_presev_ref_datExpr.RData")

dynamicColors=colorLabels
load("Data/Pan_NPD/LMSM_ASD_SCZ_BPD.RData")
sign_module_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj<0.05,]))
cerna_modules=c()
i=1
for (s in sign_module_index){
  cerna_modules[[i]]=LMSM_WGCNA_Modulegenes[[s]]
  i=i+1
}
cerna_modulegenes=unlist(cerna_modules)
colorLabels=dynamicColors[intersect(cerna_modulegenes,intersect(colnames(asd_datExpr)
                                                                ,colnames(LMSM_Modulegenes_expr)))]

nSets = 2
multiExpr = list()


multiExpr[[1]] = list(data = LMSM_Modulegenes_expr[,names(colorLabels)])
multiExpr[[2]] = list(data = asd_datExpr[,names(colorLabels)])


setLabels = c("pan_npd", "blood")
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)
multiColor = list(pan_npd = colorLabels)
system.time( { mp = modulePreservation(multiExpr, multiColor, referenceNetworks = 1,nPermutations = 200, 
                                       randomSeed = 1, quickCor = 0, verbose = 3) } )

plot_preservation(mp)


##--SCZ Module preservation

##----- SCZ input
load("Data/Validation/Bloodbased/SCZ_Blood_GSE165604_normalizedData.RData")
scz_datExpr=data.frame(t(scz_datExpr))

load("Data/Pan_NPD/mod_presev_ref_datExpr.RData")

dynamicColors=colorLabels
sign_module_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj<0.05,]))
cerna_modules=c()
i=1
for (s in sign_module_index){
  cerna_modules[[i]]=LMSM_WGCNA_Modulegenes[[s]]
  i=i+1
}
cerna_modulegenes=unlist(cerna_modules)
colorLabels=dynamicColors[intersect(cerna_modulegenes,intersect(colnames(scz_datExpr)
                                                                ,colnames(LMSM_Modulegenes_expr)))]
nSets = 2
multiExpr = list()


multiExpr[[1]] = list(data = LMSM_Modulegenes_expr[,names(colorLabels)])
multiExpr[[2]] = list(data = scz_datExpr[,names(colorLabels)])


setLabels = c("pan_npd", "blood")
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)
multiColor = list(pan_npd = colorLabels)
system.time( { mp = modulePreservation(multiExpr, multiColor, referenceNetworks = 1,nPermutations = 200, 
                                       randomSeed = 1, quickCor = 0, verbose = 3) } )

plot_preservation(mp)



##----- BD input


load("Data/Validation/Bloodbased/BPD_Blood_GSE124326_normalizedData.RData")

bpd_datExpr=t(bpd_datExpr)
load("Data/Pan_NPD/mod_presev_ref_datExpr.RData")

dynamicColors=colorLabels
sign_module_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj<0.05,]))
cerna_modules=c()
i=1
for (s in sign_module_index){
  cerna_modules[[i]]=LMSM_WGCNA_Modulegenes[[s]]
  i=i+1
}
cerna_modulegenes=unlist(cerna_modules)
colorLabels=dynamicColors[intersect(cerna_modulegenes,intersect(colnames(bpd_datExpr)
                                                                ,colnames(LMSM_Modulegenes_expr)))]
nSets = 2
multiExpr = list()


multiExpr[[1]] = list(data = LMSM_Modulegenes_expr[,names(colorLabels)])
multiExpr[[2]] = list(data = bpd_datExpr[,names(colorLabels)])


setLabels = c("pan_npd", "blood")
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)
multiColor = list(pan_npd = colorLabels)
system.time( { mp = modulePreservation(multiExpr, multiColor, referenceNetworks = 1,nPermutations = 200, 
                                       randomSeed = 1, quickCor = 0, verbose = 3) } )


plot_preservation(mp)