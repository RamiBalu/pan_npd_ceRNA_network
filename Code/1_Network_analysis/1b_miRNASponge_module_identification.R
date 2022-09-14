library(dplyr)
library(ggplot2)
library(WGCNA)
library(igraph)
library(SPONGE)
library(parallel)
library(doParallel)

########
setwd("D:/1_CeRNA_Modules/PanData/Post_COMBAT/")

## expression_matrix
load("Data/Pan_NPD/ASD_SCZ_BD_normalized.RData")
load("Data/Pan_NPD/Pan_npd_moduleData.RData")
geneProbes=datProbes[datProbes$gene_biotype == "protein_coding" | datProbes$gene_biotype == "lncRNA",]
miRNAProbes=datProbes[datProbes$gene_biotype == "miRNA" & datProbes$mirbase_id!="",]

# expression_matrix

mrna_expr=datExpr[,datProbes$external_gene_id[datProbes$gene_biotype == "protein_coding"]]
lncRNA_expr=datExpr[,datProbes$external_gene_id[datProbes$gene_biotype == "lncRNA"]]
mir_expr=t(mir_expr);
miRNAProbes$mirbase_id=gsub("mir","miR",miRNAProbes$mirbase_id)


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

lncInter=read.csv('Data/lncRNA-miRNA_interactions.csv')
mrnaInter=read.csv('Data/mRNA-miRNA_interactions.csv')
miRTarget=rbind(lncInter,mrnaInter)

putative_interactions=read.csv("Data/Putative_Interactions.csv")
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

save(LMSM_WGCNA_miRSponge,LMSM_WGCNA_miR_distribution,LMSM_WGCNA_miRSpongeTarget,LMSM_WGCNA_Filter_modules,
     LMSM_WGCNA_modules_p.values,LMSM_WGCNA_CommonmiRs,LMSM_WGCNA_Modulegenes,file="Data/Pan_NPD/LMSM_ASD_SCZ_BPD.RData")


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

write.csv(LMSM_WGCNA_miRSpongeTarget,"Results/LMSM_Pan_norms_miRSpongeTarget.csv",quote = F)
write.csv(LMSM_WGCNA_modules_p.values,"Results/LMSM_Pan_norm_miRSpongeP_value.csv",row.names = F,quote = F)

