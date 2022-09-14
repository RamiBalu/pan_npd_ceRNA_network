library(igraph)
library(SPONGE)
library(parallel)
library(doParallel)

library(SummarizedExperiment)
library(GSEABase)
library(WGCNA)
library(flashClust)

library(WGCNA)
library(PMA)
library(genefu) #
library(varhandle)
library(broom)
library(GSVA) #
library(pheatmap)
library(ggplot2)
library(survival)
library(SPONGE)
library(mldr) #
library(utiml) #
library(e1071)
library(corpcor)
library(miRspongeR)


setwd("D:/ExpData/GVEX/")
source("LMSM/LMSM.R")
load("Data/SCZ_input.RData")
mirnaProbes=read.csv("Data/miRNAProbes.csv")
geneProbes=read.csv("Data/GeneProbes.csv")


lncbase=read.csv("D:/DataSource/miRNA_interactions/LncBasev2_download.csv",sep="\t")
lncbase_hsa=lncbase[(lncbase$species=="Homo sapiens" & lncbase$mirna %in% mirnaProbes$mirbase_id),][,c('mirna','geneName')]
lncbase_hsa=lncbase_hsa[!is.na(lncbase_hsa$geneName),]


mirtarbase=read.csv("D:/DataSource/miRNA_interactions/hsa_MTI.csv")
mirtarbase_hsa=mirtarbase[mirtarbase$miRNA %in% unique(lncbase_hsa$mirna),][,c('miRNA','Target.Gene')]
colnames(mirtarbase_hsa)=c('mirna','geneName')

miRTarget=unique(rbind(mirtarbase_hsa,lncbase_hsa))

scz_mRNA_expr=t(scz_gene_expr[geneProbes$external_gene_id[geneProbes$gene_biotype=="protein_coding"],])
scz_lncRNA_expr=t(scz_gene_expr[geneProbes$external_gene_id[geneProbes$gene_biotype=="lncRNA"],])


miRNA_Exp_SummarizedExperiment=t(scz_mir_expr)
miRNA_Exp_SummarizedExperiment<- SummarizedExperiment(assays=list(miRNA_Exp_SummarizedExperiment = miRNA_Exp_SummarizedExperiment))

lncRNA_Exp_SummarizedExperiment=scz_mRNA_expr
lncRNA_Exp_SummarizedExperiment<- SummarizedExperiment(assays=list(lncRNA_Exp_SummarizedExperiment = lncRNA_Exp_SummarizedExperiment))


mRNA_Exp_SummarizedExperiment=scz_lncRNA_expr
mRNA_Exp_SummarizedExperiment<- SummarizedExperiment(assays=list(mRNA_Exp_SummarizedExperiment = mRNA_Exp_SummarizedExperiment))

miRTarget_lncR_vs_mR_SummarizedExperiment <- miRTarget
miRTarget_lncR_vs_mR_SummarizedExperiment <- SummarizedExperiment(assays=list(miRTarget_lncR_vs_mR_SummarizedExperiment = miRTarget_lncR_vs_mR_SummarizedExperiment))




num.cores=6
cl <- makeCluster(num.cores)
registerDoParallel(cl)
# LMSM method for miRNA sponge modules
CandidateModulegenes_WGCNA <- module_WGCNA(scz_lncRNA_expr, scz_mRNA_expr
                                              , RsquaredCut = 0.8)
CommonmiRs_WGCNA <- share_miRs(t(scz_mir_expr), scz_lncRNA_expr, scz_mRNA_expr, miRTarget, CandidateModulegenes_WGCNA)

LMSM_WGCNA <- LMSM(t(scz_mir_expr), scz_lncRNA_expr, scz_mRNA_expr, miRTarget, CandidateModulegenes_WGCNA)

LMSM_WGCNA_Filter_modules <- LMSM_WGCNA[which(LMSM_WGCNA[, 3] >=3 & LMSM_WGCNA[, 5] < 0.05 & 
                                                LMSM_WGCNA[, 6] > 0.8 & LMSM_WGCNA[, 10] > 0.05), ]
LMSM_WGCNA_Modulegenes <- lapply(which(LMSM_WGCNA[, 3] >=3 & LMSM_WGCNA[, 5] < 0.05 & 
                                         LMSM_WGCNA[, 6] > 0.8 & LMSM_WGCNA[, 10] > 0.05), 
                                 function(i) CandidateModulegenes_WGCNA[[i]])
LMSM_WGCNA_CommonmiRs <- lapply(which(LMSM_WGCNA[, 3] >=3 & LMSM_WGCNA[, 5] < 0.05 & 
                                        LMSM_WGCNA[, 6] > 0.8 & LMSM_WGCNA[, 10] > 0.05), #
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

expr <- t(cbind(scz_lncRNA_expr, scz_mRNA_expr))
SCZ_gsva_es <- gsva(expr, LMSM_WGCNA_Modulegenes, mx.diff = FALSE)
rownames(SCZ_gsva_es) <- paste("LMSM", seq_along(LMSM_WGCNA_Modulegenes), sep=" ")


## Extract lncRNA-related miRNA sponge interactions of each LMSM module
LMSM_WGCNA_miRSponge <- lapply(seq(LMSM_WGCNA_Modulegenes), function(i) Extract.miRSponge(scz_lncRNA_expr, scz_mRNA_expr, LMSM_WGCNA_Modulegenes[[i]]))

## Extract miRNA-target interactions of each LMSM module
LMSM_WGCNA_miRTarget <- lapply(seq(LMSM_WGCNA_CommonmiRs), function(i) Extract.miRTarget(LMSM_WGCNA_CommonmiRs[[i]], LMSM_WGCNA_Modulegenes[[i]]))

## Understand predicted lncRNA-related miRNA sponge interactions, predicted and putative miRNA-target interactions of each LMSM module
LMSM_WGCNA_understand_miRSpongeTarget <- lapply(seq(LMSM_WGCNA_CommonmiRs), function(i) 
  Understand.miRSpongeTarget(scz_lncRNA_expr, scz_mRNA_expr, miRTarget, LMSM_WGCNA_CommonmiRs[[i]], 
                             LMSM_WGCNA_Modulegenes[[i]]))
LMSM_WGCNA_miRSpongeTarget <- do.call("rbind", LMSM_WGCNA_understand_miRSpongeTarget)
rownames(LMSM_WGCNA_miRSpongeTarget) <- names(LMSM_WGCNA_CommonmiRs)


## miRNAs distribution in LMSM modules
LMSM_WGCNA_miR_distribution <- miR.distribution(LMSM_WGCNA_CommonmiRs)

save.image("LMSM/LMSM_WGCNA_0.8_SCZ_miRNA_lncRNA_mRNA.RData")
write.table(LMSM_WGCNA_Modulegenes[[1]],"LMSM/LMSM_1.txt",row.names = F,col.names = F,quote = F)
write.table(LMSM_WGCNA_Modulegenes[[2]],"LMSM/LMSM_2.txt",row.names = F,col.names = F,quote = F)
write.table(LMSM_WGCNA_Modulegenes[[6]],"LMSM/LMSM_6.txt",row.names = F,col.names = F,quote = F)

