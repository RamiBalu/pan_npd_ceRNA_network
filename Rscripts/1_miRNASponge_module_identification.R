
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

