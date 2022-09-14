setwd("D:/1_CeRNA_Modules/PanData/Post_COMBAT/")
#######----------------- Putative_Interactions CeRNA network
###
load("Data/Pan_NPD/LMSM_ASD_SCZ_BPD.RData")
sign_module_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj < 0.05,]))
putative_interactions=read.csv("Putative_Interactions.csv")

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
  all_put_int_mod=rbind(all_put_int_mod,put_int_mod[[i]])
  i=i+1
}
all_put_int_mod=dplyr::bind_rows(put_int_mod, .id="Module")
lncRNAs=unique(unlist(lncRNAs))
write.csv((all_put_int_mod),"Putative_Interactions_miRNASpongeModules.csv",row.names = F,quote = F)

write.table(lncRNAs,"LncRNAs.txt",col.names = F,row.names = F,quote = F)
