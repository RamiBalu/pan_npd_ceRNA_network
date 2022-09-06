
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
