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