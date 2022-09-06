
### RNAs KME_Table

kme=read.csv("kME_table.csv",row.names = 1)[,-c(4:7,9:11,13:23)]
kme_lncRNA=kme[kme$gene_biotype=="lncRNA",]
kme_mRNA=kme[kme$gene_biotype=="protein_coding",]
colnames(bpd_de_lncrna) <- colnames(scz_de_lncrna) <- colnames(asd_de_lncrna) <- 
  colnames(bpd_de_mrna) <- colnames(scz_de_mrna) <- colnames(asd_de_mrna) <- c("external_gene_id","Direction")

kme_lncRNA=merge(kme_lncRNA,asd_de_lncrna,by="external_gene_id",all.x = TRUE)
kme_lncRNA=merge(kme_lncRNA,scz_de_lncrna,by="external_gene_id",all.x = TRUE)
kme_lncRNA=merge(kme_lncRNA,bpd_de_lncrna,by="external_gene_id",all.x = TRUE)

kme_mRNA=merge(kme_mRNA,asd_de_mrna,by="external_gene_id",all.x = TRUE)
kme_mRNA=merge(kme_mRNA,scz_de_mrna,by="external_gene_id",all.x = TRUE)
kme_mRNA=merge(kme_mRNA,bpd_de_mrna,by="external_gene_id",all.x = TRUE)

kme_summary=rbind.data.frame(kme_mRNA, kme_lncRNA)
colnames(kme_summary) <- colnames(kme_mRNA) <- colnames(kme_lncRNA) <- c(colnames(kme_summary)[1:5],"ASD_DEG","SCZ_DEG","BD_DEG")

