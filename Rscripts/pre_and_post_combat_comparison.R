setwd("D:/CeRNA_Modules/PanData/")
load('Data_Preparation/LMSM_ASD_SCZ_BPD_noctl.RData')
sign_module_index=as.numeric(rownames(LMSM_WGCNA_modules_p.values[LMSM_WGCNA_modules_p.values$p.adj<=0.05,]))

post_cmbt_LMSM_sign_modulegenes=c()
i=1
for(s in sign_module_index){
  m=
  post_cmbt_LMSM_sign_modulegenes[[i]]=LMSM_WGCNA_Modulegenes[[s]]
  i=i+1
}

load("WGCNA/LMSM_ASD+SCZ+BPD.RData")

pre_cmbt_LMSM_sign_modulegenes=c()
i=1
for(s in sign_module_index){
  pre_cmbt_LMSM_sign_modulegenes[[i]]=LMSM_WGCNA_Modulegenes[[s]]
  i=i+1
}

intersecting_genes=c()
for (i in seq(1:10)){
  for (j in seq(1:10)){
    m=paste(c("Post-Combat-M",i,"Pre-Combat-M",j),collapse = "-")
    intersecting_genes[[m]]=length(intersect(post_cmbt_LMSM_sign_modulegenes[[i]],pre_cmbt_LMSM_sign_modulegenes[[j]]))
  }
}
intersecting_genes=data.frame(t(data.frame(intersecting_genes)))

postcmbt=unlist(post_cmbt_LMSM_sign_modulegenes)
precmbt=unlist(pre_cmbt_LMSM_sign_modulegenes)
length(intersect(precmbt,postcmbt))