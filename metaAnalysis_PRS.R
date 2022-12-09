library(metafor)
tcga = NULL
for(cancer in c("Bladder_Cancer","Breast_Carcinoma","Colorectal_Cancer", "Endometrial_Cancer",
                "Esophagogastric_Carcinoma", "Glioma", "Head_and_Neck_Carcinoma", "Melanoma",
                "Non-Small_Cell_Lung_Cancer", "Ovarian_Cancer", "Prostate_Cancer", "Renal_Cell_Carcinoma")){
  for(burden in c("AllCNB", "CNB", "TMB")){
    data = read.delim(paste0("~/Desktop/Cancer/Burden/TCGA/PRS/", cancer,"_QNORM", burden, "_prs_match_profile.assoc.txt"), stringsAsFactors = F)
    data["Cancer"] = cancer
    data["Burden"] = burden
    data = subset(data, select = c(Cancer, Burden, rs, ps, beta, se, p_wald))
    colnames(data) = c("Cancer", "Burden", "PRS", "Level", "TCGA_Beta", "TCGA_Stderr", "TCGA_Pvalue")
    tcga = rbind(tcga, data)   
  }
}
tcga[which(tcga$Burden=="AllCNB"),"Burden"] = "All CNB"
pancancer=NULL
for(burden in unique(tcga$Burden)){
  for(prs in unique(tcga$PRS)){
    for(lvl in unique(tcga$Level)){
      this = tcga[which(tcga$Burden==burden & tcga$PRS==prs & tcga$Level==lvl),]
      z=rma(yi=as.numeric(this$TCGA_Beta), sei=as.numeric(this$TCGA_Stderr), method="FE")
      pancancer = rbind(pancancer, c("Pan-Cancer", burden, prs, lvl, z$beta, z$se, z$pval))
    }
  }
}
colnames(pancancer) = colnames(tcga)
tcga = rbind(tcga, pancancer)


source('~/Desktop/Cancer/prs_lmm_analyses.R')

PRS$Pvalue = as.numeric(PRS$Pvalue)
bestPRS = NULL
for(i in unique(PRS$Cancer)){
  for(j in unique(PRS$Burden)){
    for(k in unique(PRS$PRS)){
      check = PRS[which(PRS$Cancer==i & PRS$Burden==j & PRS$PRS==k & PRS$Source=="Profile"),]
      level = check[which(check$Pvalue==min(check$Pvalue)),"Level"]
      add = PRS[which(PRS$Cancer==i & PRS$Burden==j & PRS$PRS==k & PRS$Level==level),]
      bestPRS = rbind(bestPRS, add)
    }
  }
}

profile=bestPRS[which(bestPRS$Source=="Profile"),]
profile = subset(profile, select = c(Cancer, Burden, PRS, Level, Beta, Stderr, Pvalue))
tumor=bestPRS[which(bestPRS$Source=="Tempus" & bestPRS$Sample=="Tumor"),]
tumor = subset(tumor, select = c(Cancer, Burden, PRS, Level, Beta, Stderr, Pvalue))
normal=bestPRS[which(bestPRS$Source=="Tempus" & bestPRS$Sample=="Normal"),]
normal = subset(normal, select = c(Cancer, Burden, PRS, Level, Beta, Stderr, Pvalue))
bestPRS = merge(profile, tumor, by = c("Cancer", "Burden", "PRS", "Level"), all = T)
bestPRS = merge(bestPRS, normal, by = c("Cancer", "Burden", "PRS", "Level"), all = T)
colnames(bestPRS) = c("Cancer","Burden", "PRS", "Level", "Beta_Profile", "Stderr_Profile", "Pvalue_Profile", "Beta_Tumor", "Stderr_Tumor", "Pvalue_Tumor", "Beta_Normal", "Stderr_Normal", "Pvalue_Normal")

bestPRS = merge(bestPRS, tcga, by = c("Cancer", "Burden", "PRS", "Level"), all.x = T)
bestPRS["metaAll_Beta_T"] = NA
bestPRS["metaAll_Stderr_T"] = NA
bestPRS["metaAll_Pvalue_T"] = NA

bestPRS["metaAll_Beta_N"] = NA
bestPRS["metaAll_Stderr_N"] = NA
bestPRS["metaAll_Pvalue_N"] = NA
for(line in 1:dim(bestPRS)[1]){
  z=rma(yi=as.numeric(c(bestPRS[line,"Beta_Profile"],bestPRS[line,"Beta_Tumor"], bestPRS[line,"TCGA_Beta"])),sei=as.numeric(c(bestPRS[line,"Stderr_Profile"],bestPRS[line,"Stderr_Tumor"], bestPRS[line,"TCGA_Stderr"])),method="FE")
  bestPRS[line, "metaAll_Beta_T"] = z$beta
  bestPRS[line, "metaAll_Stderr_T"] = z$se
  bestPRS[line, "metaAll_Pvalue_T"] = z$pval
  
  z=rma(yi=as.numeric(c(bestPRS[line,"Beta_Profile"],bestPRS[line,"Beta_Normal"], bestPRS[line,"TCGA_Beta"])),sei=as.numeric(c(bestPRS[line,"Stderr_Profile"],bestPRS[line,"Stderr_Normal"], bestPRS[line,"TCGA_Stderr"])),method="FE")
  bestPRS[line, "metaAll_Beta_N"] = z$beta
  bestPRS[line, "metaAll_Stderr_N"] = z$se
  bestPRS[line, "metaAll_Pvalue_N"] = z$pval
}

bestPRS$metaAll_Pvalue_T = as.numeric(bestPRS$metaAll_Pvalue_T)
bestPRS$metaAll_Pvalue_N = as.numeric(bestPRS$metaAll_Pvalue_N)

sig = bestPRS[which(bestPRS$metaAll_Pvalue_T<0.05/14),]

write.table(bestPRS, "~/Desktop/Cancer/Burden/Paper/Tables/final_allmetaprs.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(sig, "~/Desktop/Cancer/Burden/Paper/Tables/final_sigmetaprs.txt", col.names = T, row.names = F, sep = "\t", quote = F)

