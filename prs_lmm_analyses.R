library(stringr)
library(metafor)

profile_prs = NULL
for(outcome in c("DeepCNB", "CNB", "ProteinTMB")){
  for(cancer in c("Bladder_Cancer", "Breast_Carcinoma","Cancer_of_Unknown_Primary", 
           "Colorectal_Cancer","Endometrial_Cancer","Esophagogastric_Carcinoma",
           "Glioma","Head_and_Neck_Carcinoma", "Leukemia", "Melanoma",
           "Non-Hodgkin_Lymphoma","Non-Small_Cell_Lung_Cancer", "Ovarian_Cancer",
           "Pancreatic_Cancer", "Prostate_Cancer", "Renal_Cell_Carcinoma", "Soft_Tissue_Sarcoma")){
    data = read.delim(paste0("~/Desktop/Cancer/Burden/Profile/prs_results/", cancer, "_MSS", outcome, "_prs.assoc.txt"), stringsAsFactors = F)
    data["Cancer"] = cancer
    data["Burden"] = outcome
    data["rs"] = str_remove(data$rs, ".sscore")
    x=data.frame(str_split(data$rs, ".PV", simplify = T))
    data["PRS"] = x$X1
    data["Level"] = x$X2
    data = subset(data, select = c("Cancer", "Burden", "PRS", "Level", "beta", "se", "p_wald"))
    colnames(data) = c("Cancer", "Burden", "PRS", "Level", "Beta", "Stderr", "Pvalue")
    profile_prs = rbind(profile_prs, data)
  }
}
profile_prs[which(profile_prs$Burden=="CNB"),"Burden"] = "All CNB"
profile_prs[which(profile_prs$Burden=="DeepCNB"),"Burden"] = "CNB"
profile_prs[which(profile_prs$Burden=="ProteinTMB"),"Burden"] = "TMB"
profile_prs["Source"] = "Profile"
profile_prs["Sample"] = "Tumor"

tempus_prs = NULL
for(outcome in c("CNB", "ProteinTMB")){
  for(cancer in c("bladder", "breast", "colorectal", "endometrial", "esophagogastric", "glioma",
                  "head", "lung", "melanoma","ovarian", "pancreatic", "prostate", "soft", "unknown")){
    for(sample in c("tumor", "normal")){
      data = read.delim(paste0("~/Desktop/Cancer/Burden/Tempus/prs_results/", cancer, "_Qnorm", outcome, "_", sample, "_prs.assoc.txt"), stringsAsFactors = F)
      data["Cancer"] = cancer
      data["Burden"] = outcome
      data["rs"] = str_remove(data$rs, ".sscore")
      x=data.frame(str_split(data$rs, ".PV", simplify = T))
      data["PRS"] = x$X1
      data["Level"] = x$X2
      data = subset(data, select = c("Cancer", "Burden", "PRS", "Level", "beta", "se", "p_wald"))
      colnames(data) = c("Cancer", "Burden", "PRS", "Level", "Beta", "Stderr", "Pvalue")
      data["Source"] = "Tempus"
      data["Sample"] = sample
      tempus_prs = rbind(tempus_prs, data)
    }
  }
}

tempus_prs[which(tempus_prs$Sample=="tumor"), "Sample"] = "Tumor"
tempus_prs[which(tempus_prs$Sample=="normal"), "Sample"] = "Normal"
tempus_prs[which(tempus_prs$Burden=="ProteinTMB"),"Burden"] = "TMB"
tempus_prs[which(tempus_prs$Cancer=="bladder"),"Cancer"] = "Bladder_Cancer"
tempus_prs[which(tempus_prs$Cancer=="breast"),"Cancer"] = "Breast_Carcinoma"
tempus_prs[which(tempus_prs$Cancer=="colorectal"),"Cancer"] = "Colorectal_Cancer"
tempus_prs[which(tempus_prs$Cancer=="endometrial"),"Cancer"] = "Endometrial_Cancer"
tempus_prs[which(tempus_prs$Cancer=="esophagogastric"),"Cancer"] = "Esophagogastric_Carcinoma"
tempus_prs[which(tempus_prs$Cancer=="glioma"),"Cancer"] = "Glioma"
tempus_prs[which(tempus_prs$Cancer=="head"),"Cancer"] = "Head_and_Neck_Carcinoma"
tempus_prs[which(tempus_prs$Cancer=="lung"),"Cancer"] = "Non-Small_Cell_Lung_Cancer"
tempus_prs[which(tempus_prs$Cancer=="melanoma"),"Cancer"] = "Melanoma"
tempus_prs[which(tempus_prs$Cancer=="ovarian"),"Cancer"] = "Ovarian_Cancer"
tempus_prs[which(tempus_prs$Cancer=="pancreatic"),"Cancer"] = "Pancreatic_Cancer"
tempus_prs[which(tempus_prs$Cancer=="prostate"),"Cancer"] = "Prostate_Cancer"
tempus_prs[which(tempus_prs$Cancer=="soft"),"Cancer"] = "Soft_Tissue_Sarcoma"
tempus_prs[which(tempus_prs$Cancer=="unknown"),"Cancer"] = "Cancer_of_Unknown_Primary"

panCancerPRS = NULL
for(outcome in c("TMB", "CNB", "All CNB")){
  for(prs in unique(profile_prs$PRS)){
    for(level in 0:7){
      temp = profile_prs[which(profile_prs$Burden==outcome & profile_prs$PRS==prs & profile_prs$Level==level),]
      z=rma(yi=temp$Beta, sei = temp$Stderr, method = 'FE')
      panCancerPRS = rbind(panCancerPRS, c("Pan-Cancer", outcome, prs, level, z$beta, z$se, z$pval, "Profile", "Tumor"))
    }
  }
}
  
tempus_prs = tempus_prs[-grep('scaled',tempus_prs$PRS),]
tempus_prs = tempus_prs[-grep("combined",tempus_prs$PRS),]
for(outcome in c("TMB", "CNB")){
  for(sample in c("Tumor", "Normal")){
    for(prs in unique(tempus_prs$PRS)){
      for(level in 0:7){
        temp = tempus_prs[which(tempus_prs$Burden==outcome & tempus_prs$Sample==sample & tempus_prs$PRS==prs & tempus_prs$Level==level),]
        z=rma(yi=temp$Beta, sei = temp$Stderr, method = 'FE')
        panCancerPRS = rbind(panCancerPRS, c("Pan-Cancer", outcome, prs, level, z$beta, z$se, z$pval, "Tempus", sample))
      }
    }
  }
}

panCancerPRS = as.data.frame(panCancerPRS)
colnames(panCancerPRS) = colnames(profile_prs)
PRS = rbind(profile_prs, tempus_prs)
PRS = rbind(PRS, panCancerPRS)
PRS[which(PRS$PRS=="UKB_460K.cov_EDU_YEARS"),"Beta"] = -1*as.numeric(PRS[which(PRS$PRS=="UKB_460K.cov_EDU_YEARS"),"Beta"])
rm(cancer, data, level, outcome, panCancerPRS, profile_prs, prs, sample, temp, tempus_prs, x, z)
