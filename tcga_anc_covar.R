library(metafor)
cancers = c("Bladder_Cancer","Breast_Carcinoma","Colorectal_Cancer", "Endometrial_Cancer",
            "Esophagogastric_Carcinoma", "Glioma", "Head_and_Neck_Carcinoma", "Melanoma",
            "Non-Small_Cell_Lung_Cancer", "Ovarian_Cancer", "Prostate_Cancer", "Renal_Cell_Carcinoma")

results = NULL
for(c in cancers){
  data = read.table(paste0("~/Desktop/Cancer/Burden/TCGA/analyses/", c, "/", c, "_data.txt"), head = T, sep = '\t')
  data["AJnonAJ"] = 0
  data[which(data$ProjectedPC2>=1.4E-8),"AJnonAJ"] = 1
  data[which(data$AJnonAJ==1),"ProjectedPC1"] = NA
  data["TMBH"] = 0
  data[which(data$TMB>=10),"TMBH"] = 1
  if(c %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
    sex_tmb = c(NA, NA, NA)
    sex_tmbh = c(NA, NA, NA)
    sex_cnb = c(NA, NA, NA)
    sex_cnball = c(NA, NA, NA)
    
    if(dim(data[which(data$TMBH==1),])[1]>=10){
      age_tmbh = summary(glm(formula = TMBH ~ Age + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data, family='binomial'))$coefficients[2,c(1,2,4)]
      nwse_tmbh = summary(glm(formula = TMBH ~ ProjectedPC1 + Age + Tumor_Purity, data=data, family='binomial'))$coefficients[2,c(1,2,4)]
    } else {
      age_tmbh = c(NA, NA, NA)
      nwse_tmbh = c(NA, NA, NA)
    }
    
    age_tmb = summary(lm(formula = QNORMTMB ~ Age + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data))$coefficients[2,c(1,2,4)]
    age_cnb = summary(lm(formula = QNORMCNB ~ Age + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data))$coefficients[2,c(1,2,4)]
    age_cnball = summary(lm(formula = QNORMAllCNB ~ Age + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data))$coefficients[2,c(1,2,4)]
    
    nwse_tmb = summary(lm(formula = QNORMTMB ~ ProjectedPC1 + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
    nwse_cnb = summary(lm(formula = QNORMCNB ~ ProjectedPC1 + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
    nwse_cnball = summary(lm(formula = QNORMAllCNB ~ ProjectedPC1 + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]

    if(dim(data[which(data$AJnonAJ==1),])[1]>=10){
      aj_tmb = summary(lm(formula = QNORMTMB ~ AJnonAJ + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
      if(dim(data[which(data$TMBH==1),])[1]>=10){
        aj_tmbh = summary(glm(formula = TMBH ~ AJnonAJ + Age + Tumor_Purity, data=data, family='binomial'))$coefficients[2,c(1,2,4)]
      } else{
        aj_tmbh = c(NA, NA, NA)
      }
      aj_cnb = summary(lm(formula = QNORMCNB ~ AJnonAJ + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
      aj_cnball = summary(lm(formula = QNORMAllCNB ~ AJnonAJ + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
    
    } else{
      aj_tmb = c(NA, NA, NA)
      aj_tmbh = c(NA, NA, NA)
      aj_cnb = c(NA, NA, NA)
      aj_cnball = c(NA, NA, NA)
    
    }
  } else{
    sex_tmb = summary(lm(formula = QNORMTMB ~ Sex + Age + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data))$coefficients[2,c(1,2,4)]
    age_tmb = summary(lm(formula = QNORMTMB ~ Age + Sex + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data))$coefficients[2,c(1,2,4)]
    if(dim(data[which(data$TMBH==1),])[1]>=10){
      sex_tmbh = summary(glm(formula = TMBH ~ Sex + Age + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data, family='binomial'))$coefficients[2,c(1,2,4)]
      age_tmbh = summary(glm(formula = TMBH ~ Age + Sex + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data, family='binomial'))$coefficients[2,c(1,2,4)]
    } else{
      sex_tmbh = c(NA, NA, NA)
      age_tmbh = c(NA, NA, NA)
    }
    sex_cnb = summary(lm(formula = QNORMCNB ~ Sex + Age + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data))$coefficients[2,c(1,2,4)]
    age_cnb = summary(lm(formula = QNORMCNB ~ Age + Sex + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data))$coefficients[2,c(1,2,4)]

    sex_cnball = summary(lm(formula = QNORMAllCNB ~ Sex + Age + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data))$coefficients[2,c(1,2,4)]
    age_cnball = summary(lm(formula = QNORMAllCNB ~ Age + Sex + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, data=data))$coefficients[2,c(1,2,4)]
    
    nwse_tmb = summary(lm(formula = QNORMTMB ~ ProjectedPC1 + Sex + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
    if(dim(data[which(data$TMBH==1),])[1]>=10){
      nwse_tmbh = summary(glm(formula = TMBH ~ ProjectedPC1 + Sex + Age + Tumor_Purity, data=data, family='binomial'))$coefficients[2,c(1,2,4)]
    } else{
      nwse_tmbh = c(NA, NA, NA)
    }
    nwse_cnb = summary(lm(formula = QNORMCNB ~ ProjectedPC1 + Sex + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
    nwse_cnball = summary(lm(formula = QNORMAllCNB ~ ProjectedPC1 + Sex + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
    
    if(dim(data[which(data$AJnonAJ==1),])[1]>=10){
      aj_tmb = summary(lm(formula = QNORMTMB ~ AJnonAJ + Sex + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
      if(dim(data[which(data$TMBH==1),])[1]>=10){
        aj_tmbh = summary(glm(formula = TMBH ~ AJnonAJ + Sex + Age + Tumor_Purity, data=data, family='binomial'))$coefficients[2,c(1,2,4)]
      } else{
        aj_tmbh = c(NA, NA, NA)
      }
      aj_cnb = summary(lm(formula = QNORMCNB ~ AJnonAJ + Sex + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
      aj_cnball = summary(lm(formula = QNORMAllCNB ~ AJnonAJ + Sex + Age + Tumor_Purity, data=data))$coefficients[2,c(1,2,4)]
      
    } else{
      aj_tmb = c(NA, NA, NA)
      aj_tmbh = c(NA, NA, NA)
      aj_cnb = c(NA, NA, NA)
      aj_cnball = c(NA, NA, NA)
      
    }
  }
  
  final = rbind(age_tmb, sex_tmb)
  final = rbind(final, nwse_tmb)
  final = rbind(final, aj_tmb)
  
  final = rbind(final, age_tmbh)
  final = rbind(final, sex_tmbh)
  final = rbind(final, nwse_tmbh)
  final = rbind(final, aj_tmbh)
  
  final = rbind(final, age_cnb)
  final = rbind(final, sex_cnb)
  final = rbind(final, nwse_cnb)
  final = rbind(final, aj_cnb)
  
  final = rbind(final, age_cnball)
  final = rbind(final, sex_cnball)
  final = rbind(final, nwse_cnball)
  final = rbind(final, aj_cnball)
  
  final = data.frame(final)
  colnames(final) = c("Beta", "SE", "Pvalue")
  final["Burden"] = c(rep("TMB", 4), rep("TMBH", 4), rep("CNB", 4), rep("All CNB", 4))
  final["Ind_Var"] = rep(c("Age", "Sex", "NWSE", "AJnonAJ"),4)
  final["Cancer"] = c
  results = rbind(results, final)
}


results = subset(results, select = c(Cancer, Burden, Ind_Var, Beta, SE, Pvalue))
rownames(results) = 1:nrow(results)
for(i in unique(results$Burden)){
  for(j in unique(results$Ind_Var)){
    temp = results[which(results$Burden==i & results$Ind_Var==j),]
    temp$Beta = as.numeric(temp$Beta)
    temp$SE = as.numeric(temp$SE)
    z=rma(yi=temp$Beta, sei=temp$SE, method = 'FE')
    results = rbind(results, c("Pan-Cancer", i, j, z$beta, z$se, z$pval))
  }
}

results$Beta = as.numeric(results$Beta)
results$SE = as.numeric(results$SE)
results$Pvalue = as.numeric(results$Pvalue)
