
ancestry_assoc = NULL
for(i in unique(profile$Cancer)){
  data = profile[which(profile$Cancer==i),]
  if(i %in% c("Breast_Carcinoma","Endometrial_Cancer","Prostate_Cancer","Ovarian_Cancer")){
    #None
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "None", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "None", "None", res$coefficients[3,c(1,2,4)]))
    #BMI
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$BMI + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "BMI", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "BMI", "BMI", res$coefficients[3,c(1,2,4)]))
    #WBC
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$WBC + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "WBC", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "WBC", "WBC", res$coefficients[3,c(1,2,4)]))
    #Autoimmune EverSmoked DrinkNow College
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$Autoimmune + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "Autoimmune", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "Autoimmune", "Autoimmune", res$coefficients[3,c(1,2,4)]))
    #EverSmoked
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$EverSmoked + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "EverSmoked", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "EverSmoked", "EverSmoked", res$coefficients[3,c(1,2,4)]))
    #DrinkNow
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$DrinkNow + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "DrinkNow", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "DrinkNow", "DrinkNow", res$coefficients[3,c(1,2,4)]))
    #College
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$College + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "College", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "College", "College", res$coefficients[3,c(1,2,4)]))
  } else{
    #None
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic) + as.factor(data$SEX)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "None", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "None", "None", res$coefficients[3,c(1,2,4)]))
    #BMI
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$BMI + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic) + as.factor(data$SEX)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "BMI", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "BMI", "BMI", res$coefficients[3,c(1,2,4)]))
    #WBC
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$WBC + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic) + as.factor(data$SEX)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "WBC", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "WBC", "WBC", res$coefficients[3,c(1,2,4)]))
    #Autoimmune EverSmoked DrinkNow College
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$Autoimmune + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic) + as.factor(data$SEX)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "Autoimmune", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "Autoimmune", "Autoimmune", res$coefficients[3,c(1,2,4)]))
    #EverSmoked
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$EverSmoked + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic) + as.factor(data$SEX)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "EverSmoked", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "EverSmoked", "EverSmoked", res$coefficients[3,c(1,2,4)]))
    #DrinkNow
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$DrinkNow + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic) + as.factor(data$SEX)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "DrinkNow", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "DrinkNow", "DrinkNow", res$coefficients[3,c(1,2,4)]))
    #College
    res = summary(lm(data$MSSProteinTMB ~ data$AJnonAJ + data$College + data$AGE + as.factor(data$PANEL_VERSION) + data$TUMOR_PURITY + as.factor(data$Metastatic) + as.factor(data$SEX)))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "College", "AJnonAJ", res$coefficients[2,c(1,2,4)]))
    ancestry_assoc = rbind(ancestry_assoc, c(i, "College", "College", res$coefficients[3,c(1,2,4)]))
  }
}

ancestry_assoc <- read.delim("~/Desktop/ancestry_assoc.txt", header=T)
colnames(ancestry_assoc) = c("Cancer", "Covariate", "IndVar", "Beta", "SE", "Pvalue")

pancancer = NULL
library(metafor)
for(i in unique(ancestry_assoc$Covariate)){
  data=ancestry_assoc[which(ancestry_assoc$Covariate==i & ancestry_assoc$IndVar==i),]
  z=rma(yi = as.numeric(data$Beta), sei = as.numeric(data$SE), method = "FE")
  pancancer = rbind(pancancer, c("Pan-Cancer",i, i, z$beta, z$se, z$pval))
  data=ancestry_assoc[which(ancestry_assoc$Covariate==i & ancestry_assoc$IndVar=="AJnonAJ"),]
  z=rma(yi = as.numeric(data$Beta), sei = as.numeric(data$SE), method = "FE")
  pancancer = rbind(pancancer, c("Pan-Cancer",i, "AJnonAJ", z$beta, z$se, z$pval))
}
pancancer = data.frame(pancancer)
colnames(pancancer) = colnames(ancestry_assoc)
ancestry_assoc = rbind(ancestry_assoc, pancancer)
ancestry_assoc$Beta = as.numeric(ancestry_assoc$Beta)
ancestry_assoc$SE = as.numeric(ancestry_assoc$SE)
ancestry_assoc$Pvalue = as.numeric(ancestry_assoc$Pvalue)
ancestry_assoc = ancestry_assoc[-which(ancestry_assoc$Covariate=="None" & ancestry_assoc$IndVar=="None"),]




full_profile_data <- read.delim("~/Desktop/DFCI/Cancer/Burden/full_profile_data.txt", header=T)
full_profile_data = subset(full_profile_data, select = c("SAMPLE_ID", "College", "EverSmoked"))
NSCLC_surv_IO <- read.csv("~/Desktop/DFCI/Cancer/Burden/Survival/NSCLC_surv_IO.csv", stringsAsFactors = F)
nsclc = read.table(paste0("~/Desktop/DFCI/Cancer/Burden/Survival/Non-Small_Cell_Lung_Cancer_INPUT.txt"), head = T, sep = '\t')
nsclc = merge(nsclc, full_profile_data, by = "SAMPLE_ID")
nsclc = merge(nsclc, NSCLC_surv_IO, by = "SAMPLE_ID")
nsclc$EduYears_Per = -1*nsclc$EduYears_Per
nsclc = subset(nsclc, select = c(tstart, tstop, event, QNORMTMB, `TMB.H`, AGE, SEX, Metastatic, PANEL_VERSION, TUMOR_PURITY, 
                                 PC1, PC2, PC3, PC4, PC5, EduYears_Per, College, EverSmoked))
colnames(nsclc) = c("tstart", "tstop", "event", "TMB", "TMBH", "Age", "Sex", "Metastatic", "Panel_Version", "Tumor_Purity", 
                    "PC1", "PC2", "PC3", "PC4", "PC5", "EduYears", "College", "EverSmoked")
nsclc$Sex = factor(nsclc$Sex, levels = c("Male","Female"))
nsclc$Metastatic = factor(nsclc$Metastatic, levels = c(0,1), labels = c("Primary","Metastatic"))
nsclc$Panel_Version = factor(nsclc$Panel_Version, levels = c(1,2,3))
nsclc$College = factor(nsclc$College, levels = c(0,1), labels = c("No","Yes"))
nsclc$EverSmoked = factor(nsclc$EverSmoked, levels = c(0,1), labels = c("No","Yes"))
summary(coxph(Surv(tstart, tstop, event) ~ TMB*EduYears + College + Age + Sex + Metastatic + Panel_Version + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, nsclc))


full_profile_data <- read.delim("~/Desktop/DFCI/Cancer/Burden/full_profile_data.txt", header=T)
name = c("Bladder_Cancer", "Breast_Carcinoma","Cancer_of_Unknown_Primary", 
         "Colorectal_Cancer","Endometrial_Cancer","Esophagogastric_Carcinoma",
         "Glioma","Head_and_Neck_Carcinoma", "Leukemia", "Melanoma",
         "Non-Hodgkin_Lymphoma","Non-Small_Cell_Lung_Cancer", "Ovarian_Cancer",
         "Pancreatic_Cancer", "Prostate_Cancer", "Renal_Cell_Carcinoma", "Soft_Tissue_Sarcoma")

source("~/Desktop/DFCI/Cancer/Burden/random/lr_clin_anc.R")
outcomes["TMBH"] = 0
outcomes[which(outcomes$OrigMSSTMB>=10),"TMBH"] = 1
outcomes = subset(outcomes, select = c("SAMPLE_ID","TMBH"))

college = merge(outcomes, full_profile_data, by = "SAMPLE_ID")
results = NULL
for(i in unique(college$Cancer)){
  temp = college[which(college$Cancer==i),]
  if(i %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
    x = summary(lm(MSSProteinTMB ~ College + EverSmoked + AGE + factor(PANEL_VERSION) + TUMOR_PURITY + factor(Metastatic) + 
                     scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5), data=temp))
    results = rbind(results, c(i,"TMB", "College", x$coef[2,c(1,2,4)]))
    results = rbind(results, c(i,"TMB", "Smoked", x$coef[3,c(1,2,4)]))
    x = summary(glm(TMBH ~ College + EverSmoked + AGE + factor(PANEL_VERSION) + TUMOR_PURITY + factor(Metastatic) + 
                      scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5), data=temp, family = binomial(link = "logit")))
    results = rbind(results, c(i,"TMBH", "College", x$coef[2,c(1,2,4)]))
    results = rbind(results, c(i,"TMBH", "Smoked", x$coef[3,c(1,2,4)]))
  } else{
    x = summary(lm(MSSProteinTMB ~ College + EverSmoked + AGE + factor(SEX) + factor(PANEL_VERSION) + TUMOR_PURITY + factor(Metastatic) + 
             scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5), data=temp))
    results = rbind(results, c(i,"TMB", "College", x$coef[2,c(1,2,4)]))
    results = rbind(results, c(i,"TMB", "Smoked", x$coef[3,c(1,2,4)]))
    x = summary(glm(TMBH ~ College + EverSmoked + AGE + factor(SEX) + factor(PANEL_VERSION) + TUMOR_PURITY + factor(Metastatic) + 
                   scale(PC1) + scale(PC2) + scale(PC3) + scale(PC4) + scale(PC5), data=temp, family = binomial(link = "logit")))
    results = rbind(results, c(i,"TMBH", "College", x$coef[2,c(1,2,4)]))
    results = rbind(results, c(i,"TMBH", "Smoked", x$coef[3,c(1,2,4)]))
  } 
}
results = data.frame(results)
colnames(results) = c("Cancer", "Burden", "IndVar", "Beta", "SE", "Pvalue")
results[which(results$Cancer %in% c("Leukemia", "Pancreatic_Cancer", "Prostate_Cancer", "Renal_Cell_Carcinoma")), 4:6] = NA
results$Beta = as.numeric(results$Beta)
results$SE = as.numeric(results$SE)
results$Pvalue = as.numeric(results$Pvalue)
a = rma(yi = results[which(results$Burden=="TMB" & results$IndVar=="College"),"Beta"],
    sei = results[which(results$Burden=="TMB" & results$IndVar=="College"),"SE"], method = "FE")

b = rma(yi = results[which(results$Burden=="TMBH" & results$IndVar=="College"),"Beta"],
    sei = results[which(results$Burden=="TMBH" & results$IndVar=="College"),"SE"], method = "FE")

c = rma(yi = results[which(results$Burden=="TMB" & results$IndVar=="Smoked"),"Beta"],
    sei = results[which(results$Burden=="TMB" & results$IndVar=="Smoked"),"SE"], method = "FE")

d = rma(yi = results[which(results$Burden=="TMBH" & results$IndVar=="Smoked"),"Beta"],
    sei = results[which(results$Burden=="TMBH" & results$IndVar=="Smoked"),"SE"], method = "FE")

results = rbind(results, c("Pan-Cancer", "TMB", "College", a$beta, a$se, a$pval))
results = rbind(results, c("Pan-Cancer", "TMBH", "College", b$beta, b$se, b$pval))
results = rbind(results, c("Pan-Cancer", "TMB", "Smoked", c$beta, c$se, c$pval))
results = rbind(results, c("Pan-Cancer", "TMBH", "Smoked", d$beta, d$se, d$pval))
results$Beta = as.numeric(results$Beta)
results$SE = as.numeric(results$SE)
results$Pvalue = as.numeric(results$Pvalue)

library(survival)

Melanoma_surv_IO <- read.csv("~/Desktop/DFCI/Cancer/Burden/Survival/Melanoma_surv_IO.csv", stringsAsFactors = F)
Melanoma_INPUT <- read.delim("~/Desktop/DFCI/Cancer/Burden/Survival/Melanoma_INPUT.txt", stringsAsFactors = F)
melanoma = merge(Melanoma_INPUT, Melanoma_surv_IO, by = "SAMPLE_ID")

melanoma = subset(melanoma, select = c(SAMPLE_ID, tstart, tstop, event, QNORMAllCNB, AGE, SEX, Metastatic, PANEL_VERSION, TUMOR_PURITY, 
                                       PC1, PC2, PC3, PC4, PC5, Smoker_Per, CigsPerDay_Pan))
colnames(melanoma) = c("SAMPLE_ID", "tstart", "tstop", "event", "AllCNB", "Age", "Sex", "Metastatic", "Panel_Version", "Tumor_Purity", 
                       "PC1", "PC2", "PC3", "PC4", "PC5", "Smoker", "Cigs")
melanoma = merge(full_profile_data, melanoma, by = "SAMPLE_ID")
melanoma$Sex = factor(melanoma$Sex, levels = c("Female","Male"))
melanoma$Metastatic = factor(melanoma$Metastatic, levels = c(0,1), labels = c("Primary","Metastatic"))
melanoma$Panel_Version = factor(melanoma$Panel_Version, levels = c(3,1,2))

summary(coxph(Surv(tstart, tstop, event) ~ Smoker*AllCNB + Cigs + Age + Sex + Metastatic + Panel_Version + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, melanoma))

summary(coxph(Surv(tstart, tstop, event) ~ AllCNB + Age + Sex + Metastatic + Panel_Version + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, melanoma))


nonIO_patients <- read.csv("~/Desktop/DFCI/Cancer/Burden/Survival/nonIO_patients.csv")

Prostate_Cancer_INPUT <- read.delim("~/Desktop/DFCI/Cancer/Burden/Survival/Prostate_Cancer_INPUT.txt", stringsAsFactors = F)
prostate = merge(Prostate_Cancer_INPUT, nonIO_patients, by = "SAMPLE_ID")
prostate = subset(prostate, select = c(SAMPLE_ID, tstart, tstop, event, QNORMAllCNB, AGE, SEX, Metastatic, PANEL_VERSION, TUMOR_PURITY, 
                                       PC1, PC2, PC3, PC4, PC5))
colnames(prostate) = c("SAMPLE_ID", "tstart", "tstop", "event", "AllCNB", "Age", "Sex", "Metastatic", "Panel_Version", "Tumor_Purity", 
                       "PC1", "PC2", "PC3", "PC4", "PC5")
prostate$Sex = factor(prostate$Sex, levels = c("Male","Female"))
prostate$Metastatic = factor(prostate$Metastatic, levels = c(0,1), labels = c("Primary","Metastatic"))
prostate$Panel_Version = factor(prostate$Panel_Version, levels = c(1,2,3))
Prostate_Cancer_CodingSNVs_BRCA12 <- read.delim("~/Desktop/DFCI/Cancer/Burden/Prostate_Cancer_CodingSNVs_BRCA12.txt")
prostate["BRCA12_SNV"] = 0
prostate[which(prostate$SAMPLE_ID %in% Prostate_Cancer_CodingSNVs_BRCA12$SAMPLE_ID),"BRCA12_SNV"] = 1

Prostate_Cancer_CNV_BRCA12 <- read.delim("~/Desktop/DFCI/Cancer/Burden/cnv_brca12.txt")
prostate["BRCA12_CNV"] = 0
keep = Prostate_Cancer_CNV_BRCA12[,which(colnames(Prostate_Cancer_CNV_BRCA12) %in% prostate$SAMPLE_ID)]
final = data.frame("SAMPLE_ID" = colnames(keep), "BRCA12" = colSums(keep))
final = final[which(final$BRCA12!=0),]
prostate[which(prostate$SAMPLE_ID %in% final$SAMPLE_ID),"BRCA12_CNV"] = 1
prostate["ANY_MUT"] = 0
prostate[which(prostate$BRCA12_CNV==1 | prostate$BRCA12_SNV==1),"ANY_MUT"] = 1
summary(coxph(Surv(tstart, tstop, event) ~ AllCNB + ANY_MUT + Age + Metastatic + Panel_Version + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, prostate))$coef
