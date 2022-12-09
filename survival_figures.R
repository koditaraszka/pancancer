library(metafor)
library(ggplot2)
library(patchwork)
all_theme = theme_bw() + theme(strip.background = element_blank(), 
                               axis.text =element_text( size = 12, color = "black"), 
                               axis.title = element_text(size = 12, color = "black"), 
                               legend.text = element_text(size = 12, color = "black"), 
                               legend.title = element_text(size = 12, color = "black"), 
                               plot.title = element_text(hjust = 0.5, size=12, color="black"), 
                               strip.text = element_text(size = 12, color = "black"), 
                               plot.tag = element_text(size=12, color = "black", face = "bold"))

#OS ~ Burden IO then non-IO
burden = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_IO_burden.csv")
burden = burden[which(burden$analysis=="Survival ~ outcome + burden + covars" & burden$outcome=="AGE" & burden$variable=="burden"),]
burden = subset(burden, select = c(cancer, burden, coef, pval, N))
burden["zscore"] = qnorm(burden$pval/2, lower.tail = F)
burden[which(burden$coef<0),"zscore"] = -1 * burden[which(burden$coef<0),"zscore"] 
burden["stderr"] = burden$coef/burden$zscore
burden = burden[which(burden$N>=40),]
burden["treatment"] = "IO"
results = NULL
for(b in unique(burden$burden)){
  temp = burden[which(burden$burden==b),]
  z = rma(yi = temp$coef, sei = temp$stderr, method = 'FE')
  results = rbind(results, c("Pan-Cancer", b, z$beta, z$pval, NA, z$zval, z$se, "IO"))
}
results = data.frame(results)
colnames(results) = colnames(burden)
burden = rbind(burden, results)

burden2 = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_nonIO_burden.csv")
burden2 = burden2[which(burden2$analysis=="Survival ~ outcome + burden + covars" & burden2$outcome=="AGE" & burden2$variable=="burden"),]
burden2 = subset(burden2, select = c(cancer, burden, coef, pval, N))
burden2["zscore"] = qnorm(burden2$pval/2, lower.tail = F)
burden2[which(burden2$coef<0),"zscore"] = -1 * burden2[which(burden2$coef<0),"zscore"] 
burden2["stderr"] = burden2$coef/burden2$zscore
burden2 = burden2[which(burden2$N>=40),]
burden2["treatment"] = "non-IO"
results = NULL
for(b in unique(burden2$burden)){
  temp = burden2[which(burden2$burden==b),]
  z = rma(yi = temp$coef, sei = temp$stderr, method = 'FE')
  results = rbind(results, c("Pan-Cancer", b, z$beta, z$pval, NA, z$zval, z$se, "non-IO"))
}
results = data.frame(results)
colnames(results) = colnames(burden)
burden = rbind(burden, burden2)
burden = rbind(burden, results)

burden$pval = as.numeric(burden$pval)
burden$coef = as.numeric(burden$coef)
burden$stderr = as.numeric(burden$stderr)
burden["HR"] = exp(burden$coef)
burden["signif"] = ""
burden[which(burden$pval<0.05),"signif"] = "*"
burden[which(burden$pval<0.05 & burden$cancer=="Pan-Cancer"),"signif"] = "**"
burden[which(burden$pval<0.05/10 & burden$treatment=="IO" & burden$cancer!="Pan-Cancer"),"signif"] = "***"
burden[which(burden$pval<0.05/17 & burden$treatment=="non-IO" & burden$cancer!="Pan-Cancer"),"signif"] = "***"

burden[which(burden$burden=="TMBH"),"signif"] = ""

burden[which(burden$cancer%in% c("Endometrial_Cancer","Leukemia","Non-Small_Cell_Lung_Cancer", "Pancreatic_Cancer","Prostate_Cancer") 
             & burden$burden=="TMBH" & burden$treatment=="non-IO" & burden$pval<0.05),"signif"] = "*"
burden[which(burden$cancer%in% c("Endometrial_Cancer","Leukemia","Non-Small_Cell_Lung_Cancer", "Pancreatic_Cancer","Prostate_Cancer") 
             & burden$burden=="TMBH" & burden$treatment=="non-IO" & burden$pval<0.05/5),"signif"] = "***"
burden[which(burden$cancer=="Melanoma" & burden$pval<0.05 & burden$treatment=="IO"),"signif"] = "***"
write.table(burden, "~/Desktop/Cancer/Burden/Paper/Tables/final_allburdenOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(burden[which(burden$signif %in% c("**", "***")),], "~/Desktop/Cancer/Burden/Paper/Tables/final_sigburdenOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)


clinvar = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_IO_noburden.csv")
clinvar["treatment"] = "IO"
clinvar2 = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_nonIO_noburden.csv")
clinvar2["treatment"] = "non-IO"
clinvar = rbind(clinvar, clinvar2)
clinvar = subset(clinvar, select = c(treatment, cancer, outcome, coef, pval, N))
clinvar["zscore"] = qnorm(clinvar$pval/2, lower.tail = F)
clinvar[which(clinvar$coef<0),"zscore"] = -1 * clinvar[which(clinvar$coef<0),"zscore"] 
clinvar["stderr"] = clinvar$coef/clinvar$zscore
clinvar = clinvar[which(clinvar$N>=40),]
clinvar = clinvar[which(clinvar$outcome %in% c("AGE","SEX","Metastatic", "NWSE","AJnonAJ")),]
results = NULL
for(t in unique(clinvar$treatment)){
  for(o in unique(clinvar$outcome)){
    temp=clinvar[which(clinvar$treatment==t & clinvar$outcome==o),] 
    z = rma(yi = temp$coef, sei = temp$stderr, method = 'FE')
    results = rbind(results, c(t, "Pan-Cancer", o, z$beta, z$pval, NA, z$zval, z$se))
  }
}
results = data.frame(results)
colnames(results) = colnames(clinvar)
clinvar = rbind(clinvar, results)
clinvar$pval = as.numeric(clinvar$pval)
clinvar$coef = as.numeric(clinvar$coef)
clinvar$stderr = as.numeric(clinvar$stderr)
clinvar["HR"] = exp(clinvar$coef)
write.table(clinvar, "~/Desktop/Cancer/Burden/Paper/Tables/final_allclinOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)

triples = list(c("Cancer_of_Unknown_Primary", "AGE", "QNORMTMB", "cov"), c("Head_and_Neck_Carcinoma","AGE", "QNORMTMB", "cov"), 
               c("Soft_Tissue_Sarcoma", "AGE", "QNORMTMB", "cov"), c("Pan-Cancer", "AGE", "QNORMTMB", "cov"), 
               c("Glioma", "AGE", "QNORMCNB", "cov"), c("Ovarian_Cancer", "AGE", "QNORMCNB", "cov"), c("Pan-Cancer", "AGE", "QNORMCNB", "cov"),
               c("Endometrial_Cancer", "AGE", "QNORMAllCNB", "cov"), c("Glioma", "AGE", "QNORMAllCNB", "cov"), 
               c("Ovarian_Cancer", "AGE", "QNORMAllCNB", "cov"), c("Soft_Tissue_Sarcoma", "AGE", "QNORMAllCNB", "cov"), 
               c("Pan-Cancer", "AGE", "QNORMAllCNB", "cov"), c("Melanoma", "SEX", "QNORMTMB", "cov"), 
               c("Pan-Cancer", "SEX", "QNORMTMB", "cov"), c("Pan-Cancer", "SEX", "QNORMCNB", "cov"), c("Esophagogastric_Carcinoma", "SEX", "QNORMAllCNB", "cov"),  
               c("Pan-Cancer", "SEX", "QNORMAllCNB", "cov"), c("Breast_Carcinoma", "Metastatic", "QNORMTMB", "cov"), 
               c("Pan-Cancer", "Metastatic", "QNORMTMB", "cov"), c("Breast_Carcinoma", "Metastatic", "QNORMAllCNB", "cov"), 
               c("Cancer_of_Unknown_Primary", "Metastatic", "QNORMAllCNB", "cov"), c("Colorectal_Cancer", "Metastatic", "QNORMAllCNB", "cov"), 
               c("Non-Small_Cell_Lung_Cancer", "Metastatic", "QNORMAllCNB", "cov"), c("Prostate_Cancer", "Metastatic", "QNORMAllCNB", "cov"), 
               c("Soft_Tissue_Sarcoma", "Metastatic", "QNORMAllCNB", "cov"), c("Pan-Cancer", "Metastatic", "QNORMAllCNB", "cov"), 
               c("Colorectal_Cancer", "AJnonAJ", "QNORMTMB", "anc"), c("Endometrial_Cancer", "AJnonAJ", "QNORMTMB", "anc"), 
               c("Leukemia", "AJnonAJ", "QNORMTMB", "anc"), c("Non-Small_Cell_Lung_Cancer", "AJnonAJ", "QNORMTMB", "anc"), 
               c("Pan-Cancer", "AJnonAJ", "QNORMTMB", "anc"), c("Breast_Carcinoma", "NWSE", "QNORMTMB", "anc"), c("Colorectal_Cancer", "NWSE", "QNORMTMB", "anc"), 
               c("Glioma", "NWSE", "QNORMTMB", "anc"), c("Head_and_Neck_Carcinoma", "NWSE", "QNORMTMB", "anc"), 
               c("Melanoma", "NWSE", "QNORMTMB", "anc"), c("Non-Hodgkin_Lymphoma", "NWSE", "QNORMTMB", "anc"), 
               c("Ovarian_Cancer", "NWSE", "QNORMTMB", "anc"), c("Pancreatic_Cancer", "NWSE", "QNORMTMB", "anc"),
               c("Prostate_Cancer", "NWSE", "QNORMTMB", "anc"), c("Soft_Tissue_Sarcoma", "NWSE", "QNORMTMB", "anc"), 
               c("Pan-Cancer", "NWSE", "QNORMTMB", "anc"))

triples = data.frame(t(data.frame(triples)))
rownames(triples) = 1:dim(triples)[1]
colnames(triples) = c("Cancer", "Outcome", "Burden", "Source")
triples$Cancer = as.character(triples$Cancer)
triples$Outcome = as.character(triples$Outcome)
triples$Burden = as.character(triples$Burden)
triples = triples[,1:2]
triples = unique(triples)
clin_results = NULL
for(t in 1:dim(triples)[1]){
  temp = clinvar[which(clinvar$cancer==triples[t,"Cancer"] & clinvar$outcome==triples[t,"Outcome"]),]
  clin_results = rbind(clin_results, temp)
}

clin_results$pval = as.numeric(clin_results$pval)

clin_results["signif"] = ""
clin_results[which(clin_results$pval<0.05),"signif"] = "*"

clin_results[which(clin_results$treatment=="IO" & clin_results$pval<0.05/9 & clin_results$outcome=="AGE"),"signif"] = "***"
clin_results[which(clin_results$treatment=="IO" & clin_results$pval<0.05/2 & clin_results$outcome=="AJnonAJ"),"signif"] = "***"
clin_results[which(clin_results$treatment=="IO" & clin_results$pval<0.05/6 & clin_results$outcome=="Metastatic"),"signif"] = "***"
clin_results[which(clin_results$treatment=="IO" & clin_results$pval<0.05/6 & clin_results$outcome=="NWSE"),"signif"] = "***"
clin_results[which(clin_results$treatment=="IO" & clin_results$pval<0.05/7 & clin_results$outcome=="SEX"),"signif"] = "***"

clin_results[which(clin_results$treatment=="non-IO" & clin_results$pval<0.05/12 & clin_results$outcome=="AGE"),"signif"] = "***"
clin_results[which(clin_results$treatment=="non-IO" & clin_results$pval<0.05/5 & clin_results$outcome=="AJnonAJ"),"signif"] = "***"
clin_results[which(clin_results$treatment=="non-IO" & clin_results$pval<0.05/9 & clin_results$outcome=="Metastatic"),"signif"] = "***"
clin_results[which(clin_results$treatment=="non-IO" & clin_results$pval<0.05/11 & clin_results$outcome=="NWSE"),"signif"] = "***"
clin_results[which(clin_results$treatment=="non-IO" & clin_results$pval<0.05/5 & clin_results$outcome=="SEX"),"signif"] = "***"

clin_results[which(clin_results$pval<0.05 & clin_results$cancer=="Pan-Cancer"),"signif"] = "**"

sig_clin = clin_results[which(clin_results$signif%in%c("***","**")),]
write.table(sig_clin, "~/Desktop/Cancer/Burden/Paper/Tables/final_sigclinOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)


io = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_IO_burden.csv")
io = io[which(io$analysis=="Survival ~ outcome + burden + covars" & io$variable=="outcome"),]
io["treatment"] = "IO"
io = subset(io, select = c(treatment, cancer, burden, outcome, coef, pval, N))
io["zscore"] = qnorm(io$pval/2, lower.tail = F)
io[which(io$coef<0),"zscore"] = -1 * io[which(io$coef<0),"zscore"] 
io["stderr"] = io$coef/io$zscore
io = io[which(io$N>=40),]

nonio = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_nonIO_burden.csv")
nonio = nonio[which(nonio$analysis=="Survival ~ outcome + burden + covars" & nonio$variable=="outcome"),]
nonio["treatment"] = "non-IO"
nonio = subset(nonio, select = c(treatment, cancer, burden, outcome, coef, pval, N))
nonio["zscore"] = qnorm(nonio$pval/2, lower.tail = F)
nonio[which(nonio$coef<0),"zscore"] = -1 * nonio[which(nonio$coef<0),"zscore"] 
nonio["stderr"] = nonio$coef/nonio$zscore
nonio = nonio[which(nonio$N>=40),]

clinvar = rbind(io, nonio)
clinvar = clinvar[which(clinvar$outcome %in% c("AGE","SEX","Metastatic", "AJnonAJ", "NWSE")),]
results = NULL
`%notin%` <- Negate(`%in%`)
for(t in unique(clinvar$treatment)){
  for(b in unique(clinvar$burden)){
    for(o in unique(clinvar$outcome)){
      if((o %notin% c("Metastatic", "AJnonAJ", "NWSE") | b!="QNORMCNB") & (o %notin% c("AJnonAJ", "NWSE") | b!="QNORMAllCNB")){
        temp = clinvar[which(clinvar$treatment==t & clinvar$outcome==o & clinvar$burden==b),]
        z = rma(yi = temp$coef, sei = temp$stderr, method = 'FE')
        results = rbind(results, c(t, "Pan-Cancer", b, o, z$beta, z$pval, NA, z$zval, z$se))
      }
    }
  }
}
results = data.frame(results)
colnames(results) = colnames(clinvar)
clinvar = rbind(clinvar, results)
clinvar$pval = as.numeric(clinvar$pval)
clinvar$coef = as.numeric(clinvar$coef)
clinvar$stderr = as.numeric(clinvar$stderr)
clinvar["HR"] = exp(clinvar$coef)

triples = list(c("Cancer_of_Unknown_Primary", "AGE", "QNORMTMB", "cov"), c("Head_and_Neck_Carcinoma","AGE", "QNORMTMB", "cov"), 
               c("Soft_Tissue_Sarcoma", "AGE", "QNORMTMB", "cov"), c("Pan-Cancer", "AGE", "QNORMTMB", "cov"), 
               c("Glioma", "AGE", "QNORMCNB", "cov"), c("Ovarian_Cancer", "AGE", "QNORMCNB", "cov"), c("Pan-Cancer", "AGE", "QNORMCNB", "cov"),
               c("Endometrial_Cancer", "AGE", "QNORMAllCNB", "cov"), c("Glioma", "AGE", "QNORMAllCNB", "cov"), 
               c("Ovarian_Cancer", "AGE", "QNORMAllCNB", "cov"), c("Soft_Tissue_Sarcoma", "AGE", "QNORMAllCNB", "cov"), 
               c("Pan-Cancer", "AGE", "QNORMAllCNB", "cov"), c("Melanoma", "SEX", "QNORMTMB", "cov"), 
               c("Pan-Cancer", "SEX", "QNORMTMB", "cov"), c("Pan-Cancer", "SEX", "QNORMCNB", "cov"), c("Esophagogastric_Carcinoma", "SEX", "QNORMAllCNB", "cov"),  
               c("Pan-Cancer", "SEX", "QNORMAllCNB", "cov"), c("Breast_Carcinoma", "Metastatic", "QNORMTMB", "cov"), 
               c("Pan-Cancer", "Metastatic", "QNORMTMB", "cov"), c("Breast_Carcinoma", "Metastatic", "QNORMAllCNB", "cov"), 
               c("Cancer_of_Unknown_Primary", "Metastatic", "QNORMAllCNB", "cov"), c("Colorectal_Cancer", "Metastatic", "QNORMAllCNB", "cov"), 
               c("Non-Small_Cell_Lung_Cancer", "Metastatic", "QNORMAllCNB", "cov"), c("Prostate_Cancer", "Metastatic", "QNORMAllCNB", "cov"), 
               c("Soft_Tissue_Sarcoma", "Metastatic", "QNORMAllCNB", "cov"), c("Pan-Cancer", "Metastatic", "QNORMAllCNB", "cov"), 
               c("Colorectal_Cancer", "AJnonAJ", "QNORMTMB", "anc"), c("Endometrial_Cancer", "AJnonAJ", "QNORMTMB", "anc"), 
               c("Leukemia", "AJnonAJ", "QNORMTMB", "anc"), c("Non-Small_Cell_Lung_Cancer", "AJnonAJ", "QNORMTMB", "anc"), 
               c("Pan-Cancer", "AJnonAJ", "QNORMTMB", "anc"), c("Breast_Carcinoma", "NWSE", "QNORMTMB", "anc"), c("Colorectal_Cancer", "NWSE", "QNORMTMB", "anc"), 
               c("Glioma", "NWSE", "QNORMTMB", "anc"), c("Head_and_Neck_Carcinoma", "NWSE", "QNORMTMB", "anc"), 
               c("Melanoma", "NWSE", "QNORMTMB", "anc"), c("Non-Hodgkin_Lymphoma", "NWSE", "QNORMTMB", "anc"), 
               c("Ovarian_Cancer", "NWSE", "QNORMTMB", "anc"), c("Pancreatic_Cancer", "NWSE", "QNORMTMB", "anc"),
               c("Prostate_Cancer", "NWSE", "QNORMTMB", "anc"), c("Soft_Tissue_Sarcoma", "NWSE", "QNORMTMB", "anc"), 
               c("Pan-Cancer", "NWSE", "QNORMTMB", "anc"))

triples = data.frame(t(data.frame(triples)))
rownames(triples) = 1:dim(triples)[1]
colnames(triples) = c("Cancer", "Outcome", "Burden", "Source")
triples$Cancer = as.character(triples$Cancer)
triples$Outcome = as.character(triples$Outcome)
triples$Burden = as.character(triples$Burden)
triples = triples[,1:3]
triples = unique(triples)
clin_results = NULL
for(t in 1:dim(triples)[1]){
  temp = clinvar[which(clinvar$cancer==triples[t,"Cancer"] & clinvar$outcome==triples[t,"Outcome"] & clinvar$burden==triples[t,"Burden"]),]
  clin_results = rbind(clin_results, temp)
}
clinvar = clin_results
clin_results = NULL
for(t in 1:dim(sig_clin)[1]){
  temp = clinvar[which(clinvar$treatment==sig_clin[t,"treatment"] & clinvar$cancer==sig_clin[t,"cancer"] & clinvar$outcome==sig_clin[t,"outcome"]),]
  clin_results = rbind(clin_results, temp)
}
write.table(clin_results, "~/Desktop/Cancer/Burden/Paper/Tables/final_sigclinburdenOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)


triples = rbind(c("Melanoma", "Tanning_Per", "QNORMTMB", "prs"), c("Non-Small_Cell_Lung_Cancer", "CigsPerDay_Per", "QNORMTMB", "prs"),
c("Non-Small_Cell_Lung_Cancer", "EduYears_Per", "QNORMTMB", "prs"), c("Pan-Cancer", "EduYears_Pan", "QNORMTMB", "prs"),
c("Pan-Cancer", "CigsPerDay_Pan", "QNORMTMB", "prs"), c("Pan-Cancer", "WBC_Pan", "QNORMTMB", "prs"),
c("Bladder_Cancer", "PRCA_Per", "QNORMAllCNB", "prs"), c("Colorectal_Cancer", "DrinksPerWeek_Per", "QNORMAllCNB", "prs"), 
c("Endometrial_Cancer", "EduYears_Per", "QNORMAllCNB", "prs"), c("Melanoma", "Smoker_Per", "QNORMAllCNB", "prs"), 
c("Non-Hodgkin_Lymphoma", "AutoimmuneSure_Per", "QNORMAllCNB", "prs"), c("Pan-Cancer", "AutoimmuneSure_Pan", "QNORMAllCNB", "prs"))

triples = data.frame(triples)
rownames(triples) = 1:dim(triples)[1]
colnames(triples) = c("Cancer", "Outcome", "Burden", "Source")
triples$Cancer = as.character(triples$Cancer)
triples$Outcome = as.character(triples$Outcome)
triples$Burden = as.character(triples$Burden)

clinvar = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_IO_noburden.csv")
clinvar["treatment"] = "IO"
clinvar2 = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_nonIO_noburden.csv")
clinvar2["treatment"] = "non-IO"
clinvar = rbind(clinvar, clinvar2)
clinvar = subset(clinvar, select = c(treatment, cancer, outcome, coef, pval, N))
clinvar["zscore"] = qnorm(clinvar$pval/2, lower.tail = F)
clinvar[which(clinvar$coef<0),"zscore"] = -1 * clinvar[which(clinvar$coef<0),"zscore"] 
clinvar["stderr"] = clinvar$coef/clinvar$zscore
clinvar = clinvar[which(clinvar$N>=40),]
clinvar = clinvar[which(clinvar$outcome %in% c("AutoimmuneSure_Pan", "CigsPerDay_Pan", "EduYears_Pan", "WBC_Pan")),]
results = NULL
for(t in unique(clinvar$treatment)){
  for(o in unique(clinvar$outcome)){
    temp=clinvar[which(clinvar$treatment==t & clinvar$outcome==o),] 
    z = rma(yi = temp$coef, sei = temp$stderr, method = 'FE')
    results = rbind(results, c(t, "Pan-Cancer", o, z$beta, z$pval, NA, z$zval, z$se))
  }
}
results = data.frame(results)
colnames(results) = colnames(clinvar)
clinvar = rbind(clinvar, results)
clinvar$pval = as.numeric(clinvar$pval)
clinvar$coef = as.numeric(clinvar$coef)
clinvar$stderr = as.numeric(clinvar$stderr)
clinvar["HR"] = exp(clinvar$coef)
write.table(clinvar, "~/Desktop/Cancer/Burden/Paper/Tables/final_allprsOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(clinvar[which(clinvar$pval<0.05),], "~/Desktop/Cancer/Burden/Paper/Tables/final_sigprsOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)

io = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_IO_burden.csv")
io = io[which(io$analysis=="Survival ~ outcome + burden + covars" & io$variable=="outcome"),]
io["treatment"] = "IO"
io = subset(io, select = c(treatment, cancer, burden, outcome, coef, pval, N))
io["zscore"] = qnorm(io$pval/2, lower.tail = F)
io[which(io$coef<0),"zscore"] = -1 * io[which(io$coef<0),"zscore"] 
io["stderr"] = io$coef/io$zscore
io = io[which(io$N>=40),]

nonio = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_nonIO_burden.csv")
nonio = nonio[which(nonio$analysis=="Survival ~ outcome + burden + covars" & nonio$variable=="outcome"),]
nonio["treatment"] = "non-IO"
nonio = subset(nonio, select = c(treatment, cancer, burden, outcome, coef, pval, N))
nonio["zscore"] = qnorm(nonio$pval/2, lower.tail = F)
nonio[which(nonio$coef<0),"zscore"] = -1 * nonio[which(nonio$coef<0),"zscore"] 
nonio["stderr"] = nonio$coef/nonio$zscore
nonio = nonio[which(nonio$N>=40),]

clinvar = rbind(io, nonio)
clinvar = clinvar[which(clinvar$outcome %in% c("AutoimmuneSure_Pan", "CigsPerDay_Pan", "EduYears_Pan", "WBC_Pan")),]
results = NULL
`%notin%` <- Negate(`%in%`)
for(t in unique(clinvar$treatment)){
  for(o in c("EduYears_Pan", "CigsPerDay_Pan", "WBC_Pan")){
    temp = clinvar[which(clinvar$treatment==t & clinvar$outcome==o & clinvar$burden=="QNORMTMB"),]
    z = rma(yi = temp$coef, sei = temp$stderr, method = 'FE')
    results = rbind(results, c(t, "Pan-Cancer", b, o, z$beta, z$pval, NA, z$zval, z$se))
  }
  o="AutoimmuneSure_Pan"
  temp = clinvar[which(clinvar$treatment==t & clinvar$outcome==o & clinvar$burden=="QNORMAllCNB"),]
  z = rma(yi = temp$coef, sei = temp$stderr, method = 'FE')
  results = rbind(results, c(t, "Pan-Cancer", b, o, z$beta, z$pval, NA, z$zval, z$se))
}
results = data.frame(results)
colnames(results) = colnames(clinvar)
clinvar = rbind(clinvar, results)
clinvar$pval = as.numeric(clinvar$pval)
clinvar$coef = as.numeric(clinvar$coef)
clinvar$stderr = as.numeric(clinvar$stderr)
clinvar["HR"] = exp(clinvar$coef)
write.table(clinvar, "~/Desktop/Cancer/Burden/Paper/Tables/final_allprsburdenOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(clinvar[which(clinvar$pval<0.05),], "~/Desktop/Cancer/Burden/Paper/Tables/final_sigprsburdenOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)

result_pancancer_AJ_IO <- read.csv("~/Downloads/result_pancancer_AJ_IO.csv")
result_pancancer_AJ_IO = result_pancancer_AJ_IO[which(result_pancancer_AJ_IO$analysis=="Survival ~ outcome + burden + outcome:burden + covars" & result_pancancer_AJ_IO$variable=="outcome:burden"),]
result_pancancer_AJ_IO = subset(result_pancancer_AJ_IO, select = -c(AJ))
io = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_IO_burden.csv")
io = io[which(io$analysis=="Survival ~ outcome + burden + outcome:burden + covars" & io$variable=="outcome:burden"),]
io = rbind(io, result_pancancer_AJ_IO)
io["treatment"] = "IO"
io = subset(io, select = c(treatment, cancer, burden, outcome, coef, pval, N))
io["zscore"] = qnorm(io$pval/2, lower.tail = F)
io[which(io$coef<0),"zscore"] = -1 * io[which(io$coef<0),"zscore"] 
io["stderr"] = io$coef/io$zscore
io = io[which(io$N>=40),]

result_pancancer_AJ_nonIO <- read.csv("~/Downloads/result_pancancer_AJ_nonIO.csv")
result_pancancer_AJ_nonIO = result_pancancer_AJ_nonIO[which(result_pancancer_AJ_nonIO$analysis=="Survival ~ outcome + burden + outcome:burden + covars" & result_pancancer_AJ_nonIO$variable=="outcome:burden"),]
result_pancancer_AJ_nonIO = subset(result_pancancer_AJ_nonIO, select = -c(AJ))
nonio = read.csv("~/Desktop/Cancer/survival_analyses/survival_pancancer_nonIO_burden.csv")
nonio = nonio[which(nonio$analysis=="Survival ~ outcome + burden + outcome:burden + covars" & nonio$variable=="outcome:burden"),]
nonio = rbind(nonio, result_pancancer_AJ_nonIO)
nonio["treatment"] = "non-IO"
nonio = subset(nonio, select = c(treatment, cancer, burden, outcome, coef, pval, N))
nonio["zscore"] = qnorm(nonio$pval/2, lower.tail = F)
nonio[which(nonio$coef<0),"zscore"] = -1 * nonio[which(nonio$coef<0),"zscore"] 
nonio["stderr"] = nonio$coef/nonio$zscore
nonio = nonio[which(nonio$N>=40),]


clinvar = rbind(io, nonio)
clinvar = clinvar[which(clinvar$outcome %in% c("NWSE", "AJnonAJ", "AutoimmuneSure_Pan", "CigsPerDay_Pan", "EduYears_Pan", "WBC_Pan")),]

results = NULL
`%notin%` <- Negate(`%in%`)
for(t in unique(clinvar$treatment)){
  for(o in c("EduYears_Pan", "CigsPerDay_Pan", "WBC_Pan", "NWSE", "AJnonAJ")){
    temp = clinvar[which(clinvar$treatment==t & clinvar$outcome==o & clinvar$burden=="QNORMTMB"),]
    z = rma(yi = temp$coef, sei = temp$stderr, method = 'FE')
    results = rbind(results, c(t, "Pan-Cancer", "QNORMTMB", o, z$beta, z$pval, NA, z$zval, z$se))
  }
  o="AutoimmuneSure_Pan"
  temp = clinvar[which(clinvar$treatment==t & clinvar$outcome==o & clinvar$burden=="QNORMAllCNB"),]
  z = rma(yi = temp$coef, sei = temp$stderr, method = 'FE')
  results = rbind(results, c(t, "Pan-Cancer", "QNORMAllCNB", o, z$beta, z$pval, NA, z$zval, z$se))
}
results = data.frame(results)
colnames(results) = colnames(clinvar)
clinvar = rbind(clinvar, results)
clinvar$pval = as.numeric(clinvar$pval)
clinvar$coef = as.numeric(clinvar$coef)
clinvar$stderr = as.numeric(clinvar$stderr)
clinvar["HR"] = exp(clinvar$coef)
write.table(clinvar, "~/Desktop/Cancer/Burden/Paper/Tables/final_allsinteractOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(clinvar[which(clinvar$pval<0.05),], "~/Desktop/Cancer/Burden/Paper/Tables/final_sigsinteractOS.txt", col.names = T, row.names = F, sep = "\t", quote = F)




survival_percancer_IO_noburden <- read.csv("~/Desktop/Cancer/survival_analyses/survival_percancer_IO_noburden.csv")
survival_percancer_IO_noburden = survival_percancer_IO_noburden[which(survival_percancer_IO_noburden$outcome %in% c("Tanning_Per","Smoker_Per","DrinksPerWeek_Per","CigsPerDay_Per","PRCA_Per", "EduYears_Per", "AutoimmuneSure_Per")),]
survival_percancer_IO_noburden["treatment"] = "IO"

survival_percancer_nonIO_noburden <- read.csv("~/Desktop/Cancer/survival_analyses/survival_percancer_nonIO_noburden.csv")
survival_percancer_nonIO_noburden = survival_percancer_nonIO_noburden[which(survival_percancer_nonIO_noburden$outcome %in% c("Tanning_Per","Smoker_Per","DrinksPerWeek_Per","CigsPerDay_Per","PRCA_Per", "EduYears_Per", "AutoimmuneSure_Per")),]
survival_percancer_nonIO_noburden["treatment"] = "nonIO"

percancer_noburden = rbind(survival_percancer_IO_noburden, survival_percancer_nonIO_noburden)

survival_percancer_IO_burden <- read.csv("~/Desktop/Cancer/survival_analyses/survival_percancer_IO_burden.csv")
survival_percancer_IO_burden = survival_percancer_IO_burden[which(survival_percancer_IO_burden$outcome %in% c("Tanning_Per","Smoker_Per","DrinksPerWeek_Per","CigsPerDay_Per","PRCA_Per", "EduYears_Per", "AutoimmuneSure_Per")),]
survival_percancer_IO_burden["treatment"] = "IO"

survival_percancer_nonIO_burden <- read.csv("~/Desktop/Cancer/survival_analyses/survival_percancer_nonIO_burden.csv")
survival_percancer_nonIO_burden = survival_percancer_nonIO_burden[which(survival_percancer_nonIO_burden$outcome %in% c("Tanning_Per","Smoker_Per","DrinksPerWeek_Per","CigsPerDay_Per","PRCA_Per", "EduYears_Per", "AutoimmuneSure_Per")),]
survival_percancer_nonIO_burden["treatment"] = "nonIO"

percancer_burden = rbind(survival_percancer_IO_burden, survival_percancer_nonIO_burden)
percancer_burden = percancer_burden[which(percancer_burden$variable=="outcome:burden"),]

NSCLC_surv_IO <- read.csv("~/Downloads/NSCLC_surv_IO (1).csv", stringsAsFactors = F)
NSCLC_INPUT <- read.delim("~/Desktop/Cancer/Burden/Survival/Non-Small_Cell_Lung_Cancer_INPUT.txt", stringsAsFactors = F)

demographics <- read.csv("~/Desktop/Cancer/Burden/Profile/INPUT/demographics.csv", stringsAsFactors = F)
colnames(demographics)[1] = "SAMPLE_ID"
colnames(demographics)[3] = "RELIGION"
demographics["College"] = NA
demographics[which(demographics$EDUCATION_LEVEL_NM %in% c("8TH GRADE OR LESS", "SOME HIGH SCHOOL",  "OBTAINED GED",
                                                          "SOME TECHNICAL PROGRAM", "SOME VOCATIONAL PROGRAM", "SOME COLLEGE")),"College"] = 0
demographics[which(demographics$EDUCATION_LEVEL_NM %in% c("GRADUATED - COLLEGE", "GRADUATED - GRAD SCHOOL", "GRADUATED - POST GRADUATE")),"College"] = 1
demographics$SAMPLE_ID = paste(demographics$SAMPLE_ID,"S1", sep = "_")
demographics = subset(demographics, select = c(SAMPLE_ID, College))

nsclc = merge(NSCLC_surv_IO, NSCLC_INPUT, by = "SAMPLE_ID")
nsclc = merge(nsclc, demographics, by = "SAMPLE_ID")
rm(demographics, NSCLC_INPUT, NSCLC_surv_IO)
nsclc = subset(nsclc, select = c(tstart, tstop, event, QNORMTMB, `TMB.H`, AGE, SEX, Metastatic, PANEL_VERSION, TUMOR_PURITY, 
                                 PC1, PC2, PC3, PC4, PC5, EduYears_Per, College))
colnames(nsclc) = c("tstart", "tstop", "event", "TMB", "TMBH", "Age", "Sex", "Metastatic", "Panel_Version", "Tumor_Purity", 
                    "PC1", "PC2", "PC3", "PC4", "PC5", "EduYears", "College")

nsclc$Sex = factor(nsclc$Sex, levels = c("Male","Female"))
nsclc$Metastatic = factor(nsclc$Metastatic, levels = c(0,1), labels = c("Primary","Metastatic"))
nsclc$Panel_Version = factor(nsclc$Panel_Version, levels = c(1,2,3))
nsclc$College = factor(nsclc$College, levels = c(0,1), labels = c("No","Yes"))
library(survival)

summary(coxph(Surv(tstart, tstop, event) ~ TMB*EduYears + College + Age + Sex + Metastatic + Panel_Version + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, nsclc))

Melanoma_surv_IO <- read.csv("~/Downloads/Melanoma_surv_IO.csv", stringsAsFactors = F)
Melanoma_INPUT <- read.delim("~/Desktop/Cancer/Burden/Survival/Melanoma_INPUT.txt", stringsAsFactors = F)
melanoma = merge(Melanoma_INPUT, Melanoma_surv_IO, by = "SAMPLE_ID")

melanoma = subset(melanoma, select = c(tstart, tstop, event, QNORMAllCNB, AGE, SEX, Metastatic, PANEL_VERSION, TUMOR_PURITY, 
                                 PC1, PC2, PC3, PC4, PC5, Smoker_Per, CigsPerDay_Pan))
colnames(melanoma) = c("tstart", "tstop", "event", "AllCNB", "Age", "Sex", "Metastatic", "Panel_Version", "Tumor_Purity", 
                    "PC1", "PC2", "PC3", "PC4", "PC5", "Smoker", "Cigs")
melanoma$Sex = factor(melanoma$Sex, levels = c("Male","Female"))
melanoma$Metastatic = factor(melanoma$Metastatic, levels = c(0,1), labels = c("Primary","Metastatic"))
melanoma$Panel_Version = factor(melanoma$Panel_Version, levels = c(1,2,3))


summary(coxph(Surv(tstart, tstop, event) ~ Smoker*AllCNB + Cigs + Age + Sex + Metastatic + Panel_Version + Tumor_Purity + PC1 + PC2 + PC3 + PC4 + PC5, melanoma))

