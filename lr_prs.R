setwd("~/Desktop/Cancer/Burden/Survival/")
source("setup_input.R")

all_cancer = c("Bladder_Cancer", "Breast_Carcinoma", "Cancer_of_Unknown_Primary", "Colorectal_Cancer",
               "Endometrial_Cancer", "Esophagogastric_Carcinoma", "Glioma", "Head_and_Neck_Carcinoma",
               "Leukemia", "Melanoma", "Non-Hodgkin_Lymphoma", "Non-Small_Cell_Lung_Cancer", "Ovarian_Cancer",
               "Pancreatic_Cancer", "Prostate_Cancer", "Renal_Cell_Carcinoma", "Soft_Tissue_Sarcoma")

triples = triples[which(triples$Source=="prs"),]

process_covars = function(data, cancer, outcome){
  covars = c("AGE", "SEX", "Metastatic", "PANEL_VERSION", 
             "TUMOR_PURITY", "PC1", "PC2", "PC3", "PC4", "PC5")
  data = subset(data, select = covars)
  if(cancer %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
    data = subset(data, select = -c(SEX))
  }
  return(data)
}

process = function(data, this_row){
  x = this_row$Outcome
  x_num = which(colnames(data)==x)
  indvar = data[,x_num]
  y="TMB.H"
  y_num = which(colnames(data)==y)
  burden = data[,y_num]
  covars = process_covars(data, cancer, x)
  gen = data.frame(cbind(burden, indvar, covars))
  gen = gen[complete.cases(gen),]
  return(gen)
}

results = NULL
for(cancer in all_cancer){
  data = read.table(paste0(cancer, "_INPUT.txt"), head = T, sep = '\t')
  if(dim(data[which(data$TMB.H==1),])[1]>10){
    data$SEX = factor(data$SEX, levels = c("Male", "Female"))
    data$PANEL_VERSION = factor(data$PANEL_VERSION, levels = c("1","2","3"))
    data$Metastatic = factor(data$Metastatic, levels = c("0","1"), labels = c("Primary", "Metastatic"))
    subset=triples[which(triples$Cancer %in% c(cancer, "Pan-Cancer")),]
    for(line in 1:dim(subset)[1]){
        if(subset[line,"Outcome"] == "SEX" & cancer %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
        } else{
          final = process(data, subset[line,])
          logreg = summary(glm(data=final, burden ~ indvar + ., family = binomial(link = "logit")))$coef[2,c(1,2,4)]
          results = rbind(results, c(subset[line,"Cancer"], cancer, subset[line,"Outcome"], logreg, dim(final)[1]))
          if(subset[line,"Outcome"]=="AGE"){
            final["BinAge"] = NA
            final[which(final$indvar>=18 & final$indvar<=50),"BinAge"] = "18-50"
            final[which(final$indvar>=72 & final$indvar<=97),"BinAge"] = "72-97"
            final$BinAge = factor(final$BinAge, levels = c("18-50", "72-97"))
            final = subset(final, select = -c(indvar))
            logreg = summary(glm(data=final, burden ~ BinAge + ., family = binomial(link = "logit")))$coef[1,c(1,2,4)]
            baseline = logreg[1]
            logreg = summary(glm(data=final, burden ~ BinAge + ., family = binomial(link = "logit")))$coef[2,c(1,2,4)]
            change=logreg[1]
            newprob=exp(-1*(baseline+change))
            oldprob=exp(-1*(baseline))
            newprob = 1/(1+newprob)
            oldprob = 1/(1+oldprob)
            change=exp(change)
            baseline=exp(baseline)
            #print(head(final))
            results = rbind(results, c(subset[line,"Cancer"], cancer, "Baseline", baseline, NA, NA, dim(final)[1]))
            results = rbind(results, c(subset[line,"Cancer"], cancer, "Change", change, NA, NA, dim(final)[1]))
            results = rbind(results, c(subset[line,"Cancer"], cancer, "RatioChange", oldprob, newprob, (newprob-oldprob)/oldprob, dim(final)[1]))
          }
        }
      }
  }
}
results = data.frame(results)
pancancer = results[which(results$V1=="Pan-Cancer"),2:7]
colnames(pancancer) = c("Cancer", "IndVar", "Beta", "Stderr", "Pvalue", "SampleSize")
results = results[which(results$V1!="Pan-Cancer"),2:7]
colnames(results) = c("Cancer", "IndVar", "Beta", "Stderr", "Pvalue", "SampleSize")

pancancer$Beta=as.numeric(pancancer$Beta)
pancancer$Stderr=as.numeric(pancancer$Stderr)
pancancer$Pvalue=as.numeric(pancancer$Pvalue)

meta = NULL
for(i in unique(pancancer$IndVar)){
  temp = pancancer[which(pancancer$IndVar==i),]
  z=rma(yi=temp$Beta, sei=temp$Stderr, method = 'FE')
  meta = rbind(meta, c("Pan-Cancer", i, z$beta, z$se, z$pval,NA))
}

meta= data.frame(meta)
colnames(meta) = colnames(results)
results = rbind(results, meta)
results$Beta=as.numeric(results$Beta)
results$Stderr=as.numeric(results$Stderr)
results$Pvalue=as.numeric(results$Pvalue)

tmbh_prs = results
tmbh_prs = tmbh_prs[-which(tmbh_prs$IndVar %in% c('EduYears_Pan', 'WBC_Pan')),]
### NSCLC -- Edu 
#data = read.table(paste0("Non-Small_Cell_Lung_Cancer_INPUT.txt"), head = T, sep = '\t')
#data$SEX = factor(data$SEX, levels = c("Male", "Female"))
#data$PANEL_VERSION = factor(data$PANEL_VERSION, levels = c("1","2","3"))
#data$Metastatic = factor(data$Metastatic, levels = c("0","1"), labels = c("Primary", "Metastatic"))

#demographics <- read.csv("~/Desktop/Cancer/Burden/Profile/INPUT/demographics.csv", stringsAsFactors = F)
#colnames(demographics)[1] = "SAMPLE_ID"
#colnames(demographics)[3] = "RELIGION"
#demographics["College"] = NA
#demographics[which(demographics$EDUCATION_LEVEL_NM %in% c("8TH GRADE OR LESS", "SOME HIGH SCHOOL",  "OBTAINED GED", "GRADUATED - HIGH SCHOOL",
#                                                          "SOME TECHNICAL PROGRAM", "SOME VOCATIONAL PROGRAM", "SOME COLLEGE")),"College"] = 0
#demographics[which(demographics$EDUCATION_LEVEL_NM %in% c("GRADUATED - COLLEGE", "GRADUATED - GRAD SCHOOL", "GRADUATED - POST GRADUATE")),"College"] = 1
#demographics$SAMPLE_ID = paste(demographics$SAMPLE_ID,"S1", sep = "_")
#demographics = subset(demographics, select = c("SAMPLE_ID", "EDUCATION_LEVEL_NM", "College"))

#data = merge(data, demographics, by = "SAMPLE_ID")

#summary(lm(data = data, QNORMTMB ~ EduYears_Per + AGE + SEX + Metastatic + PANEL_VERSION + TUMOR_PURITY + PC1 + PC2 + PC3 + PC4 + PC5))

#summary(lm(data = data, QNORMTMB ~ College + AGE + SEX + Metastatic + PANEL_VERSION + TUMOR_PURITY + PC1 + PC2 + PC3 + PC4 + PC5))

#summary(glm(data = data, TMB.H ~ EduYears_Per + AGE + SEX + Metastatic + PANEL_VERSION + TUMOR_PURITY + PC1 + PC2 + PC3 + PC4 + PC5, family=binomial(link = "logit")))

#summary(glm(data = data, TMB.H ~ College + AGE + SEX + Metastatic + PANEL_VERSION + TUMOR_PURITY + PC1 + PC2 + PC3 + PC4 + PC5, family=binomial(link = "logit")))

#summary(glm(data = data, TMB.H ~ College + EduYears_Per + AGE + SEX + Metastatic + PANEL_VERSION + TUMOR_PURITY + PC1 + PC2 + PC3 + PC4 + PC5, family=binomial(link = "logit")))


