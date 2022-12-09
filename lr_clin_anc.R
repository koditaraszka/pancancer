library(metafor)
outcomes = NULL
for(i in name){
  x = read.table(paste0("~/Desktop/Cancer/Burden/Profile/INPUT/",i,"/",i,"_phenos.txt"), head = T)
  x["Cancer"] = i
  outcomes = rbind(outcomes, x)
}
nw = read.delim("~/Desktop/Cancer/Burden/Profile/INPUT/PROF_2018.EA.PC1.total", stringsAsFactors = F)
nw = subset(nw, select = c("IID","SCORE1_AVG"))
colnames(nw) = c("SAMPLE_ID","ProjectedPC1")
nw = subset(nw, select = c("SAMPLE_ID", "ProjectedPC1"))
outcomes = merge(outcomes, nw, by = "SAMPLE_ID", all.x = T)
aj = read.delim("~/Desktop/Cancer/Burden/Profile/INPUT/PROF_2018.EA.PC2.total", stringsAsFactors = F)
aj = subset(aj, select = c("IID","SCORE1_AVG"))
colnames(aj) = c("SAMPLE_ID","ProjectedPC2")
aj = subset(aj, select = c("SAMPLE_ID", "ProjectedPC2"))
outcomes = merge(outcomes, aj, by = "SAMPLE_ID", all.x = T)
source('~/Desktop/Cancer/Burden/process_prs.R')
outcomes = merge(outcomes, profile_prs, by = c("SAMPLE_ID","Cancer"), all.x = T)
outcomes = outcomes[complete.cases(outcomes),]
outcomes$SEX = as.character(outcomes$SEX)
outcomes[which(outcomes$SEX=="Male"),"SEX"] = 0
outcomes[which(outcomes$SEX=="Female"),"SEX"] = 1
outcomes$CigsPerDay_Pan = -1*outcomes$CigsPerDay_Pan
orig_values = subset(outcomes, select = c("Cancer", "ProjectedPC1", "ProjectedPC2",
                                          "OrigMSSProteinTMB", "OrigMSSDeepCNB", "OrigMSSCNB", 
                                          "AGE", "SEX", "Metastatic", "PANEL_VERSION", "TUMOR_PURITY", 
                                          "PC1", "PC2", "PC3", "PC4", "PC5",
                                          "CigsPerDay_Pan", "EduYears_Pan", "AutoimmuneSure_Pan"))
colnames(orig_values) = c("Cancer", "ProjectedPC1", "ProjectedPC2", 
                          "TMB", "CNB", "CNB_All", "AGE", "SEX", "Metastatic", 
                          "Panel", "Tumor_Purity", "PC1", "PC2", "PC3", "PC4", "PC5",
                          "CigsPerDay", "EduYears", "AutoimmuneSure")
orig_values["MATCH"] = NA
orig_values["TMB-H"] = 0
orig_values[which(orig_values$TMB>=10),"TMB-H"] = 1
orig_values["Source"] = "Profile"


outcomes = subset(outcomes, select = c("SAMPLE_ID", "Cancer", "ProjectedPC1", "ProjectedPC2", "MSSProteinTMB", "MSSDeepCNB", "MSSCNB", "AGE", "PANEL_VERSION", "TUMOR_PURITY", "SEX", "Metastatic", "PC1", "PC2", "PC3", "PC4", "PC5", "CigsPerDay_Pan", "EduYears_Pan", "AutoimmuneSure_Pan", "OrigMSSTMB"))

tempus_details <- read.delim("~/Desktop/Cancer/Burden/Tempus/panCancer_details_TUMOR.txt", stringsAsFactors = F)
old_details = read.delim("~/Desktop/Cancer/Burden/Tempus/tempus_details.txt", stringsAsFactors = F)
old_details = subset(old_details, select = c(MSSCNB, CNB_All, MATCH, ProjectedPC1_tumor, ProjectedPC2_tumor, 
                                             sex, age_at_dx, metastatic, tumor_purity, assay, 
                                             PC1, PC2, PC3, PC4, PC5))
tempus_details = merge(tempus_details, old_details, by = c("sex", "age_at_dx", "metastatic", "tumor_purity", "assay", 
                                                           "PC1", "PC2", "PC3", "PC4", "PC5"))
tempus_details = subset(tempus_details, select = c("Cancer", "ProjectedPC1_tumor", "ProjectedPC2_tumor",
                                                   "MSSTMB", "MSSCNB", "CNB_All",
                                                   "age_at_dx", "sex", "metastatic", "assay", "tumor_purity", 
                                                   "PC1", "PC2", "PC3", "PC4", "PC5", "MATCH",
                                                   "GSCAN.CigarettesPerDay.PV1.sscore", "UKB_460K.cov_EDU_YEARS.PV0.sscore", 
                                                   "UKB_460K.disease_AID_SURE.PV1.sscore"))
colnames(tempus_details) = c("Cancer", "ProjectedPC1", "ProjectedPC2", 
                             "TMB", "CNB", "CNB_All", "AGE", "SEX", "Metastatic", 
                             "Panel", "Tumor_Purity", "PC1", "PC2", "PC3", "PC4", "PC5", "MATCH",
                             "CigsPerDay", "EduYears", "AutoimmuneSure")
tempus_details[which(tempus_details$Cancer=="bladder"),"Cancer"] = "Bladder_Cancer"
tempus_details[which(tempus_details$Cancer=="breast"),"Cancer"] = "Breast_Carcinoma"
tempus_details[which(tempus_details$Cancer=="colorectal"),"Cancer"] = "Colorectal_Cancer"
tempus_details[which(tempus_details$Cancer=="endometrial"),"Cancer"] = "Endometrial_Cancer"
tempus_details[which(tempus_details$Cancer=="esophagogastric"),"Cancer"] = "Esophagogastric_Carcinoma"
tempus_details[which(tempus_details$Cancer=="glioma"),"Cancer"] = "Glioma"
tempus_details[which(tempus_details$Cancer=="head"),"Cancer"] = "Head_and_Neck_Carcinoma"
tempus_details[which(tempus_details$Cancer=="lung"),"Cancer"] = "Non-Small_Cell_Lung_Cancer"
tempus_details[which(tempus_details$Cancer=="lymphoma"),"Cancer"] = "Non-Hodgkin_Lymphoma"
tempus_details[which(tempus_details$Cancer=="melanoma"),"Cancer"] = "Melanoma"
tempus_details[which(tempus_details$Cancer=="ovarian"),"Cancer"] = "Ovarian_Cancer"
tempus_details[which(tempus_details$Cancer=="pancreatic"),"Cancer"] = "Pancreatic_Cancer"
tempus_details[which(tempus_details$Cancer=="prostate"),"Cancer"] = "Prostate_Cancer"
tempus_details[which(tempus_details$Cancer=="renal"),"Cancer"] = "Renal_Cell_Carcinoma"
tempus_details[which(tempus_details$Cancer=="soft"),"Cancer"] = "Soft_Tissue_Sarcoma"
tempus_details[which(tempus_details$Cancer=="unknown"),"Cancer"] = "Cancer_of_Unknown_Primary"
tempus_details["TMB-H"] = 0
tempus_details[which(tempus_details$TMB>=10),"TMB-H"] = 1
tempus_details["Source"] = "Tempus"

data = rbind(orig_values, tempus_details)

data$SEX = factor(data$SEX, levels = c(0,1), labels = c("Male", "Female"))
data$Metastatic = factor(data$Metastatic, levels = c(0,1), labels = c("Primary", "Metastatic"))
data["AgeQuintile"] = "18-50"
data[which(data$AGE>=51),"AgeQuintile"] = "51-58"
data[which(data$AGE>=59),"AgeQuintile"] = "59-65"
data[which(data$AGE>=66),"AgeQuintile"] = "66-71"
data[which(data$AGE>=72),"AgeQuintile"] = "72-97"
data$AgeQuintile = factor(data$AgeQuintile, levels = c("18-50", "51-58", "59-65", "66-71", "72-97"))

data["AJ - non-AJ"] = "non-AJ"
data["NW - SE"] = "NW"
data[which(data$ProjectedPC2>=1e-8), "AJ - non-AJ"] = "AJ"
data[which(data$ProjectedPC2>=1e-8), "NW - SE"] = NA
data[which(data$ProjectedPC2<1e-8 & data$ProjectedPC1>=1.3e-8), "NW - SE"] = "SE"
data$`AJ - non-AJ` = factor(data$`AJ - non-AJ`, levels = c("non-AJ", "AJ"))
data$`NW - SE` = factor(data$`NW - SE`, levels = c("NW", "SE"))

data[which(data$Panel=="xT"),"Panel"] = "1"
data[which(data$Panel=="xT.v2"),"Panel"] = "2"
data[which(data$Panel=="xT.v3"),"Panel"] = "3"
data$Panel = factor(data$Panel, levels = c("1", "2", "3"))

returnDetails = function(temp){
  if(dim(temp)[1]==0){
    x = data.frame(c("SampleSize", "18-50", "51-58", "59-65", "66-71", "72-97",
                     "Metastatic", "Primary", "Female", "Male", "NW", "SE", "AJ", "non-AJ"))
    colnames(x) = "Outcome"
    x["Category"] = c("Sample",rep("Age5",5), rep("Met",2), rep("Sex",2), rep("Ancestry",4))
    x["Count"] = rep(NA, 24) 
    x["Carrier"] = rep(NA, 24) 
    x["Prop"] = rep(NA, 24) 
    x["SE"] = rep(NA, 24) 
    return (x)
  }
  total = dim(temp)[1]
  carriers = dim(temp[which(temp$`TMB-H`==1),])[1]
  #TOTAL
  #Age Quintile
  age18_50 = dim(temp[which(temp$AgeQuintile=="18-50"),])[1]
  age51_58 = dim(temp[which(temp$AgeQuintile=="51-58"),])[1]
  age59_65 = dim(temp[which(temp$AgeQuintile=="59-65"),])[1]
  age66_71 = dim(temp[which(temp$AgeQuintile=="66-71"),])[1]
  age72_97 = dim(temp[which(temp$AgeQuintile=="72-97"),])[1]
  #Metastatic
  met = dim(temp[which(temp$Metastatic=="Metastatic"),])[1]
  nomet = dim(temp[which(temp$Metastatic=="Primary"),])[1]
  #female
  female = dim(temp[which(temp$SEX=="Female"),])[1]
  male = dim(temp[which(temp$SEX=="Male"),])[1]
  # ancestry
  nw = dim(temp[which(temp$`NW - SE`=="NW"),])[1]
  se = dim(temp[which(temp$`NW - SE`=="SE"),])[1]
  aj = dim(temp[which(temp$`AJ - non-AJ`=="AJ"),])[1]
  nonaj = dim(temp[which(temp$`AJ - non-AJ`=="non-AJ"),])[1]
  #CARRIERS
  #Age Quintile
  age18_50_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$AgeQuintile=="18-50"),])[1]
  age51_58_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$AgeQuintile=="51-58"),])[1]
  age59_65_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$AgeQuintile=="59-65"),])[1]
  age66_71_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$AgeQuintile=="66-71"),])[1]
  age72_97_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$AgeQuintile=="72-97"),])[1]
  #Metastatic
  met_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$Metastatic=="Metastatic"),])[1]
  nomet_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$Metastatic=="Primary"),])[1]
  #female
  female_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$SEX=="Female"),])[1]
  male_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$SEX=="Male"),])[1]
  # ancestry
  nw_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$`NW - SE`=="NW"),])[1]
  se_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$`NW - SE`=="SE"),])[1]
  aj_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$`AJ - non-AJ`=="AJ"),])[1]
  nonaj_carrier = dim(temp[which(temp$`TMB-H`==1 & temp$`AJ - non-AJ`=="non-AJ"),])[1]
  
  x = data.frame(c("SampleSize", "18-50", "51-58", "59-65", "66-71", "72-97",
                   "Metastatic", "Primary", "Female", "Male", "NW", "SE", "AJ", "non-AJ"))
  colnames(x) = "Variable"
  x["Category"] = c("Sample",rep("Age5",5), rep("Met",2), rep("Sex",2), rep("Ancestry",4))
  x["Count"] = c(total, age18_50, age51_58, age59_65, age66_71, age72_97, met, nomet, female, male, nw, se, aj, nonaj)
  
  x["Carrier"] =c(carriers, age18_50_carrier, age51_58_carrier, age59_65_carrier, age66_71_carrier, age72_97_carrier,
                  met_carrier, nomet_carrier, female_carrier, male_carrier, nw_carrier, se_carrier, aj_carrier, nonaj_carrier)
  x["Prop"] = x$Carrier/x$Count
  x["SE"] = sqrt(x$Prop*(1-x$Prop)/x$Count)
  return(x)
}

returnGLM = function(temp){
  if(dim(temp)[1]==0){
    return(data.frame(-1)) ###HERE!!!###
  }
  if(unique(temp$Source)=="Tempus"){
    if(unique(temp$Cancer)=="Glioma"){
      glm = data.frame(summary(glm(data=temp, `TMB-H` ~ AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[,c(1,2,4)])
      glm = rbind(glm, c(NA, NA, NA))
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = c("Intercept", "Age", "Sex", "PC1", "PC2", "PC3", "PC4", "PC5", "Panel2", "Panel3", "Tumor_Purity", "MATCH", "Metastatic")
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
     } else if(unique(temp$Cancer) %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
       glm = data.frame(summary(glm(data=temp, `TMB-H` ~ AGE + Metastatic + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[,c(1,2,4)])
       glm = rbind(glm, c(NA, NA, NA))
       colnames(glm) = c("Beta", "SE", "P")
       glm["Variable"] = c("Intercept", "Age", "Metastatic", "PC1", "PC2", "PC3", "PC4", "PC5", "Panel2", "Panel3", "Tumor_Purity", "MATCH", "Sex")
       glm["Source"] = unique(temp$Source)
       glm["Cancer"] = unique(temp$Cancer)
       rownames(glm) = 1:nrow(glm)
     } else{
       glm = data.frame(summary(glm(data=temp, `TMB-H` ~ AGE + SEX + Metastatic + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[,c(1,2,4)])
       colnames(glm) = c("Beta", "SE", "P")
       glm["Variable"] = c("Intercept", "Age", "Sex", "Metastatic", "PC1", "PC2", "PC3", "PC4", "PC5", "Panel2", "Panel3", "Tumor_Purity", "MATCH")
       glm["Source"] = unique(temp$Source)
       glm["Cancer"] = unique(temp$Cancer)
       rownames(glm) = 1:nrow(glm)
     }
  } else{
    if(unique(temp$Cancer) %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
      glm = data.frame(summary(glm(data=temp, `TMB-H` ~ AGE + Metastatic + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity, family = 'binomial'))$coef[,c(1,2,4)])
      glm = rbind(glm, c(NA, NA, NA))
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = c("Intercept", "Age", "Metastatic", "PC1", "PC2", "PC3", "PC4", "PC5", "Panel2", "Panel3","Tumor_Purity", "Sex")
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    } else {
      glm = data.frame(summary(glm(data=temp, `TMB-H` ~ AGE + SEX + Metastatic + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity, family = 'binomial'))$coef[,c(1,2,4)])
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = c("Intercept", "Age", "Sex", "Metastatic", "PC1", "PC2", "PC3", "PC4", "PC5", "Panel2", "Panel3", "Tumor_Purity")
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    }
  }
  return (glm)
}    

returnGLM2 = function(temp, colname){
  if(dim(temp)[1]==0){
    return(data.frame(-1)) ###HERE!!!###
  }
  if(unique(temp$Source)=="Tempus"){
    if(unique(temp$Cancer)=="Glioma"){
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[2,c(1,2,4)]))
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = colname
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    } else if(unique(temp$Cancer) %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ AGE + Metastatic + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[2,c(1,2,4)]))
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = colname
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    } else{
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ AGE + SEX + Metastatic + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[2,c(1,2,4)]))
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = colname
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    }
  } else{
    if(unique(temp$Cancer) %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ AGE + Metastatic + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity, family = 'binomial'))$coef[2,c(1,2,4)]))
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = colname
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    } else {
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ AGE + SEX + Metastatic + PC1 + PC2 + PC3 + PC4 + PC5 + Panel + Tumor_Purity, family = 'binomial'))$coef[2,c(1,2,4)]))
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = colname
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    }
  }
  return (glm)
}  

returnGLM3 = function(temp){
  if(dim(temp)[1]==0){
    return(data.frame(-1)) ###HERE!!!###
  }
  temp$ProjectedPC1 = scale(temp$ProjectedPC1)
  temp["Z"] = temp$`NW - SE`
  temp$`NW - SE` = temp$ProjectedPC1
  temp[which(is.na(temp$Z)),"NW - SE"] = NA
  if(unique(temp$Source)=="Tempus"){
    if(unique(temp$Cancer)=="Glioma"){
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `AJ - non-AJ` + AGE + SEX + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm2 = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `NW - SE` + AGE + SEX + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm = rbind(glm, glm2)
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = c("AJ-nonAJ","NW-SE")
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    } else if(unique(temp$Cancer) %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `AJ - non-AJ` + AGE + Metastatic + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm2 = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `NW - SE` + AGE + Metastatic + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm = rbind(glm, glm2)
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = c("AJ-nonAJ","NW-SE")
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    } else{
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `AJ - non-AJ` + AGE + SEX + Metastatic + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm2 = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `NW - SE` + AGE + SEX + Metastatic + Panel + Tumor_Purity + MATCH, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm = rbind(glm, glm2)
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = c("AJ-nonAJ","NW-SE")
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    }
  } else{
    if(unique(temp$Cancer) %in% c("Breast_Carcinoma", "Endometrial_Cancer", "Ovarian_Cancer", "Prostate_Cancer")){
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `AJ - non-AJ` + AGE + Metastatic + Panel + Tumor_Purity, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm2 = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `NW - SE` + AGE + Metastatic + Panel + Tumor_Purity, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm = rbind(glm, glm2)
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = c("AJ-nonAJ","NW-SE")
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    } else {
      glm = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `AJ - non-AJ` + AGE + SEX + Metastatic + Panel + Tumor_Purity, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm2 = data.frame(t(summary(glm(data=temp, `TMB-H` ~ `NW - SE` + AGE + SEX + Metastatic + Panel + Tumor_Purity, family = 'binomial'))$coef[2,c(1,2,4)]))
      glm = rbind(glm,glm2)
      colnames(glm) = c("Beta", "SE", "P")
      glm["Variable"] = c("AJ-nonAJ","NW-SE")
      glm["Source"] = unique(temp$Source)
      glm["Cancer"] = unique(temp$Cancer)
      rownames(glm) = 1:nrow(glm)
    }
  }
  return (glm)
}

glm_details = NULL
final_glm = NULL
for(source in c("Profile", "Tempus")){
  for(cancer in name){
    temp=data[which(data$Source==source & data$Cancer==cancer),]
    if(dim(temp[which(temp$`TMB-H`==1),])[1]>10){
      dets = data.frame(returnDetails(temp))
      dets["Source"] = source
      dets["Cancer"] = cancer
      glm_details = rbind(glm_details, dets)
      anc = returnGLM3(temp)
      one = returnGLM(temp)
      temp["18-50_51-58"] = NA
      temp["51-58_59-65"] = NA
      temp["59-65_66-71"] = NA
      temp["66-71_72-97"] = NA
      temp["18-50_72-97"] = NA
      temp[which(temp$AgeQuintile=="18-50"),"18-50_51-58"] = "18-50"
      temp[which(temp$AgeQuintile=="51-58"),"18-50_51-58"] = "51-58"
      temp[which(temp$AgeQuintile=="51-58"),"51-58_59-65"] = "51-58"
      temp[which(temp$AgeQuintile=="59-65"),"51-58_59-65"] = "59-65"
      temp[which(temp$AgeQuintile=="59-65"),"59-65_66-71"] = "59-65"
      temp[which(temp$AgeQuintile=="66-71"),"59-65_66-71"] = "66-71"
      temp[which(temp$AgeQuintile=="66-71"),"66-71_72-97"] = "66-71"
      temp[which(temp$AgeQuintile=="72-97"),"66-71_72-97"] = "72-97"
      temp[which(temp$AgeQuintile=="18-50"),"18-50_72-97"] = "18-50"
      temp[which(temp$AgeQuintile=="72-97"),"18-50_72-97"] = "72-97"
      temp$`18-50_51-58` = factor(temp$`18-50_51-58`, levels = c("18-50", "51-58"))
      temp$`51-58_59-65` = factor(temp$`51-58_59-65`, levels = c("51-58", "59-65"))
      temp$`59-65_66-71` = factor(temp$`59-65_66-71`, levels = c("59-65", "66-71"))
      temp$`66-71_72-97` = factor(temp$`66-71_72-97`, levels = c("66-71", "72-97"))
      temp$`18-50_72-97` = factor(temp$`18-50_72-97`, levels = c("18-50","72-97"))
      temp$AGE = temp$`18-50_51-58`
      two = returnGLM2(temp, "18-50_51-58")
      temp$AGE = temp$`51-58_59-65`
      three = returnGLM2(temp,"51-58_59-65")
      temp$AGE = temp$`59-65_66-71`
      four = returnGLM2(temp, "59-65_66-71")
      temp$AGE = temp$`66-71_72-97`
      five = returnGLM2(temp,"66-71_72-97")
      temp$AGE = temp$`18-50_72-97`
      six = returnGLM2(temp,"18-50_72-97")
      if (length(one) != 1){
        final_glm = rbind(final_glm, one, two, three, four, five, six, anc)
      }
    }
  }
}
final_glm["Ratio"] = NA


final_glm = final_glm[-which(final_glm$Cancer=="Glioma" & final_glm$Variable%in% c("18-50_51-58","51-58_59-65", "59-65_66-71", "66-71_72-97", "18-50_72-97") & final_glm$Source=="Profile"),]
final_glm = final_glm[-which(final_glm$Cancer=="Endometrial_Cancer" & final_glm$Variable=="Metastatic" & final_glm$Source=="Profile"),]
final_glm = final_glm[-which(final_glm$Cancer %in% c("Non-Hodgkin_Lymphoma", "Colorectal_Cancer") & final_glm$Variable=="AJ-nonAJ" & final_glm$Source=="Profile"),]

pandets = NULL
for(source in unique(glm_details$Source)){
  for(outcome in unique(glm_details$Variable)){
    temp = glm_details[which(glm_details$Source==source & glm_details$Variable==outcome),]
    count = sum(temp$Count,na.rm=T)
    carrier = sum(temp$Carrier,na.rm=T)
    sqrt(x$Prop*(1-x$Prop)/x$Count)
    pandets = rbind(pandets, c(outcome, unique(temp$Category), count, carrier, carrier/count, sqrt(((carrier/count)*(1-(carrier/count)))/count), source, "Pan-Cancer"))
  }
}
colnames(pandets) = colnames(glm_details)
glm_details = rbind(glm_details, pandets)
glm_details$Count = as.numeric(glm_details$Count)
glm_details$Carrier = as.numeric(glm_details$Carrier)
glm_details$Prop = as.numeric(glm_details$Prop)
glm_details$SE = as.numeric(glm_details$SE)

pancancer = NULL
for(source in unique(final_glm$Source)){
  for(outcome in unique(final_glm$Variable)){
    if(source=="Profile" & outcome=="MATCH"){
      print("skip")
    } else{
      temp = final_glm[which(final_glm$Source==source & final_glm$Variable==outcome),]
      ratio = -1
      if(dim(temp)[1] > 0){
        ratio = max(temp$SE, na.rm = T)/min(temp$SE, na.rm = T)
      }
      z = rma(yi = temp$Beta, sei = temp$SE, method = 'FE')
      pancancer = rbind(pancancer, c(z$beta, z$se, z$pval, outcome, source, "Pan-Cancer", ratio))
    }
  }
}
colnames(pancancer) = colnames(final_glm)
final_glm = rbind(final_glm, pancancer)
final_glm$P = as.numeric(final_glm$P)
final_glm$SE = as.numeric(final_glm$SE)
final_glm$Ratio = as.numeric(final_glm$Ratio)
final_glm = subset(final_glm, select = c(Source, Cancer, Variable, Beta, SE, P))


sig_glm = final_glm[which(final_glm$P<0.05 & final_glm$Variable %in% c("Age", "18-50_51-58", "18-50_72-97", "51-58_59_65","59-65_66-71","66-71_72-97", "Metastatic", "Sex", "AJ-nonAJ", "NW-SE")),]
cancerGLM = sig_glm[which(sig_glm$Cancer=="Pan-Cancer" & sig_glm$Variable!="Age"),]
cancerGLM = rbind(cancerGLM, sig_glm[which(sig_glm$Variable=="AJ-nonAJ"|sig_glm$Variable=="NW-SE"),])

cancerGLM = cancerGLM[order(cancerGLM$Variable),]
rownames(cancerGLM) = 1:nrow(cancerGLM)
sigBinary = data.frame(c(rep("18-50",2),"66-71", rep("AJ", 2), "Primary", rep("NW",5), rep("Male",2)))
colnames(sigBinary) = "group1"
sigBinary["group2"] = c(rep("72-97",2),"72-97", rep("non-AJ", 2), "Metastatic", rep("SE",5), rep("Female",2))
sigBinary["Beta"] = round(as.numeric(cancerGLM$Beta),3)
sigBinary["OR"] = round(exp(as.numeric(cancerGLM$Beta)),2)
sigBinary["Pvalue"] = scientific(as.numeric(cancerGLM$P),3)
sigBinary["Label"] = paste0("OR=",sigBinary$OR, "; P=",sigBinary$Pvalue)
sigBinary["Category"] = c(rep("Age5",3), rep("AJnonAJ",2),"Met", rep("NWSE",5), rep("Sex",2))
sigBinary["Source"] = c("Profile","Tempus", rep("Profile",2),"Tempus","Profile","Tempus","Profile",rep("Tempus",3),"Profile","Tempus")


profile = NULL
tempus = NULL
for(i in name){
  p1 = final_glm[which(final_glm$Source=="Profile" & final_glm$Cancer==i & final_glm$Variable=="Intercept"),"Beta"]
  p2 = final_glm[which(final_glm$Source=="Profile" & final_glm$Cancer==i & final_glm$Variable=="Age"),"Beta"]
  p3 = final_glm[which(final_glm$Source=="Profile" & final_glm$Cancer==i & final_glm$Variable=="Sex"),"Beta"]
  p4 = final_glm[which(final_glm$Source=="Profile" & final_glm$Cancer==i & final_glm$Variable=="Sex"),"SampleSize"]
  p5 = dim(data[which(data$Cancer==i & data$Source=="Profile"),])[1]
  t1 = final_glm[which(final_glm$Source=="Tempus" & final_glm$Cancer==i & final_glm$Variable=="Intercept"),"Beta"]
  t2 = final_glm[which(final_glm$Source=="Tempus" & final_glm$Cancer==i & final_glm$Variable=="Age"),"Beta"]
  t3 = final_glm[which(final_glm$Source=="Tempus" & final_glm$Cancer==i & final_glm$Variable=="Sex"),"Beta"]
  t4 = final_glm[which(final_glm$Source=="Tempus" & final_glm$Cancer==i & final_glm$Variable=="Sex"),"SampleSize"]
  t5=dim(data[which(data$Cancer==i & data$Source=="Tempus"),])[1]
  profile = rbind(profile, c(i,p1,p2,p3,p4,p5))
  tempus = rbind(tempus, c(i,t1,t2,t3,t4,t5))
}
profile = as.data.frame(profile)
tempus = as.data.frame(tempus)
colnames(profile) = c("Cancer", "Intercept", "Age", "Sex", "SampleSize")
colnames(tempus) = c("Cancer", "Intercept", "Age", "Sex", "SampleSize")

profile = profile[-which(profile$Cancer %in%c("Glioma", "Esophagogastric_Carcinoma", "Leukemia", "Pancreatic_Cancer", "Prostate_Cancer", "Renal_Cell_Carcinoma")),]
profile$Intercept = as.numeric(profile$Intercept)
profile$SampleSize = as.numeric(profile$SampleSize)
profile$Sex = as.numeric(profile$Sex)
profile$Age = as.numeric(profile$Age)
profile["YoungMan"] = profile$SampleSize*(1/(1+exp(-1*(profile$Intercept + profile$Sex))))
profile["OldWoman"] = profile$SampleSize*(1/(1+exp(-1*(profile$Intercept + profile$Age))))
#sum(profile$YoungMan,na.rm=T)/sum(profile$SampleSize,na.rm=T)
#sum(profile$OldWoman,na.rm=T)/sum(profile$SampleSize,na.rm=T)
#(0.1160456-0.05236009)/0.1160456
tempus = tempus[-which(tempus$Cancer %in% c( "Leukemia", "Non-Hodgkin_Lymphoma", "Renal_Cell_Carcinoma")),]
tempus$Intercept = as.numeric(tempus$Intercept)
tempus$SampleSize = as.numeric(tempus$SampleSize)
tempus$Sex = as.numeric(tempus$Sex)
tempus$Age = as.numeric(tempus$Age)
tempus["YoungMan"] = tempus$SampleSize*(1/(1+exp(-1*(tempus$Intercept + tempus$Sex))))
tempus["OldWoman"] = tempus$SampleSize*(1/(1+exp(-1*(tempus$Intercept + tempus$Age))))
#sum(tempus$YoungMan,na.rm=T)/sum(tempus$SampleSize,na.rm=T)
#sum(tempus$OldWoman,na.rm=T)/sum(tempus$SampleSize,na.rm=T)
#(0.173409-0.1447516)/0.173409

rm(aj,anc,cancer,dets,five,four,i,nw,old_details,one,outcome,p1,p2,p3,p4,p5,t1,t2,t3,t4,t5,pancancer,pandets,
   profile_prs,ratio,returnDetails,returnGLM,returnGLM2,returnGLM3,six,source,temp,three,two,x,z)