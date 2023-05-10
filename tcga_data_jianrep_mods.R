library(ggplot2)
library(readxl)
library(stringr)
#rm(list=ls())
all_theme = theme_bw() + theme(strip.background = element_blank(), 
                         axis.text =element_text( size = 12, color = "black"), 
                         axis.title = element_text(size = 12, color = "black"),
                         legend.text = element_text(size = 12, color = "black"), 
                         legend.title = element_text(size = 12, color = "black"), 
                         plot.title = element_text(hjust = 0.5, size=12, color="black"), 
                         strip.text = element_text(size = 12, color = "black"), 
                         plot.tag = element_text(size=12, color = "black", face = "bold"))

tempus_profile_panels <- read.delim("~/Desktop/Cancer/tempus_profile_panels.txt")
profile_genes = tempus_profile_panels[which(tempus_profile_panels$Source %in% c("Profile", "Both")),"gene"]
TCGA <- read.delim("~/Desktop/Cancer/Burden/TCGA/ANCESTRY/TCGA.CO.PC1.sscore", header=FALSE, comment.char="#")
TCGA = subset(TCGA, select = c(V1, V5))
colnames(TCGA) =c("FID","PC1")
TCGA_pc2 <- read.delim("~/Desktop/Cancer/Burden/TCGA/ANCESTRY/TCGA.CO.PC2.sscore", header=FALSE, comment.char="#")
TCGA_pc2 = subset(TCGA_pc2, select = c(V1, V5))
colnames(TCGA_pc2) =c("FID","PC2")
TCGA = merge(TCGA, TCGA_pc2, by = "FID")

cancers = c("Bladder_Cancer","Breast_Carcinoma","Colorectal_Cancer", "Endometrial_Cancer",
            "Esophagogastric_Carcinoma", "Glioma", "Head_and_Neck_Carcinoma", "Leukemia", "Melanoma",
            "Non-Small_Cell_Lung_Cancer", "Ovarian_Cancer", "Pancreatic_Cancer",
            "Prostate_Cancer", "Renal_Cell_Carcinoma", "Soft_Tissue_Sarcoma")

names = NULL
for(c in cancers){
  temp = read.delim(paste0("~/Desktop/Cancer/Burden/TCGA/", c, "/", c, "_final.txt"), header=FALSE)
  colnames(temp) = c("FID", "PATIENT_ID")
  temp["Cancer"] = c
  temp2 = read.delim(paste0("~/Desktop/Cancer/Burden/TCGA/", c, "/clean_patient.txt"), header=TRUE)
  temp2 = subset(temp2, select = c(PATIENT_ID, SEX, AGE, RACE, ETHNICITY))
  temp = merge(temp, temp2, by = "PATIENT_ID")
  names = rbind(names, temp)
}

TCGA = merge(TCGA, names, by = "FID")
TCGA["Metastatic"] = 0
clinical_metastatic <- read.delim("~/Desktop/Cancer/Burden/TCGA/clinical_metastatic.tsv")
clinical_metastatic = subset(clinical_metastatic, select = c(case_submitter_id))
colnames(clinical_metastatic) = "PATIENT_ID"
TCGA[which(TCGA$PATIENT_ID %in% clinical_metastatic$PATIENT_ID),"Metastatic"] = 1
TCGA[which(is.na(TCGA$RACE)),"RACE"] = "Unknown"
TCGA[which(TCGA$RACE=="WHITE"),"RACE"] = "White"
TCGA[which(TCGA$RACE=="ASIAN"),"RACE"] = "Asian"
TCGA[which(TCGA$RACE=="BLACK"),"RACE"] = "Black"
TCGA[which(TCGA$RACE=="OTHER"),"RACE"] = "Other"
TCGA[which(is.na(TCGA$ETHNICTY)),"ETHNICITY"] = "Unknown"

colnames(TCGA)[6:9] = c("Sex", "Age", "Race", "Ethnicity")
TCGA["Race/Ethnicity"] = "Unknown"
TCGA[which(TCGA$Ethnicity=="YES"),"Race/Ethnicity"] = "Hispanic"
TCGA[which(TCGA$Race=="White"),"Race/Ethnicity"] = "White"
TCGA[which(TCGA$Race=="Asian"),"Race/Ethnicity"] = "Asian"
TCGA[which(TCGA$Race=="Black"),"Race/Ethnicity"] = "Black"
TCGA[which(TCGA$Race=="Other"),"Race/Ethnicity"] = "Other"
TCGA["Alpha"] = 1
TCGA[which(TCGA$`Race/Ethnicity`=="Unknown"), "Alpha"] = 0.5
TCGA[which(TCGA$`Race/Ethnicity`=="White"), "Alpha"] = 0.75
mean_pc1 = mean(TCGA[which(TCGA$Race=="White"),"PC1"])
std_pc1 = sd(TCGA[which(TCGA$Race=="White"),"PC2"])
mean_pc2 = mean(TCGA[which(TCGA$Race=="White"),"PC2"])
std_pc2 = sd(TCGA[which(TCGA$Race=="White"),"PC2"])
#ggplot(TCGA, aes(x=PC1,y=PC2,color=`Race/Ethnicity`)) + geom_point(alpha=TCGA$Alpha) + guides(alpha='none') + geom_rect(xmin=mean_pc1-2*std_pc1, xmax=mean_pc1+2*std_pc1, ymin=mean_pc2-2*std_pc2, ymax=mean_pc2+2*std_pc2, fill = NA, color = 'black') + all_theme + xlim(c(1.58e-8,2.23e-08))
mmc5 <- read_excel("Desktop/Cancer/Burden/TCGA/mmc5.xlsx",col_types = c("text", "text", "text"), skip = 2)
msi = data.frame(mmc5)
msi = msi[,1:2]
colnames(msi) = c("PATIENT_ID", "MSIsensor_score")
msi = msi[which(msi$MSIsensor_score>=4),]
TCGA["MSI"] = 0
TCGA[which(TCGA$PATIENT_ID %in% msi$PATIENT_ID),"MSI"] = 1
TCGA["Drop"] = 0
TCGA[which(TCGA$PC1<mean_pc1-2*std_pc1 | TCGA$PC1>mean_pc1+2*std_pc1 |
           TCGA$PC2<mean_pc2-2*std_pc2 | TCGA$PC2>mean_pc2+2*std_pc2),"Drop"] = 1
#TCGA[which(TCGA$Metastatic==1 & TCGA$Cancer!="Melanoma"),"Drop"] = 1
#TCGA = TCGA[-which(TCGA$Drop==1 | TCGA$MSI==1 | TCGA$Cancer %in% c("Leukemia", "Pancreatic_Cancer", "Soft_Tissue_Sarcoma")),]
#write.table(TCGA$FID, "~/Desktop/Cancer/Burden/TCGA/eur_samples.txt", col.names = T, quote = F, row.names = F)
burden = NULL
for(c in unique(TCGA$Cancer)){
  print(c)
  somatic = read.csv(paste0("~/Desktop/Cancer/Burden/TCGA/", c, "/somatic_mutations.txt"), sep="")
  barcode = subset(somatic, select = c("PATIENT_ID","Tumor_Sample_Barcode"))
  barcode = unique(barcode)
  if(c=="Melanoma"){
    barcode = barcode[-which(barcode$Tumor_Sample_Barcode %in% c("TCGA-ER-A19T-06", "TCGA-ER-A2NF-06")),]
  }
  tmb=data.frame(table(somatic$Tumor_Sample_Barcode))
  colnames(tmb) = c("Tumor_Sample_Barcode", "TMB")
  tmb$TMB = tmb$TMB/38
  tmb["QNORMTMB"] = qnorm(rank(tmb$TMB, ties.method = "random")/(length(tmb$TMB)+1))
  
  somatic_profile = somatic[which(somatic$Hugo_Symbol %in% profile_genes),]
  tmb_profile=data.frame(table(somatic_profile$Tumor_Sample_Barcode))
  colnames(tmb_profile) = c("Tumor_Sample_Barcode", "TMB_Profile")
  tmb_profile$TMB_Profile = tmb_profile$TMB_Profile/3
  tmb = merge(tmb, tmb_profile, by = "Tumor_Sample_Barcode", all.x = T)
  tmb[which(is.na(tmb$TMB_Profile)),"TMB_Profile"] = 0
  tmb["QNORMTMB_Profile"] = qnorm(rank(tmb$TMB_Profile, ties.method = "random")/(length(tmb$TMB_Profile)+1))
  
  somatic_funct = somatic[grep("Frame_Shift_Del|Frame_Shift_Ins|Nonsense_Mutation|Nonstop_Mutation|Splice_Region|Splice_Site", somatic$Variant_Classification),]
  tmb_funct=data.frame(table(somatic_funct$Tumor_Sample_Barcode))
  colnames(tmb_funct) = c("Tumor_Sample_Barcode", "TMB_Funct")
  tmb_funct$TMB_Funct = tmb_funct$TMB_Funct/3
  tmb = merge(tmb, tmb_funct, by = "Tumor_Sample_Barcode", all.x = T)
  tmb = merge(tmb, barcode, by = "Tumor_Sample_Barcode")
  tmb[which(is.na(tmb$TMB_Funct)),"TMB_Funct"] = 0
  tmb["QNORMTMB_Funct"] = qnorm(rank(tmb$TMB_Funct, ties.method = "random")/(length(tmb$TMB_Funct)+1))
  
  data_CNA <- read.delim(paste0("~/Desktop/Cancer/Burden/TCGA/", c, "/data_CNA.txt"), check.names = F)
  profile_CNA = data_CNA[which(data_CNA$Hugo_Symbol %in% profile_genes),]
  data_CNA=data_CNA[,which(colnames(data_CNA) %in% barcode$Tumor_Sample_Barcode)]
  data_CNA = t(data_CNA)
  profile_CNA = profile_CNA[,which(colnames(profile_CNA) %in% barcode$Tumor_Sample_Barcode)]
  profile_CNA = t(profile_CNA)
  cnb = data.frame(barcode$Tumor_Sample_Barcode)
  colnames(cnb) = "Tumor_Sample_Barcode"
  cnb["AllCNB"] = NA
  cnb["QNORMAllCNB"] = NA
  cnb["CNB"] = NA
  cnb["QNORMCNB"] = NA
  cnb["AllCNB_Profile"] = NA
  cnb["QNORMAllCNB_Profile"] = NA
  cnb["CNB_Profile"] = NA
  cnb["QNORMCNB_Profile"] = NA
  for(n in cnb$Tumor_Sample_Barcode){
    values=data_CNA[which(rownames(data_CNA)==n),]
    values = abs(values)
    cnb[which(cnb$Tumor_Sample_Barcode==n),"AllCNB"] = length(which(values==1 | values==2))
    cnb[which(cnb$Tumor_Sample_Barcode==n),"CNB"] = length(which(values==2))
    values=profile_CNA[which(rownames(profile_CNA)==n),]
    values = abs(values)
    cnb[which(cnb$Tumor_Sample_Barcode==n),"AllCNB_Profile"] = length(which(values==1 | values==2))
    cnb[which(cnb$Tumor_Sample_Barcode==n),"CNB_Profile"] = length(which(values==2))
  }
  cnb$AllCNB/38
  cnb["QNORMAllCNB"] = qnorm(rank(cnb$AllCNB, ties.method = "random")/(length(cnb$AllCNB)+1))
  cnb$CNB/38
  cnb["QNORMCNB"] = qnorm(rank(cnb$CNB, ties.method = "random")/(length(cnb$CNB)+1))
  cnb$AllCNB_Profile/3
  cnb["QNORMAllCNB_Profile"] = qnorm(rank(cnb$AllCNB_Profile, ties.method = "random")/(length(cnb$AllCNB_Profile)+1))
  cnb$CNB_Profile/38
  cnb["QNORMCNB_Profile"] = qnorm(rank(cnb$CNB_Profile, ties.method = "random")/(length(cnb$CNB_Profile)+1))
  x = merge(tmb, cnb, by = "Tumor_Sample_Barcode", all = T)
  burden = rbind(burden, x)
}

TCGA = subset(TCGA, select = c(FID, PATIENT_ID, Cancer, Sex, Age, Metastatic, MSI))
TCGA = merge(TCGA, burden, by = "PATIENT_ID", all.x = T)
rm(TCGA_pc2, barcode, clinical_metastatic, mmc5, msi, names, somatic, temp, temp2, 
    tmb, x, c, mean_pc1, mean_pc2, std_pc1, std_pc2, burden, cnb, data_CNA, n, values)

purity <- read.delim("~/Desktop/Cancer/Burden/TCGA/purity.txt")
#purity["AVE"] = rowMeans(purity[,3:6], na.rm=T)
#purity[which(purity$AVE==0),"AVE"] = NA
#purity[which(is.na(purity$CPE)),"CPE"] = purity[which(is.na(purity$CPE)),"AVE"]
x=data.frame(t(as.data.frame(str_split(purity$Sample.ID,"-"))))
y=data.frame(paste(x$X1, x$X2, x$X3, sep = '-'))
colnames(y) = "PATIENT_ID"
y["SAMPLE"] = x$X4
y["Sample.ID"] = data.frame(paste(x$X1, x$X2, x$X3, x$X4, sep = '-'))
y["NUM"] = NA
y["LTR"] = NA
y[which(y$SAMPLE %in% c("01A", "01B", "01C", "01D", "01R")),"NUM"] = "01"
y[which(y$SAMPLE %in% c("02A", "02B")),"NUM"] = "02"
y[which(y$SAMPLE %in% c("05A")),"NUM"] = "05"
y[which(y$SAMPLE %in% c("06A", "06B")),"NUM"] ="06"
y[which(y$SAMPLE %in% c("07A")),"NUM"] = "07"
y[which(y$SAMPLE %in% c("01A", "02A", "05A", "06A", "07A")),"LTR"] = "A"
y[which(y$SAMPLE %in% c("01B", "02B", "06B")),"LTR"] = "B"
y[which(y$SAMPLE %in% c("01C")),"LTR"] = "C"
y[which(y$SAMPLE %in% c("01D")),"LTR"] = "D"
y[which(y$SAMPLE %in% c("01R")),"LTR"] = "R"
y["Tumor_Sample_Barcode"] = paste(y$PATIENT_ID, y$NUM, sep = "-")

purity = merge(y, purity, by = "Sample.ID")
purity = purity[-is.na(purity$CPE),]
final_purity = NULL
for(i in unique(purity$Tumor_Sample_Barcode)){
  temp = purity[which(purity$Tumor_Sample_Barcode==i),]
  mean_purity = mean(temp$CPE, na.rm = T)
  final_purity = rbind(final_purity, c(i, unique(temp$PATIENT_ID), unique(temp$NUM), mean_purity))
}
final_purity = data.frame(final_purity)
colnames(final_purity) = c("Tumor_Sample_Barcode", "PATIENT_ID","SampleNumber", "Tumor_Purity")
TCGA = merge(TCGA, final_purity, by = "PATIENT_ID", all.x = T)
fix = TCGA$PATIENT_ID[duplicated(TCGA$PATIENT_ID)]
figure_out = TCGA[which(TCGA$PATIENT_ID %in% fix),]
TCGA = TCGA[-which(TCGA$PATIENT_ID %in% fix),]

keep = figure_out[which(figure_out$Tumor_Sample_Barcode.x==figure_out$Tumor_Sample_Barcode.y),]
TCGA = rbind(TCGA, keep)
figure_out = figure_out[-which(figure_out$PATIENT_ID %in% keep$PATIENT_ID),]
figure_out = figure_out[which(figure_out$SampleNumber=="01"),]
figure_out$Tumor_Sample_Barcode.y = figure_out$Tumor_Sample_Barcode.x
TCGA = rbind(TCGA, figure_out)
esoph_purity <- read.delim("~/Desktop/Cancer/Burden/TCGA/esoph_purity.txt")
esoph_purity = subset(esoph_purity, select = c(barcode, Absolute.extract.purity))
colnames(esoph_purity) = c("PATIENT_ID", "Tumor_Purity")
for(line in 1:dim(esoph_purity)[1]){
  temp=esoph_purity[line,]
  TCGA[which(TCGA$PATIENT_ID==temp$PATIENT_ID),"Tumor_Purity"] = temp$Tumor_Purity
}

TCGA = subset(TCGA, select = c(Cancer, FID, PATIENT_ID, Sex, Age, Tumor_Purity, Metastatic, MSI,
                               TMB, CNB, AllCNB, TMB_Profile, TMB_Funct, CNB_Profile, AllCNB_Profile, 
                               QNORMTMB, QNORMCNB, QNORMAllCNB, QNORMTMB_Profile, QNORMTMB_Funct, QNORMCNB_Profile, QNORMAllCNB_Profile))
TCGA[which(TCGA$Cancer %in% c("Breast_Carcinoma", "Ovarian_Cancer", "Endometrial_Cancer") & TCGA$Sex=="Male"),"Sex"] = NA
TCGA[which(TCGA$Cancer=="Prostate_Cancer" & TCGA$Sex=="Female"), "Sex"] = NA
TCGA = TCGA[complete.cases(TCGA),]
TCGA[which(TCGA$Sex=="Male"),"Sex"] = 0
TCGA[which(TCGA$Sex=="Female"),"Sex"] = 1
#TCGA$Age = floor(TCGA$Age)
TCGA = TCGA[which(TCGA$Age>=18),]
TCGA_PCS <- read.delim("~/Desktop/Cancer/Burden/TCGA/top20.eigenvec", header=TRUE, comment.char="")
#TCGA_PCS = subset(TCGA_PCS, select = c(IID, PC1, PC2, PC3, PC4, PC5))
TCGA_PCS[,-1]
colnames(TCGA_PCS)[1] = "FID"
TCGA = merge(TCGA, TCGA_PCS, by = "FID")
TCGA$PC1 = scale(TCGA$PC1)
TCGA$PC2 = scale(TCGA$PC2)
TCGA$PC3 = scale(TCGA$PC3)
TCGA$PC4 = scale(TCGA$PC4)
TCGA$PC5 = scale(TCGA$PC5)
TCGA$PC6 = scale(TCGA$PC6)
TCGA$PC7 = scale(TCGA$PC7)
TCGA$PC8 = scale(TCGA$PC8)
TCGA$PC9 = scale(TCGA$PC9)
TCGA$PC10 = scale(TCGA$PC10)
TCGA$PC11 = scale(TCGA$PC11)
TCGA$PC12 = scale(TCGA$PC12)
TCGA$PC13 = scale(TCGA$PC13)
TCGA$PC14 = scale(TCGA$PC14)
TCGA$PC15 = scale(TCGA$PC15)
TCGA$PC16 = scale(TCGA$PC16)
TCGA$PC17 = scale(TCGA$PC17)
TCGA$PC18 = scale(TCGA$PC18)
TCGA$PC19 = scale(TCGA$PC19)
TCGA$PC20 = scale(TCGA$PC20)

TCGA_EUR <- read.delim("~/Desktop/Cancer/Burden/TCGA/ANCESTRY/TCGA.EA.PC1.sscore")
TCGA_EUR = subset(TCGA_EUR, select = c("IID", "SCORE1_AVG"))
colnames(TCGA_EUR) = c("FID", "ProjectedPC1")
TCGA_EUR_PC2 <- read.delim("~/Desktop/Cancer/Burden/TCGA/ANCESTRY/TCGA.EA.PC2.sscore")
TCGA_EUR_PC2 = subset(TCGA_EUR_PC2, select = c("IID", "SCORE1_AVG"))
colnames(TCGA_EUR_PC2) = c("FID", "ProjectedPC2")
TCGA_EUR = merge(TCGA_EUR, TCGA_EUR_PC2, by = "FID")
TCGA = merge(TCGA, TCGA_EUR, by = "FID")
for(c in unique(TCGA$Cancer)){
  temp = TCGA[which(TCGA$Cancer==c),]
  ordered = read.table(paste0("~/Desktop/Cancer/Burden/TCGA/", c, "_snps_order.txt"), head = F)
  colnames(ordered) = "FID"
  temp = merge(temp, ordered, by = "FID", all.y = T)
  temp = temp[match(ordered$FID, temp$FID),]
  if(c=="Melanoma"){
    temp = subset(temp, select = -c(Cancer, ProjectedPC1, ProjectedPC2))
  } else{
    #else if(c %in% c("Breast_Carcinoma", "Ovarian_Cancer", "Endometrial_Cancer", "Prostate_Cancer")){
    #temp = subset(temp, select = -c(Cancer, Metastatic, Sex, ProjectedPC1, ProjectedPC2))
  #}else{
    temp = subset(temp, select = -c(Cancer, Metastatic, ProjectedPC1, ProjectedPC2))
  }
  covars = subset(temp, select = -c(PATIENT_ID, TMB, CNB, AllCNB, TMB_Profile, CNB_Profile, AllCNB_Profile, QNORMTMB, QNORMCNB, QNORMAllCNB, QNORMTMB_Profile, QNORMCNB_Profile, QNORMAllCNB_Profile))
  phenos = subset(temp, select = c(FID, TMB, CNB, AllCNB, TMB_Profile, CNB_Profile, AllCNB_Profile, QNORMTMB, QNORMCNB, QNORMAllCNB, QNORMTMB_Profile, QNORMCNB_Profile, QNORMAllCNB_Profile))
  write.table(temp, paste0("~/Desktop/Cancer/Burden/TCGA/analyses/", c, "/",c, "_data.txt"), col.names = T, row.names = F, sep='\t', quote = F)
  write.table(phenos, paste0("~/Desktop/Cancer/Burden/TCGA/analyses/", c, "/",c, "_phenos.txt"), col.names = T, row.names = F, sep='\t', quote = F)
  write.table(covars, paste0("~/Desktop/Cancer/Burden/TCGA/analyses/", c, "/",c, "_covars.txt"), col.names = T, row.names = F, sep='\t', quote = F)
  covars$FID=1 #replace FID with intercept
  write.table(phenos, paste0("~/Desktop/Cancer/Burden/TCGA/analyses/", c, "/",c, "_phenos_gemma.txt"), col.names = F, row.names = F, sep='\t', quote = F)
  write.table(covars, paste0("~/Desktop/Cancer/Burden/TCGA/analyses/", c, "/",c, "_covars_gemma.txt"), col.names = F, row.names = F, sep='\t', quote = F)
}
