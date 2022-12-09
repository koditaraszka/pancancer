MetaCNB_GWAS <- read.table("~/Desktop/Cancer/Burden/Profile/MetaCNB_GWAS.txt", quote="\"", comment.char="", stringsAsFactors = F)
MetaCNB_GWAS = MetaCNB_GWAS[which(MetaCNB_GWAS$V1=="rs4636983"),]
MetaCNB_GWAS = t(MetaCNB_GWAS[,2:ncol(MetaCNB_GWAS)])
MetaDeepCNB_GWAS <- read.table("~/Desktop/Cancer/Burden/Profile/MetaDeepCNB_GWAS.txt", quote="\"", comment.char="", stringsAsFactors = F)
MetaDeepCNB_GWAS = MetaDeepCNB_GWAS[which(MetaDeepCNB_GWAS$V1=="rs4636983"),]
MetaDeepCNB_GWAS = t(MetaDeepCNB_GWAS[,2:ncol(MetaDeepCNB_GWAS)])
name = c("Bladder_Cancer", "Breast_Carcinoma","Cancer_of_Unknown_Primary", 
         "Colorectal_Cancer","Endometrial_Cancer","Esophagogastric_Carcinoma",
         "Glioma","Head_and_Neck_Carcinoma", "Leukemia", "Melanoma",
         "Non-Hodgkin_Lymphoma","Non-Small_Cell_Lung_Cancer", "Ovarian_Cancer",
         "Pancreatic_Cancer", "Prostate_Cancer", "Renal_Cell_Carcinoma", "Soft_Tissue_Sarcoma")

row_names = NULL
values = NULL
for(i in name){
  row_names = c(row_names, i, i)
  values = c(values, "Beta", "Stderr")
}
CNB = data.frame(row_names)
colnames(CNB) = "Cancer"
CNB["Value"] = values
CNB["CNB"] = c(MetaDeepCNB_GWAS)
CNB["CNB_AllCNVs"] = c(MetaCNB_GWAS)

rs4636983 = CNB[which(CNB$Value=="Beta"),c(1,3,4)]
colnames(rs4636983) = c("Cancer","CNB_Beta", "CNB_AllCNVs_Beta")
temp = CNB[which(CNB$Value=="Stderr"),c(1,3,4)]
colnames(temp) = c("Cancer","CNB_Stderr", "CNB_AllCNVs_Stderr")
rs4636983 = merge(rs4636983, temp, by = "Cancer")
rs4636983$CNB_Beta = as.numeric(as.character(rs4636983$CNB_Beta))
rs4636983$CNB_AllCNVs_Beta = as.numeric(as.character(rs4636983$CNB_AllCNVs_Beta))
rs4636983$CNB_Stderr = as.numeric(as.character(rs4636983$CNB_Stderr))
rs4636983$CNB_AllCNVs_Stderr = as.numeric(as.character(rs4636983$CNB_AllCNVs_Stderr))
rs4636983["CNB_Pvalue"] = 2*pnorm(abs(rs4636983$CNB_Beta/rs4636983$CNB_Stderr), lower.tail = F)
rs4636983["CNB_AllCNVs_Pvalue"] = 2*pnorm(abs(rs4636983$CNB_AllCNVs_Beta/rs4636983$CNB_AllCNVs_Stderr), lower.tail = F)
rs4636983$Cancer = as.character(rs4636983$Cancer)
z_cnb = rma(yi=rs4636983$CNB_Beta, sei=rs4636983$CNB_Stderr, method = 'FE')
z_cnb_allcnvs = rma(yi=rs4636983$CNB_AllCNVs_Beta, sei=rs4636983$CNB_AllCNVs_Stderr, method = 'FE')
rs4636983 = rbind(rs4636983, c("Pan-Cancer", z_cnb$beta, z_cnb_allcnvs$beta, z_cnb$se, z_cnb_allcnvs$se, z_cnb$pval, z_cnb_allcnvs$pval))
rs4636983$CNB_Pvalue = as.numeric(rs4636983$CNB_Pvalue)
rs4636983$CNB_AllCNVs_Pvalue = as.numeric(rs4636983$CNB_AllCNVs_Pvalue)
rs4636983$Cancer = as.character(rs4636983$Cancer)
name = c("Bladder_Cancer", "Breast_Carcinoma","Cancer_of_Unknown_Primary", 
         "Colorectal_Cancer","Endometrial_Cancer","Esophagogastric_Carcinoma",
         "Glioma","Head_and_Neck_Carcinoma", #"Leukemia", 
         "Melanoma",
         #"Non-Hodgkin_Lymphoma",
         "Non-Small_Cell_Lung_Cancer", "Ovarian_Cancer",
         "Pancreatic_Cancer", "Prostate_Cancer", #"Renal_Cell_Carcinoma", 
         "Soft_Tissue_Sarcoma")
col_names = "rsid"
for(i in name){
  col_names = c(col_names, paste0(i,"_Beta"), paste0(i, "_Stderr"))
}
col_names = c(col_names, "orig")
significant_gwas <- read.table("~/Desktop/Cancer/Burden/Tempus/significant_gwas.txt", quote="\"", comment.char="")
colnames(significant_gwas) = col_names
significant_gwas = data.frame(t(significant_gwas))
colnames(significant_gwas) = c("Tumor", "Normal")
significant_gwas = significant_gwas[-1,]
significant_gwas = significant_gwas[-dim(significant_gwas)[1],]
temp = significant_gwas[grep("Stderr", rownames(significant_gwas)),]
temp["Cancer"] = name
colnames(temp) = c("Tumor_Stderr", "Normal_Stderr", "Cancer")
rownames(temp) = 1:dim(temp)[1]
significant_gwas = significant_gwas[grep("Beta", rownames(significant_gwas)),]
significant_gwas["Cancer"] = name
colnames(significant_gwas) = c("Tumor_Beta", "Normal_Beta", "Cancer")
rownames(significant_gwas) = 1:dim(significant_gwas)[1]
significant_gwas = merge(significant_gwas, temp, by = "Cancer", all = T)
significant_gwas$Tumor_Beta = as.numeric(as.character(significant_gwas$Tumor_Beta))
significant_gwas$Normal_Beta = as.numeric(as.character(significant_gwas$Normal_Beta))
significant_gwas$Tumor_Stderr = as.numeric(as.character(significant_gwas$Tumor_Stderr))
significant_gwas$Normal_Stderr = as.numeric(as.character(significant_gwas$Normal_Stderr))
significant_gwas["Tumor_Pvalue"] = 2*pnorm(abs((significant_gwas$Tumor_Beta/significant_gwas$Tumor_Stderr)), lower.tail=F)
significant_gwas["Normal_Pvalue"] = 2*pnorm(abs((significant_gwas$Normal_Beta/significant_gwas$Normal_Stderr)), lower.tail=F)
z_normal = rma(yi=significant_gwas$Normal_Beta, sei=significant_gwas$Normal_Stderr, method = 'FE')
z_tumor = rma(yi=significant_gwas$Tumor_Beta, sei=significant_gwas$Tumor_Stderr, method = 'FE')
significant_gwas = rbind(significant_gwas, c("Pan-Cancer", z_tumor$beta, z_normal$beta, z_tumor$se, z_normal$se, z_tumor$pval, z_normal$pval))
significant_gwas$Normal_Pvalue = as.numeric(significant_gwas$Normal_Pvalue)
significant_gwas$Tumor_Pvalue = as.numeric(significant_gwas$Tumor_Pvalue)
significant_gwas$Cancer = as.character(significant_gwas$Cancer)
significant_gwas = merge(rs4636983, significant_gwas, by = "Cancer", all=T)
rs4636983 = data.frame(c(rep(significant_gwas$Cancer,4)))
colnames(rs4636983) = "Cancer"
rs4636983["Beta"] = c(significant_gwas$CNB_Beta, significant_gwas$CNB_AllCNVs_Beta, 
                      significant_gwas$Tumor_Beta, significant_gwas$Normal_Beta)
rs4636983["Stderr"] = c(significant_gwas$CNB_Stderr, significant_gwas$CNB_AllCNVs_Stderr, 
                        significant_gwas$Tumor_Stderr, significant_gwas$Normal_Stderr)
rs4636983["Pvalue"] = c(significant_gwas$CNB_Pvalue, significant_gwas$CNB_AllCNVs_Pvalue, 
                        significant_gwas$Tumor_Pvalue, significant_gwas$Normal_Pvalue)
rs4636983["Source"] = c(rep("Profile: CNB", dim(significant_gwas)[1]), rep("Profile: All CNB", dim(significant_gwas)[1]),
                        rep("Tempus: CNB", dim(significant_gwas)[1]), rep("Tempus: CNB", dim(significant_gwas)[1]))
rs4636983["Sample"] = c(rep("Tumor", 3*dim(significant_gwas)[1]), rep("Normal", dim(significant_gwas)[1]))

rs4636983$Beta = as.numeric(rs4636983$Beta)
rs4636983$Stderr = as.numeric(rs4636983$Stderr)
rs4636983["CI_Upper"] = rs4636983$Beta + rs4636983$Stderr * qnorm(0.975)
rs4636983["CI_Lower"] = rs4636983$Beta - rs4636983$Stderr * qnorm(0.975)


rs4636983["Significant"] = ""
rs4636983[which(rs4636983$Pvalue<0.05),"Significant"] = "*"
rs4636983[which(rs4636983$Pvalue<0.05/32),"Significant"] = "***"
rs4636983[which(rs4636983$Pvalue<5e-8 & rs4636983$Cancer=="Pan-Cancer"),"Significant"] = "" #4.2E-13"
rs4636983[which(rs4636983$Pvalue<0.05/32 & rs4636983$Cancer=="Pan-Cancer" 
                & rs4636983$Source=="Profile: CNB"),"Significant"] = "" #7.8E-5"
rs4636983[which(rs4636983$Pvalue<0.05 & rs4636983$Cancer=="Pan-Cancer" 
                & rs4636983$Source=="Tempus: CNB"),"Significant"] = "" #3.7E-3"
rs4636983$Cancer = str_replace_all(rs4636983$Cancer, "_", " ")
rs4636983[which(rs4636983$Cancer=="Soft\nTissue\nSarcoma"),"Cancer"] = "Soft Tissue\nSarcoma"
rs4636983[which(rs4636983$Cancer=="Cancer\nof\nUnknown\nPrimary"),"Cancer"] = "Cancer of\nUnknown Primary"
rs4636983[which(rs4636983$Cancer=="Non-Small\nCell\nLung\nCancer"),"Cancer"] = "Non-Small Cell\nLung Cancer"
rs4636983[which(rs4636983$Cancer=="Head\nand\nNeck\nCarcinoma"),"Cancer"] = "Head and Neck\nCarcinoma"
rs4636983[which(rs4636983$Cancer=="Renal\nCell\nCarcinoma"),"Cancer"] = "Renal Cell\nCarcinoma"
rs4636983$Cancer = factor(rs4636983$Cancer,levels = ordered_cancers)
#rs4636983$Source[which(rs4636983$Source=="Profile: All CNB")] = "Profile: All CNB"
rs4636983$Source = factor(rs4636983$Source, levels = c("Profile: All CNB", "Profile: CNB", "Tempus: CNB"))
rs4636983_fig = ggplot(rs4636983, aes(x=Beta, y=Cancer, xmin = CI_Lower, xmax = CI_Upper, color = Sample, 
                                      label = Significant)) + all_theme + geom_point(position = position_dodge(0.5)) + 
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
  geom_errorbarh(position = position_dodge(0.5), height=0.0) + facet_grid(.~Source, scales= "free") + 
  ggtitle("rs4636983") + guides(colour = guide_legend(override.aes = list(size=1))) +
  theme(axis.title.x = element_blank(), panel.spacing = unit(0.5, "lines")) + 
  scale_color_manual(values = c("red","blue")) + scale_y_discrete(limits=rev) + 
  geom_text(position = position_dodge(0.5), size = 10, hjust = 0)


MetaCNB_GWAS_results <- read.delim("~/Desktop/Cancer/Burden/Profile/MetaCNB_GWAS_results.txt", header=FALSE)


colnames(MetaCNB_GWAS_results) = c("RSID", "NUM_STUDY", "PVALUE_FE", "BETA_FE", "STD_FE","PVALUE_RE", "BETA_RE", "STD_RE","PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE")
MetaCNB_GWAS_results = subset(MetaCNB_GWAS_results, select = c("RSID", "NUM_STUDY", "PVALUE_FE", "BETA_FE", "STD_FE","PVALUE_RE", "BETA_RE", "STD_RE","PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE"))

MetaCNB_GWAS_results["I2_Group"] = "0-25%"
MetaCNB_GWAS_results[which(MetaCNB_GWAS_results$I_SQUARE>=25),"I2_Group"] = "25-50%"
MetaCNB_GWAS_results[which(MetaCNB_GWAS_results$I_SQUARE>=50),"I2_Group"] = "50-75%"
MetaCNB_GWAS_results[which(MetaCNB_GWAS_results$I_SQUARE>=75),"I2_Group"] = "75-100%"

gwas_details = read.table("~/Desktop/Cancer/Burden/Profile/gwas_details.txt", head = T, stringsAsFactors = F)
gwas_details["CHR"] = as.data.frame(str_split_fixed(gwas_details$LOC,":",2))[[1]]
gwas_details["BP"] = as.data.frame(str_split_fixed(gwas_details$LOC,":",2))[[2]]
gwas_details$CHR = as.numeric(as.character(gwas_details$CHR))
gwas_details$BP = as.numeric(as.character(gwas_details$BP))
MetaCNB_GWAS_results = merge(gwas_details, MetaCNB_GWAS_results, by = "RSID", all.y = T)

man_fig = subset(MetaCNB_GWAS_results, select = c("RSID", "CHR", "BP", "BETA_FE", "STD_FE", "PVALUE_FE"))
colnames(man_fig) = c("RSID", "chr", "pos", "BETA", "STDERR", "pvalue")
man_fig = man_fig[-which(is.na(man_fig$pvalue)),]
man_fig$pos = as.numeric(as.character(man_fig$pos))
man_fig["color"] = NA
man_fig[which(man_fig$RSID=="rs4636983"), "color"] = "blue"

rm(gwas_details, MetaCNB_GWAS_results)

manplot = fast_manhattan(man_fig, speed = "slow", pointsize = 3) +
  geom_hline(yintercept = -log10(5e-08), color ="red") + all_theme + ylab("")

#(rs4636983_fig/manplot) + plot_annotation(tag_levels = 'A')
#ggsave("~/Desktop/Cancer/Burden/Paper/Figures/cnb_rs4636983_gwas.png", 
#    grid.arrange(rs4636983_fig),
#    width = 6, height = 4, units = "in", scale = 3)

ggsave("~/Desktop/Cancer/Burden/Paper/Figures/main_panCancerAncestry_Dec1.png", 
       main_anc_plot, width = 5, height = 2.5, units = "in", scale = 4)
ggsave("~/Desktop/Cancer/Burden/Paper/Figures/supp_panCancerAncestry_Dec1.png", 
       supp_anc_plot, width = 6, height = 3, units = "in", scale = 4)
ggsave("~/Desktop/Cancer/Burden/Paper/Figures/main_panCancerSexAgeMet_Dec1.png", 
       main_covar_plot, width = 5, height = 2.5, units = "in", scale = 4)
ggsave("~/Desktop/Cancer/Burden/Paper/Figures/main_panCancerSexAgeMet_table.png", 
       grid.arrange(main_covar_table), width = 6, height = 3, units = "in", scale = 4)
ggsave("~/Desktop/Cancer/Burden/Paper/Figures/supp_panCancerSexAgeMet_table.png", 
       grid.arrange(cnball_covar_table), width = 6, height = 3, units = "in", scale = 4)
ggsave("~/Desktop/Cancer/Burden/Paper/Figures/supp_panCancerSexAgeMet_Dec1.png", 
       supp_covar_plot, width = 6, height = 3, units = "in", scale = 4)



bothcohorts$Cancer = factor(bothcohorts$Cancer,levels = c("Bladder\nCancer", 
                                                          "Breast\nCarcinoma", "Cancer of\nUnknown Primary",
                                                          "Colorectal\nCancer", "Endometrial\nCancer",
                                                          "Esophagogastric\nCarcinoma", "Head and Neck\nCarcinoma",
                                                          "Glioma", "Leukemia", "Melanoma", "Non-Hodgkin\nLymphoma",
                                                          "Non-Small Cell\nLung Cancer", "Ovarian\nCancer",
                                                          "Pancreatic\nCancer", "Prostate\nCancer", "Renal Cell\nCarcinoma",
                                                          "Soft Tissue\nSarcoma"))      



tmb_all = ggplot(bothcohorts, aes(x=Cancer, y = TMB, fill = Source)) + geom_violin() + all_theme + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank())

tmb_zoom = ggplot(bothcohorts[which(bothcohorts$TMB<25),], aes(x=Cancer, y = TMB, fill = Source)) + geom_violin() + all_theme + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank())

cnb_all = ggplot(bothcohorts, aes(x=Cancer, y = CNB, fill = Source)) + geom_violin() + all_theme + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank())

cnb_zoom = ggplot(bothcohorts[which(bothcohorts$CNB<5),], aes(x=Cancer, y = CNB, fill = Source)) + geom_violin() + all_theme + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank())

cnball_all = ggplot(bothcohorts, aes(x=Cancer, y = CNB_All, fill = Source)) + geom_violin() + all_theme + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank()) + guides(fill = 'none')

cnball_zoom = ggplot(bothcohorts[which(bothcohorts$CNB_All<20),], aes(x=Cancer, y = CNB_All, fill = Source)) + geom_violin() + all_theme + theme(axis.title.y = element_blank()) + guides(fill = 'none')

alldist = wrap_plots((tmb_all /tmb_zoom / cnb_all / cnb_zoom / cnball_all / cnball_zoom), guides='collect', tag_level = 'new') & theme(legend.position='bottom')

zoomdist = wrap_plots((tmb_zoom / cnb_zoom / cnball_zoom + guides(fill='none')), guides='collect', tag_level = 'new') & theme(legend.position='bottom')
