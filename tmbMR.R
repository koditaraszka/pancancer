library(stringr)
prune = function(df){
 x = data.frame(t(data.frame(strsplit(df$loc, ":"))))
 rownames(x) = 1:nrow(x)
 colnames(x) = c("CHR", "BP")
 df["CHR"] = as.numeric(x$CHR)
 df["BP"] = as.numeric(x$BP)
 df$z = abs(df$z)
 keep = df[which(df$z==max(df$z, na.rm=T)),][1,]
 chr = keep$CHR
 bp = keep$BP
 minBP = max(0, bp-500000)
 maxBP = bp+500000
 df = df[which(df$CHR!=chr | df$BP>maxBP | df$BP < minBP),]
 while(dim(df)[1]>0){
   add = df[which(df$z==max(df$z, na.rm=T)),][1,]
   chr = add$CHR
   bp = add$BP
   minBP = max(0, bp-500000)
   maxBP = bp+500000
   df = df[which(df$CHR!=chr | df$BP>maxBP | df$BP < minBP),]
   keep = rbind(keep, add)
 } 
 return(keep)
}

PROF_2018.prune <- read.table("~/Desktop/Cancer/Burden/Profile/PROF_2018.prune.in", quote="\"", comment.char="", stringsAsFactors = F)
colnames(PROF_2018.prune) = "loc"

orig <- read.table("~/Desktop/Cancer/Burden/Profile/GSCAN.CigarettesPerDay.PV7.prs", quote="\"", comment.char="")
colnames(orig) = c("rsid", "loc", "alt", "z", "pvalue")
cpd7 = orig
cpd7 = cpd7[which(cpd7$loc %in% PROF_2018.prune$loc),]
pruned = prune(orig)
if(dim(pruned)[1] > dim(cpd7)[1]){
  print(dim(cpd7)[1])
  print(dim(pruned)[1])
  print("CIG PER DAY PRUNED")
  cpd7 = pruned
}
cpdgwas <- read.delim("~/Desktop/Cancer/Burden/Profile/CigarettesPerDay.txt", stringsAsFactors = F)
cpdgwas["loc"] = paste(cpdgwas$CHROM, cpdgwas$POS, sep = ":")
colnames(cpdgwas)[5] = "alt"
cpdgwas = merge(cpdgwas, cpd7, by = c("loc","alt"))
cpdgwas = subset(cpdgwas, select = c("loc", "alt","BETA"))
colnames(cpdgwas) = c("loc", "alt", "orig_beta")

nsclc_tmb_gwas <- read.delim("~/Desktop/Cancer/Burden/Profile/Non-Small_Cell_Lung_Cancer_OrigMSSProteinTMB.txt",stringsAsFactors = F)
colnames(nsclc_tmb_gwas)[5] = "alt"
nsclc_tmb_gwas = merge(nsclc_tmb_gwas, cpd7, by = c("loc","alt"))
nsclc_tmb_gwas = subset(nsclc_tmb_gwas, select = c("loc", "alt", "beta"))
colnames(nsclc_tmb_gwas) = c("loc", "alt", "tmb_beta")
final_cpd = merge(nsclc_tmb_gwas, cpdgwas, by = c('loc','alt'))
print(summary(lm(final_cpd$tmb_beta ~ final_cpd$orig_beta)))

orig <- read.table("~/Desktop/Cancer/Burden/Profile/UKB_460K.pigment_TANNING.PV7.prs", quote="\"", comment.char="")
colnames(orig) = c("rsid", "loc", "alt", "z", "pvalue")
tan7 = orig
tan7 = tan7[which(tan7$loc %in% PROF_2018.prune$loc),]
pruned = prune(orig)
if(dim(pruned)[1] > dim(tan7)[1]){
  print(dim(tan7)[1])
  print(dim(pruned)[1])
  print("TANNING PRUNED")
  tan7 = pruned
}

tangwas <- read.delim("~/Desktop/Cancer/Burden/Profile/pigment_TANNING.sumstats", stringsAsFactors = F)
tangwas["loc"] = paste(tangwas$CHR, tangwas$POS, sep = ":")
colnames(tangwas)[5] = "alt"
tangwas = merge(tangwas, tan7, by = c('loc', 'alt')) 
tangwas = subset(tangwas, select = c("loc","alt", "Beta"))
colnames(tangwas) = c("loc", "alt", "orig_beta")

melanoma_tmb_gwas <- read.delim("~/Desktop/Cancer/Burden/Profile/Melanoma_OrigMSSProteinTMB.txt",stringsAsFactors = F)
colnames(melanoma_tmb_gwas)[5] = 'alt'
melanoma_tmb_gwas = merge(melanoma_tmb_gwas, tan7, by = c('loc', 'alt')) 
melanoma_tmb_gwas = subset(melanoma_tmb_gwas, select = c("loc","alt", "beta"))
colnames(melanoma_tmb_gwas) = c("loc", "alt", "tmb_beta")
final_tan = merge(melanoma_tmb_gwas, tangwas, by = c('loc', 'alt'))
print(summary(lm(final_tan$tmb_beta ~ final_tan$orig_beta)))


#RSID	#STUDY	PVALUE_FE	BETA_FE	STD_FE	PVALUE_RE
pancancer_tmb_gwas <- read.delim("~/Desktop/Cancer/Burden/Profile/pancancer_meta_originaltmb.txt",stringsAsFactors = F, header = F)
colnames(pancancer_tmb_gwas)[1:6] = c('RSID', 'NUM_STUDY', 'PVALUE_FE', 'BETA_FE', 'STD_FE', 'PVALUE_RE')
pancancer_tmb_gwas = subset(pancancer_tmb_gwas, select = c("RSID", "BETA_FE"))
colnames(pancancer_tmb_gwas) = c('rsid', 'tmb_beta')
pancancer_tmb_gwas = merge(pancancer_tmb_gwas, cpd7, by = "rsid")
final_pancpd = merge(pancancer_tmb_gwas, cpdgwas, by = c('loc','alt'))
print(summary(lm(final_pancpd$tmb_beta ~ final_pancpd$orig_beta)))

