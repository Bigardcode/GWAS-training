rm(list = ls())
library(qqman)

setwd("F:/My_Projects/workshops/My_workshops/GWAS/Motamed/3.GWAS_In_R/1.Our_codes")
load("GWAS_Data.RData")

# Example for SNP 1
ind=phenotype
ind$Trait=ind$Trait-1
ind$SNP=SNP[,1]
fit <- glm(Trait~SNP,data=ind,family=binomial())
summary(fit)

# GWAS
P_val=as.data.frame(matrix(NA,nrow = ncol(SNP),ncol = 4))
colnames(P_val)=c("SNP", "CHR", "BP", "P")
P_val$SNP=colnames(SNP)
P_val$CHR=map$Chr
P_val$BP=map$BPPos
for (i in 1:ncol(SNP)) {
  ind$SNP=SNP[,i]
  fit <- glm(Trait~SNP,data=ind,family=binomial())
  P_val[i,"P"]=summary(fit)$coefficients["SNP","Pr(>|z|)"]
  if(i%%100==0) cat(paste0(i,"\n"))
}
# defult
manhattan(P_val, chr="CHR", bp="BP", snp="SNP", p="P" )
# Bonferroni
png("GWAS_Manhattan_plot.png", height = 500,width = 1000)
manhattan(P_val, chr="CHR", bp="BP", snp="SNP", p="P",suggestiveline=-log10(0.05/nrow(P_val)),genomewideline= -log10(0.01/nrow(P_val)) )
dev.off()
rm(list = ls())
