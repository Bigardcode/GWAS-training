# Genotypic and phenotypic data of rice accessions can be downloaded from the RDP website:
# http://www.ricediversity.org/data/
rm(list = ls())
library(SNPRelate)
library(gdsfmt)
library(scales)
library("rrBLUP")
library(qqman)
# Set  the directory
setwd("F:/My_Projects/workshops/My_workshops/GWAS/Motamed/3.GWAS_In_R/3.MLM_GWAS")
# load the data file for rice
load("SNP_pheno_rice.RData")
# load the map file
load("map.RData")
#
head(map)
SNP[1:6,1:6]
# Making GDS file for SNPRelate package
data(hapmap_geno)
SNPgds=hapmap_geno
rm(hapmap_geno)
SNPgds$sample.id=as.character(colnames(SNP))
SNPgds$snp.id=as.character(rownames(SNP))
SNPgds$snp.position=map$Pos_bp
SNPgds$snp.chromosome=map$Chr
SNPgds$snp.allele=rep("A/G",nrow(map))
SNPgds$genotype=SNP
snpgdsCreateGeno("SNPgds.gds", genmat = SNPgds$genotype,
                 sample.id = SNPgds$sample.id, snp.id = SNPgds$snp.id,
                 snp.chromosome = SNPgds$snp.chromosome,
                 snp.position = SNPgds$snp.position,
                 snp.allele = SNPgds$snp.allele, snpfirstdim=TRUE)
snpgdsSummary("SNPgds.gds")
genofile <- snpgdsOpen("SNPgds.gds")
# PCA analyses
pca <- snpgdsPCA(genofile,num.thread=4)
# create a function that make first letter to upper case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
# plotting the PC results
# PC1 and PC2
temp=as.data.frame(matrix(NA,nrow=nrow(pheno),ncol=3))
temp[,1]=pca$eigenvect[,1]
temp[,2]=pca$eigenvect[,2]
temp[,3]=firstup(pheno$Subpop)
temp[,3]=as.factor(temp[,3])
colnames(temp)=c("PCA1","PCA2","Subpop")
attach(temp)
jpeg('Scatter_polt_of_PC1_PC2.jpg',width = 600,height=600)
par(mar=c(6,7,4,2))
plot(PCA1,PCA2,xlab= "",ylab = "",
     col=alpha(c("yellow2","red","gray19","seagreen1","dodgerblue4","steelblue1")[Subpop],0.7),pch=16,cex.lab=2,cex=2,
     main ="",
     cex.main=1.5,cex.axis=2,font.axis=1)
detach(temp)
par(mar=c(5,5.5,4,2),cex.lab = 2.5)
title(ylab =paste0("PC2 (",round(pca$varprop[2]*100,2),"%)"),xlab= paste0("PC1 (",round(pca$varprop[1]*100,2),"%)"))
legend("topleft",legend = levels(temp$Subpop),
       col=c("yellow2","red","gray19","seagreen1","dodgerblue4","steelblue1"),pch=16)
dev.off()
# PC3 and PC4
temp=as.data.frame(matrix(NA,nrow=nrow(pheno),ncol=3))
temp[,1]=pca$eigenvect[,3]
temp[,2]=pca$eigenvect[,4]
temp[,3]=firstup(pheno$Subpop)
temp[,3]=as.factor(temp[,3])
colnames(temp)=c("PCA3","PCA4","Subpop")
attach(temp)
jpeg('Scatter_polt_of_PC3_PC4.jpg',width = 600,height=600)
par(mar=c(6,7,4,2))
plot(PCA3,PCA4,xlab= "",ylab = "",col=alpha(c("yellow2","red","gray19","seagreen1","dodgerblue4","steelblue1")[Subpop],0.7),pch=16,cex.lab=2,cex=2,
     main ="",
     cex.main=1.5,cex.axis=2,font.axis=1)
detach(temp)
par(mar=c(5,5.5,4,2),cex.lab = 2.5)
title(ylab =paste0("PC3 (",round(pca$varprop[3]*100,2),"%)"),xlab= paste0("PC4 (",round(pca$varprop[4]*100,2),"%)"))
legend("topleft",legend = levels(temp$Subpop),
       col=c("yellow2","red","gray19","seagreen1","dodgerblue4","steelblue1"),pch=16)
dev.off()

# PC5 and PC6
temp=as.data.frame(matrix(NA,nrow=nrow(pheno),ncol=3))
temp[,1]=pca$eigenvect[,5]
temp[,2]=pca$eigenvect[,6]
temp[,3]=firstup(pheno$Subpop)
temp[,3]=as.factor(temp[,3])
colnames(temp)=c("PCA5","PCA6","Subpop")
attach(temp)
jpeg('Scatter_polt_of_PC5_PC6.jpg',width = 600,height=600)
par(mar=c(6,7,4,2))
plot(PCA5,PCA6,xlab= "",ylab = "",col=alpha(c("yellow2","red","gray19","seagreen1","dodgerblue4","steelblue1")[Subpop],0.7),pch=16,cex.lab=2,cex=2,
     main ="",
     cex.main=1.5,cex.axis=2,font.axis=1)
detach(temp)
par(mar=c(5,5.5,4,2),cex.lab = 2.5)
title(ylab =paste0("PC6 (",round(pca$varprop[6]*100,2),"%)"),xlab= paste0("PC5 (",round(pca$varprop[5]*100,2),"%)"))
legend("topright",legend = levels(temp$Subpop),
       col=c("yellow2","red","gray19","seagreen1","dodgerblue4","steelblue1"),pch=16)
dev.off()

# PC4 and PC6
temp=as.data.frame(matrix(NA,nrow=nrow(pheno),ncol=3))
temp[,1]=pca$eigenvect[,4]
temp[,2]=pca$eigenvect[,6]
temp[,3]=firstup(pheno$Subpop)
temp[,3]=as.factor(temp[,3])
colnames(temp)=c("PCA4","PCA6","Subpop")
attach(temp)
jpeg('Scatter_polt_of_PC4_PC6.jpg',width = 600,height=600)
par(mar=c(6,7,4,2))
plot(PCA4,PCA6,xlab= "",ylab = "",col=alpha(c("yellow2","red","gray19","seagreen1","dodgerblue4","steelblue1")[Subpop],0.7),pch=16,cex.lab=2,cex=2,
     main ="",
     cex.main=1.5,cex.axis=2,font.axis=1)
detach(temp)
par(mar=c(5,5.5,4,2),cex.lab = 2.5)
title(ylab =paste0("PC6 (",round(pca$varprop[6]*100,2),"%)"),xlab= paste0("PC4 (",round(pca$varprop[4]*100,2),"%)"))
legend("topleft",legend = levels(temp$Subpop),
       col=c("yellow2","red","gray19","seagreen1","dodgerblue4","steelblue1"),pch=16)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
# GWAS using rrBLUP (MLM GWAS)
#----------------------------------------------------------------------------------------------------------------------
# kinship matrix
SNP=t(SNP)
X <- scale(SNP, center = TRUE, scale = TRUE)
G_snp <- tcrossprod(X)/ncol(X)
# Geno file prepration
SNP=t(SNP)
SNP[1:6,1:6]
SNP=SNP-1
SNP[1:6,1:6]
SNP=as.data.frame(SNP)
SNP$name=rownames(SNP)
SNP$chr=map$Chr
SNP$pos=map$Pos_bp
SNP=SNP[,c((ncol(SNP)-2):ncol(SNP),1:(ncol(SNP)-3))]
SNP[1:6,1:6]
SNP$chr=as.numeric(SNP$chr)
class(SNP$pos)
class(SNP$chr)

# pheno prepration
pheno=as.data.frame(pheno)
pheno=pheno[colnames(SNP)[4:ncol(SNP)],]
pheno$line=rownames(pheno)
# reducing the number of SNPs
set.seed(123456)
selected=sample(1:nrow(SNP),round(nrow(SNP)/2,0))
SNP=SNP[-selected,]
map=map[-selected,]
# select a trait
Trait="Flowering.time.at.Arkansas"
#-------------------------------------------
# Without population stratification
P_val=as.data.frame(matrix(NA,nrow = nrow(SNP),ncol = 4))
colnames(P_val)=c("SNP", "CHR", "BP", "P")
P_val$SNP=SNP$name
P_val$CHR=SNP$chr
P_val$BP=SNP$pos
SNP2=as.matrix(SNP[,4:ncol(SNP)])
for (i in 1:nrow(SNP)) {
  pheno1=pheno[,c("line","Subpop",Trait)]
  pheno1$SNP=SNP2[i,]
  pheno1=pheno1[-which(is.na(pheno1[,Trait])),]
  fit <- lm(pheno1[,3]~SNP+Subpop,data=pheno1)
  P_val[i,"P"]=summary(fit)$coefficients["SNP","Pr(>|t|)"]
  if(i%%100==0) cat(paste0(i,"\n"))
}
# Bonferroni
png(paste0("Without_PS_",Trait,".png"), height = 500,width = 1000)
manhattan(P_val, chr="CHR", bp="BP", snp="SNP", p="P",suggestiveline=-log10(0.05/nrow(P_val)),genomewideline= -log10(0.01/nrow(P_val)) )
dev.off()
#-------------------------------------------
# Without PCs 
P_val=as.data.frame(matrix(NA,nrow = nrow(SNP),ncol = 4))
colnames(P_val)=c("SNP", "CHR", "BP", "P")
P_val$SNP=SNP$name
P_val$CHR=SNP$chr
P_val$BP=SNP$pos
SNP2=as.matrix(SNP[,4:ncol(SNP)])

for (i in 1:nrow(SNP)) {
  pheno1=pheno[,c("line","Subpop",Trait, "PCA1","PCA2","PCA3","PCA4")]
  pheno1$SNP=SNP2[i,]
  pheno1=pheno1[-which(is.na(pheno1[,Trait])),]
  fit <- lm(pheno1[,3]~SNP+Subpop+PCA1+PCA2+PCA3+PCA4,data=pheno1)
  P_val[i,"P"]=summary(fit)$coefficients["SNP","Pr(>|t|)"]
  if(i%%100==0) cat(paste0(i,"\n"))
}
# Bonferroni
png(paste0("With_PCs_",Trait,".png"), height = 500,width = 1000)
manhattan(P_val, chr="CHR", bp="BP", snp="SNP", p="P",suggestiveline=-log10(0.05/nrow(P_val)),genomewideline= -log10(0.01/nrow(P_val)) )
dev.off()
#-------------------------------------------
# With kinship
pheno1=pheno[,c("line","Subpop",Trait)]
# GWAS
scores <- GWAS(pheno=pheno1, geno=SNP, fixed=c("Subpop"), K=G_snp, n.PC=0,
               min.MAF=0.01, n.core=1, P3D=TRUE, plot=T)
dev.print(pdf, paste0("./With_kinship_",Trait,'.pdf'),width=10, height=5)



