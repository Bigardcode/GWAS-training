rm(list = ls())
library(snpStats)
library(doParallel)
library(SNPRelate)
library(dplyr)
library(NAM)
library(qqman)
# set directory
setwd("F:/My_Projects/workshops/My_workshops/GWAS/Motamed/3.GWAS_In_R/2.snpStats")
#-------------------------------------------------------------------------------------------------------
# read data in plink format
sample <- read.plink("data2.bed", "data2.bim", "data2.fam")

sample$genotypes
SNP=sample$genotypes@.Data
col.summary(sample$genotypes)$MAF[1:10]
min(col.summary(sample$genotypes)$MAF)
head(sample$fam)
head(sample$map)

# subset reading
subset <- read.plink("data2.bed", "data2.bim", "data2.fam", select.snps=1:10)
subset$genotypes
#-------------------------------------------------------------------------------------------------------
# Use SNP call rate of 100%, MAF of 0.1 (very stringent)
maf <- 0.1
callRate <- 1
SNPstats <- col.summary(sample$genotypes)
maf_call <- with(SNPstats, MAF > maf & Call.rate == callRate)
sample$genotypes <- sample$genotypes[,maf_call]
sample$map <- sample$map[maf_call,]
SNPstats <- SNPstats[maf_call,]
#
#-------------------------------------------------------------------------------------------------------
# LD and kinship coeff
ld <- .2
kin <- .1
snpgdsBED2GDS(bed.fn = "data2.bed", bim.fn = "data2.bim",
              fam.fn = "data2.fam", out.gdsfn = "myGDS", cvt.chr = "char")
genofile <- snpgdsOpen("myGDS", readonly = F)
#gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
#add.gdsn(genofile, "sample.id", gds.ids, replace = T)
geno.sample.ids <- rownames(sample$genotypes)
# First filter for LD
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld,
                          sample.id = geno.sample.ids,
                          snp.id = colnames(sample$genotypes))
snpset.ibd <- unlist(snpSUB, use.names = F)
# And now filter 
ibd <- snpgdsIBDMoM(genofile, kinship = T,
                    sample.id = geno.sample.ids,
                    snp.id = snpset.ibd,
                    num.thread = 1)
ibdcoef <- snpgdsIBDSelection(ibd)
ibdcoef <- ibdcoef[ibdcoef$kinship >= kin,]

# Filter samples out
related.samples <- NULL
ind=table(c(ibdcoef[,1],ibdcoef[,2]))
related.samples=names(ind[which(ind>1)])
ibdcoef[which(ibdcoef[,1]%in%related.samples),1]=NA
ibdcoef[which(ibdcoef[,2]%in%related.samples),2]=NA
ibdcoef=ibdcoef[complete.cases(ibdcoef),]
related.samples=c(related.samples,ibdcoef[,1])
sample$genotypes <- sample$genotypes[!(rownames(sample$genotypes) %in% related.samples),]
sample$genotypes
sample$fam <- sample$fam[!(rownames(sample$fam) %in% related.samples),]

#-------------------------------------------------------------------------------------------------------
# PCA
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids, snp.id = snpset.ibd, num.thread = 1)
pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],
                    PC2 = pca$eigenvect[,2],
                    stringsAsFactors = F)



pcaCol <- rep(rgb(0,0,0,.3), length(pca$sample.id)) # Set black for chinese
pcaCol[1:45] <- rgb(1,0,0,.3) # red for indian
pcaCol[46:90] <- rgb(0,.7,0,.3) # green for malay

png("PCApopulation.png", width = 500, height = 500)
plot(pctab$PC1, pctab$PC2, xlab = "PC1", ylab = "PC2", col = pcaCol, pch = 16)
abline(h = 0, v = 0, lty = 2, col = "grey")
legend("top", legend = c("Case", "Control"), col = pcaCol[c(1,90)], pch = 16, bty = "n")
dev.off()

#-------------------------------------------------------------------------------------------------------
# GWAS
phenodata <- data.frame("id" = rownames(sample$genotypes),
                        "phenotype" = sample$fam$affected, stringsAsFactors = F)
y=sample$fam$affected
gen=sample$genotypes@.Data
gen[1:6,1:6]
gen <- apply(gen,1 ,as.numeric)
gen=gen-1
gen=t(gen)
colnames(gen)=sample$map$snp.name
test=gwas(y=y,gen=gen)

# plotting the results
P_val=as.data.frame(matrix(NA,nrow = ncol(gen),ncol = 4))
colnames(P_val)=c("SNP", "CHR", "BP", "P")
P_val$SNP=colnames(gen)
P_val$CHR=sample$map$chromosome
P_val$BP=sample$map$position
P_val$P=test$PolyTest$pval
png("GWAS_Manhattan_plot.png", height = 500,width = 1000)
manhattan(P_val, chr="CHR", bp="BP", snp="SNP", p="P",suggestiveline=-log10(0.05/nrow(P_val)),genomewideline= -log10(0.01/nrow(P_val)),logp =F )
dev.off()
