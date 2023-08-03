# With subpopulation fitted as fixed
rm(list = ls())
library("BGLR")
library(qqman)
# set the directory
setwd("F:/My_Projects/workshops/My_workshops/GWAS/Motamed/3.GWAS_In_R/4.BayesC")
#load data
load("SNP_pheno_rice.RData")
# prepare the data for BGLR
SNP=t(SNP)
pheno=pheno[rownames(SNP),]
name=colnames(pheno)[4:13]
# Specify the regression function (or linear predictor)
set.seed(123456)
selected=sample(1:ncol(SNP),round(ncol(SNP)/2,0))
SNP=SNP[,-c(selected)]
ETA<-list(list(~factor(Subpop),data=pheno,model='FIXED'),
          list(X=SNP,model="BayesC"))
# Running BGLR function (BayesC)
start_time <- Sys.time()
fm=BGLR(y=as.numeric(pheno[,name[1]]), response_type = "gaussian",ETA = ETA,
        nIter = 25000, burnIn = 5000, thin = 5, saveAt = "",
        verbose = T, rmExistingFiles = TRUE, groups=NULL)
end_time <- Sys.time()

end_time - start_time
#
load("map.RData")
P_val=as.data.frame(matrix(NA,nrow = ncol(SNP),ncol = 4))
colnames(P_val)=c("SNP", "CHR", "BP", "PIP")
P_val$SNP=colnames(SNP)
rownames(P_val)=P_val$SNP
map=map[rownames(P_val),]
P_val$CHR=map$Chr
P_val$BP=map$Pos_bp
P_val$PIP=fm$ETA[[2]]$d
png(paste0("BayesC_PIP_",name[1],".png"), height = 500,width = 1000)
manhattan(P_val, chr="CHR", bp="BP", snp="SNP", p="PIP",suggestiveline=F,genomewideline= F,ylim=c((min(P_val$PIP)+max(P_val$PIP))/2,max(P_val$PIP)*1.05),logp = F )
dev.off()

