rm(list = ls())
library(qqman)
library("biomaRt")
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(topGO)
library(Rgraphviz)
setwd("F:/My_Projects/workshops/My_workshops/GWAS/Motamed/4.post-GWAS")
# load the genotype data
load("chr6.RData")

#--------------------------------------------------------------------------------------------------
### GWAS
#--------------------------------------------------------------------------------------------------
P_val=as.data.frame(matrix(NA,nrow = ncol(SNP),ncol = 4))
colnames(P_val)=c("SNP", "CHR", "BP", "P")
P_val$SNP=colnames(SNP)
P_val$CHR=map$Chr
P_val$BP=map$BP
ind=y
for (i in 1:ncol(SNP)) {
  ind$SNP=SNP[,i]
  fit <- glm(trait~SNP,data=ind,family=binomial())
  P_val[i,"P"]=summary(fit)$coefficients["SNP","Pr(>|z|)"]
  if(i%%100==0) cat(paste0(i,"\n"))
}
P_val=P_val[order(P_val$P),]
head(P_val)
# Manhattan plot
manhattan(P_val, chr="CHR", bp="BP", snp="SNP", p="P",suggestiveline=-log10(0.05/nrow(P_val)),genomewideline= -log10(0.01/nrow(P_val)) )
#
# Exclude the significant SNP information
#
Signif_SNP=P_val[which(P_val$P<0.01/nrow(P_val)),]
#
paste0(min(Signif_SNP$BP),":",max(Signif_SNP$BP))
# using LocusZoom
Locus_data=P_val
Locus_data=Locus_data[,c(1,4)]
colnames(Locus_data)=c("MarkerName","P.value")
write.table(Locus_data,"Locus_data.txt",sep = "\t",row.names = F,quote = F)
# go to: http://locuszoom.sph.umich.edu/genform.php?type=yourdata

#####
####### Finding genes located in this area using ensemble biomart database
# Go to: ensembl biomart 
# Then Go to: https://david.ncifcrf.gov/ for annotation

#----------------------------------------------------------------
# Using R to extract biomart information
#----------------------------------------------------------------

listMarts(host="www.ensembl.org")
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "www.ensembl.org")
filters = listFilters(ensembl)
Attributes_list=listAttributes(mart = ensembl)
Gen_list=getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region"),values = list(chromosomal_region="6:32390011:32892654"), mart = ensembl)

#----------------------------------------------------------------
# GO analysis  
#----------------------------------------------------------------
library(org.Hs.eg.db)
selected_gens=Gen_list$entrezgene[which(!is.na(Gen_list$entrezgene))]
#------------------------------------
# CC: cellular component
egoCC <- enrichGO(gene          =selected_gens ,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(egoCC))

# ensembl
ego2 <- enrichGO(gene         = Gen_list$ensembl_gene_id[which(!is.na(Gen_list$ensembl_gene_id))],
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))
# Symbol
ego3 <- enrichGO(gene         = Gen_list$hgnc_symbol[which(!is.na(Gen_list$hgnc_symbol))],
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego3))
#------------------------------------
# MF: molecular function
egoMF <- enrichGO(gene          =selected_gens ,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
head(summary(egoMF))


#------------------------------------
# BP: biological process
egoBP <- enrichGO(gene          =selected_gens ,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
head(summary(egoBP))
#------------------------------------
# Visualizating the results
data(geneList)
dotplot(egoCC, showCategory=30)
dotplot(egoBP, showCategory=30)

cnetplot(egoCC, foldChange=geneList)

plotGOgraph(egoCC)





