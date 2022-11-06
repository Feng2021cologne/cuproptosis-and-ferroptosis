library(dplyr)
library(tibble)
library(stringr)
library(data.table)
library(limma)

tcga = data.frame(fread("TPM.txt",header = T))
tcga = tcga[,which(str_sub(colnames(tcga),start = 14L,end = 15L)!="11")]
tcga[,-1] = log2(tcga[,-1] + 1)
colnames(tcga)[1] = "row.names"
GSE30219 = read.table("GSE30219.txt",header = T,row.names = NULL)
GSE31210 = read.table("GSE31210.txt",header = T,row.names = NULL)
GSE3141 = read.table("GSE3141.txt",header = T,row.names = NULL)
GSE37745 = read.table("GSE37745.txt",header = T,row.names = NULL)
GSE81089 = read.table("GSE81089.txt",header = T,row.names = NULL,sep = "\t")

data = inner_join(tcga,GSE30219,by=c("row.names"="gene"))
data = inner_join(data,GSE31210,by=c("row.names"="gene"))
data = inner_join(data,GSE3141,by=c("row.names"="gene"))
data = inner_join(data,GSE37745,by=c("row.names"="gene"))
data = inner_join(data,GSE81089,by=c("row.names"="gene"))
data = column_to_rownames(data,var = "row.names")
colnames(data) = gsub("\\.","-",colnames(data))

batch = data.frame(id = colnames(data),
                   batch = c(rep(1,ncol(tcga)-1),
                             rep(2,ncol(GSE30219)-1),
                             rep(3,ncol(GSE31210)-1),
                             rep(4,ncol(GSE3141)-1),
                             rep(5,ncol(GSE37745)-1),
                             rep(6,ncol(GSE81089)-1)))

sameSample= Reduce(intersect,list(colnames(data),batch$id))
data=data[,sameSample]
batch = batch[sameSample,]

score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)
high = rownames(score)[score$MEscore >= median(score$MEscore)]
low = rownames(score)[score$MEscore < median(score$MEscore)]

data = cbind(data[,low],data[,high])
batch = rbind(batch[low,],batch[high,])

design=model.matrix(~0+factor(Type) + batch$batch)
colnames(design)=c(levels(factor(Type)),"batch")

fit=lmFit(data, design)
cont.matrix=makeContrasts(high - low, levels=design)
fit2=contrasts.fit(fit, cont.matrix)
fit2=eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)

diffSig=allDiff[with(allDiff, (abs(logFC)>0.5 & adj.P.Val < 0.05)),]
diffSigOut=rbind(id=colnames(diffSig),diffSig)

#################GO、kegg
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
genes=rownames(diffSig)
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       

kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG = KEGG[KEGG$pvalue < 0.05,]

path = read.table("selected kegg pathways.txt",sep = "\t")
load("kk.RData")
kk2 = kk
kk2@result =  kk2@result[kk2@result$Description %in% path$V1,]

pdf("",width = 6,height = 8)
barplot(kk2) 
dev.off()  
  
  
  
  