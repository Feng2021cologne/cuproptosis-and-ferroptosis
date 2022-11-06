library(limma) 
library(VennDiagram)
library(data.table)
library(stringr)
library(dplyr)
library(tibble)

cluster=read.table("Cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
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
                             rep(6,ncol(GSE81089)-1)),
                   row.names = colnames(data))

sameSample= Reduce(intersect,list(colnames(data), row.names(cluster),batch$id))
data=data[,sameSample]
cluster=cluster[sameSample,,drop=F]
batch = batch[sameSample,]

geneList=list()
design=model.matrix(~0+factor(Type) + batch$batch)
colnames(design)=c(levels(factor(Type)),"batch")
comp=combn(levels(factor(Type)), 2)
allDiffGenes=c()
for(i in 1:ncol(comp)){
  fit=lmFit(data, design)
  contrast=paste0(comp[2,i], "-", comp[1,i])
  
  cont.matrix=makeContrasts(contrast, levels=design)
  fit2=contrasts.fit(fit, cont.matrix)
  fit2=eBayes(fit2)
  
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiffOut=rbind(id=colnames(allDiff),allDiff)
  write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
  
  diffSig=allDiff[with(allDiff, (abs(logFC)>0.5 & adj.P.Val < 0.05 )), ]
  diffSigOut=rbind(id=colnames(diffSig),diffSig)
  write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
  geneList[[contrast]]=row.names(diffSig)
}

venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file="", width=5, height=5)
grid.draw(venn.plot)
dev.off()
