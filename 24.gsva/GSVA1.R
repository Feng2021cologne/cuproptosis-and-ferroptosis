
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(stringr)

rt=read.table("data.txt", header=T, sep="\t",row.names = 1)
data=as.matrix(rt)
genelist = read.table("genelist.txt",sep = "\t")
geneSets = strsplit(genelist$V2,", ")
names(geneSets) = genelist$V1
gsvaResult=gsva(data, 
                geneSets)
colnames(gsvaResult) = gsub("\\.","-",colnames(gsvaResult))
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)

cluster=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster) = gsub("\\.","-",rownames(cluster))

gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
colnames(cluster) = "CuFescore"
gsvaCluster=cbind(gsvaResult, cluster)

library(corrplot)
data = cor(gsvaCluster)
testRes = cor.mtest(gsvaCluster, conf.level = 0.95)
cor.plot<-corrplot(corr =data,type="upper",tl.pos="tp",
                   tl.col="black",p.mat = testRes$p, insig = "label_sig", 
                   sig.level = c(.01, .05),
                   pch.cex=1,pch.col = "black",order = "original",
                   col=rev(COL2('RdBu', 200)))
cor.plot<-corrplot(corr = data,type="lower",add=TRUE,method="number",
                   tl.pos="n",tl.col="black",
                   col="black",tl.cex=1.2,diag=FALSE, cl.pos="n",pch.col = "black",
                   number.cex = 0.8,order = "original")
