
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
                geneSets, 
                min.sz=0, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=6)
colnames(gsvaResult) = gsub("\\.","-",colnames(gsvaResult))
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)

cluster=read.table("Cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster) = gsub("\\.","-",rownames(cluster))

gsvaResult=t(gsvaResult)
sameSample=intersect(row.names(gsvaResult), row.names(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)

library(tidyr)
library(ggplot2)
library(ggpubr)
data0 = gather(gsvaCluster,genesets,value,-ncol(gsvaCluster))
colnames(data0)[1] = "Cup-Ferroptosis cluster"
ggplot(data0,aes(genesets,value,fill=`Cup-Ferroptosis cluster`)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#519dcd","#ddac29","#e23121")) +
  labs(x="",y="Enrichment Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "top")
ggsave("",width = 6,height = 4)

