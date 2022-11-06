
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
library(stringr)
library(data.table)
library(tibble)

rt=data.frame(fread("data.txt", header=T, sep="\t"))
rt = column_to_rownames(rt,var = "V1")
data=as.matrix(rt)

geneSets=getGmt("immune.gmt", geneIdType=SymbolIdentifier()) 
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
ssgseaScore=normalize(ssgseaScore)
colnames(ssgseaScore) = gsub("\\.","-",colnames(ssgseaScore))
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)

cluster=read.table('Cluster.txt', header=T, sep="\t", check.names=F, row.names=1)

ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

data=melt(scoreCluster, id.vars="cluster")
colnames(data)=c("Cup-Ferroptosis cluster", "Immune", "Fraction")

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"Cup-Ferroptosis cluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="Cup-Ferroptosis cluster", 
     ylab="Immune infiltration",
     xlab="",
     legend.title="Cup-Ferroptosis cluster",
     palette=bioCol)
pdf(file="", width=8, height=6.5)                        
p+stat_compare_means(aes(group=`Cup-Ferroptosis cluster`),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()

