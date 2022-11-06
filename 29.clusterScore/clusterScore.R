
library(limma)
library(ggpubr)
m6aClu=read.table("Cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table("geneCluster.txt", header=T, sep="\t", check.names=F, row.names=1)
rownames(geneClu) = gsub("\\.","-",rownames(geneClu))
score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)
rownames(score) = gsub("\\.","-",rownames(score))

geneClu = geneClu[rownames(m6aClu),,drop=F]
twoCluster=cbind(m6aClu, geneClu)
sameSample=intersect(row.names(twoCluster), row.names(score))
data=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])

colnames(data) = c("CuFescore","CuFecluster","geneCluster")
data$CuFecluster=factor(data$CuFecluster, levels=levels(factor(data$CuFecluster)))
group=levels(factor(data$CuFecluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$CuFecluster)))]
	
boxplot=ggboxplot(data, x="CuFecluster", y="CuFescore", color="CuFecluster",
			      xlab="CuFecluster",
			      ylab="CuFescore",
			      legend.title="CuFecluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)

pdf(file="", width=5, height=4.5)
print(boxplot)
dev.off()

data$geneCluster=factor(data$geneCluster, levels=levels(factor(data$geneCluster)))
group=levels(factor(data$geneCluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$geneCluster)))]
	
boxplot=ggboxplot(data, x="geneCluster", y="CuFescore", color="geneCluster",
			      xlab="geneCluster",
			      ylab="CuFescore",
			      legend.title="geneCluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
pdf(file="", width=5, height=4.5)
print(boxplot)
dev.off()

