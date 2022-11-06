## PCA
library(limma)
library(ggplot2)
library(tibble)
    
rt=read.table("../9.cluster/uniSigExp.txt", header=T, sep="\t",row.names = 1)
rt = rt[,-c(1,2)]
rt=as.matrix(rt)

data.pca=prcomp(rt, scale. = TRUE)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="newTab.txt", quote=F, sep="\t")

cluster=read.table("../9.cluster/Cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
m6Acluster=as.vector(cluster[,1])

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
m6aCluCol=bioCol[1:length(levels(factor(m6Acluster)))]

PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], m6Acluster=m6Acluster)
PCA.mean=aggregate(PCA[,1:2], list(m6Acluster=PCA$m6Acluster), mean)
pdf(file="PCA.pdf", height=5, width=6.5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = m6Acluster)) +
	scale_colour_manual(name="Cup-Ferroptosis cluster", values =m6aCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$m6Acluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

