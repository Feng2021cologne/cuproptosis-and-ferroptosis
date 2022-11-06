##genecluster
library(limma)
library(ConsensusClusterPlus)

data=read.table("../18.uniCox/uniSigGeneExp.txt", header=T, sep="\t", row.names = 1)
data = t(data)

maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=".",
              clusterAlg="pam",
              distance="euclidean",
              seed=123456,
              plot="pdf")
icl <- calcICL(results,plot = 'pdf')
clusterNum=3      
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("geneCluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$geneCluster))
cluster$geneCluster=letter[match(cluster$geneCluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="geneCluster.txt", sep="\t", quote=F, col.names=F)

