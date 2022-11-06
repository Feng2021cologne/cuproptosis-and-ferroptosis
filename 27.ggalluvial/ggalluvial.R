
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(data.table)
library(tibble)

m6aClu=read.table("Cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table("geneCluster.txt", header=T, sep="\t", check.names=F, row.names=1)
rownames(geneClu) = gsub("\\.","-",rownames(geneClu))
score=read.table("Cup-FerroptosisScore.group.txt", header=T, sep="\t", check.names=F, row.names=1)
rownames(score) = NULL
score = column_to_rownames(score,var = "id")

sur1 = data.frame(fread("survival.tsv"))
sur1 = sur1[,c(1,2,4)]
colnames(sur1) = c("id","fustat","futime")
sur2 = data.frame(fread("survival2.txt"))
sur2 = sur2[,c(1,3,2)]
sur3 = data.frame(fread("survival2.txt"))
sur4 = data.frame(fread("survival2.txt"))
sur5 = data.frame(fread("survival2.txt"))
sur6 = data.frame(fread(".survival2.txt"))

cli = Reduce(rbind.data.frame,list(sur1,sur2,sur3,sur4,sur5,sur6))
cli = column_to_rownames(cli,var = "id")

geneClu = geneClu[rownames(m6aClu),,drop=F]
twoCluster=cbind(m6aClu, geneClu)
sameSample=intersect(row.names(twoCluster), row.names(score))
scoreClu=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])
sameSample=intersect(row.names(scoreClu), row.names(cli))
rt=cbind(scoreClu[sameSample,], cli[sameSample,])

rt=rt[,c("cluster", "geneCluster", "group", "fustat")]
colnames(rt)=c("CuFeCluster", "geneCluster", "CuFescore", "fustat")
rt$fustat = as.character(rt$fustat)
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

pdf(file="", width=8, height=8)
mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
   scale_x_discrete(expand = c(0, 0)) +  
   geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + 
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()

