
library(maftools)
library(stringr)

laml <- read.maf(maf ="")
x = tmb(maf = laml) 
data = data.frame(x[,c(1,3)])
colnames(data) = c("id",'TMB')
data$id = str_sub(data$id,start = 1L,end = 16L)

library(ggpubr)
library(reshape2)
tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)      
score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)    
clu=read.table("geneCluster.txt", header=T, sep="\t", check.names=F, row.names=1)    

sameSample=intersect(row.names(tmb), row.names(score))
tmb=tmb[sameSample,,drop=F]
score=score[sameSample,,drop=F]
rownames(clu)=gsub("\\.", "-", rownames(clu))
clu=clu[sameSample,,drop=F]
data=cbind(score, tmb, clu)
data$group = ifelse(data$MEscore >= median(data$MEscore),"High","Low")
colnames(data)[1] = "CuFescore"

data=data[,c("CuFescore", "group", "geneCluster", "TMB")]

data$group=factor(data$group, levels=c("Low", "High"))
group=levels(factor(data$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#519dcd","#e23121","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]

boxplot=ggboxplot(data, x="group", y="TMB", fill="group",
		          xlab="",
		          ylab="Tumor Burden Mutation",
		          legend.title="CuFescore",
		          palette = bioCol )+ 
	    stat_compare_means(comparisons = my_comparisons)
pdf(file="",width=5,height=4.5)
print(boxplot)
dev.off()

length=length(levels(factor(data$geneCluster)))
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
p1=ggplot(data, aes(CuFescore, TMB)) + 
		  xlab("CuFescore")+ylab("Tumor Burden Mutation")+
		  geom_point(aes(colour=geneCluster))+
		  scale_color_manual(values=bioCol[1:length])+ 
		  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =CuFescore, y =TMB))
pdf(file="", width=6, height=4.5)
print(p1)
dev.off()


