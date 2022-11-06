
library(ggpubr)                  
library(stringr)

ips=read.table("TCIA-ClinicalData.tsv", header=T, sep="\t", check.names=F, row.names=1) 
ips = ips[,c(23:26)]

score=read.table("Cup-FerroptosisScore.group.txt", header=T, sep="\t", check.names=F, row.names=1)
score = score[grep("TCGA",score$id),]
score$id = str_sub(score$id,start =1L,end = 12L)
rownames(score) = score$id
score = score[,-ncol(score)]
score$group = ifelse(score$MEscore >= median(score$MEscore),"High","Low")

sameSample=intersect(row.names(ips), row.names(score))
ips=ips[sameSample, , drop=F]
score=score[sameSample, "group", drop=F]
data=cbind(ips, score)

data$group=factor(data$group, levels=c("Low", "High"))
group=levels(factor(data$group))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

for(i in colnames(data)[1:(ncol(data)-1)]){
	rt=data[,c(i, "group")]
	colnames(rt)=c("IPS", "group")
	gg1=ggviolin(rt, x="group", y="IPS", fill = "group", 
	         xlab="", ylab=i,
	         legend.title="CuFescore",
	         palette=c("#0066FF", "#FF0000"),
	         add = "boxplot", add.params = list(fill="white"))+ 
	         stat_compare_means(comparisons = my_comparisons)
	        
	pdf(file=paste0(i, ".pdf"), width=6, height=5)
	print(gg1)
	dev.off()
}


