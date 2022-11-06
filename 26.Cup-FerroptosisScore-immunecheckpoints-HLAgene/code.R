library(stringr)
library(limma)
library(tidyr)
library(ggplot2)
library(ggpubr)
score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)
score$group = ifelse(score$MEscore >= median(score$MEscore),"High","Low")
score = score[,-1,drop=F]

rt=read.table("data.txt", header=T, sep="\t",row.names = 1)
data=as.matrix(rt)
colnames(data) = gsub("\\.","-",colnames(data))

immune = read.table("immune checkpoints.txt",header=F)
data = data.frame(t(data))
sameSample=intersect(row.names(data), row.names(score))
data=data[sameSample,,drop=F]
score=score[sameSample,,drop=F]
scoreCluster=cbind(data, score)

scoreCluster=gather(scoreCluster,gene,value,-ncol(scoreCluster))
scoreCluster$`Group of Cup-FerroptosisScore` = factor(scoreCluster$`Group of Cup-FerroptosisScore`,levels = c("Low","High"))
ggplot(scoreCluster,aes(gene,value,fill=`Group of Cup-FerroptosisScore`)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#41a2d3","#bd3225")) +
  labs(x="",y="log2(TPM)") +
  stat_compare_means(aes(group=`Group of Cup-FerroptosisScore`),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("",width = 12,height = 6)

###############################################lncRNAscore和HLA
score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)
score$group = ifelse(score$MEscore >= median(score$MEscore),"High","Low")
score = score[,-1,drop=F]

rt=read.table("data.txt", header=T, sep="\t",row.names = 1)
data=as.matrix(rt)
colnames(data) = gsub("\\.","-",colnames(data))

HLA = read.table("HLA.txt",header=F)
data = data[rownames(data) %in% HLA$V1,]
data = data.frame(t(data))

sameSample=intersect(row.names(data), row.names(score))
data=data[sameSample,,drop=F]
score=score[sameSample,,drop=F]
scoreCluster=cbind(data, score)

scoreCluster=gather(scoreCluster,gene,value,-ncol(scoreCluster))
colnames(scoreCluster)[1] = "Group of Cup-FerroptosisScore"
scoreCluster$`Group of Cup-FerroptosisScore` = factor(scoreCluster$`Group of Cup-FerroptosisScore`,levels = c("Low","High"))
scoreCluster$gene = gsub("\\.","-",scoreCluster$gene)
ggplot(scoreCluster,aes(gene,value,fill=`Group of Cup-FerroptosisScore`)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#41a2d3","#bd3225")) +
  labs(x="",y="log2(TPM)") +
  stat_compare_means(aes(group=`Group of Cup-FerroptosisScore`),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("",width = 12,height = 6)

#############lncRNAscore and immune checkpoints
score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)
rt=read.table("data.txt", header=T, sep="\t",row.names = 1)
data=as.matrix(rt)
colnames(data) = gsub("\\.","-",colnames(data))
immune = read.table("immune checkpoints.txt",header=F)
data = data.frame(t(data))
sameSample=intersect(row.names(data), row.names(score))
data=data[sameSample,,drop=F]
score=score[sameSample,,drop=F]
scoreCluster=cbind(data, score)
immune = scoreCluster

score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)
rt=read.table("data.txt", header=T, sep="\t",row.names = 1)
data=as.matrix(rt)
colnames(data) = gsub("\\.","-",colnames(data))

HLA = read.table("HLA.txt",header=F)
data = data[rownames(data) %in% HLA$V1,]
dim(data)
data = data.frame(t(data))

sameSample=intersect(row.names(data), row.names(score))
data=data[sameSample,,drop=F]
score=score[sameSample,,drop=F]
scoreCluster=cbind(data, score)
HLA = scoreCluster

sameSample=intersect(row.names(immune), row.names(HLA))
immune=immune[sameSample,,drop=F]
immune=immune[,-ncol(immune)]
HLA=HLA[sameSample,,drop=F]
scoreCluster=cbind(immune, HLA)
data0 = data.frame()
for (i in 1:(ncol(scoreCluster)-1)) {
  corr = cor.test(scoreCluster[,ncol(scoreCluster)],scoreCluster[,i])$estimate
  pvalue = cor.test(scoreCluster[,ncol(scoreCluster)],scoreCluster[,i])$p.value
  data = data.frame(gene = colnames(scoreCluster)[i],
                    corr = corr,
                    pvalue = pvalue)
  data0 = rbind.data.frame(data0,data)
}

data0$gene = gsub("\\.","-",data0$gene)
data0$gene = factor(data0$gene,levels =data0$gene)
data0$corr = round(data0$corr,2)
library(ggplot2)
library(dplyr)
library(patchwork)
p1 = ggplot(data0,aes("1",gene)) +
  geom_tile(aes(fill=corr),) +
  scale_fill_gradient2(low = "#2f4191",mid = "white",high = "#e52324",name="correlation coefficient") +
  geom_text(aes(label = corr)) +
  labs(x="",y="") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

data01 = mutate(data0,pvalue=ifelse(pvalue>0.05,"ns",
                                    ifelse(pvalue<=0.05 & pvalue>0.01,"*",
                                                  ifelse(pvalue<=0.01 & pvalue>0.001,"**","***"))))
data01 = mutate(data01,corr=ifelse(pvalue=="ns",0,corr))
p2 = ggplot(data01,aes("1",gene)) +
  geom_tile(aes(fill=corr)) +
  scale_fill_gradient2(low = "#2f4191",mid = "white",high = "#e52324",name="correlation coefficient") +
  geom_text(aes(label = pvalue)) +
  labs(x="",y="") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) 

data02 = cbind(data01,group = c(rep("Immune checkpoints",35),rep("HLA family",17)))
p3 = ggplot(data02,aes("1",gene)) +
  geom_tile(aes(fill=group)) +
  scale_fill_manual(values = c("#f4a620","#2e7cbd")) +
  labs(x="",y="") +
  guides(fill=FALSE) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
ggarrange(p3, p1, p2,ncol  = 3,common.legend = TRUE,legend = "right")  
ggsave("",width = 5.5,height = 12)

