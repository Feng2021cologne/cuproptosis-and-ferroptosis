
library(maftools)     
library(dplyr)
library(stringr)

score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F,row.names = NULL)
colnames(score) = c("Tumor_Sample_Barcode", "Cup-FerroptosisScore")
score = score[grep("TCGA",score$Tumor_Sample_Barcode),]
score$group = ifelse(score$`Cup-FerroptosisScore` >= median(score$`Cup-FerroptosisScore`),"High","Low")

maf=read.maf(maf="")
high = maf@data[str_sub(maf@data$Tumor_Sample_Barcode,start = 1L,end=16L) %in% score$Tumor_Sample_Barcode[score$group=="High"],]
low = maf@data[str_sub(maf@data$Tumor_Sample_Barcode,start = 1L,end=16L) %in% score$Tumor_Sample_Barcode[score$group=="Low"],]

score$Tumor_Sample_Barcode = str_sub(score$Tumor_Sample_Barcode,start = 1L,end = 12L)
colnames(score)[3] = "CuFescore"
score = score[,-2]

gene = read.table("selected gene.txt",header = F)$V1
ann_colors=list()
col=c("#0066FF","#FF0000")
names(col)=c("Low", "High")
ann_colors[["CuFescore"]]=col

pdf(file="", width=6, height=6)
maf=read.maf(maf="low.maf", clinicalData="ann.txt")
oncoplot(maf=maf, 
         genes=gene, 
         keepGeneOrder=T)
dev.off()

pdf(file="", width=6, height=6)
maf=read.maf(maf="high.maf", clinicalData="ann.txt")
oncoplot(maf=maf, 
         genes=gene,
         keepGeneOrder=T)
dev.off()

library(ggplot2)
low=read.maf(maf="low.maf")
high=read.maf(maf="high.maf")

pdf("",width=8,height=3)
lollipopPlot2(m1=low, m2=high, m1_name="Low", m2_name="High", gene="TP53")
dev.off()

pdf("",width=7,height=7)
par(oma=c(2,2,2,2))
somaticInteractions(maf = low, genes = gene, pvalue = c(0.05, 0.1),showSum =F)
dev.off()

pdf("",width=7,height=7)
par(oma=c(2,2,2,2))
somaticInteractions(maf = high, genes = gene, pvalue = c(0.05, 0.1),showSum =F)
dev.off()

pt.vs.rt <- mafCompare(m1 = low, m2 = high, m1Name = 'low', m2Name = 'high', minMut = 5)
print(pt.vs.rt)
pt.vs.rt$results = pt.vs.rt$results[pt.vs.rt$results$adjPval < 0.05,]

gene2 = gene
pt.vs.rt$results = pt.vs.rt$results[pt.vs.rt$results$Hugo_Symbol %in% gene2,]
pdf("",width = 7,height = 5)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 1)
dev.off()

maf=read.maf(maf="")
pdf("",width=7,height=7)
par(oma=c(2,2,2,2))
somaticInteractions(maf = maf, genes = gene, pvalue = c(0.05, 0.1)) 
dev.off()

data = maf@variant.classification.summary
data$`Non-synonymous mutation` = data$Missense_Mutation + data$Nonsense_Mutation + data$Nonstop_Mutation
data = data.frame(data)
data = data[,c(1,(ncol(data)-1),ncol(data))]
data$Tumor_Sample_Barcode = str_sub(data$Tumor_Sample_Barcode,start = 1L,end = 16L)

score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F,row.names = NULL)
colnames(score) = c("Tumor_Sample_Barcode", "Cup-FerroptosisScore")
score = score[grep("TCGA",score$Tumor_Sample_Barcode),]
score$group = ifelse(score$`Cup-FerroptosisScore` >= median(score$`Cup-FerroptosisScore`),"High","Low")
data2 = inner_join(data,score,by="Tumor_Sample_Barcode")

library(ggpubr)
library(patchwork)
p1 = ggplot(data2,aes(`Cup-FerroptosisScore`,Non.synonymous.mutation)) +
  geom_point(aes(color=group)) +
  geom_smooth(method = "lm") +
  scale_color_manual(values=c("#cdcacb","#7289b7")) +
  stat_cor() +
  labs(x="CuFescore",y="Non-synonymous mutation counts") +
  guides(color=FALSE) +
  theme_classic()
  
p2 = ggplot(data2,aes(group,Non.synonymous.mutation)) +
  geom_boxplot(aes(fill=group)) +
  scale_fill_manual(values=c("#cdcacb","#7289b7")) +
  theme_void()

p1 + p2 + plot_layout(widths = c(5, 1))
ggsave("",width = 6,height = 4)

p1 = ggplot(data2,aes(`Cup-FerroptosisScore`,total)) +
  geom_point(aes(color=group)) +
  geom_smooth(method = "lm") +
  scale_color_manual(values=c("#cdcacb","#7289b7")) +
  stat_cor() +
  xlab("CuFescore") +
  ylab("All mutation counts") +
  guides(color=FALSE) +
  theme_classic()

p2 = ggplot(data2,aes(group,total)) +
  geom_boxplot(aes(fill=group)) +
  scale_fill_manual(values=c("#cdcacb","#7289b7")) +
  theme_void()

p1 + p2 + plot_layout(widths = c(5, 1))
ggsave("",width = 6,height = 4)

