library(dplyr)
library(ggsci)
library(cowplot)

lncRNAcluster = read.table("Cluster.txt",header = T)
lncRNAscore = read.table("Cup-FerroptosisScore.txt",header = T)
lncRNAscore$group = ifelse(lncRNAscore$MEscore < median(lncRNAscore$MEscore),"low","high")

lncRNAscore_high = rownames(lncRNAscore)[lncRNAscore$group=="high"]
lncRNAscore_low = rownames(lncRNAscore)[lncRNAscore$group=="low"]

lncRNAcluster_A = lncRNAcluster$ID[lncRNAcluster$cluster =="A"]
lncRNAcluster_B = lncRNAcluster$ID[lncRNAcluster$cluster =="B"]
lncRNAcluster_C = lncRNAcluster$ID[lncRNAcluster$cluster =="C"]

library(VennDiagram)
venn.plot <- venn.diagram(
  list(`CuFescore_high`=lncRNAscore_high,`CuFescore_low`=lncRNAscore_low,`CuFecluster_A`=lncRNAcluster_A,`CuFecluster_B`=lncRNAcluster_B,`CuFecluster_C`=lncRNAcluster_C),
  filename = NULL,
  lty = 1,
  lwd = 0.75,
  col = "black",  
  fill = c("#7178ba", "#f5f085", "#92c976","#ed7c7c","#6ac2c5"),
  alpha = 0.60,
  cat.col = "black",
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 0.8
)

pdf("",width = 10,height = 10)
grid.draw(venn.plot)
dev.off()

lncRNAcluster = read.table("geneCluster.txt",header = T)
lncRNAcluster$ID = gsub("\\.","-",lncRNAcluster$ID)
lncRNAscore = read.table("Cup-FerroptosisScore.txt",header = T)
lncRNAscore$group = ifelse(lncRNAscore$MEscore < median(lncRNAscore$MEscore),"low","high")
rownames(lncRNAscore) = gsub("\\.","-",rownames(lncRNAscore))

lncRNAscore_high = rownames(lncRNAscore)[lncRNAscore$group=="high"]
lncRNAscore_low = rownames(lncRNAscore)[lncRNAscore$group=="low"]

lncRNAcluster_A = lncRNAcluster$ID[lncRNAcluster$geneCluster=="A"]
lncRNAcluster_B = lncRNAcluster$ID[lncRNAcluster$geneCluster=="B"]
lncRNAcluster_C = lncRNAcluster$ID[lncRNAcluster$geneCluster=="C"]

library(VennDiagram)
venn.plot <- venn.diagram(
  list(`CuFescore_high`=lncRNAscore_high,`CuFescore_low`=lncRNAscore_low,genecluster_A=lncRNAcluster_A,genecluster_B=lncRNAcluster_B,genecluster_C=lncRNAcluster_C),
  filename = NULL,
  lty = 1,
  lwd = 0.75,
  col = "black", 
  fill = c("#7178ba", "#f5f085", "#92c976","#ed7c7c","#6ac2c5"),
  alpha = 0.60,
  cat.col = "black",
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 0.8
)


pdf("",width = 10,height = 10)
grid.draw(venn.plot)
dev.off()

lncRNAcluster = read.table("Cluster.txt",header = T)

lncRNAscore = read.table("Cup-FerroptosisScore.txt",header = T,row.names = NULL)
lncRNAscore$group = ifelse(lncRNAscore$MEscore < median(lncRNAscore$MEscore),"low","high")

data = inner_join(lncRNAscore,lncRNAcluster,by=c("row.names"="ID"))
data = mutate(data,group=ifelse(group=="low","CuFescore_low","CuFescore_high"))
library(ggplot2)
library(ggpubr)
chisq.test(table(data$group,data$cluster))   

count <- group_by(data, group, cluster) %>%
  summarise(count = n()) %>%
  arrange(group, desc(cluster)) %>%
  mutate(cumsum = cumsum(count),
         prop = count / sum(count),
         cumprop = cumsum(count) / sum(count))

data =  data.frame(count)
ggplot(data = data) +
  geom_col(aes(x = group, y = count, fill = cluster),
           position = 'fill',
           color = 'black',
           width = .7) +
  scale_fill_lancet() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="",y="") +
  guides(fill=guide_legend(title = "CuFeCluster")) +
  theme_half_open() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("",width = 8,height = 6)

lncRNAcluster = read.table("geneCluster.txt",header = T)
lncRNAcluster$ID = gsub("\\.","-",lncRNAcluster$ID)

lncRNAscore = read.table("Cup-FerroptosisScore.txt",header = T,row.names = NULL)
lncRNAscore$group = ifelse(lncRNAscore$MEscore < median(lncRNAscore$MEscore),"low","high")
lncRNAscore$row.names = gsub("\\.","-",lncRNAscore$row.names)

data = inner_join(lncRNAscore,lncRNAcluster,by=c("row.names"="ID"))
data = mutate(data,group=ifelse(group=="low","CuFescore_low","CuFescore_high"))
library(ggplot2)
library(ggpubr)
chisq.test(table(data$group,data$geneCluster)) 

count <- group_by(data, group, geneCluster) %>%
  summarise(count = n()) %>%
  arrange(group, desc(geneCluster)) %>%
  mutate(cumsum = cumsum(count),
         prop = count / sum(count),
         cumprop = cumsum(count) / sum(count))

data =  data.frame(count)
ggplot(data = data) +
  geom_col(aes(x = group, y = count, fill = geneCluster),
           position = 'fill',
           color = 'black',
           width = .7) +
  scale_fill_lancet() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="",y="") +
  theme_half_open() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("",width = 8,height = 6)




