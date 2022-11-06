library(data.table)
library(tibble)
library(ggplot2)
library(tidyr)
library(dplyr)

rt=read.table("data.txt", header=T, sep="\t",row.names = 1)

gene1 = read.table("gene1.txt",encoding = "UTF-8")
data1 = rt[rownames(rt) %in% gene1$V1,]
data1 = data.frame(t(data1))
data1 = rownames_to_column(data1,var = "id")

cluster = read.table("geneCluster.txt",header = T)
data11 = inner_join(cluster,data1,by=c("ID"="id"))
data11 = data11[,-1]

data111 = gather(data11,gene,value,-1)

library(ggpubr)
ggplot(data111,aes(gene,value,fill=geneCluster)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#3c8fbe","#e42f2e","#d4ab3d")) +
  labs(x="",y="Relative Expression") +
  theme_classic()
ggsave("",width = 8,height = 6)
  
 
gene1 = read.table("gene2.txt",encoding = "UTF-8")
data1 = rt[rownames(rt) %in% gene1$V1,]
data1 = data.frame(t(data1))
data1 = rownames_to_column(data1,var = "id")
data1$id = gsub("\\.","-",data1$id)

cluster = read.table("geneCluster.txt",header = T)
data11 = inner_join(cluster,data1,by=c("ID"="id"))
data11 = data11[,-1]

data111 = gather(data11,gene,value,-1)

library(ggpubr)
ggplot(data111,aes(gene,value,fill=geneCluster)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#3c8fbe","#e42f2e","#d4ab3d")) +
  labs(x="",y="Relative Expression") +
  theme_classic()
ggsave("",width = 8,height = 6)
  
gene1 = read.table("gene3.txt",encoding = "UTF-8")
data1 = rt[rownames(rt) %in% gene1$V1,]
data1 = data.frame(t(data1))
data1 = rownames_to_column(data1,var = "id")
data1$id = gsub("\\.","-",data1$id)

cluster = read.table("geneCluster.txt",header = T)
data11 = inner_join(cluster,data1,by=c("ID"="id"))
data11 = data11[,-1]

data111 = gather(data11,gene,value,-1)

library(ggpubr)
ggplot(data111,aes(gene,value,fill=geneCluster)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#3c8fbe","#e42f2e","#d4ab3d")) +
  labs(x="",y="Relative Expression") +
  theme_classic()
ggsave("",width = 8,height = 6)


gene1 = read.table("gene3.txt",encoding = "UTF-8")
data1 = rt[rownames(rt) %in% gene1$V1,]
data1 = data.frame(t(data1))
data1 = rownames_to_column(data1,var = "id")

cluster = read.table("Cluster.txt",header = T)

data11 = inner_join(cluster,data1,by=c("ID"="id"))
data11 = data11[,-1]

data111 = gather(data11,gene,value,-1)

library(ggpubr)
ggplot(data111,aes(gene,value,fill=cluster)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#3c8fbe","#e42f2e","#d4ab3d")) +
  labs(x="",y="Relative Expression") +
  theme_classic()
ggsave("",width = 8,height = 6)


