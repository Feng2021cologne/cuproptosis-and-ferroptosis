
library(limma)
library(estimate)
library(tibble)
library(dplyr)

filterCommonGenes(input.f=".data.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds="commonGenes.gct",
              output.ds="estimateScore.gct",
              platform = "illumina")

scores=read.table("estimateScore.gct", skip=2, header=T, check.names=F)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
out=rbind(ID=colnames(scores), scores)

scores2 = data.frame(scores)
scores2 = rownames_to_column(scores2,var = "row.names")
cufe = read.table("Cup-FerroptosisScore.txt",header = T,row.names = NULL)

scores3 = inner_join(cufe,scores2,by="row.names")

library(ggplot2)
library(ggpubr)
ggplot(scores3,aes(StromalScore,MEscore)) +
  geom_point(color="#a05ba9") +
  geom_smooth(method = "lm",color="#ebc8a7") +
  geom_rug(color="#81ca81") +
  labs(y="CuFescore") +
  theme_classic()
ggsave("",width = 6,height = 4)

ggplot(scores3,aes(ImmuneScore,MEscore)) +
  geom_point(color="#a05ba9") +
  geom_smooth(method = "lm",color="#ebc8a7") +
  geom_rug(color="#81ca81") +
  labs(y="CuFescore") +
  theme_classic()
ggsave("",width = 6,height = 4)

ggplot(scores3,aes(ESTIMATEScore,MEscore)) +
  geom_point(color="#a05ba9") +
  geom_smooth(method = "lm",color="#ebc8a7") +
  geom_rug(color="#81ca81") +
  labs(y="CuFescore") +
  theme_classic()
ggsave("",width = 6,height = 4)
