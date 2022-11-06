
library(readxl)
library(tibble)
library(WGCNA)
library(dplyr)
rt = data.frame(read_excel("GSE78220_FPKM.xlsx"))
rt = column_to_rownames(rt,var = "Gene")

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm+1) - log(sum(fpkm+1)) + log(1e6))
}
rt <- apply(rt,2,fpkmToTpm)

data=read.table("uniSigGeneExp.txt", header=T, sep="\t", check.names=F, row.names=1)
rt = rt[rownames(rt) %in% colnames(data),]
rt = data.frame(t(rt))
score = moduleEigengenes(rt, rep("score",ncol(rt)))$eigengenes
score = arrange(score,MEscore)
score = rownames_to_column(score,var = "id")

cli = read.table("clinical.txt",header = T,row.names = NULL,sep = "\t")
cli = mutate(cli,group=ifelse(group=="Progressive Disease","PD",
                              ifelse(group=="Partial Response","PR","CR")))

score2 = inner_join(score,cli,by="id")
score2$id = factor(score2$id,levels = score2$id)
library(ggplot2)
ggplot(score2,aes(id,MEscore,fill=group)) +
  geom_bar(stat = "identity",color="black",size=1) +
  scale_fill_manual(values = c("#ce293c","#91dde4","#cc79ce")) +
  labs(x="",y="CuFescore") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))
ggsave("",width = 6,height = 4)

