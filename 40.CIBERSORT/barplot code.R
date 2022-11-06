
library(limma)
immuneFile="CIBERSORT-Results.txt"     
cluFile="../23.Cup-FerroptosisScore/Cup-FerroptosisScore.txt"           
pFilter=0.05       

immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
dim(immune) 
immune=immune[immune[,"P-value"]<pFilter,]
dim(immune) 
immune=as.matrix(immune[,1:(ncol(immune)-3)])
rownames(immune) = gsub("\\.","-",rownames(immune))

cluster=read.table(cluFile, header=T, sep="\t", row.names=1, check.names=F)
dim(cluster)
cluster = cluster[rownames(cluster) %in% rownames(immune),,drop=F]
cluster$group = ifelse(cluster$MEscore >= median(cluster$MEscore),"High","Low")

lowName=row.names(cluster)[cluster[,ncol(cluster)]=="Low"]
highName=row.names(cluster)[cluster[,ncol(cluster)]=="High"]

lowImm=intersect(row.names(immune), lowName)
highImm=intersect(row.names(immune), highName)
rt=rbind(immune[lowImm,], immune[highImm,])
lowNum=length(lowImm) 
highNum=length(highImm) 


######################################
exp <- t(rt)
logFoldChange=0.5
adjustP=0.05

class <- c(rep("normal",lowNum),rep("tumor",highNum))
class<- factor(class,levels=c("normal","tumor"),labels=c("normal","tumor"))
design <- model.matrix(~0+class)
colnames(design) <- c("normal","tumor")  
fit <- lmFit(exp,design)  
cont.matrix<-makeContrasts(tumor - normal,levels=design)
fit1 <- contrasts.fit(fit, cont.matrix)   
fit1 <- eBayes(fit1)   
allDiff<-topTable(fit1,adjust='fdr',number=100000)
diff = allDiff[allDiff$adj.P.Val <= 0.05,]


library(ggplot2)
library(dplyr)
library(tibble)
diff2 = rownames_to_column(diff2,var = "cell")
diff2$cell = factor(diff2$cell,levels = diff2$cell)
ggplot(diff2,aes(logFC,cell,fill=-log10(P.Value))) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#844697",
                      high = "#e61b1e") +
  labs(x="Difference",y="") +
  theme_classic()
ggsave("",width = 6,height = 8)
