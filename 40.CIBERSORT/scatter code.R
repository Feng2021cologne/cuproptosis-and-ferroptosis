
library(limma)
immuneFile="datasets/GSE30219/CIBERSORT-Results.txt"      
cluFile="../23.Cup-FerroptosisScore/Cup-FerroptosisScore.txt"                  

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

exp <- t(rt)
class <- c(rep("normal",lowNum),rep("tumor",highNum))
class<- factor(class,levels=c("normal","tumor"),labels=c("normal","tumor"))
design <- model.matrix(~0+class)
colnames(design) <- c("normal","tumor")  
fit <- lmFit(exp,design)  
cont.matrix<-makeContrasts(tumor - normal,levels=design)
fit1 <- contrasts.fit(fit, cont.matrix)
fit1 <- eBayes(fit1)  
allDiff<-topTable(fit1,adjust='fdr',number=100000)
diff = allDiff[rownames(allDiff) %in% rownames(diff),]

#########################################GSE31210
library(limma)
immuneFile="datasets/GSE31210/CIBERSORT-Results.txt"     
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

exp <- t(rt)
class <- c(rep("normal",lowNum),rep("tumor",highNum))
class<- factor(class,levels=c("normal","tumor"),labels=c("normal","tumor"))
design <- model.matrix(~0+class)
colnames(design) <- c("normal","tumor")   
fit <- lmFit(exp,design)  
cont.matrix<-makeContrasts(tumor - normal,levels=design)
fit1 <- contrasts.fit(fit, cont.matrix)   
fit1 <- eBayes(fit1)   
allDiff<-topTable(fit1,adjust='fdr',number=100000)
diff = allDiff[rownames(allDiff) %in% rownames(diff),]


#########################################GSE3141

library(limma)
immuneFile="datasets/GSE3141/CIBERSORT-Results.txt"     
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

exp <- t(rt)
class <- c(rep("normal",lowNum),rep("tumor",highNum))
class<- factor(class,levels=c("normal","tumor"),labels=c("normal","tumor"))
design <- model.matrix(~0+class)
colnames(design) <- c("normal","tumor")   
fit <- lmFit(exp,design) 
cont.matrix<-makeContrasts(tumor - normal,levels=design)
fit1 <- contrasts.fit(fit, cont.matrix)  
fit1 <- eBayes(fit1)   
allDiff<-topTable(fit1,adjust='fdr',number=100000)
dim(allDiff)
head(allDiff)
diff = allDiff[rownames(allDiff) %in% rownames(diff),]
dim(diff)
write.table(diff,"GSE3141 immune.txt",sep="\t",quote=F)


#########################################GSE37745

library(limma)
immuneFile="datasets/GSE37745/CIBERSORT-Results.txt"      
cluFile="../23.Cup-FerroptosisScore/Cup-FerroptosisScore.txt"                
     
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

exp <- t(rt)
class <- c(rep("normal",lowNum),rep("tumor",highNum))
class<- factor(class,levels=c("normal","tumor"),labels=c("normal","tumor"))
design <- model.matrix(~0+class)
colnames(design) <- c("normal","tumor")  
fit <- lmFit(exp,design)  
cont.matrix<-makeContrasts(tumor - normal,levels=design)
fit1 <- contrasts.fit(fit, cont.matrix)  
fit1 <- eBayes(fit1) 
allDiff<-topTable(fit1,adjust='fdr',number=100000)
dim(allDiff)
head(allDiff)
diff = allDiff[rownames(allDiff) %in% rownames(diff),]
dim(diff)

#########################################GSE81089
library(limma)
immuneFile="datasets/GSE81089/CIBERSORT-Results.txt"     
cluFile="../23.Cup-FerroptosisScore/Cup-FerroptosisScore.txt"                      
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

exp <- t(rt)
class <- c(rep("normal",lowNum),rep("tumor",highNum))
class<- factor(class,levels=c("normal","tumor"),labels=c("normal","tumor"))
design <- model.matrix(~0+class)
colnames(design) <- c("normal","tumor")   
fit <- lmFit(exp,design)
cont.matrix<-makeContrasts(tumor - normal,levels=design)
fit1 <- contrasts.fit(fit, cont.matrix)   
fit1 <- eBayes(fit1)  
allDiff<-topTable(fit1,adjust='fdr',number=100000)
diff = allDiff[rownames(allDiff) %in% rownames(diff),]

GSE30219 = read.table("GSE30219 immune.txt",header = T,sep = "\t",row.names = NULL)
GSE31210 = read.table("GSE31210 immune.txt",header = T,sep = "\t",row.names = NULL)
GSE3141 = read.table("GSE3141 immune.txt",header = T,sep = "\t",row.names = NULL)
GSE37745 = read.table("GSE37745 immune.txt",header = T,sep = "\t",row.names = NULL)
GSE81089 = read.table("GSE81089 immune.txt",header = T,sep = "\t",row.names = NULL)

data = Reduce(rbind,list(GSE30219,GSE31210,GSE3141,GSE37745,GSE81089))
data$dataset = rep(c("GSE30219","GSE31210","GSE3141","GSE37745","GSE81089"),each = 14)
data = data[data$adj.P.Val <= 0.05,]

library(dplyr)
cell = data %>% 
  group_by(row.names) %>%
  summarise(mean = mean(logFC)) %>%
  arrange(desc(mean)) %>%
  mutate(group=ifelse(mean>0,"up","down")) %>%
  as.data.frame()

data$row.names = factor(data$row.names,levels = cell$row.names)

library(ggplot2)
library(dplyr)
library(tibble)
data0 = full_join(data,cell,by="row.names")
data0$row.names = factor(data0$row.names,levels = cell$row.names)
p1 = ggplot(data0,aes(row.names,fill=group)) +
  geom_bar(stat = "count") +
  labs(x="",y="#Datasets") +
  theme_classic() +
  theme(
  axis.text.x =  element_blank(),
    axis.ticks.x =  element_blank())

p2 = ggplot(data,aes(row.names,dataset,color=logFC)) +
  geom_point(aes(size= -log10(P.Value))) +
  scale_size(range = c(3, 9)) +
  labs(x="",y="") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

library(patchwork)
p1 / p2
ggsave("",width = 10,height = 8)




