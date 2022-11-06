library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(timeROC)

rt = read.table("Cup-FerroptosisScore.txt",header = T,row.names = NULL)

tcga = data.frame(fread("survival.tsv"))

GSE30219 = data.frame(fread("survival.txt"))
GSE30219 = GSE30219[,c(1,3,2)]
GSE31210 = data.frame(fread("survival.txt"))
GSE3141 = data.frame(fread("survival.txt"))
GSE37745 = data.frame(fread("survival.txt"))
GSE81089 = data.frame(fread("survival.txt"))

tcga = inner_join(rt,tcga,by=c("row.names"="id"))
GSE30219 = inner_join(rt,GSE30219,by=c("row.names"="id"))
GSE31210 = inner_join(rt,GSE31210,by=c("row.names"="id"))
GSE3141 = inner_join(rt,GSE3141,by=c("row.names"="id"))
GSE37745 = inner_join(rt,GSE37745,by=c("row.names"="id"))
GSE81089 = inner_join(rt,GSE81089,by=c("row.names"="id"))

data = Reduce(rbind,list(tcga,GSE30219,GSE31210,GSE3141,GSE37745,GSE81089))
data$dataset = c(rep("tcga",nrow(tcga)),
                 rep("GSE30219",nrow(GSE30219)),
                 rep("GSE31210",nrow(GSE31210)),
                 rep("GSE3141",nrow(GSE3141)),
                 rep("GSE37745",nrow(GSE37745)),
                 rep("GSE81089",nrow(GSE81089)))

data2 = data[data$dataset=="tcga",]
ROC1 <- timeROC(T=data2$futime, 
                 delta=data2$fustat, marker=data2$MEscore,
                 cause=1,
                 weighting="marginal",
                 times=seq(2,10,2),
                 iid=TRUE)
data2 = data[data$dataset=="GSE30219",]
ROC2 <- timeROC(T=data2$futime, 
                 delta=data2$fustat, marker=data2$MEscore,
                 cause=1,
                 weighting="marginal",
                 times=seq(2,10,2),
                 iid=TRUE)
data2 = data[data$dataset=="GSE31210",]
ROC3 <- timeROC(T=data2$futime, 
                 delta=data2$fustat, marker=data2$MEscore,
                 cause=1,
                 weighting="marginal",
                 times=seq(2,9,2),
                 iid=TRUE)
data2 = data[data$dataset=="GSE3141",]
ROC4 <- timeROC(T=data2$futime, 
                 delta=data2$fustat, marker=data2$MEscore,
                 cause=1,
                 weighting="marginal",
                 times=seq(2,10,2),
                 iid=TRUE)
data2 = data[data$dataset=="GSE37745",]
ROC5 <- timeROC(T=data2$futime, 
                 delta=data2$fustat, marker=data2$MEscore,
                 cause=1,
                 weighting="marginal",
                 times=seq(2,10,2),
                 iid=TRUE)
data2 = data[data$dataset=="GSE81089",]
ROC6 <- timeROC(T=data2$futime, 
                 delta=data2$fustat, marker=data2$MEscore,
                 cause=1,
                 weighting="marginal",
                 times=seq(2,10,2),
                 iid=TRUE)


pdf("timeROC.pdf", 6, 5)
plotAUCcurve(ROC1, conf.int=FALSE, col="red")
plotAUCcurve(ROC2, conf.int=FALSE, col="darkblue", add=TRUE)
plotAUCcurve(ROC3, conf.int=FALSE, col="darkgreen", add=TRUE)
plotAUCcurve(ROC5, conf.int=FALSE, col="#f1c744", add=TRUE)

legend("topright", c("TCGA","GSE30219","GSE31210","GSE37745"),
       col=c("red","darkblue","darkgreen","#f1c744"),
       bty='n', lty=1, lwd=2, cex=0.8)
dev.off()

