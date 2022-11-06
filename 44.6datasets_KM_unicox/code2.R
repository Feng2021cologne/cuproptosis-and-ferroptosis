library(data.table)
library(dplyr)
library(survival)
library(survminer)
rt = read.table("Cup-FerroptosisScore.txt",header = T,row.names = NULL)

tcga = data.frame(fread("survival.tsv"))

GSE30219 = data.frame(fread("survival.txt"))
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
data = data %>%
  group_by(dataset) %>%
  mutate(group=ifelse(MEscore >= median(MEscore),"High","Low"))

i= unique(data$dataset)[2]
data2 = data[data$dataset==i,]
length=length(levels(factor(data2$group)))
diff=survdiff(Surv(futime, fustat) ~ group, data = data2)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data2)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data2=data2,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="CuFescore",
                   legend.labs=levels(factor(data2[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 5,
                   palette = bioCol,
                   #surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file=paste(i,"survival.pdf",sep = " "),onefile = FALSE,width=10,height=8)
print(surPlot)
dev.off()


i= unique(data$dataset)[3]
data2 = data[data$dataset==i,]
length=length(levels(factor(data2$group)))
diff=survdiff(Surv(futime, fustat) ~ group, data = data2)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data2)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data2=data2,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="CuFescore",
                   legend.labs=levels(factor(data2[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 2,
                   palette = bioCol,
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25) +
  theme_survminer(legend.position = "left")
pdf(file=paste(i,"survival.pdf",sep = " "),onefile = FALSE,width=10,height=8)
print(surPlot)
dev.off()

i= unique(data$dataset)[4]
data2 = data[data$dataset==i,]
length=length(levels(factor(data2$group)))
diff=survdiff(Surv(futime, fustat) ~ group, data = data2)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data2)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data2=data2,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="CuFescore",
                   legend.labs=levels(factor(data2[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 2,
                   palette = bioCol,
                   #surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file=paste(i,"survival.pdf",sep = " "),onefile = FALSE,width=10,height=8)
print(surPlot)
dev.off()

i= unique(data$dataset)[5]
data2 = data[data$dataset==i,]
length=length(levels(factor(data2$group)))
diff=survdiff(Surv(futime, fustat) ~ group, data = data2)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data2)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data2=data2,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="CuFescore",
                   legend.labs=levels(factor(data2[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 5,
                   palette = bioCol,
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file=paste(i,"survival.pdf",sep = " "),onefile = FALSE,width=10,height=8)
print(surPlot)
dev.off()


i= unique(data$dataset)[1]
data2 = data[data$dataset==i,]
length=length(levels(factor(data2$group)))
diff=survdiff(Surv(futime, fustat) ~ group, data = data2)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data2)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data2=data2,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="CuFescore",
                   legend.labs=levels(factor(data2[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 5,
                   palette = bioCol,
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file=paste(i,"survival.pdf",sep = " "),onefile = FALSE,width=10,height=8)
print(surPlot)
dev.off()


i= unique(data$dataset)[6]
data2 = data[data$dataset==i,]
length=length(levels(factor(data2$group)))
diff=survdiff(Surv(futime, fustat) ~ group, data = data2)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data2)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data2=data2,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="CuFescore",
                   legend.labs=levels(factor(data2[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 2,
                   palette = bioCol,
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file=paste(i,"survival.pdf",sep = " "),onefile = FALSE,width=10,height=8)
print(surPlot)
dev.off()

