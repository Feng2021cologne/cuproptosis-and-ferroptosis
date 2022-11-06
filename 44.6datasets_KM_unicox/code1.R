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
outTab=data.frame()
for(i in unique(data$dataset)){
  rt = data[data$dataset==i,]
  cox <- coxph(Surv(futime, fustat) ~ rt[,2], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}


