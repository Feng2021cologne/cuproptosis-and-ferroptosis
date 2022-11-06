##cluster survival
library(survival)
library(survminer)
library(data.table)        
library(dplyr)

cluster=read.table("../9.cluster/Cluster.txt", header=T, sep="\t", check.names=F, row.names=NULL)

sur1 = data.frame(fread("../5.cuprotosisgene_network/TCGA-LUAD.survival.tsv"))
sur1 = sur1[,c(1,2,4)]
sur1$OS.time = sur1$OS.time / 365
colnames(sur1) = c("id","fustat","futime")
sur2 = data.frame(fread("../6.merge/prognosis_data/GSE30219/survival2.txt"))
sur2 = sur2[,c(1,3,2)]
sur3 = data.frame(fread("../6.merge/prognosis_data/GSE31210/survival2.txt"))
sur4 = data.frame(fread("../6.merge/prognosis_data/GSE3141/survival2.txt"))
sur5 = data.frame(fread("../6.merge/prognosis_data/GSE37745/survival2.txt"))
sur6 = data.frame(fread("../6.merge/prognosis_data/GSE81089/survival2.txt"))

sur = Reduce(rbind.data.frame,list(sur1,sur2,sur3,sur4,sur5,sur6))

data = inner_join(sur,cluster,by=c("id"="ID"))


length=length(levels(factor(data$cluster)))
diff=survdiff(Surv(futime, fustat) ~ cluster, data = data)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ cluster, data = data)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=data,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="CuFecluster",
		           legend.labs=levels(factor(data[,"cluster"])),
		           legend = c(0.8, 0.8),
		           font.legend=10,
		           xlab="Time(years)",
		           break.time.by = 5,
		           palette = bioCol,
		           #surv.median.line = "hv",
		           risk.table=T,
		           cumevents=F,
		           risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width=10,height=8)
print(surPlot)
dev.off()

