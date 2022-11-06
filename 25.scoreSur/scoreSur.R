##score survival
library(survival)
library(survminer)
library(data.table)
library(dplyr)

data = inner_join(sur,score,by=c("id"="row.names"))

res.cut=surv_cutpoint(data, time="futime", event="fustat", variables=c("MEscore"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,"MEscore"]<=cutoff, "Low", "High")
data$group=Type
outTab=rbind(id=colnames(data), data)
write.table(outTab, file="Cup-FerroptosisScore.group.txt", sep="\t", quote=F, col.names=F)

data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data)

bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
			       data=data,
			       conf.int=F,
			       pval=pValue,
			       pval.size=6,
			       legend.title="CuFeScore",
			       legend.labs=levels(factor(data[,"group"])),
			       legend = c(0.8, 0.8),
			       font.legend=12,
			       xlab="Time(years)",
			       break.time.by = 5,
			       palette = bioCol,
			       #surv.median.line = "hv",
			       risk.table=T,
			       cumevents=F,
			       risk.table.height=.25)

pdf(file="scoresurvival.pdf", onefile = FALSE, width=10, height=8)
print(surPlot)
dev.off()


