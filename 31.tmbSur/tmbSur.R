
library(survival)
library(survminer)
library(tibble)

score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)    #??ȡm6A???ֵķ????ļ?
tmb=read.table("TMB.txt", header=T, sep="\t", check.names=F, row.names=1)        #??ȡTMB?????ļ?
sur = read.table("Cup-FerroptosisScore.group.txt",header = T,row.names = NULL)
sur = sur[,c(2,3,4)]
sur = column_to_rownames(sur,var = "id")

sameSample=Reduce(intersect,list(row.names(tmb), row.names(score),rownames(sur))) 
tmb=tmb[sameSample,,drop=F]
score=score[sameSample,,drop=F]
sur = sur[sameSample,]
data= Reduce(cbind.data.frame,list(score, tmb,sur)) 
res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")

res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("MEscore"))
cutoff=as.numeric(res.cut$cutpoint[1])
scoreType=ifelse(data[,"MEscore"]<=cutoff, "L-CuFescore", "H-CuFescore")

mergeType=paste0(tmbType, "+", scoreType)

bioSurvival=function(surData=null, outFile=null,width=NULL,height=NULL){
	diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
	length=length(levels(factor(surData[,"group"])))
	pValue=1-pchisq(diff$chisq, df=length-1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
	bioCol=c("#FF0000","#0066FF","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
	bioCol=bioCol[1:length]
	surPlot=ggsurvplot(fit, 
			           data=surData,
			           conf.int=F,
			           pval=pValue,
			           pval.size=6,
			           legend.title="",
			           legend.labs=levels(factor(surData[,"group"])),
			           font.legend=10,
			           legend = c(0.8, 0.8),
			           xlab="Time(years)",
			           break.time.by = 5,
			           palette = bioCol,
			           risk.table=T,
			           cumevents=F,
			           risk.table.height=.25)
	pdf(file=outFile, onefile = FALSE, width=width, height=height)
	print(surPlot)
	dev.off()
}

data$group=tmbType
bioSurvival(surData=data, outFile="",width = 10,height = 8)

data$group=mergeType
bioSurvival(surData=data, outFile="",width = 14,height = 10)

