
library(survival)
library(survminer)
library(data.table)
library(tibble)

score=read.table("Cup-FerroptosisScore.group.txt", header=T, sep="\t",check.names=F, row.names=1)
score = score[grep("TCGA",score$id),]
score$group = ifelse(score$MEscore >= median(score$MEscore),"High","Low")
rownames(score) = NULL
score = column_to_rownames(score,var = "id")

cli=fread("phenotype.tsv", header=T)
cli = cli[,c(1,44)]
colnames(cli) = c("id","T")
cli = cli[cli$T != "TX" & cli$T != "",]
cli = mutate(cli,T = ifelse(T=="T1" | T=="T1a" | T=="T1b" | T=="T2" | T=="T2a" | T=="T2b","T1-2","T3-4"))
cli = column_to_rownames(cli,var = "id")

sameSample=intersect(rownames(cli), rownames(score))
score=score[sameSample,]
cli=cli[sameSample,,drop=F]
data=cbind(futime=score[,2], fustat=score[,1], cli, group=score[,"group"])
colnames(data)[3] = "clinical"
rt = data
tab=table(data[,"clinical"])

trait = "T"
for(j in names(tab)){
	rt1=rt[(rt[,"clinical"]==j),]
	tab1=table(rt1[,"group"])
	labels=names(tab1)
	diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
	pValue=1-pchisq(diff$chisq, df=1)
	if(pValue<0.001){
	  pValue="p<0.001"
	}else{
	  pValue=paste0("p=",sprintf("%.03f", pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
	bioCol=c("#FF0000","#0066FF","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
	bioCol=bioCol[1:length(levels(factor(rt1[,"group"])))]
	surPlot=ggsurvplot(fit, 
	                   data=rt1,
	                   conf.int=F,
	                   pval=pValue,
	                   pval.size=6,
	                   title=paste0("Patients with ",j),
	                   legend.title="CuFescore",
	                   legend.labs=labels,
	                   font.legend=12,
	                   xlab="Time(years)",
	                   break.time.by = 5,
	                   palette=bioCol,
	                   risk.table=TRUE,
	                   risk.table.title="",
	                   risk.table.col = "strata",
	                   risk.table.height=.25)
	pdf(file=paste0(trait,"_",j,".pdf"), onefile = FALSE, width = 10, height =8)
	print(surPlot)
	dev.off()
}

