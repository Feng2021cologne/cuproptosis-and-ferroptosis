## cuprotosisgene survival
library(limma)
library(survival)
library(survminer)
library(data.table)
library(dplyr)
library(tibble)
library(stringr)

library(survival)
colnames(rt)[c(1,2)] = c("fustat","futime")

pFilter=1                                                    
rt$futime <- rt$futime/365
rt = data.frame(rt)
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

##forest
rt2 <- read.table("uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt2)
hr <- sprintf("%.3f",rt2$"HR")
hrLow  <- sprintf("%.3f",rt2$"HR.95L")
hrHigh <- sprintf("%.3f",rt2$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt2$pvalue<0.001, "<0.001", sprintf("%.3f", rt2$pvalue))

pdf(file="forest.pdf", width = 8,height = 6.5)
n <- nrow(rt2)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)


par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
#boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
boxcolor = ifelse(as.numeric(hr) > 1, '#2556a6', '#2556a6')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()


dir.create("kmplot")
for (i in colnames(data0)[3:ncol(data0)]) {
  data=data0[,c("futime", "fustat", i)]
  colnames(data)=c("futime", "fustat", "gene")
  res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
  res.cat=surv_categorize(res.cut)
  fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)
  diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
  pValue=1-pchisq(diff$chisq, df=1)
  if(pValue<0.05){
    if(pValue<0.001){
      pValue="p<0.001"
    }else{
      pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    surPlot=ggsurvplot(fit,
                       data=res.cat,
                       pval=pValue,
                       pval.size=6,
                       legend.title=i,
                       legend.labs=c("high","low"),
                       xlab="Time(years)",
                       ylab="Overall survival",
                       palette=c("#e0b759", "#84b5d8"),
                       break.time.by=5,
                       conf.int=F,
                       risk.table=TRUE,
                       risk.table.title="",
                       risk.table.height=.25)
    pdf(file=paste0("kmplot/sur.", i, ".pdf"),onefile = FALSE,
        width = 8,
        height =6)
    print(surPlot)
    dev.off()
  }
}


