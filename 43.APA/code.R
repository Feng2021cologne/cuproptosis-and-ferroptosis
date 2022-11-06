library(tibble)
library(stringr)

rt = read.table("PDUIs.txt",header = T,check.names = T)  
rt = rt[,-c(2,3)]
rt$event_id = sapply(strsplit(rt$event_id,"\\|"),"[",2)

rownames(rt) = NULL
rt = column_to_rownames(rt,var = "event_id")
rt = rt[,!duplicated(str_sub(colnames(rt),start = 1L,end = 16L))] 
colnames(rt) = str_sub(colnames(rt),start = 1L,end = 16L)
colnames(rt) = gsub("\\.","-",colnames(rt))

data = rt
data = data[rowSums(is.na(data)) < ncol(data) / 2, colSums(is.na(data)) < nrow(data) / 2]

rt<-as.matrix(data) 
library("impute")
library(limma)
mat=impute.knn(rt)
rt=mat$data     
rt=avereps(rt)

group = read.table("Cup-FerroptosisScore.txt",header = T)
group = group[grep("TCGA",rownames(group)),,drop=F]
group$group = ifelse(group$MEscore >= median(group$MEscore),"High","Low")
low = rt[,colnames(rt) %in% rownames(group)[group$group =="Low"]]
high = rt[,colnames(rt) %in% rownames(group)[group$group=="High"]]

data = cbind(low,high)
pudi = data
######################################
exp <- data

class <- c(rep("normal",ncol(low)),rep("tumor",ncol(high)))
class<- factor(class,levels=c("normal","tumor"),labels=c("normal","tumor"))
design <- model.matrix(~0+class)
colnames(design) <- c("normal","tumor")   
fit <- lmFit(exp,design)  
cont.matrix<-makeContrasts(tumor - normal,levels=design)
fit1 <- contrasts.fit(fit, cont.matrix)   
fit1 <- eBayes(fit1) 
allDiff<-topTable(fit1,adjust='fdr',number=100000)

diff = allDiff[allDiff$adj.P.Val < 0.05,]
diff = allDiff[abs(allDiff$logFC) > 0.05 & allDiff$adj.P.Val < 0.05,]

library(data.table)
library(dplyr)
library(survival)
library(survminer)
rt = pudi
rt = rt[rownames(rt) %in% rownames(diff),]
rt =data.frame(t(rt))
rt = rownames_to_column(rt,var = "id")

cli=fread("survival.tsv", header=T, sep="\t", check.names=F) 
cli=cli[,-3]
colnames(cli)=c("id","fustat","futime")
cli$futime=cli$futime/365
rt = data.frame(inner_join(cli,rt,by="id"))

dir.create("kmplot")
data0 = rt
res = data.frame()
for (i in colnames(data0)[4:ncol(data0)]) {
  data=data0[,c("futime", "fustat", i)]
  colnames(data)=c("futime", "fustat", "gene")

  res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
  res.cat=surv_categorize(res.cut)
  fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)

  diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
  pValue=1-pchisq(diff$chisq, df=1)
  res0 = data.frame(gene = i,
                    pvalue = pValue)
  res = rbind.data.frame(res,res0)

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
                       legend.labs=c("Lengthening","Shortening"),
                       xlab="Time(years)",
                       ylab="Overall survival",
                       palette=c("#fb0608", "#89cef0"),
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
res2 = res[res$pvalue <= 0.05,]

data1 = data.frame()
for (i in colnames(data0)[4:ncol(data0)]) {
  data=data0[,c("futime", "fustat", i)]
  colnames(data)=c("futime", "fustat", "gene")

  res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
  rownames(res.cut$cutpoint) = i
  data1 = rbind.data.frame(data1,res.cut$cutpoint)
}

#############################
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
data <- dplyr::mutate(allDiff,direction = ifelse(logFC >= 0.05 &  adj.P.Val <= 0.05,"up",
                                                 ifelse(logFC <= -0.05 & adj.P.Val <= 0.05,"down","none")))
data21 = data[data$direction == "down",]
data22 = data[data$direction == "up",]

data2 = rbind(data22[1:10,],data21[1:10,])

data2$gene = rownames(data2)
ggplot(data = data, aes(x = logFC, 
                        y = -log10(adj.P.Val))) +
  geom_point(size = 4,
             aes(color = direction),
             show.legend = F) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dotdash', 
             color = 'grey30') + 
  geom_vline(xintercept = c(-0.05, 0.05),
             linetype = 'dotdash', color = 'grey30') +
  scale_color_manual(
    values = c('#1500FF', '#A9A9A9', '#FF0102')) +
  labs(x = 'Log2(fold change)', 
       y = '-log10(adj.P.Val)') +
  theme_half_open() +
  theme(plot.title = element_text(hjust = 0.5)) 

rt = data0
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,4:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.01){
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

library(patchwork)
diff2 = allDiff[abs(allDiff$logFC) > 0.05 & allDiff$adj.P.Val < 0.05,]
diff2  = diff2[outTab$id,]  
diff2$gene = rownames(diff2)
diff2 = mutate(diff2,sig=ifelse(adj.P.Val<0.05 & adj.P.Val >0.01,"*",
                                ifelse(adj.P.Val<0.01 & adj.P.Val >0.001,"**","***")))
ggplot(diff2,aes(logFC,gene,fill = adj.P.Val)) +
  geom_bar(stat ="identity") +
  scale_fill_gradient(low ="#9817e7" ,high = "#ff2936") +
  scale_y_discrete(position = "right") +
  scale_x_reverse(position = "top") +
  labs(x="APA_PDUI difference",y="") +
  theme_classic() +
  theme(legend.position = "left") -> p1

rt <- read.table("uniCox.txt",header=T,sep="\t",row.names = 1, check.names=F)
rt$HR <- sprintf("%.3f",rt$"HR")
rt$HR.95L <- sprintf("%.3f",rt$"HR.95L")
rt$HR.95H <- sprintf("%.3f",rt$"HR.95H")
rt$`HR of APA_PDUI` =  paste0(rt$HR,"(",rt$HR.95L,"-",rt$HR.95H,")")
rt$pvalue <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
rt = mutate(rt,sig= ifelse(pvalue <= 0.05 & pvalue > 0.01,'*',
                           ifelse(pvalue <= 0.01 & pvalue > 0.001,'**',"***")))
rt$HR = as.numeric(rt$HR)
rt$HR.95L = as.numeric(rt$HR.95L)
rt$HR.95H = as.numeric(rt$HR.95H)
rt = rt[rownames(diff2),]
rt = rownames_to_column(rt,var = "id")
rt$num  = seq(1,12,1)
ggplot(rt,aes(HR,id)) +
  geom_point(shape=15,size = 2.5) +
  geom_segment(aes(x=HR.95L,xend=HR.95H,y=id,yend=id)) +
  geom_segment(aes(x=HR.95L,xend=HR.95L,y=num-0.07,yend=num+0.07)) +
  geom_segment(aes(x=HR.95H,xend=HR.95H,y=num-0.07,yend=num+0.07)) +
  scale_x_continuous(limits = c(-3,5),breaks = seq(-3,3,2),labels = seq(-3,3,2)) +
  geom_vline(aes(xintercept=1),linetype="longdash") +
  geom_text(aes(x=-2,label=`HR of APA_PDUI`)) +
  geom_text(aes(x=3,label=pvalue)) +
  geom_text(aes(x=4,label=sig)) +
  geom_text(aes(x=6,y=12.5,label="p value")) +
  geom_text(aes(x=-2,y=12.5,label="HR of APA_PDUI")) +
  labs(y="") +
  theme_classic() -> p2

p1 + p2
ggsave("",width = 14,height = 8)
dev.off()

dir.create("kmplot2")
data00 = cbind(data0[,1:3],data0[colnames(data0) %in% rt$id])
res = data.frame()
for (i in colnames(data00)[4:ncol(data00)]) {
  data=data00[,c("futime", "fustat", i)]
  colnames(data)=c("futime", "fustat", "gene")

  res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("gene"))
  res.cat=surv_categorize(res.cut)
  fit=survfit(Surv(futime, fustat) ~gene, data = res.cat)

  diff=survdiff(Surv(futime, fustat) ~gene,data =res.cat)
  pValue=1-pchisq(diff$chisq, df=1)
  res0 = data.frame(gene = i,
                    pvalue = pValue)
  res = rbind.data.frame(res,res0)

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
                       legend.labs=c("Lengthening","Shortening"),
                       xlab="Time(years)",
                       ylab="Overall survival",
                       palette=c("#fb0608", "#89cef0"),
                       break.time.by=5,
                       conf.int=F,
                       risk.table=TRUE,
                       risk.table.title="",
                       risk.table.height=.25)
    pdf(file=paste0("kmplot2/sur.", i, ".pdf"),onefile = FALSE,
        width = 8,
        height =6)
    print(surPlot)
    dev.off()
  }
}

