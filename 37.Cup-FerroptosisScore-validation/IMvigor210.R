
library(IOBR)
library(IMvigor210CoreBiologies)   
library(readxl)
library(tibble)

data(cds)  
expMatrix <- counts(cds)
eff_length2 <- fData(cds)[,c("entrez_id","length","symbol")]
rownames(eff_length2) <- eff_length2$entrez_id
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix),rownames(eff_length2))
eff_length2 <- eff_length2[mm,]

x <- expMatrix/eff_length2$length
eset <- t(t(x)/colSums(x))*1e6
summary(duplicated(rownames(eset)))

eset <- anno_eset(eset = eset,
                  annotation = eff_length2,
                  symbol = "symbol",
                  probe = "entrez_id",
                  method = "mean")
if(max(eset)>100) eset <- log2(eset+1)

pdata <- pData(cds)
colnames(pdata) <- gsub(colnames(pdata),pattern = " ",replacement = "_")
pdata <- rownames_to_column(pdata[,c("Best_Confirmed_Overall_Response",
                                     "binaryResponse",
                                     "Neoantigen_burden_per_MB",
                                     "Immune_phenotype",
                                     "censOS","os")],var = "ID")
pdata<-pdata[!is.na(pdata$binaryResponse),]
rownames(pdata) = NULL
pdata = column_to_rownames(pdata,var = "ID")
sam = intersect(colnames(eset),rownames(pdata))
eset2 = eset[,sam]
pdata = pdata[sam,]

load("expcli_IMvigor210.Rdata",verbose = T)
rt0=read.table("uniSigGeneExp.txt", header=T, sep="\t", check.names=F, row.names=1)
rt = eset2
data = rt[rownames(rt) %in% colnames(rt0),]
data=t(data)

library(WGCNA)
score = moduleEigengenes(data, rep("IMvigor210score",ncol(data)))$eigengenes

library(survival)
library(survminer)
score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)
cli=pdata
colnames(cli)[c(5,6)]=c("fustat","futime")

sameSample=intersect(row.names(score), row.names(cli))
data=cbind(cli[sameSample,], `Cup-FerroptosisScore` = score[sameSample,])

res.cut=surv_cutpoint(data, time="futime", event="fustat", variables=c("Cup-FerroptosisScore"))
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(data[,"Cup-FerroptosisScore"]<=cutoff, "Low", "High")
data$group=Type

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

bioCol=c("#d9ab3a","#549ace")
surPlot=ggsurvplot(fit, 
                   data=data,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="CuFescore",
                   legend.labs=levels(factor(data[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=12,
                   xlab="Time(months)",
                   break.time.by = 5,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)

pdf(file="", onefile = FALSE, width=8, height=8)
print(surPlot)
dev.off()

#################################
library(ggplot2)
library(cowplot)
library(ggpubr)

rt = read.table("Cup-FerroptosisScore.group.txt",header = T)

ggplot(rt, aes(x = group)) +
  geom_bar(aes(fill = response),
           position = 'fill',
           color = 'black',
           width = .7) +
  scale_fill_manual(values = c("#ea4c1c","#49b4bc")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") +
  ylab("Relative Percent(x100%)") +
  theme_half_open() 
ggsave("",width = 6,height = 6)

my_comparisons <- list( c("CR", "PD"), c("CR", "PR"), c("CR", "SD"),c("PD","PR"),c("PD","SD"),c("PR","SD"))
ggplot(rt,aes(Best_Confirmed_Overall_Response,Cup.FerroptosisScore)) +
  geom_boxplot(aes(color=Best_Confirmed_Overall_Response)) +
  geom_jitter(aes(color=Best_Confirmed_Overall_Response)) +
  scale_color_manual(values = c("#82cbd4","#e1bf6f","#d9572b","#a3cc5a")) +
  stat_compare_means(label.y = 0.4) +
  stat_compare_means(comparisons = my_comparisons) +
  ylab("CuFescore") +
  theme_half_open()  +
  theme(legend.position = 'none') 
ggsave("",width = 6,height = 6)

tpm = eset2
gene = data.frame(t(tpm[rownames(tpm) == "CD274",]))

gene = gene[rownames(data),,drop=F]

data2 = cbind(gene,data)
ggplot(data2,aes(group,CD274)) +
  geom_boxplot(aes(fill=group)) +
  scale_fill_manual(values = c("#ea4c1c","#49b4bc")) +
  stat_compare_means() +
  ylab("CD274 Expression") +
  theme_half_open() 
ggsave("",width = 6,height = 6)

data2 = data2[!is.na(data2$Immune_phenotype),]
my_comparisons = list(c("desert","excluded"),c("desert","inflamed"),c("excluded","inflamed"))
ggplot(data2,aes(Immune_phenotype,`Cup-FerroptosisScore`)) +
  geom_boxplot(aes(color=Immune_phenotype)) +
  geom_jitter(aes(color=Immune_phenotype)) +
  scale_color_manual(values = c("#82cbd4","#e1bf6f","#d9572b")) +
  stat_compare_means(label.y = 0.25) +
  stat_compare_means(comparisons = my_comparisons) +
  xlab("Immune Phenotype") +
  ylab("CuFescore") +
  theme_half_open()  +
  theme(legend.position = 'none') 
ggsave("",width = 6,height = 6)

data3 = data2[!is.na(data2$Neoantigen_burden_per_MB),]
data3 = mutate(data3,group = ifelse(Neoantigen_burden_per_MB > median(Neoantigen_burden_per_MB) & `Cup-FerroptosisScore` > median(`Cup-FerroptosisScore`),"CuFescore-H+NEO-H",
                                    ifelse(Neoantigen_burden_per_MB < median(Neoantigen_burden_per_MB) & `Cup-FerroptosisScore` > median(`Cup-FerroptosisScore`),"CuFescore-L+NEO-H",
                                           ifelse(Neoantigen_burden_per_MB > median(Neoantigen_burden_per_MB) & `Cup-FerroptosisScore` < median(`Cup-FerroptosisScore`),"CuFescore-H+NEO-L","CuFescore-L+NEO-L"))))

length=length(levels(factor(data3$group)))
diff=survdiff(Surv(futime, fustat) ~ group, data = data3)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data3)

bioCol=c("#7aabcc","#e6c769","#52566b","#d94e40")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.title="CuFescore",
                   legend.labs=levels(factor(data3[,"group"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(months)",
                   break.time.by = 5,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=F,
                   cumevents=F,
                   risk.table.height=.25)
pdf(file="",onefile = FALSE,width=8,height=6)
print(surPlot)
dev.off()









