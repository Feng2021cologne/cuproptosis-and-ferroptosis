##unicox
library(limma)
library(survival)
library(data.table)
library(tibble)
library(dplyr)
data=read.table("../16.clusterDiff/interGeneExp.txt", header=T, sep="\t", row.names = 1)
merged = data.frame(fread("../6.merge/prognosis_data/merge/corrected_data.txt",header = T))
merged = column_to_rownames(merged,var = "V1")
data = merged[rownames(merged) %in% rownames(data),]
data=data.frame(t(data))

rownames(data)=gsub("\\.", "-", rownames(data))
data = rownames_to_column(data,var = "id")

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

data = inner_join(sur,data,by="id")
rt = column_to_rownames(data,var = "id")

outTab=data.frame()
sigGenes=c()
for(i in colnames(rt[,3:ncol(rt)])){
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<=0.1){
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
sigGeneExp=rt[,outTab$id]
sigGeneExp=rbind(id=colnames(sigGeneExp), sigGeneExp)
write.table(sigGeneExp, file="uniSigGeneExp.txt", sep="\t", quote=F, col.names=F)


