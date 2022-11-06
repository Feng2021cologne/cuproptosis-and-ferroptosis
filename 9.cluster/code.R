##ConsensusClusterPlus
library(data.table)
library(tibble)
library(dplyr)
library(stringr)
library(survival)
library(WGCNA)

rt = data.frame(fread("../6.merge/prognosis_data/merge/corrected_data.txt"))
rt = rt[rt$V1 %in% c(cu$V1,fe),]
rownames(rt)  = NULL
rt = column_to_rownames(rt,var = "V1")
rt = data.frame(t(rt))
rt = rownames_to_column(rt,var = "id")
rt$id = gsub("\\.","-",rt$id)

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

rt = inner_join(sur,rt,by="id")
rt = column_to_rownames(rt,var = "id")

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<0.1){
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
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

selected_cu = outTab$id[which(outTab$id %in% cu$V1)]
selected_fe = outTab$id[which(outTab$id %in% fe)]


fe1 = read.csv("../8.Cu-Fe_interaction/driver.csv",header = T)
fe2 = read.csv("../8.Cu-Fe_interaction/suppressor.csv",header = T)
fe3 = read.csv("../8.Cu-Fe_interaction/marker.csv",header = T)
cu_exp = rt[,selected_cu]

fe_driver_exp = rt[,colnames(rt) %in% fe1$symbol]
fe_suppressor_exp = rt[,colnames(rt) %in% fe2$symbol]
fe_marker_exp = rt[,colnames(rt) %in% fe3$symbol]

ME_cu = moduleEigengenes(cu_exp, rep("cu",ncol(cu_exp)))$eigengenes
ME_fe_driver = moduleEigengenes(fe_driver_exp, rep("fe_driver",ncol(fe_driver_exp)))$eigengenes
ME_fe_suppressor = moduleEigengenes(fe_suppressor_exp, rep("fe_suppressor",ncol(fe_suppressor_exp)))$eigengenes
ME_fe_marker = moduleEigengenes(fe_marker_exp, rep("fe_marker",ncol(fe_marker_exp)))$eigengenes

ME = Reduce(cbind.data.frame,list(ME_cu,ME_fe_driver,ME_fe_suppressor,ME_fe_marker))
data = t(ME)
library(ConsensusClusterPlus)
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             title=".",
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123456,
                             plot="pdf")
icl <- calcICL(results,plot = 'pdf')
clusterNum=3
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
letter=c("A","B","C","D")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="Cluster.txt", sep="\t", quote=F, col.names=F)





