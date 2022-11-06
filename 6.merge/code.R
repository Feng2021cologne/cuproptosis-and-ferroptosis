
library(dplyr)
library(tibble)
library(data.table)
library(stringr)

tcga = data.frame(fread("TPM.txt",header = T))
tcga = tcga[,which(str_sub(colnames(tcga),start = 14L,end = 15L)!="11")]
tcga[,-1] = log2(tcga[,-1] + 1)
colnames(tcga)[1] = "row.names"

GSE12667 = read.table("../GSE12667/GSE12667.txt",header = T,row.names = NULL)
GSE43458 = read.table("../GSE43458/GSE43458.txt",header = T,row.names = NULL)
GSE62949 = read.table("../GSE62949/GSE62949.txt",header = T,row.names = NULL)
GSE68465 = read.table("../GSE68465/GSE68465.txt",header = T,row.names = NULL)
GSE115002 = read.table("../GSE115002/GSE115002.txt",header = T,row.names = NULL)
GSE116959 = read.table("../GSE116959/GSE116959.txt",header = T,row.names = NULL)

data = inner_join(tcga,GSE12667,by="row.names")
data = inner_join(data,GSE43458,by="row.names")
data = inner_join(data,GSE62949,by="row.names")
data = inner_join(data,GSE68465,by="row.names")
data = inner_join(data,GSE115002,by="row.names")
data = inner_join(data,GSE116959,by="row.names")
data = column_to_rownames(data,var = "row.names")


library(sva)
library(limma)
batchType=c(rep(1,ncol(tcga)-1),
            rep(2,ncol(GSE12667)-1),rep(3,ncol(GSE43458)-1),rep(4,ncol(GSE62949)-1),
            rep(5,ncol(GSE68465)-1),rep(6,ncol(GSE115002)-1),rep(7,ncol(GSE116959)-1))
outTab=ComBat(data, batchType)


