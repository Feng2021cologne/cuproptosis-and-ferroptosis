
library(dplyr)
library(tibble)
library(stringr)
library(data.table)
library(limma)
library(tidyr)
library(ggpubr)

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm+1) - log(sum(fpkm+1)) + log(1e6))
}

rt = fread("fpkm.tsv",header = T)
rt<-as.matrix(rt)  
rownames(rt)=rt[,1]    
exp<-rt[,2:ncol(rt)] 
dimnames<-list(rownames(exp),colnames(exp))   
exp<-matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)  
rt=avereps(exp)
rt = data.frame(2^rt - 1) 
rt <- apply(rt,2,fpkmToTpm)
rt2 = data.frame(rt)


ann = read.table("",header = T)
ann = ann[,c(1,2)]

data = inner_join(ann,rt2,by = c("id"="gene"))
data = data[,-1]

data$median=apply(data[,-1],1,median)
data=data[order(data$gene,data$median,decreasing = T),]
data=data[!duplicated(data$gene),]
rownames(data) <- NULL
data <- tibble::column_to_rownames(data,var = "gene")
data <- data[,-ncol(data)]

sam <- str_sub(colnames(data),start = 14L,end=15L)
con = data[,sam=="11"]
nN = ncol(con)
tumor = data[,sam != "11"]
nT = ncol(tumor)
data2 = cbind(con,tumor)
fwrite(data2,"LUAD TPM.txt", row.names = T, sep = "\t", quote = F)

gene = read.table("gene.txt")
data3 = data2[gene$V1,]
dim(data3)

exp=log2(data3)
exp=as.data.frame(t(exp))
exp=cbind(exp, Type=c(rep("Normal",nN),rep("tumor",nT)))
data=gather(exp,gene,value,-ncol(exp))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            ylab="Gene expression",
            xlab="",
            legend.title="Type",
            palette = c("blue", "red"),
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="", width=7, height=5)
print(p1)
dev.off()

res=data.frame()
for (i in unique(data$Gene)) {
  data2 = data[data$Gene == i,]
  pvalue = wilcox.test(Expression~Type,data = data2)$p.value
  res0 = data.frame(gene = i,
                   pvalue = pvalue)
  if(pvalue<=0.05){
    res = rbind.data.frame(res,res0)
  }
}
write.table(res,"",quote = F,sep = "\t",row.names = F)







