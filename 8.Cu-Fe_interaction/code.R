library(data.table)
library(tibble)
library(stringr)
cu = read.table("gene.txt")

fe1 = read.csv("driver.csv",header = T)
fe2 = read.csv("suppressor.csv",header = T)
fe3 = read.csv("marker.csv",header = T)
fe = unique(c(fe1$symbol,fe2$symbol,fe3$symbol))

rt =  fread("TPM.txt",header = T)
rt = column_to_rownames(rt,var = "V1")
rt = rt[,str_sub(colnames(rt),start = 14L,end = 15L) != 11]
rt = log2(rt + 1)
rt_cu = rt[rownames(rt) %in% cu$V1,]
dim(rt_cu) 
rt_fe = rt[rownames(rt) %in% fe,]
dim(rt_fe) 
rt2 = rbind.data.frame(rt_cu,rt_fe)
rt2 = data.frame(t(rt2))

data = data.frame(cu = "",
                   fe = "",
                   cor = "",
                   pvalue = "")
for (i in 1:nrow(rt_cu)) {
  for (j in (nrow(rt_cu)+1):ncol(rt2)) {
    corr = cor.test(rt2[,i],rt2[,j])$estimate
    pval = cor.test(rt2[,i],rt2[,j])$p.value
    data0 = data.frame(cu = colnames(rt2)[i],
                      fe = colnames(rt2)[j],
                      cor = corr,
                      pvalue = pval)
    data = rbind.data.frame(data,data0)
  }
}

data2 = data[-1,]
data2$cor = as.numeric(data2$cor)
data2$pvalue = as.numeric(data2$pvalue)
write.csv(data2,"",quote = F,row.names = F)
data2 = data2[abs(data2$cor)> 0.5 & data2$pvalue<0.05,]
dim(data2) 
length(unique(data2$cu))
length(unique(data2$fe))

x = read.table("genemania-interactions.txt",header = T,sep = "\t") 
x = x[x$Gene.1 %in% c(unique(data2$cu),unique(data2$fe)),]
x = x[x$Gene.2 %in% c(unique(data2$cu),unique(data2$fe)),]

data3 = cbind.data.frame(rt2[,unique(data2$cu)],rt2[,unique(data2$fe)])
library(corrplot)
data4 = cor(data3)
testRes = cor.mtest(data3, conf.level = 0.95)
pdf("",width = 13,height = 13,onefile = F)
cor.plot<-corrplot(corr =data4,type="upper",tl.pos="tp",
                   tl.col="black",p.mat = testRes$p, insig = "label_sig", 
                   sig.level = c(.01, .05),
                   pch.cex=1,pch.col = "black",order = "original",
                   col=rev(COL2('RdBu', 200)))
cor.plot<-corrplot(corr = data4,type="lower",add=TRUE,method="number",
                   tl.pos="n",tl.col="black",
                   col="black",tl.cex=1.2,diag=FALSE, cl.pos="n",pch.col = "black",
                   number.cex = 0.8,order = "original")
dev.off()

