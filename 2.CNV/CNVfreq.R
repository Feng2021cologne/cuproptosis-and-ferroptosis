library(data.table)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggpubr)

rt = data.frame(fread("",header = T))
rt = rt[,!duplicated(colnames(rt))]
ann = fread("gencode.v22.annotation.gene.probeMap",header = T)
ann = ann[,c(1,2)]
rt2 = inner_join(ann,rt,by=c("id"="Gene.Symbol"))
rt2 = rt2[,-1]
rt2 = data.frame(rt2)
gene = read.table("gene.txt")
rt3 = rt2[rt2$gene %in% gene$V1,]
rownames(rt3) = NULL
fwrite(rt3,"cnvMatrix.txt",sep = "\t",row.names = F,quote = F)

rt=column_to_rownames(rt3,var = "gene")
GAIN=rowSums(rt> 0)       
LOSS=rowSums(rt< 0)     
GAIN=GAIN/ncol(rt)*100     
LOSS=LOSS/ncol(rt)*100     
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

data.max = apply(data, 1, max)
pdf(file="", width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col=2, cex=3)
points(bar,data[,"LOSS"], pch=20, col=3, cex=3)
legend("top", legend=c('GAIN','LOSS'), col=2:3, pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1)
dev.off()

data = column_to_rownames(rt3,var = "gene")
colnames(data) = str_sub(colnames(data),start = 1L,end = 16L)

tpm = data.frame(fread("LUAD TPM.txt",header = T))
tpm = column_to_rownames(tpm,var = "V1")
tpm = log2(tpm + 1)
normal = colnames(tpm)[str_sub(colnames(tpm),start = 14L,end = 15L)=="11"]

for (i in 1:nrow(data)) {
  cnv = data[i,]
  CNV_gain = colnames(cnv)[as.character(cnv) == "1"]
  CNV_loss = colnames(cnv)[as.character(cnv) == "-1"]
  none_CNV = colnames(cnv)[as.character(cnv) == "0"]
  sam = data.frame(id = c(CNV_gain,CNV_loss,none_CNV,normal),
                   group = c(rep("CNV_gain",length(CNV_gain)),
                             rep("CNV_loss",length(CNV_loss)),
                             rep("none_CNV",length(none_CNV)),
                             rep("normal",length(normal))))
  
  tpm2 = data.frame(t(tpm[rownames(cnv)[1],]))
  tpm2 = rownames_to_column(tpm2,var = "id")
  tpm2 = inner_join(tpm2,sam,by="id")
  my_comparisons <- list( c("CNV_gain", "CNV_loss"), 
                          c("CNV_gain", "none_CNV"), 
                          c("CNV_gain", "normal"),
                          c("CNV_loss","none_CNV"),
                          c("CNV_loss","normal"),
                          c("none_CNV","normal"))
  ggplot(tpm2,aes_string("group",colnames(tpm2)[2])) +
    geom_boxplot(aes(fill=group)) +
    scale_fill_manual(values = c("#346388","#8e4a97","#489746","#d3221f")) +
    stat_compare_means(aes(label = ..p.signif..)) +
    labs(x="",y=paste("log2(TPM) of",colnames(tpm2)[2],sep = " ")) +
    theme_classic()
  ggsave(paste(colnames(tpm2)[2],".pdf",sep = ""),width = 5,height = 4)
}









