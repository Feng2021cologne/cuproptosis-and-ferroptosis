
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

rt=read.table("interGene.txt", header=F, sep="\t", check.names=F)     #??ȡ?????ļ?

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       

GO=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=0.05, qvalueCutoff=0.05, ont="all", readable=T)
GO2=as.data.frame(GO)
GO2=GO2[GO2$p.adjust<0.05,]

pdf(file="", width=10, height=15)
bar=barplot(GO, drop=TRUE,showCategory = showNum,  color="p.adjust")
print(bar)
dev.off()

showNum=20
if(nrow(GO)<showNum){
  showNum=nrow(GO)
}

kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=0.05, qvalueCutoff=0.05)
KEGG=as.data.frame(kk)

KEGG=KEGG[KEGG$p.adjust < 0.05,]

showNum=20
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="", width=7, height=9)
barplot(kk, drop=TRUE, showCategory=showNum, color="p.adjust")
dev.off()

