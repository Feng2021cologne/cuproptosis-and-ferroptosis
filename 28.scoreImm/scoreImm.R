
library(corrplot)

score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=1)
rownames(score) = gsub("\\.","-",rownames(score))
immune=read.table("ssGSEA.result.txt", header=T, sep="\t", check.names=F, row.names=1)
immune=t(immune)

sameSample=intersect(row.names(score), row.names(immune))
data=cbind(score[sameSample,,drop=F], immune[sameSample,,drop=F])
colnames(data)[1] = "CuFescore"
colnames(data) = gsub("\\."," ",colnames(data))

M=cor(data)
res1=cor.mtest(data, conf.level = 0.95)

pdf(file="", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()

