##Cup-FerroptosisScore
library(WGCNA)
data=read.table("", header=T, sep="\t", check.names=F, row.names=1)
score = moduleEigengenes(data, rep("score",ncol(data)))$eigengenes
write.table(score, file="Cup-FerroptosisScore.txt", sep="\t", quote=F,row.names = T)

