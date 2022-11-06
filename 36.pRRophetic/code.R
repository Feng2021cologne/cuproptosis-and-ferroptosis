library(pRRophetic)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggpubr)
library(dplyr)
library(tibble)
library(data.table)


jco <- c("#EABF00", "#2874C5", "red")、
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表
plotp <- list()

for (drug in GCP.drug) {
  set.seed(1248103) 
  cat(drug," starts!\n") 
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(dat[,rownames(ann)]),
                                              drug = drug,
                                              tissueType = "allSolidTumors",
                                              selection = 1) 
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")}
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "ImmClust"=ann$group, 
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust) # 把类改成因子变量

}

p <- vector()
for (drug in GCP.drug) {
  tmp <- kruskal.test(est.ic50 ~ ImmClust,
                     data = predictedBoxdat[[drug]])$p.value
  p <- append(p,tmp) 
}
names(p) <- GCP.drug

p2 = data.frame(p)
p2$FDR = p.adjust(p2$p)
write.table(p2,"output_pvalue.txt", quote = F, sep = "\t")

sig = p2[p2$p <= 0.05,,drop=F]
dim(sig)
write.table(sig,"output_pvalue_sig.txt", quote = F, sep = "\t")

data = data.frame()
for (i in 1:length(predictedBoxdat)) {
 lowmed =  median(predictedBoxdat[[i]]$est.ic50[predictedBoxdat[[i]]$ImmClust == "Low"])
 highmed =  median(predictedBoxdat[[i]]$est.ic50[predictedBoxdat[[i]]$ImmClust == "High"])
 data0 = data.frame(drug = names(predictedBoxdat[i]),
                   lowmed = lowmed,
                   highmed = highmed)
 data = rbind.data.frame(data,data0)
}

p3 = rownames_to_column(p2,var = "drug")
p3  = mutate(p3,Statistical_test = ifelse(p < 0.05 & FDR < 0.05,"P < 0.05 and FDR < 0.05",
                                          ifelse(p < 0.05 & FDR > 0.05,"P < 0.05 and FDR > 0.05","P > 0.05 and FDR > 0.05")))

data2  = inner_join(p3[,-c(2,3)],data,by="drug")
data2 = mutate(data2,value = (exp(highmed) / exp(lowmed)) - 1) %>%
  arrange(desc(value))
data2$drug = factor(data2$drug,levels = data2$drug)
ggplot(data2,aes(drug,value)) +
  geom_bar(stat="identity",aes(fill = Statistical_test)) +
  scale_fill_manual(values = c("#61adad","#f4d240","#caccce")) +
  ylab("Exp(Median estimated IC50 (H))\n —————————————————— - 1\n Exp(Median estimated IC50 (L))") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 5),
        panel.grid = element_blank(),
        legend.position = c(0.8,0.8))
ggsave("drug-barplot.pdf",width = 14,height = 6)
