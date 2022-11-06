library(data.table)
library(RCircos)  
library(tibble)
library(dplyr)
ann = fread("",header = T)
ann = ann[,c(2,3,4,5)]
colnames(ann) = c("Gene","Chr","Start","End")

gene = read.table("gene.txt")
ann2 = ann[ann$Gene %in% gene$V1,]

RCircos.Gene.Label.Data = data.frame(Chromosome=ann$Chr,
                                     chromStart =ann$Start,
                                     chromEnd = ann$End,
                                     Gene = ann$Gene)
RCircos.Gene.Label.Data = RCircos.Gene.Label.Data[RCircos.Gene.Label.Data$Gene %in% gene$V1,]

RCircos.Scatter.Data = data.frame(Chromosome=ann2$Chr,
                                  start=ann2$Start,
                                  stop=ann2$End,
                                  gene=ann2$Gene)

rt = fread("cnvMatrix.txt",header = T)
rt = column_to_rownames(rt,var = "gene")
data = data.frame()
for (i in RCircos.Scatter.Data$gene) {
  loss = length(which(as.numeric(rt[i,])== -1))
  gain= length(which(as.numeric(rt[i,])== 1))
  data0 = data.frame(gene=i,
                    seg.mean = ifelse(gain>loss,1,-1))
  data = rbind.data.frame(data,data0)
}

RCircos.Scatter.Data = inner_join(RCircos.Scatter.Data,data,by="gene")
RCircos.Scatter.Data = RCircos.Scatter.Data[,-4]

cytoBandIdeogram=read.table("refer.txt", header=T, sep="\t") 
chr.exclude <- NULL
cyto.info <- cytoBandIdeogram
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

rcircos.params$text.size=1
rcircos.params$point.size=5
RCircos.Reset.Plot.Parameters(rcircos.params)

pdf(file="", width=8, height=8)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

data.col <- 4
track.num <- 1
side <- "in"
RCircos.Scatter.Plot(RCircos.Scatter.Data, data.col, track.num, side, by.fold=0.1)

name.col <- 4
side <- "in"
track.num <- 2
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num, side)
track.num <- 3
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col, track.num, side)
dev.off()


