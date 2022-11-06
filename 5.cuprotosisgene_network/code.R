
library(reshape2)
library(corrplot)
library(plyr)
library(igraph) 
library(data.table)
library(tibble)
library(stringr)
library(dplyr)

poscol <- "#FB9A99" 
negcol <- "#C6DBEF" 

mycol <- c("#FDBF6F", "#1F78B4", "#E31A1C", "#8C510A") 

rt = fread("TPM.txt",header = T)
rt = column_to_rownames(rt,var = "V1")
cu = read.table("gene.txt")
input_data = rt[cu$V1,str_sub(colnames(rt),start = 14L,end = 15L) !=11]
input_data = log2(input_data)
input_data = data.frame(t(input_data))
dim(input_data)
input_data[1:3,1:3]

input_data0 = data.frame(input_data)
input_data0 = rownames_to_column(input_data0,var = "id")
input_data0$id = gsub("\\.","-",input_data0$id)
sur = fread("survival.tsv",header = T)
input_data0 = inner_join(sur,input_data0,by=c("sample"="id"))
input_data0 = input_data0[,-c(1,3)]

library(survival)
colnames(input_data0)[c(1,2)] = c("fustat","futime")
pFilter=1                                                    
input_data0$futime <- input_data0$futime/365
input_data0 = data.frame(input_data0)
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(input_data0[,3:ncol(input_data0)])){
  cox <- coxph(Surv(futime, fustat) ~ input_data0[,i], data = input_data0)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
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

bb <- outTab
head(bb)
colnames(bb)[1] <- c("ID")

bb$weight <- abs(log10(as.numeric(bb$pvalue)))
bb$weight_HR <- (as.numeric(bb$HR)-1)*100
bb$colr <- ifelse(bb$weight_HR<0, "green", "black")
head(bb)

corr <- cor(input_data, method = "spearman")
corrplot(corr,title = "", 
         method = "pie", 
         outline = T, addgrid.col = "darkgray", 
         order="hclust", addrect = 4, 
         mar = c(4,0,4,0), 
         rect.col = "black", rect.lwd = 5, cl.pos = "b", 
         tl.col = "black", tl.cex = 1.08, cl.cex = 1.5, tl.srt=60)

p.corr <- cor.mtest(input_data)
head(p.corr[, 1:5])

rr <- as.data.frame(corr);
rr$ID <- rownames(rr)
cor <- melt(rr,"ID",value.name = "cor"); 

pp <- as.data.frame(p.corr);
pp$ID <- rownames(pp)
pvalue <- melt(pp,"ID",value.name = "pvalue");
colnames(pvalue) <- c("from","to","pvalue")

corpvlue <- cbind(pvalue, cor)
head(corpvlue)
corpvlue <- corpvlue[, -c(4:5)]
head(corpvlue)
dim(corpvlue)


corpvlue <- corpvlue[corpvlue$pvalue < 0.05,] 
dim(corpvlue)
corpvlue$weight <- corpvlue$pvalue
corpvlue$weight <- -log10(corpvlue$weight)
head(corpvlue)

corpvlue <- corpvlue[!corpvlue$cor==1,]
summary(duplicated(corpvlue$weight))
corpvlue <- corpvlue[!duplicated(corpvlue$weight),]

corpvlue$color <- ifelse(corpvlue$cor<0, negcol, poscol)

cellcluster <- as.data.frame(t(input_data))

hc <- hclust(dist((cellcluster)))
hcd <- as.dendrogram(hc)
(clus4 <- cutree(hc, 2)) 

A <- as.character(rownames(as.data.frame(subset(clus4,clus4==1))))
B <- as.character(rownames(as.data.frame(subset(clus4,clus4==2))))
cls <- list(A,B)

nodes <- as.data.frame(unlist(cls))
nodes$type <- c(rep("A",7),rep("B",7))
names(nodes) <- c("media","type.label")

nodes <- as.data.frame(nodes)
nodes$media <- as.character(nodes$media)
nodes

summary(nodes$media %in% bb$ID) 
nodes <- merge(nodes, bb, by.x = "media", "ID", all.x = T, all.y = T) 

nodes$Fraction <- abs(nodes$weight_HR)
nodes <- nodes[order(nodes$type.label),]
nodes <- nodes[,c(ncol(nodes),1:ncol(nodes)-1)]
nodes <- nodes[order(nodes$type.label),]
nodes

paste0("'",nodes$media,"'","=","'",nodes$id,"'",collapse = ",")
corpvlue$from <- revalue(corpvlue$from,c('ATP7A'='S1','ATP7B'='S2','DBT'='S3',
                                         'DLAT'='S4','DLD'='S5','DLST'='S6','FDX1'='S7',
                                         'GCSH'='S8','LIAS'='S9','LIPT1'='S10',
                                         'LIPT2'='S11','PDHA1'='S12','PDHB'='S13',
                                         'SLC31A1'='S14'))
corpvlue$to <- revalue(corpvlue$to,c('ATP7A'='S1','ATP7B'='S2','DBT'='S3',
                                     'DLAT'='S4','DLD'='S5','DLST'='S6','FDX1'='S7',
                                     'GCSH'='S8','LIAS'='S9','LIPT1'='S10',
                                     'LIPT2'='S11','PDHA1'='S12','PDHB'='S13',
                                     'SLC31A1'='S14'))
(links <- corpvlue)

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
V(net)$color <- revalue(nodes$type.label,c("A"=mycol[2],"B"=mycol[3]))
V(net)$size <- (1 + V(net)$weight)*7 
V(net)$label <- V(net)$media 
E(net)$arrow.mode <- 0 
E(net)$edge.color <- "tomato" 
E(net)$width <- 1+E(net)$weight/6  

pdf("", width = 9.75, height = 8.78 )
plot(net,
     layout=layout_in_circle, 
     edge.curved=.2, 
     vertex.label.color=V(net)$color, 
     vertex.label.dist= -2, 
     edge.color=links$color)

legend("topright", 
       c("Gene cluster-A", "Gene cluster-B"),
       pch=21, col="black", pt.bg=mycol, pt.cex=3,
       cex=1.3, bty="n", ncol=1)

f <- c(0.05, 0.001, 0.00001, 0.00000001)
s <- sqrt(abs(log10(f)))*3
legend("bottomright", 
       inset=c(0,-.1), 
       legend=f, text.width = .2, 
       title = "P value", title.adj = -.3,
       pch=21, pt.cex=s, bty='n',
       horiz = TRUE, 
       col = "black",pt.bg=rep("grey",4))

legend("bottomright",
       c("Positive correlation with P < 0.05", 
         "Negative correlation with P < 0.05"),
       col = c(poscol, negcol), bty="n", 
       cex = 1, lty = 1, lwd = 5)

dev.off()


