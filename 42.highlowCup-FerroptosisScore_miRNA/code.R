
library(dplyr)
library(tibble)
library(stringr)
library(pacman)
p_load(TCGAbiolinks)
p_load(stringi)

project_id <- "TCGA-LUAD"
TCGAbiolinks:::getProjectSummary(project_id)
data_category <- "Transcriptome Profiling"

data_type <- "Isoform Expression Quantification"

query <- GDCquery(project = project_id,
                  data.category = data_category,
                  data.type = data_type)

GDCdownload(query)

dataAssy <- GDCprepare(query = query,
                       summarizedExperiment=F)
sample_id <- as.data.frame(table(dataAssy$barcode))
miRNA_id <- as.data.frame(table(dataAssy$miRNA_region))

miRNA_matrue_RPM <- matrix(NA,ncol = nrow(sample_id),nrow = nrow(miRNA_id))
colnames(miRNA_matrue_RPM) <- sample_id$Var1
rownames(miRNA_matrue_RPM) <- as.character(miRNA_id$Var1)
for(i in 1:nrow(sample_id)){
  temp1 <- dataAssy[which(dataAssy$barcode==as.character(sample_id[i,1])),]
  
  for(j in 1:nrow(miRNA_id)){
    loc <- which(temp1$miRNA_region==as.character(miRNA_id[j,1]))
    if(length(loc)>0){
      miRNA_matrue_RPM[j,i] <- sum(temp1[loc,4])
    }else{
      miRNA_matrue_RPM[j,i] <- 0
    }
    
  }
  
  print(i)
}

rownames(miRNA_matrue_RPM1) <- substr(rownames(miRNA_matrue_RPM1),8,nchar(rownames(miRNA_matrue_RPM1)))

p_load(miRBaseVersions.db)
items <- select(miRBaseVersions.db,
                keys = rownames(miRNA_matrue_RPM1),
                keytype = "MIMAT",
                columns = c("ACCESSION","NAME","VERSION"))
id_name <- items[items$VERSION == 21.0, c("ACCESSION","NAME")]

miRNA_matrue_RPM1 = data.frame(miRNA_matrue_RPM1)
miRNA_matrue_RPM2 = rownames_to_column(miRNA_matrue_RPM1,var = "id")
miRNA_matrue_RPM2 <- inner_join(id_name,miRNA_matrue_RPM2,by = c("ACCESSION"="id"))
miRNA_matrue_RPM2 = miRNA_matrue_RPM2[,-1]
save(miRNA_matrue_RPM2,file = paste(project_id,"_miRNA_matrue_RPM.RData",sep=""))

load(file = "TCGA-LUAD_miRNA_matrue_RPM.RData",verbose = T)
data = miRNA_matrue_RPM2
data = column_to_rownames(data,var = "NAME")
colnames(data) = str_sub(colnames(data),start = 1L,end = 16L)
colnames(data) = gsub("\\.","-",colnames(data))

sam = read.table('Cup-FerroptosisScore.txt',header = T,row.names = NULL)
sam = sam[sam$row.names %in% colnames(data),]
sam$group = ifelse(sam$MEscore >= median(sam$MEscore),"High","Low")
data2 = cbind(data[,sam$row.names[sam$group=="Low"]],data[,sam$row.names[sam$group=="High"]])

nl = length(sam$row.names[sam$group=="Low"])
nh = length(sam$row.names[sam$group=="High"])

outTab=data.frame()
grade=c(rep(1,nl),rep(2,nh))
data = data2

for(i in row.names(data)){
  geneName=unlist(strsplit(i,"\\|",))[1]
  geneName=gsub("\\/", "_", geneName)
  rt=rbind(expression=data[i,],grade=grade)
  rt=as.matrix(t(rt))
  wilcoxTest=wilcox.test(expression ~ grade, data=rt)
  conGeneMeans=mean(as.numeric(data[i,1:nl]))
  treatGeneMeans=mean(as.numeric(data[i,(nl+1):ncol(data)]))
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  pvalue=wilcoxTest$p.value
  conMed=median(as.numeric(data[i,1:nl]))
  treatMed=median(as.numeric(data[i,(nl+1):ncol(data)]))
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>1 & as.numeric(as.vector(outTab$fdr))<0.05),]

rt = read.table("target gene.txt",header = T,sep = "\t")
gene = unique(unlist(strsplit(rt$Target.genes,"; ")))

#################GO、kegg
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

entrezIDs=mget(gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]       

kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=0.05, qvalueCutoff=0.05)

KEGG=as.data.frame(kk)
KEGG = KEGG[KEGG$p.adjust < 0.05,]

selected = read.table("selected miRNA-targeted pathways.txt",sep = "\t")
KEGG = read.csv("target gene KEGG_diff.csv",header = T)   
kk2 = KEGG[KEGG$Description %in% selected$V1,]
miRNA = read.table("target gene.txt",header = T,sep = "\t")
miRNA = miRNA$miRNA.Product

target = strsplit(kk2$geneID,'/')
names(target) = kk2$Description
target0 = data.frame(key="",
                     path="")
for (i in 1:length(target)) {
  target1 = data.frame(key=target[[i]],
                       path = names(target[i]))
  target0 = rbind.data.frame(target0,target1)
}
target0 = target0[-1,]

eg = bitr(target0$key,fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
target0 = inner_join(eg,target0,by=c("ENTREZID"="key"))
target0 = target0[,-1]

miRNA = data.frame(key = miRNA)
diff = read.table("tcga LUAD miRNA diff.txt",header = T)
diff = mutate(diff,type = ifelse(logFC > 0,"src_up","src_dn"))
diff = diff[,c(1,7)]
diff$path = "source"
miRNA = inner_join(miRNA,diff,by=c("key"="gene"))
miRNA2 = read.table("target gene.txt",header = T,sep = "\t")
miRNA = miRNA[miRNA$key %in% miRNA2$miRNA.Product[miRNA2$Target.genes!="-"],]

gene = read.table("high-low DEA.txt",header = T)
gene = mutate(gene,type = ifelse(logFC > 0,"tar_up","tar_dn"))
gene = gene[,c(1,8)]
target = inner_join(gene,target0,by=c("id"="SYMBOL"))
colnames(target)[1] = "key"

node = rbind.data.frame(miRNA,target)
node$signif = ""


link = data.frame(src="",
                  tar="")
for (i in 1:nrow(miRNA2)) {
  link0 = data.frame(src=miRNA2$miRNA.Product[i],
                     tar = unlist(strsplit(miRNA2$Target.genes[i],"; ")))
  link = rbind.data.frame(link,link0)
}
link = link[-1,]
link = inner_join(link,miRNA,by=c("src"="key"))
link = link[,-4]
colnames(link)[3] = "source_type"


library(magrittr)
library(tidyverse)
library(ggplot2)
library(crosslink) 


links <- link
nodes <- node

paths <- unique(nodes[nodes$path != "source", ]$path) 
nodes$path <- factor(nodes$path, levels = c("source", paths))

src_up_col <- "red"
src_dn_col <- "blue"

tar_up_col <- "red"
tar_dn_col <- "blue"

nodes$key2 = c(nodes$key[1:14],paste("target",1:491))
data = nodes[15:nrow(nodes),c(1,5)]
links = inner_join(links,data,by=c("tar"="key"))
links = links[,c(1,4,3)]
colnames(links)[2] = "key"
nodes$key = c(nodes$key[1:14],paste("target",1:491))
nodes = nodes[,-5]

toy <- crosslink(
  nodes = nodes, 
  edges = links,
  cross.by = "path", 
  xrange = c(-5, 15),
  yrange = c(-5, 5),
  spaces = "partition")

cl_plot(toy)

toCircle <- function(x, y, rx = 1, ry =1, intensity = 2){
  mapTo2pi <- function(x) {scales::rescale(c(0, x), to = c(0, 2*pi))[-1]}
  data.frame(x, y) %>%
    mutate(group = paste0("group", x)) %>%
    mutate(yy = scales::rescale(-x, to = range(y))) %>%
    mutate(xx = mean(x) + intensity * sin(yy %>% mapTo2pi),) %>%
    group_by(group) %>%
    mutate(tri = rank(y, ties.method = "first") %>% mapTo2pi)  %>%
    ungroup() %$%
    data.frame(
      x = xx + rx*sin(tri),
      y = yy + ry*cos(tri))
}

toy_circle <- toy %>% tf_fun(
  crosses = paths, 
  along = "xy",
  fun = toCircle,
  rx = 0.2, ry = 0.2)

toy_circle %>% cl_plot(label = NA)

toy_final <- toy_circle %>% 
  tf_rotate(angle = -90) %>% 
  tf_shift(y = -6, crosses = paths, relative = F) %>%
  set_header()

toy_final %>% cl_plot(label = NA) %>% cl_void()
show_aes(toy_final)

ggplot() +
  ggforce::geom_circle(
    mapping = aes(x0 = x0, y0 = y0, r = r),
    data = get_cross(toy_final) %>% filter(cross != "source") %>% 
      group_by(path) %>%
      transmute(
        x0 = mean(x),
        y0 = mean(y),
        r = 0.2
      ) %>% unique(),
    show.legend = F
  ) +

  geom_segment(
    mapping = aes(x, y, xend = xend, yend = yend, color = source_type),
    data = get_link(toy_final),
    alpha = 0.3 
  ) + 
  
  geom_point(
    mapping = aes(x, y, 
                  color = type),
    data = get_cross(toy_final) %>% filter(cross != "source")
  ) +
  
  ggrepel::geom_text_repel(
    mapping = aes(x, y, label = header), nudge_y = 0.3, 
    data = get_header(toy_final) %>% filter(cross != "source"),
    segment.color = NA
  ) +

  geom_text(
    mapping = aes(x, y, label = key), angle = 90, hjust = 1, nudge_y = -0.1,
    data = get_cross(toy_final) %>% filter(cross == "source")
  ) +

  geom_text(
    mapping = aes(x, y, label = num),
    data = get_cross(toy_final) %>% filter(cross != "source") %>% 
      group_by(path) %>%
      transmute(
        x = mean(x),
        y = mean(y),
        num = n()
      ) %>% unique()
  ) +
  

  scale_color_manual(values = c(
    src_up = src_up_col, src_dn = src_dn_col,
    tar_up = tar_up_col, tar_dn = tar_dn_col)) + 
  labs(x = NULL, y = "Target_Pathway") +
  scale_y_continuous(expand = expansion(mult = c(0.25,0.1))) -> p

p

p <- p + geom_point(
  mapping = aes(x, y,color=type),
  data = get_cross(toy_final) %>% filter(cross == "source")
)
p
ggsave("",width = 15,height = 8)


