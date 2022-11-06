
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(stringr)
library(dplyr)
library(tibble)
library(data.table)

tcga = data.frame(fread("TPM.txt",header = T))
tcga = tcga[,which(str_sub(colnames(tcga),start = 14L,end = 15L)!="11")]
colnames(tcga) = gsub("\\.","-",colnames(tcga))

GSE30219 = read.table("../GSE30219/GSE30219.txt",header = T,row.names = NULL)
GSE31210 = read.table("../GSE31210/GSE31210.txt",header = T,row.names = NULL)
GSE3141 = read.table("../GSE3141/GSE3141.txt",header = T,row.names = NULL)
GSE37745 = read.table("../GSE37745/GSE37745.txt",header = T,row.names = NULL)
GSE81089 = read.table("../GSE81089/GSE81089_log2(TPM).txt",header = T,row.names = NULL,sep = "\t")

data = data.frame(fread("data.txt",header = T))
data = column_to_rownames(data,var = "V1")
data=as.matrix(data)
geneSets=getGmt("symbols.gmt", geneIdType=SymbolIdentifier())
gsvaResult=gsva(data, 
                geneSets, 
                min.sz=10, 
                max.sz=500, 
                verbose=TRUE,
                parallel.sz=6)
colnames(gsvaResult) = gsub("\\.","-",colnames(gsvaResult))
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)

gsvaResult=data.frame(fread("gsvaOut.txt",header = T))
gsvaResult=column_to_rownames(gsvaResult,var = "id")
cluster=read.table("Cluster.txt", header=T, sep="\t", check.names=F, row.names=1)
gsvaResult=t(gsvaResult)
row.names(gsvaResult) = gsub("\\.","-",row.names(gsvaResult))
sameSample=intersect(row.names(gsvaResult),rownames(cluster))
gsvaResult=gsvaResult[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
gsvaCluster=cbind(gsvaResult, cluster)

tcga = tcga[,colnames(tcga) %in% sameSample]
GSE30219 = GSE30219[,colnames(GSE30219) %in% sameSample]
GSE31210 = GSE31210[,colnames(GSE31210) %in% sameSample]
GSE3141 = GSE3141[,colnames(GSE3141) %in% sameSample]
GSE37745 = GSE37745[,colnames(GSE37745) %in% sameSample]
GSE81089 = GSE81089[,colnames(GSE81089) %in% sameSample]

all = data.frame(id = rownames(gsvaCluster),
                 `Cup-Ferroptosis cluster` = cluster$cluster,
                 Project = c(rep("TCGA",ncol(tcga)),rep("GSE30219",ncol(GSE30219)),
                             rep("GSE31210",ncol(GSE31210)),rep("GSE3141",ncol(GSE3141)),
                             rep("GSE37745",ncol(GSE37745)),rep("GSE81089",ncol(GSE81089))))

adj.P.Val.Filter=0.05
allType=as.vector(gsvaCluster$cluster)
comp=combn(levels(factor(allType)), 2)

i=1
treat=gsvaCluster[gsvaCluster$cluster==comp[2,i],]
con=gsvaCluster[gsvaCluster$cluster==comp[1,i],]
data=rbind(con, treat)

Type=as.vector(data$cluster)
ann=data[,ncol(data),drop=F]
data=t(data[,-ncol(data)])
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
fit=lmFit(data, design)
contrast=paste0(comp[2,i], "-", comp[1,i])
cont.matrix=makeContrasts(contrast, levels=design)
fit2=contrasts.fit(fit, cont.matrix)
fit2=eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)

diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
diffTermName=as.vector(rownames(diffSig))
diffLength=length(diffTermName)
selected = read.table("selected B-A pathways.txt",header = F)
hmExp=data[selected$V1,]

annotation_row = data.frame(GSVA = rownames(hmExp),
                            row.names = rownames(hmExp))
all2 = column_to_rownames(all,var = "id")
annotation_col = all2[colnames(hmExp),,drop=F]

ann_colors=list(
  `Cup-Ferroptosis cluster` = c(A="#304ca7",C="#e12f1f"),
  Project = c(TCGA="#7b73b8",
              GSE30219 = "#d43788",
              GSE31210 = "#45a27f",
              GSE3141 = "#cf571a",
              GSE37745 = "#c55734",
              GSE81089 = "#ccc662"
              ),
  GSVA = c(
    KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS = "#333b64",
    KEGG_BASAL_TRANSCRIPTION_FACTORS = "#c12922",
    KEGG_MISMATCH_REPAIR = "#f7c63d",
    KEGG_CELL_CYCLE ="#377e47",
    KEGG_COMPLEMENT_AND_COAGULATION_CASCADES ="#42aa47",
    KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION ="#80d9de",
    KEGG_DNA_REPLICATION ="#701f29",
    KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION ="#608cb1"
  ))
  
pdf(file=paste0(contrast,".heatmap.pdf"),height=6,width=15)
pheatmap(hmExp, 
         #annotation=ann,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#4a5eaf",2), "black", rep("#eee66d",2)))(50),
         cluster_cols =F,
         cluster_rows = F,
         show_colnames = F,
         gaps_col=as.vector(cumsum(table(Type))),
         scale="row",
         fontsize = 10,
         fontsize_row=7,
         fontsize_col=10)
dev.off()	


i=3
treat=gsvaCluster[gsvaCluster$cluster==comp[2,i],]
con=gsvaCluster[gsvaCluster$cluster==comp[1,i],]
data=rbind(con, treat)

Type=as.vector(data$cluster)
ann=data[,ncol(data),drop=F]
data=t(data[,-ncol(data)])
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
fit=lmFit(data, design)
contrast=paste0(comp[2,i], "-", comp[1,i])
cont.matrix=makeContrasts(contrast, levels=design)
fit2=contrasts.fit(fit, cont.matrix)
fit2=eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)

diffSig=allDiff[with(allDiff, (abs(logFC)>0.1 & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")

selected = read.table("selected C-B pathways.txt",header = F)
hmExp=data[selected$V1,]

annotation_row = data.frame(GSVA = rownames(hmExp),
                            row.names = rownames(hmExp))
all2 = column_to_rownames(all,var = "id")
annotation_col = all2[colnames(hmExp),,drop=F]

ann_colors=list(
  `Cup-Ferroptosis cluster` = c(B="#304ca7",C="#e12f1f"),
  Project = c(TCGA="#7b73b8",
              GSE30219 = "#d43788",
              GSE31210 = "#45a27f",
              GSE3141 = "#cf571a",
              GSE37745 = "#c55734",
              GSE81089 = "#ccc662"
  ),
  GSVA = c(
    KEGG_CELL_CYCLE = "#333b64",
    KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS = "#c12922",
    KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION = "#f7c63d",
    KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION ="#377e47",
    KEGG_MAPK_SIGNALING_PATHWAY ="#42aa47",
    KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION ="#80d9de",
    KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY ="#701f29",
    KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY ="#608cb1",
    KEGG_CHEMOKINE_SIGNALING_PATHWAY ="#389049",
    KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION ="#44203a",
    KEGG_REGULATION_OF_AUTOPHAGY ="#71578b"
  ))

pdf(file=paste0(contrast,".heatmap.pdf"),height=6,width=15)
pheatmap(hmExp, 
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#4a5eaf",2), "black", rep("#eee66d",2)))(50),
         cluster_cols =F,
         cluster_rows = F,
         show_colnames = F,
         gaps_col=as.vector(cumsum(table(Type))),
         scale="row",
         fontsize = 10,
         fontsize_row=7,
         fontsize_col=10)
dev.off()	

