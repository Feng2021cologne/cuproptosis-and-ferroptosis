library(data.table)
library(magrittr)
library(dplyr)
library(tibble)
library(Seurat)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggpubr)
library(ggsci)
library(cowplot)

rt = fread("rt.txt",header=T) %>%
  column_to_rownames(var = "V1")
ann = read.table("GSE131907_Lung_Cancer_cell_annotation.txt",header = T,sep = "\t")
ann = ann$Index[ann$Sample_Origin %in% c("nLung","tLung")]
length(table(ann[ann$Sample_Origin %in% c("nLung","tLung"),]$Sample))#22

merged <- CreateSeuratObject(counts = rt, project = "Lung", min.cells = 0, min.features = 0)
dim(merged)

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
VizDimLoadings(merged, dims = 1:2, reduction = "pca")
DimHeatmap(merged, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(merged, ndims = 50)
dim.use = 1:30
merged <- FindNeighbors(merged, dims = dim.use)
merged <- FindClusters(merged, resolution = 0.5)

pdf("",width = 10,height = 10)
DotPlot(merged,features = c("MDK","SOX4","EPCAM", #cancer
                            "AGR3","FOLR1","SFTPD",#"PEBP4","AQP4",#"EPCAM", #Alveolar
                            #"AGER","SFTPC","LAMP3","SCGB1A1","FOXJ1","RFX2",  #Epithelial cells
                            "C1QB","LYZ","CD68", #Myeloid
                            "CLDN5","FCN3","RAMP2", #Endothelial
                            "C1R","COL1A1","DCN",#Fibroblasts
                            "CPA3","TPSAB1","TPSB2", #Mast cells
                            "CD79A","IGHG3","IGKC", #B cells
                            "CD3D","TRAC","TRBC2" #T cells
)) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

genes =c("MDK","SOX4","EPCAM", #cancer
         "AGR3","FOLR1","SFTPD",#"EPCAM", #ALveolar
         "AGER","SFTPC","LAMP3","SCGB1A1","FOXJ1","RFX2",  #Epithelial cells
         "C1QB","LYZ","CD68", #Myeloid
         "CLDN5","FCN3","RAMP2", #Endothelial
         "C1R","COL1A1","DCN",#Fibroblasts
         "CPA3","TPSAB1","TPSB2", #Mast cells
         "CD79A","IGHG3","IGKC", #B cells
         "CD3D","TRAC","TRBC2" #T cells
)
vln.df = data.frame(merged[["RNA"]]@data[rownames(merged[["RNA"]]@data) %in% genes,])
vln.df$gene = rownames(vln.df)
vln.df = melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)] = c("CB","exp")

anno = data.frame(CB = rownames(merged@meta.data),
                  celltype = merged@active.ident)
vln.df = inner_join(vln.df,anno,by="CB")
vln.df$gene = factor(vln.df$gene,levels = genes)

pdf("",width = 7,height = 9)
vln.df %>%
  ggplot(aes(gene,exp)) +
  geom_violin(aes(fill=celltype),scale = "width") +
  facet_grid(vln.df$celltype~.,scales = "free_y") +
  scale_fill_brewer(palette = "Set3",direction = 1) +
  scale_x_discrete("")+
  scale_y_continuous("")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
dev.off()

merged <- RunUMAP(merged, dims = dim.use)

new.cluster.ids <- c("0"="T cells",
                     "1"="T cells",
                     "2"="Myeloid cells",
                     "3"="T cells",
                     "4"="Myeloid cells",
                     "5"="B cells",
                     "6"="Myeloid cells",
                     "7"="Alveolar cells",
                     "8"="Fibroblasts",
                     "9"="Mast cells",
                     "10"="Cancer cells",
                     "11"="Myeloid cells",
                     "12"="Myeloid cells",
                     "13"="T cells",
                     "14"="Endothelial cells",
                     "15"="B cells",
                     "16"="Cancer cells",
                     "17"="Cancer cells",
                     "18"="T cells",
                     "19"="Alveolar cells", 
                     "20"="Cancer cells", 
                     "21"="Myeloid cells",
                     "22"="Alveolar cells",
                     "23"="Alveolar cells",
                     "24"="Cancer cells",
                     "25"="Endothelial cells",
                     "26"="unknow",
                     "27"="Cancer cells",
                     "28"="Fibroblasts"
)
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
merged = subset(merged,idents = c("T cells","Myeloid cells","B cells", "Alveolar cells","Fibroblasts","Mast cells","Cancer cells","Endothelial cells"))

mycols <- c("T cells"= "#d7301f",
            "Myeloid cells" = "#31a354",
            "B cells" = "#ffeda0",
            "Alveolar cells" = "#b4e28b" ,
            "Fibroblasts" = "#1f78b4",
            "Mast cells"  = "#984ea3", 
            "Cancer cells" = "#006d2c",
            "Endothelial cells"  = "#68c2cb"
)

pdf("", width = 6,height = 6)
DimPlot(merged, 
        reduction = "umap",
        pt.size = 0.5,
        label = TRUE,
        cols =  mycols) + NoLegend()
dev.off()

merged@meta.data$group = str_sub(merged@meta.data$patient_id,start=1L,end=1L)
pdf("", width = 7,height = 6)
DimPlot(merged, 
        reduction = "umap",
        pt.size = 0.5,
        group.by = 'group',
        label = F,
        cols =  c("#b0557e","#cccbc7")) 
dev.off()

##
library(tibble)
cell1 = data.frame(Idents(merged))
cell2 = merged@meta.data
cell1 = rownames_to_column(cell1,var = "cell")
cell2 = rownames_to_column(cell2,var = "cell")
cell = inner_join(cell1,cell2,by="cell")

num = cell %>%
  group_by(Idents.merged.,group) %>%
  summarise(sum = n()) %>%
  as.data.frame() 
num2 = cell %>%
  group_by(Idents.merged.) %>%
  summarise(total = n())
num = inner_join(num,num2,by="Idents.merged.")
num$percent = num$sum / num$total
num = num[num$group=="T",]
num = arrange(num,percent)
cell$Idents.merged. = factor(cell$Idents.merged.,levels = num$Idents.merged.)
pdf("",width = 6,height = 6)
ggplot(cell,aes(group)) +
  geom_bar(stat = "count",aes(fill=Idents.merged.),position = "fill") +
  scale_fill_manual(values = c("#4a9776","#bc683b","#67679b","#c64483","#71a251","#d2a549","#8b6d39","#cfd0ca")) +
  labs(x="",y="Fraction of cell") +
  guides(fill=guide_legend(title="Cell Type")) +
  theme_classic() 
dev.off()