
library(maftools)     
library(TCGAbiolinks)
library(dplyr)
library(data.table)
library(stringr)

clinical <- GDCquery(project = "TCGA-LUAD", 
                     data.category = "Clinical", 
                     file.type = "xml")
GDCdownload(clinical)
cliquery <- GDCprepare_clinic(clinical, clinical.info = "patient")
colnames(cliquery)[1] <- "Tumor_Sample_Barcode"
str(cliquery)

cliquery = mutate(cliquery,stage_event_pathologic_stage=ifelse(stage_event_pathologic_stage=="Stage I" | stage_event_pathologic_stage=="Stage IA" | stage_event_pathologic_stage=="Stage IB","Stage I",
                                                               ifelse(stage_event_pathologic_stage=="Stage II" | stage_event_pathologic_stage=="Stage IIA" | stage_event_pathologic_stage=="Stage IIB","Stage II",
                                                                      ifelse(stage_event_pathologic_stage=="Stage IIIA" | stage_event_pathologic_stage=="Stage IIIB","Stage III","Stage IV"))),
                  age = ifelse(age_at_initial_pathologic_diagnosis >= 65,">=65",
                               ifelse(age_at_initial_pathologic_diagnosis < 65,"<65","NA")))

maf <- read.maf(maf = "", clinicalData = cliquery, isTCGA = T)

gene = read.table("gene.txt",header = F)


col = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit')

gender = c("#ea4c1c","#49b4ba")
names(gender) = c("FEMALE","MALE")
age = c("#fb2afd","#191781")
names(age) = c(">=65","<65")
vital_status = c("#404040","white")
names(vital_status) = c("Alive","Dead")
has_new_tumor_events_information = c("#973228","#b1e98f")
names(has_new_tumor_events_information) = c("NO","YES")
stage_event_pathologic_stage = c("#95332a","#aeeb93","#2f51a9","#dbae29")
names(stage_event_pathologic_stage) = c("Stage I" ,"Stage II" ,"Stage III"  ,"Stage IV" )


annocolors = list(age = age, 
                  gender = gender, 
                  vital_status = vital_status,
                  has_new_tumor_events_information = has_new_tumor_events_information,
                  stage_event_pathologic_stage = stage_event_pathologic_stage)

pdf("",width = 12,height = 6)
oncoplot(maf = maf,
         genes = gene$V1,
         colors = col,
         annotationColor = annocolors,
         top = 20,
         clinicalFeatures = c("age","gender","stage_event_pathologic_stage", "vital_status","has_new_tumor_events_information"),
         writeMatrix =T,
         anno_height = 2.2)
dev.off()


maf <- read.maf(maf = "")
pdf("",width = 6,height = 6)               
Interact <- somaticInteractions(maf = maf, genes = gene$V1, pvalue = c(0.05, 0.1))
Interact$pValue
dev.off()

data = data.frame(maf@maf.silent)[,c(1,16,9,10)]
data = data[data$Hugo_Symbol %in% gene$V1,]
mut_sample = str_sub(data$Tumor_Sample_Barcode,start=1L,end=16L) 

sur = data.frame(fread("",header = T))
sur = sur[,c(1,2,4)]
sur = sur[str_sub(sur$sample,start = 14L,end = 15L) != "11",]
colnames(sur) = c("id","fustat","futime")
sur$group = ifelse(sur$id %in% mut_sample,"Mut_type=Mutation","Mut_type=NonMutation")
sur$futime = sur$futime / 365

library(survival)
library(survminer)
fit=survfit(Surv(futime, fustat) ~group, data = sur)
diff=survdiff(Surv(futime, fustat) ~group, data = sur)
pValue=1-pchisq(diff$chisq, df=1)
surPlot=ggsurvplot(fit,
                   data=sur,
                   pval=pValue,
                   pval.size=6,
                   xlab="Time(years)",
                   ylab="Overall survival",
                   palette=c("#e0b759", "#84b5d8"),
                   break.time.by=5,
                   conf.int=F,
                   risk.table=TRUE,
                   risk.table.title="",
                   risk.table.height=.25)
pdf(file="",onefile = FALSE,
    width = 10,
    height =8)
print(surPlot)
dev.off()



