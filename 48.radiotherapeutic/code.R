
library(survival)
library(survminer)
library(data.table)
library(dplyr)
score=read.table("Cup-FerroptosisScore.group.txt", header=T, sep="\t", check.names=F, row.names=1)
cli = data.frame(fread("phenotype.tsv",header = T))

score = inner_join(score,cli,by=c("id"="submitter_id.samples"))
score = mutate(score,group2 = ifelse(group=="Low" & radiation_therapy=="NO","score=Low,rad=NO",
                                     ifelse(group=="Low" & radiation_therapy=="YES","score=Low,rad=YES",
                                            ifelse(group=="High" & radiation_therapy=="NO","score=High,rad=NO","score=High,rad=YES"))))

length=length(levels(factor(score$group2)))
diff=survdiff(Surv(futime, fustat) ~ group2, data = score)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group2, data = score)

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(score[,"group2"])))]
surPlot=ggsurvplot(fit, 
                   data=score,
                   conf.int=F,
                   pval=pValue,
                   pval.size=6,
                   legend.labs=levels(factor(score[,"group2"])),
                   legend.title="group",
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 2,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)

pdf(file="", onefile = FALSE, width=8, height=8)
print(surPlot)
dev.off()

#################################
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggsci)
score=read.table("Cup-FerroptosisScore.group.txt", header=T, sep="\t", check.names=F, row.names=1)
cli = data.frame(fread("phenotype.tsv",header = T))
score = inner_join(score,cli,by=c("id"="submitter_id.samples"))
score = mutate(score,response = ifelse(Response=="Complete Remission/Response","CR",
                                 ifelse(Response=="Partial Remission/Response","PR",
                                        ifelse(Response=="Progressive Disease","PD","SD"))))
       
num = score %>%
  group_by(group,response) %>%
  summarise(count = n()) %>%
  arrange(group, desc(response)) %>%
  mutate(cumsum = cumsum(count),
         prop = count / sum(count),
         cumprop = cumsum(count) / sum(count))

M <- as.table(rbind(num$count[1:4], num$count[5:8]))
dimnames(M) <- list(group = c("High", "Low"),
                    response = num$response[1:4])
pvalue =  chisq.test(M)$p.value  
ggplot(data = num) +
  geom_col(aes(x = group, y = count, fill = response),
           position = 'fill',
           color = 'black',
           width = .7) +
  geom_text(x=1.5,y=1.1,label=pvalue) +
  labs(x="",y="") +
  scale_fill_lancet() +
  scale_y_continuous(expand = c(0, 0)) +
  theme_half_open()
ggsave("",width = 6,height = 6)

my_comparisons <- list( c("CR", "PD"), c("CR", "PR"), c("CR", "SD"),c("PD","PR"),c("PD","SD"),c("PR","SD"))
ggplot(score,aes(response,MEscore)) +
  geom_boxplot(aes(color=response)) +
  geom_jitter(aes(color=response)) +
  scale_color_manual(values = c("#82cbd4","#e1bf6f","#d9572b","#a3cc5a")) +
  stat_compare_means(label.y = 0.2) +
  stat_compare_means(comparisons = my_comparisons) +
  labs(x="",y="CuFescore") +
  theme_half_open()  +
  theme(legend.position = 'none') 
ggsave("",width = 6,height = 6)

score=read.table("Cup-FerroptosisScore.group.txt", header=T, sep="\t", check.names=F, row.names=1)
cli = data.frame(fread("phenotype.tsv",header = T))
cli = cli[,c(1,28,55)]
cli = cli[cli$radiation_therapy != "",]
colnames(cli)[2] = "Response"
score = inner_join(score,cli,by=c("id"="submitter_id.samples"))
score = mutate(score,response = ifelse(Response=="Complete Remission/Response","CR",
                                       ifelse(Response=="Partial Remission/Response","PR",
                                              ifelse(Response=="Progressive Disease","PD","SD"))))


num = score %>%
  group_by(group,radiation_therapy) %>%
  summarise(count = n()) %>%
  arrange(group, desc(radiation_therapy)) %>%
  mutate(cumsum = cumsum(count),
         prop = count / sum(count),
         cumprop = cumsum(count) / sum(count))

M <- as.table(rbind(num$count[1:2], num$count[3:4]))
pvalue =  chisq.test(M)$p.value  
ggplot(data = num) +
  geom_col(aes(x = group, y = count, fill = radiation_therapy),
           position = 'fill',
           color = 'black',
           width = .7) +
  geom_text(x=1.5,y=1.1,label=pvalue) +
  labs(x="",y="") +
  scale_fill_lancet() +
  scale_y_continuous(expand = c(0, 0)) +
  theme_half_open()
ggsave("",width = 6,height = 6)

ggplot(score,aes(radiation_therapy,MEscore)) +
  geom_boxplot(aes(color=radiation_therapy)) +
  geom_jitter(aes(color=radiation_therapy)) +
  scale_color_manual(values = c("#82cbd4","#e1bf6f","#d9572b","#a3cc5a")) +
  stat_compare_means(label.y = 0.1) +
  stat_compare_means(comparisons = my_comparisons) +
  labs(x="radiation_therapy",y="CuFescore") +
  theme_half_open()  +
  theme(legend.position = 'none') 
ggsave("",width = 6,height = 6)

library(pROC)
roc <- roc(score$radiation_therapy, score$MEscore,
                 ci=TRUE, print.auc=TRUE) 
pdf("",width = 6,height = 6)
plot(roc,
     legacy.axes = TRUE,
     print.auc=TRUE,
     col="#b6c8d8") 
dev.off()

