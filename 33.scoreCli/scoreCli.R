
library(plyr)
library(ggplot2)
library(ggpubr)

rt=read.table("Cup-FerroptosisScore.group.txt", header=T, sep="\t", check.names=F, row.names=1)

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,"fustat"]))]

rt1=rt[,c("fustat", "group")]
colnames(rt1)=c("trait", "group")

df=as.data.frame(table(rt1))

df=ddply(df, .(group), transform, percent = Freq/sum(Freq) * 100)
df=ddply(df, .(group), transform, pos = (cumsum(Freq) - 0.5 * Freq))

df$label=paste0(sprintf("%.0f", df$percent), "%")
df$group=factor(df$group, levels=c("Low", "High"))

p=ggplot(df, aes(x = factor(group), y = percent, fill = trait)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  scale_fill_manual(values=bioCol)+
  xlab("Cup-FerroptosisScore")+ ylab("Percent weight")+  guides(fill=guide_legend(title='fustat'))+
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_bw()
pdf(file="", width=4, height=5)
print(p)
dev.off()

rt2=rt[,c("fustat", "MEscore")]
colnames(rt2)=c("trait", "CuFescore")
type=levels(factor(rt2[,"trait"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

rt2$trait = factor(rt2$trait,levels = c("Alive","Dead"))
boxplot=ggboxplot(rt2, x="trait", y="CuFescore", fill="trait",
                  xlab="fustat",
                  ylab="CuFescore",
                  legend.title="Fustat",
                  palette=bioCol) + 
  stat_compare_means(comparisons=my_comparisons)
pdf(file="",width=4,height=4.5)
print(boxplot)
dev.off()

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
cli = data.table::fread("phenotype.tsv",header = T)
cli = cli[,c(1,6,78,80,98,44,42,43)]
colnames(cli) = c("id","age","gender","status","stage","T","M","N")
cli = cli[complete.cases(cli),]
cli = cli %>%
  filter(stage != "not reported"  & M != "" & N != "")  %>%
  mutate(stage = ifelse(stage == "stage i" | stage == "stage ia" | stage == "stage ib","stage I",
                        ifelse(stage == "stage ii" | stage == "stage iia" | stage == "stage iib","stage II",
                               ifelse(stage == "stage iiia" | stage == "stage iiib","stage III","stage IV"))),
         T = ifelse(T == "T1" | T == "T1a" | T == "T1b","T1",
                    ifelse(T == "T2" | T == "T2a" | T == "T2b","T2",
                           ifelse(T == "T3","T3",
                                  ifelse(T=="T4","T4","TX")))),
         M = ifelse(M == "M0","M0",
                    ifelse(M == "M1" | M == "M1a" | M == "M1b","M1","MX")))
colnames(cli) = c("id","Age","Gender","Status","Stage","T","M","N")

score=read.table("Cup-FerroptosisScore.txt", header=T, sep="\t", check.names=F, row.names=NULL)

data = data.frame(inner_join(cli,score,by=c("id"="row.names")))
data = gather(data,clinical,value,-c(1,ncol(data)))

data$value = factor(data$value,levels = c(
  "age>65"  ,  "<=65"    ,  "female"  ,  "male" ,"Alive" ,"Dead" ,"stage I"  ,"stage II" ,"stage III", "stage IV" ,"T1" ,"T2", "T3" ,"T4", "TX" ,"M0" ,"M1" ,"MX", "N0" ,"N1", "N2" , "N3" ,"NX") 
)
data$clinical = factor(data$clinical,levels = c("Age","Gender","T","M","N","Stage","Status"))
color = c("#a76ba0","#bd8eb5",
          "#eea4ab","#bc2c45",
          "#9db2c5","#69c2c5","#2f7673",
          "#66b877","#99a396","#536034","#6e7c88","#40b7b1",
          "#66823e","#faf3c3","#e2d889","#f1c744",
          "#44629e","#2e2e51",
          "#7e478a","#3c4287","#3fa56e","#79a34a","#963b3e"
         )
ggplot(data,aes(clinical,MEscore,fill=value)) +
  geom_boxplot() +
  stat_compare_means(aes(label = ..p.signif..)) +
  scale_fill_manual(values = color) +
  labs(x="",y="CuFescore") +
  theme_classic()
ggsave("",width = 6,height = 5)

