#setwd("~/MemComjournal/result")
setwd("/home/hesn/MemComJournal/result")
library(ggplot2)
library(reshape2)
library(tikzDevice)
rm(list=ls())

#figures in the journal paper
#-----------------makespan of stage 1, figure 6----------------------
makespan.long<-read.table("./result_split_1",header=TRUE)

long.copy<-makespan.long
makespan.wide<-dcast(makespan.long,TreeName+NPR+CCR+AmountProcessors~Stage1,value.var = "Makespan")
makespan.wide$ASAP<-makespan.wide$ASAP/makespan.wide$Sequence
makespan.wide$AvoidChain<-makespan.wide$AvoidChain/makespan.wide$Sequence
makespan.wide$ImprovedSplit<-makespan.wide$ImprovedSplit/makespan.wide$Sequence
makespan.wide$SplitSubtrees<-makespan.wide$SplitSubtrees/makespan.wide$Sequence

temp.wide<-makespan.wide[makespan.wide$NPR==10000,]
temp.wide<-temp.wide[temp.wide$CCR%in%c(0.1,1,10),]
sum(temp.wide$ImprovedSplit<=1)/nrow(temp.wide)
#on 56% cases, ImprovedSplit is better than or equal to SplitSubtrees

long<-melt(makespan.wide,id.vars = c("TreeName","NPR","CCR","AmountProcessors"),variable.name = "Heuristics",value.name = "makespan" )
long<-long[long$NPR%in%c(100,1000,10000),]
long<-long[long$CCR%in%c(0.1,1,10),]
long$NPR<-1/long$NPR
long$NPR<-as.factor(long$NPR)
long$CCR<-as.factor(long$CCR)

levels(long$NPR)[levels(long$NPR)=="1e-04"]<-"PNR = 1e-04"
levels(long$NPR)[levels(long$NPR)=="0.001"]<-"PNR = 0.001"
levels(long$NPR)[levels(long$NPR)=="0.01"]<-"PNR = 0.01"

long<-long[which(long$Heuristics!="Sequence"),]

aggregate(long$makespan,list(long$NPR,long$Heuristics),quantile,c(0.5))

hues = seq(15, 375, length = length(levels(long$Heuristics)) + 1)
hcl(h = hues, l = 65, c = 100)[1:length(levels(long$Heuristics))]
levels(long$Heuristics)
cb_palette <- c(ASAP="#F8766D", AvoidChain="#B79F00",
                ImprovedSplit="#00BA38", SplitSubtrees="#619CFF")

#tikz("~/ChangjiangGou/MemCom/journal/figure_makespan1.tex",width = 7.2,height = 3)
ggplot(long,aes(x=CCR,y=makespan,fill=Heuristics))+geom_boxplot(outlier.size = 0.05)+
  scale_x_discrete(name="CCR")+scale_y_continuous(name="Makespan normalized to \\textbf{Sequence}")+
  guides(fill=FALSE)+
  scale_fill_manual(labels=c("ASAP","ASAPnochain","ImprovedSplit","SplitSubtrees"),values = cb_palette)+
  labs(fill="")+theme(legend.position = c(0.799,0.85),legend.background = element_blank(),legend.key = element_blank())+
  facet_grid(.~NPR)#+
  #geom_hline(yintercept=0.55)
dev.off()

long<-long.copy
long<-long[long$CCR%in%c(0.1,1,10),]
long<-long[long$NPR%in%c(100,1000,10000),]
long$NPR<-1/long$NPR
long$NPR<-as.factor(long$NPR)
long$CCR<-as.factor(long$CCR)
long$AmountProcessors[long$AmountProcessors<=2]=3
long$perc<-long$AmountSubtrees/long$AmountProcessors
long<-long[long$Stage1!="Sequence",]

levels(long$NPR)[levels(long$NPR)=="1e-04"]<-"PNR = 1e-04"
levels(long$NPR)[levels(long$NPR)=="0.001"]<-"PNR = 0.001"
levels(long$NPR)[levels(long$NPR)=="0.01"]<-"PNR = 0.01"

temp<-aggregate(long$perc,list(long$NPR,long$Stage1),quantile,c(0.5))
View(temp)

temp<-aggregate(long$perc,list(long$NPR,long$Stage1),mean)
View(temp)

wide<-dcast(long,TreeName+NPR+CCR~Stage1,value.var = "perc")
wide$AvoidChain<-wide$AvoidChain-wide$ASAP
mean(wide$AvoidChain)

levels(long$Stage1)[levels(long$Stage1)=="AvoidChain"]<-"ASAPnochain"

cb_palette <- c(ASAP="#F8766D", ASAPnochain="#B79F00",
                ImprovedSplit="#00BA38", SplitSubtrees="#619CFF")#,
                #FirstFit="#FFFFFF")
