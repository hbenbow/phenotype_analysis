library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(multcompView)

# read in all files
setwd("~/Documents/Sobia/Watkins/Data/")
for(file in c(7, 14, 21, 28)){
  path=paste(getwd(), "/tilling_", file, ".csv", sep="")
  data<-read.csv(path)
  data$Factor<-paste(data$Genotype, data$Treatment, data$Plant.no, data$Replication)
  assign(paste("day", file, sep=""), data)
}
rm(file)
rm(path)
rm(data)

day7.2<-day7[,c(8, 4, 5)]
day14.2<-day14[,c(8, 4, 5)]
day21.2<-day21[,c(8, 4, 5)]
day28.2<-day28[,c(8, 4, 5)]
colnames(day7.2)<-c("Factor", "Necrosis.7", "Pycnidia.7")
colnames(day14.2)<-c("Factor", "Necrosis.14", "Pycnidia.14")
colnames(day21.2)<-c("Factor", "Necrosis.21", "Pycnidia.21")
colnames(day28.2)<-c("Factor", "Necrosis.28", "Pycnidia.28")
meta<-day7[,c(8, 1, 2, 3, 7)]
start<-merge(day7.2, day14.2, by="Factor")
mid<-merge(start, day21.2, by="Factor")
last<-merge(mid, day28.2, by="Factor")
all<-merge(last, meta, by="Factor")


all$AUDPCnecrosis<-(((all$Necrosis.7 + all$Necrosis.14) / 2) * 7) + 
  (((all$Necrosis.14 + all$Necrosis.21) / 2) * 7) +
  (((all$Necrosis.21 + all$Necrosis.28) / 2) * 7)

all$AUDPCPycnidia<-(((all$Pycnidia.7 + all$Pycnidia.14) / 2) * 7) + 
  (((all$Pycnidia.14 + all$Pycnidia.21) / 2) * 7) +
  (((all$Pycnidia.21 + all$Pycnidia.28) / 2) * 7)


 # define function that will calculate letters for pairwise comparisons
tri.to.squ<-function(x)
{
  rn<-row.names(x)
  cn<-colnames(x)
  an<-unique(c(cn,rn))
  myval<-x[!is.na(x)]
  mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
    
  }
  return(mymat)
}
# Read data in and wrangle
# ============================================================================
AUDPC<-na.omit(all)[,c(1, 10:15)]
AUDPC$Factor<-  paste(AUDPC$Genotype, AUDPC$Treatment)
AUDPC$Factor<- as.factor(AUDPC$Factor)
# ============================================================================

necrosis<-
  aggregate(AUDPC$AUDPCnecrosis, 
            by=list(AUDPC$Genotype, AUDPC$Treatment), FUN=mean)
necrosis$SD<- 
  aggregate(AUDPC$AUDPCnecrosis, 
            by=list(AUDPC$Genotype, AUDPC$Treatment), FUN=sd)$x
necrosis$n<- 
  aggregate(AUDPC$AUDPCnecrosis, 
            by=list(AUDPC$Genotype, AUDPC$Treatment), FUN=length)$x
necrosis$SE<-necrosis$SD/sqrt(necrosis$n)
colnames(necrosis)<-c("Genotype", "Treatment", "AUDPC", "SD", "n", "SE")
necrosis$Genotype<-
  factor(necrosis$Genotype, levels=c("Cadenza", "CAD451-BC2" ,"CAD499-BC2"))
necrosis$Treatment<-
  factor(necrosis$Treatment  , levels=c("T20","IPO323" ,  "IPO88004" ,"IPO89011", "IPO90012" ,"IPO94269"), 
         labels=c("Tween 20","IPO323" ,  "IPO88004" ,"IPO89011", "IPO90012" ,"IPO94269"))


pycnidia<-
  aggregate(AUDPC$AUDPCPycnidia, 
            by=list(AUDPC$Genotype, AUDPC$Treatment), FUN=mean)
pycnidia$SD<- 
  aggregate(AUDPC$AUDPCPycnidia, 
            by=list(AUDPC$Genotype, AUDPC$Treatment), FUN=sd)$x
pycnidia$n<- 
  aggregate(AUDPC$AUDPCPycnidia, 
            by=list(AUDPC$Genotype, AUDPC$Treatment), FUN=length)$x
pycnidia$SE<-pycnidia$SD/sqrt(pycnidia$n)
colnames(pycnidia)<-c("Genotype", "Treatment", "AUDPC", "SD", "n", "SE")
pycnidia$Genotype<-
  factor(pycnidia$Genotype, levels=c("Cadenza", "CAD451-BC2" ,"CAD499-BC2"))
pycnidia$Treatment<-
  factor(pycnidia$Treatment  , levels=c("T20","IPO323" ,  "IPO88004" ,"IPO89011", "IPO90012" ,"IPO94269"), 
         labels=c("Tween 20","IPO323" ,  "IPO88004" ,"IPO89011", "IPO90012" ,"IPO94269"))

AUDPC$Genotype<-
  factor(AUDPC$Genotype, levels=c("Cadenza", "CAD451-BC2" ,"CAD499-BC2"))
AUDPC$Treatment<-
  factor(AUDPC$Treatment  , levels=c("T20","IPO323" ,  "IPO88004" ,"IPO89011", "IPO90012" ,"IPO94269"), 
         labels=c("Tween 20","IPO323" ,  "IPO88004" ,"IPO89011", "IPO90012" ,"IPO94269"))




# test for normality
shapiro.test(AUDPC$AUDPCnecrosis)
shapiro.test(AUDPC$AUDPCPycnidia)

# do pairwise stats on each Variety
list<-list()
for(i in AUDPC$Treatment){
  data<-AUDPC[(AUDPC$Treatment==i),]
  necro<-necrosis[(necrosis$Treatment==i),]
  pp_necro<-pairwise.wilcox.test(data$AUDPCnecrosis, data$Genotype,
                                 p.adjust.method = "BH")
  mymat_necro<-tri.to.squ(pp_necro$p.value)
  myletters_necro<-as.data.frame(multcompLetters(mymat_necro,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_necro)<-"Letters"
  myletters_necro$Genotype<-row.names(myletters_necro)
  ne<-merge(necro, myletters_necro, by="Genotype")
  list[[length(list)+1]]<-ne
}

all_necrosis<-as.data.frame(do.call(rbind, list))
all_necrosis$Factor<-paste(all_necrosis$Treatment, all_necrosis$Genotype)
all_necrosis<-all_necrosis[!(duplicated(all_necrosis$Factor)),]
all_necrosis$Treatment<-factor(all_necrosis$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269"))

ggplot(all_necrosis, aes(x=Genotype, y=AUDPC)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.9))+ 
  theme_bw() + theme(legend.position = "right") + 
  facet_wrap(~Treatment, ncol=2) +
  geom_errorbar(aes(ymin=AUDPC - 2*SE, ymax= AUDPC+2*SE),position=position_dodge(width=0.9), width=0.3) + 
  ylab("AUDPC necrosis") +
  scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=20, hjust=1)) +
  geom_text(aes(x=Genotype, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0,1600))


list<-list()
for(i in AUDPC$Treatment){
  data<-AUDPC[(AUDPC$Treatment==i),]
  necro<-pycnidia[(pycnidia$Treatment==i),]
  pp_necro<-pairwise.wilcox.test(data$AUDPCPycnidia, data$Genotype,
                                 p.adjust.method = "BH")
  mymat_necro<-tri.to.squ(pp_necro$p.value)
  myletters_necro<-as.data.frame(multcompLetters(mymat_necro,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_necro)<-"Letters"
  myletters_necro$Genotype<-row.names(myletters_necro)
  ne<-merge(necro, myletters_necro, by="Genotype")
  list[[length(list)+1]]<-ne
}

all_pycnidia<-as.data.frame(do.call(rbind, list))
all_pycnidia$Factor<-paste(all_pycnidia$Treatment, all_pycnidia$Genotype)
all_pycnidia<-all_pycnidia[!(duplicated(all_pycnidia$Factor)),]
all_pycnidia$Treatment<-factor(all_pycnidia$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269"))

ggplot(all_pycnidia, aes(x=Genotype, y=AUDPC)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.9))+ 
  theme_bw() + theme(legend.position = "right") + 
  facet_wrap(~Treatment, ncol=2) +
  geom_errorbar(aes(ymin=AUDPC - 2*SE, ymax= AUDPC+2*SE),position=position_dodge(width=0.9), width=0.3) + 
  ylab("AUDPC pycnidia") +
  scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=20, hjust=1)) +
  geom_text(aes(x=Genotype, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0,600))
