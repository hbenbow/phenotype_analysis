library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(multcompView)
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

AUDPC_Pycnidia <- read.csv("~/Documents/Sobia/Watkins/AUDPC_Pycnidia.csv")
AUDPC_Chlorosis <- read.csv("~/Documents/Sobia/Watkins/AUDPC_Chlorosis.csv")
AUDPC_Chlorosis<-na.omit(AUDPC_Chlorosis)
AUDPC_Pycnidia<-na.omit(AUDPC_Pycnidia)

# necrosis analysis general
# ============================================================================
# Create a new column of "factor" which describes genotype and treatment
AUDPC_Chlorosis$Factor<-
  paste(AUDPC_Chlorosis$Variety, AUDPC_Chlorosis$Treatment)
# convert to factor
AUDPC_Chlorosis$Factor<-as.factor(AUDPC_Chlorosis$Factor)

# aggregate data for plttting
AUDPC_Chlorosis<-na.omit(AUDPC_Chlorosis)
Chlorosis<-
  aggregate(AUDPC_Chlorosis$AUDPC, 
            by=list(AUDPC_Chlorosis$Variety, AUDPC_Chlorosis$Treatment), FUN=mean)
Chlorosis$SD<- 
  aggregate(AUDPC_Chlorosis$AUDPC, 
            by=list(AUDPC_Chlorosis$Variety, AUDPC_Chlorosis$Treatment), FUN=sd)$x
Chlorosis$n<- 
  aggregate(AUDPC_Chlorosis$AUDPC, 
            by=list(AUDPC_Chlorosis$Variety, AUDPC_Chlorosis$Treatment), FUN=length)$x
Chlorosis$SE<-Chlorosis$SD/sqrt(Chlorosis$n)
colnames(Chlorosis)<-c("Genotype", "Treatment", "AUDPC", "SD", "n", "SE")
Chlorosis$Genotype<-
  factor(Chlorosis$Genotype, levels=c("Longbow","Stigg","WAT1190182", "WAT1190337","WAT1190912","WBCDB0009", "T1010004","T2220033"))

# test for normality
ks.test(AUDPC_Chlorosis$AUDPC, "pnorm", 
        alternative = "two.sided")
shapiro.test(AUDPC_Chlorosis$AUDPC)

# do pairwise stats on each Variety
list<-list()
for(i in AUDPC_Chlorosis$Variety){
  data<-AUDPC_Chlorosis[(AUDPC_Chlorosis$Variety==i),]
  chlor<-Chlorosis[(Chlorosis$Genotype==i),]
  pp_chlor<-pairwise.wilcox.test(data$AUDPC,data$Treatment,
                                 p.adjust.method = "BH")
  mymat_chlor<-tri.to.squ(pp_chlor$p.value)
  myletters_chlor<-as.data.frame(multcompLetters(mymat_chlor,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_chlor)<-"Letters"
  myletters_chlor$Treatment<-row.names(myletters_chlor)
  ch<-merge(chlor, myletters_chlor, by="Treatment")
  list[[length(list)+1]]<-ch
}

all_chlorosis<-as.data.frame(do.call(rbind, list))
all_chlorosis$Factor<-paste(all_chlorosis$Treatment, all_chlorosis$Genotype)
all_chlorosis<-all_chlorosis[!(duplicated(all_chlorosis$Factor)),]
all_chlorosis$Treatment<-factor(all_chlorosis$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))
# ============================================================================

# chlorosis stigg and longbow only
# ============================================================================
# Stigg and Longbow
data<-AUDPC_Chlorosis[(AUDPC_Chlorosis$Variety=="Longbow" | AUDPC_Chlorosis$Variety == "Stigg"),]
chlor<-Chlorosis[(Chlorosis$Genotype=="Longbow" | Chlorosis$Genotype=="Stigg"),]
chlor$Factor<-paste(chlor$Genotype,chlor$Treatment)
data$Factor<-paste(data$Variety, data$Treatment)
pp_chlor<-pairwise.wilcox.test(data$AUDPC,data$Factor,
                               p.adjust.method = "BH")
mymat_chlor<-tri.to.squ(pp_chlor$p.value)
myletters_chlor<-as.data.frame(multcompLetters(mymat_chlor,compare="<=",threshold=0.05,Letters=letters)$Letters)
colnames(myletters_chlor)<-"Letters"
myletters_chlor$Factor<-row.names(myletters_chlor)
ch<-merge(chlor, myletters_chlor, by="Factor")
ch$Treatment<-factor(ch$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))

ggplot(ch, aes(x=Treatment, y=AUDPC, group=Genotype)) + 
  geom_bar(aes(fill=Genotype), stat="identity", position=position_dodge(width=0.9))+ 
  theme_bw() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE, group=Genotype),position=position_dodge(width=0.9), width=0.3) + 
  ylab("AUDPC Necrosis") +
  scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=20, hjust=1),
        axis.text.y = element_text(colour="black", size=15)) +
  geom_text(aes(x=Treatment, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0,1500))
ggsave("~/Documents/Sobia/Watkins/S_L_necrosis.pdf")

# ============================================================================

# Compare necrosis to stigg and longbow
# ============================================================================
list<-list()
for(i in AUDPC_Chlorosis$Treatment){
  data<-AUDPC_Chlorosis[(AUDPC_Chlorosis$Treatment==i),]
  Pyc<-Chlorosis[(Chlorosis$Treatment==i),]
  pp_Pyc<-pairwise.wilcox.test(data$AUDPC,data$Variety,
                               p.adjust.method = "BH")
  mymat_Pyc<-tri.to.squ(pp_Pyc$p.value)
  myletters_Pyc<-as.data.frame(multcompLetters(mymat_Pyc,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_Pyc)<-"Letters"
  myletters_Pyc$Variety<-row.names(myletters_Pyc)
  py<-merge(Pyc, myletters_Pyc,by.x="Genotype", by.y="Variety")
  list[[length(list)+1]]<-py
}

all_Chlorosis_var<-as.data.frame(do.call(rbind, list))
all_Chlorosis_var$Factor<-paste(all_Chlorosis_var$Treatment, all_Chlorosis_var$Genotype)
all_Chlorosis_var<-all_Chlorosis_var[!(duplicated(all_Chlorosis_var$Factor)),]
all_Chlorosis_var$Treatment<-factor(all_Chlorosis_var$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))

ggplot(all_Chlorosis_var[!(all_Chlorosis_var$Treatment=="Tween 20"),], aes(x=Genotype, y=AUDPC)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.9))+ 
  theme_bw() + theme(legend.position = "right") + 
  facet_wrap(~Treatment, ncol=2) +
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE),position=position_dodge(width=0.9), width=0.3) + 
  ylab("AUDPC Chlorosis") +
  scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=20, hjust=1)) +
  geom_text(aes(x=Genotype, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0,1500))

# Pycnidia analysis general

# ============================================================================
# Create a new column of "factor" which describes genotype and treatment
AUDPC_Pycnidia$Factor<-
  paste(AUDPC_Pycnidia$Variety, AUDPC_Pycnidia$Treatment)
# convert to factor
AUDPC_Pycnidia$Factor<-as.factor(AUDPC_Pycnidia$Factor)

# aggregate data for plttting
AUDPC_Pycnidia<-na.omit(AUDPC_Pycnidia)
Pycnidia<-
  aggregate(AUDPC_Pycnidia$AUDPC, 
            by=list(AUDPC_Pycnidia$Variety, AUDPC_Pycnidia$Treatment), FUN=mean)
Pycnidia$SD<- 
  aggregate(AUDPC_Pycnidia$AUDPC, 
            by=list(AUDPC_Pycnidia$Variety, AUDPC_Pycnidia$Treatment), FUN=sd)$x
Pycnidia$n<- 
  aggregate(AUDPC_Pycnidia$AUDPC, 
            by=list(AUDPC_Pycnidia$Variety, AUDPC_Pycnidia$Treatment), FUN=length)$x
Pycnidia$SE<-Pycnidia$SD/sqrt(Pycnidia$n)
colnames(Pycnidia)<-c("Genotype", "Treatment", "AUDPC", "SD", "n", "SE")
Pycnidia$Genotype<-
  factor(Pycnidia$Genotype, levels=c("Longbow","Stigg","WAT1190182", "WAT1190337","WAT1190912","WBCDB0009", "T1010004","T2220033"))

# test for normality
ks.test(AUDPC_Pycnidia$AUDPC, "pnorm", 
        alternative = "two.sided")
shapiro.test(AUDPC_Pycnidia$AUDPC)

# do pairwise stats on each Variety
list<-list()
for(i in AUDPC_Pycnidia$Variety){
  data<-AUDPC_Pycnidia[(AUDPC_Pycnidia$Variety==i),]
  Pyc<-Pycnidia[(Pycnidia$Genotype==i),]
  pp_Pyc<-pairwise.wilcox.test(data$AUDPC,data$Treatment,
                                 p.adjust.method = "BH")
  mymat_Pyc<-tri.to.squ(pp_Pyc$p.value)
  myletters_Pyc<-as.data.frame(multcompLetters(mymat_Pyc,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_Pyc)<-"Letters"
  myletters_Pyc$Treatment<-row.names(myletters_Pyc)
  py<-merge(Pyc, myletters_Pyc, by="Treatment")
  list[[length(list)+1]]<-py
}

all_Pycnidia<-as.data.frame(do.call(rbind, list))
all_Pycnidia$Factor<-paste(all_Pycnidia$Treatment, all_Pycnidia$Genotype)
all_Pycnidia<-all_Pycnidia[!(duplicated(all_Pycnidia$Factor)),]
all_Pycnidia$Treatment<-factor(all_Pycnidia$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))

ggplot(all_Pycnidia, aes(x=Treatment, y=AUDPC)) + 
  geom_bar(stat="identity", width=.6, fill="grey60")+ 
  theme_bw() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE, group=Treatment),position=position_dodge(width=0.9), width=0.3) + 
  facet_wrap(Genotype~., ncol=2)+
  ylab("AUDPC Pycnidia") +
  theme(text = element_text(size=15, colour="black"), axis.text.x = element_text(colour="black", size=12, angle=45, hjust=1)) +
  geom_text(aes(x=Treatment, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0, 1))+
  scale_y_continuous(expand = c(0,3))+
  geom_hline(aes(yintercept=0))
ggsave("~/Documents/Sobia/Watkins/Pycnidia.pdf")
# ============================================================================

# Compare pycnidia AUDPC to Stigg and Longbow
# ===========================================================
list<-list()
for(i in AUDPC_Pycnidia$Treatment){
  data<-AUDPC_Pycnidia[(AUDPC_Pycnidia$Treatment==i),]
  Pyc<-Pycnidia[(Pycnidia$Treatment==i),]
  pp_Pyc<-pairwise.wilcox.test(data$AUDPC,data$Variety,
                               p.adjust.method = "BH")
  mymat_Pyc<-tri.to.squ(pp_Pyc$p.value)
  myletters_Pyc<-as.data.frame(multcompLetters(mymat_Pyc,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_Pyc)<-"Letters"
  myletters_Pyc$Variety<-row.names(myletters_Pyc)
  py<-merge(Pyc, myletters_Pyc,by.x="Genotype", by.y="Variety")
  list[[length(list)+1]]<-py
}

all_Pycnidia_var<-as.data.frame(do.call(rbind, list))
all_Pycnidia_var$Factor<-paste(all_Pycnidia_var$Treatment, all_Pycnidia_var$Genotype)
all_Pycnidia_var<-all_Pycnidia_var[!(duplicated(all_Pycnidia_var$Factor)),]
all_Pycnidia_var$Treatment<-factor(all_Pycnidia_var$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))

ggplot(all_Pycnidia_var[!(all_Pycnidia_var$Treatment=="Tween 20"),], aes(x=Genotype, y=AUDPC)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.9))+ 
  theme_bw() + theme(legend.position = "right") + 
  facet_wrap(~Treatment, ncol=2) +
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE),position=position_dodge(width=0.9), width=0.3) + 
  ylab("AUDPC Pycnidia") +
  scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=20, hjust=1)) +
  geom_text(aes(x=Genotype, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0,100))
# ============================================================================


# Stigg and Longbow only pycnidia
# ===========================================================
data<-AUDPC_Pycnidia[(AUDPC_Pycnidia$Variety=="Longbow" | AUDPC_Pycnidia$Variety == "Stigg"),]
chlor<-Pycnidia[(Pycnidia$Genotype=="Longbow" | Pycnidia$Genotype=="Stigg"),]
chlor$Factor<-paste(chlor$Genotype,chlor$Treatment)
data$Factor<-paste(data$Variety, data$Treatment)
pp_chlor<-pairwise.wilcox.test(data$AUDPC,data$Factor,
                               p.adjust.method = "BH")
mymat_chlor<-tri.to.squ(pp_chlor$p.value)
myletters_chlor<-as.data.frame(multcompLetters(mymat_chlor,compare="<=",threshold=0.05,Letters=letters)$Letters)
colnames(myletters_chlor)<-"Letters"
myletters_chlor$Factor<-row.names(myletters_chlor)
py<-merge(chlor, myletters_chlor, by="Factor")
py$Treatment<-factor(py$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))

ggplot(py, aes(x=Treatment, y=AUDPC, group=Genotype)) + 
  geom_bar(aes(fill=Genotype), stat="identity", position=position_dodge(width=0.9))+ 
  theme_bw() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE, group=Genotype),position=position_dodge(width=0.9), width=0.3) + 
  ylab("AUDPC Pycnidia") +
  scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=20, hjust=1)) +
  geom_text(aes(x=Treatment, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0,450))
ggsave("~/Documents/Sobia/Watkins/S_L_pycnidia.pdf")

# ============================================================================




