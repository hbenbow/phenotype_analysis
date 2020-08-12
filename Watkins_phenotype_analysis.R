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

AUDPC_pycnidia <- read.csv("~/Documents/Sobia/Watkins/Data/AUDPC_Pycnidia.csv")
AUDPC_necrosis <- read.csv("~/Documents/Sobia/Watkins/Data/AUDPC_necrosis.csv")
AUDPC_necrosis<-na.omit(AUDPC_necrosis)
AUDPC_pycnidia<-na.omit(AUDPC_pycnidia)




# ============================================================================
# necrosis analysis general
# Create a new column of "factor" which describes genotype and treatment
AUDPC_necrosis$Factor<-
  paste(AUDPC_necrosis$Variety, AUDPC_necrosis$Treatment)
# convert to factor
AUDPC_necrosis$Factor<-as.factor(AUDPC_necrosis$Factor)

# aggregate data for plotting
AUDPC_necrosis<-na.omit(AUDPC_necrosis)
necrosis<-
  aggregate(AUDPC_necrosis$AUDPC, 
            by=list(AUDPC_necrosis$Variety, AUDPC_necrosis$Treatment), FUN=mean)
necrosis$SD<- 
  aggregate(AUDPC_necrosis$AUDPC, 
            by=list(AUDPC_necrosis$Variety, AUDPC_necrosis$Treatment), FUN=sd)$x
necrosis$n<- 
  aggregate(AUDPC_necrosis$AUDPC, 
            by=list(AUDPC_necrosis$Variety, AUDPC_necrosis$Treatment), FUN=length)$x
necrosis$SE<-necrosis$SD/sqrt(necrosis$n)
colnames(necrosis)<-c("Genotype", "Treatment", "AUDPC", "SD", "n", "SE")
necrosis$Genotype<-
  factor(necrosis$Genotype, levels=c("Longbow","Stigg","WAT1190182", "WAT1190337","WAT1190912","WBCDB0009", "T1010004","T2220033"))

# test for normality
ks.test(AUDPC_necrosis$AUDPC, "pnorm", 
        alternative = "two.sided")
shapiro.test(AUDPC_necrosis$AUDPC)

# do pairwise stats on each Variety
list<-list()
for(i in AUDPC_necrosis$Variety){
  data<-AUDPC_necrosis[(AUDPC_necrosis$Variety==i),]
  necro<-necrosis[(necrosis$Genotype==i),]
  pp_necro<-pairwise.wilcox.test(data$AUDPC,data$Treatment,
                                 p.adjust.method = "BH")
  mymat_necro<-tri.to.squ(pp_necro$p.value)
  myletters_necro<-as.data.frame(multcompLetters(mymat_necro,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_necro)<-"Letters"
  myletters_necro$Treatment<-row.names(myletters_necro)
  ne<-merge(necro, myletters_necro, by="Treatment")
  list[[length(list)+1]]<-ne
}

all_necrosis<-as.data.frame(do.call(rbind, list))
all_necrosis$Factor<-paste(all_necrosis$Treatment, all_necrosis$Genotype)
all_necrosis<-all_necrosis[!(duplicated(all_necrosis$Factor)),]
all_necrosis$Treatment<-factor(all_necrosis$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))
# necrosis analysis general
# ============================================================================

# pycnidia analysis general
# ============================================================================
# Create a new column of "factor" which describes genotype and treatment
AUDPC_pycnidia$Factor<-
  paste(AUDPC_pycnidia$Variety, AUDPC_pycnidia$Treatment)
# convert to factor
AUDPC_pycnidia$Factor<-as.factor(AUDPC_pycnidia$Factor)

# aggregate data for plotting
AUDPC_pycnidia<-na.omit(AUDPC_pycnidia)
pycnidia<-
  aggregate(AUDPC_pycnidia$AUDPC, 
            by=list(AUDPC_pycnidia$Variety, AUDPC_pycnidia$Treatment), FUN=mean)
pycnidia$SD<- 
  aggregate(AUDPC_pycnidia$AUDPC, 
            by=list(AUDPC_pycnidia$Variety, AUDPC_pycnidia$Treatment), FUN=sd)$x
pycnidia$n<- 
  aggregate(AUDPC_pycnidia$AUDPC, 
            by=list(AUDPC_pycnidia$Variety, AUDPC_pycnidia$Treatment), FUN=length)$x
pycnidia$SE<-pycnidia$SD/sqrt(pycnidia$n)
colnames(pycnidia)<-c("Genotype", "Treatment", "AUDPC", "SD", "n", "SE")
pycnidia$Genotype<-
  factor(pycnidia$Genotype, levels=c("Longbow","Stigg","WAT1190182", "WAT1190337","WAT1190912","WBCDB0009", "T1010004","T2220033"))

# test for normality
ks.test(AUDPC_pycnidia$AUDPC, "pnorm", 
        alternative = "two.sided")
shapiro.test(AUDPC_pycnidia$AUDPC)

# do pairwise stats on each Variety
list<-list()
for(i in AUDPC_pycnidia$Variety){
  data<-AUDPC_pycnidia[(AUDPC_pycnidia$Variety==i),]
  pyc<-pycnidia[(pycnidia$Genotype==i),]
  pp_pyc<-pairwise.wilcox.test(data$AUDPC,data$Treatment,
                                 p.adjust.method = "BH")
  mymat_pyc<-tri.to.squ(pp_pyc$p.value)
  myletters_pyc<-as.data.frame(multcompLetters(mymat_pyc,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_pyc)<-"Letters"
  myletters_pyc$Treatment<-row.names(myletters_pyc)
  py<-merge(pyc, myletters_pyc, by="Treatment")
  list[[length(list)+1]]<-py
}

all_pycnidia<-as.data.frame(do.call(rbind, list))
all_pycnidia$Factor<-paste(all_pycnidia$Treatment, all_pycnidia$Genotype)
all_pycnidia<-all_pycnidia[!(duplicated(all_pycnidia$Factor)),]
all_pycnidia$Treatment<-factor(all_pycnidia$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))
# ============================================================================


# necrosis stigg and longbow only
# ============================================================================
# Stigg and Longbow
data<-AUDPC_necrosis[(AUDPC_necrosis$Variety=="Longbow" | AUDPC_necrosis$Variety == "Stigg"),]
chlor<-necrosis[(necrosis$Genotype=="Longbow" | necrosis$Genotype=="Stigg"),]
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
ggsave("~/Documents/Sobia/Watkins/Raw figures/S_L_necrosis.pdf")

# ============================================================================

# Compare necrosis to stigg and longbow
# ============================================================================
list<-list()
for(i in AUDPC_necrosis$Treatment){
  data<-AUDPC_necrosis[(AUDPC_necrosis$Treatment==i),]
  Chlor<-necrosis[(necrosis$Treatment==i),]
  pp_Chlor<-pairwise.wilcox.test(data$AUDPC,data$Variety,
                               p.adjust.method = "BH")
  mymat_Chlor<-tri.to.squ(pp_Chlor$p.value)
  myletters_Chlor<-as.data.frame(multcompLetters(mymat_Chlor,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_Chlor)<-"Letters"
  myletters_Chlor$Variety<-row.names(myletters_Chlor)
  py<-merge(Chlor, myletters_Chlor,by.x="Genotype", by.y="Variety")
  list[[length(list)+1]]<-py
}

all_necrosis_var<-as.data.frame(do.call(rbind, list))
all_necrosis_var$Factor<-paste(all_necrosis_var$Treatment, all_necrosis_var$Genotype)
all_necrosis_var<-all_necrosis_var[!(duplicated(all_necrosis_var$Factor)),]
all_necrosis_var$Treatment<-factor(all_necrosis_var$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))

ggplot(all_necrosis_var[!(all_necrosis_var$Treatment=="Isolate average"),], aes(x=Genotype, y=AUDPC)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.9))+ 
  theme_bw() + theme(legend.position = "right") + 
  facet_wrap(~Treatment, ncol=2) +
  geom_errorbar(aes(ymin=AUDPC - 2*SE, ymax= AUDPC+2*SE),position=position_dodge(width=0.9), width=0.3) + 
  ylab("AUDPC necrosis") +
  scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=20, hjust=1)) +
  geom_text(aes(x=Genotype, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0,1500))
ggsave("~/Documents/Sobia/Watkins/Raw figures/necrosis_by_isolate.pdf")
# ============================================================================

# Compare pycnidia AUDPC to Stigg and Longbow

list<-list()
for(i in AUDPC_pycnidia$Treatment){
  data<-AUDPC_pycnidia[(AUDPC_pycnidia$Treatment==i),]
  Chlor<-pycnidia[(pycnidia$Treatment==i),]
  pp_Chlor<-pairwise.wilcox.test(data$AUDPC,data$Variety,
                                 p.adjust.method = "BH")
  mymat_Chlor<-tri.to.squ(pp_Chlor$p.value)
  myletters_Chlor<-as.data.frame(multcompLetters(mymat_Chlor,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_Chlor)<-"Letters"
  myletters_Chlor$Variety<-row.names(myletters_Chlor)
  py<-merge(Chlor, myletters_Chlor,by.x="Genotype", by.y="Variety")
  list[[length(list)+1]]<-py
}

all_pycnidia_var<-as.data.frame(do.call(rbind, list))
all_pycnidia_var$Factor<-paste(all_pycnidia_var$Treatment, all_pycnidia_var$Genotype)
all_pycnidia_var<-all_pycnidia_var[!(duplicated(all_pycnidia_var$Factor)),]
all_pycnidia_var$Treatment<-factor(all_pycnidia_var$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))

ggplot(all_pycnidia_var[!(all_pycnidia_var$Treatment=="Isolate average"),], aes(x=Genotype, y=AUDPC)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.9))+ 
  theme_bw() + theme(legend.position = "right") + 
  facet_wrap(~Treatment, ncol=2) +
  geom_errorbar(aes(ymin=AUDPC -SE, ymax= AUDPC+SE),position=position_dodge(width=0.9), width=0.3) + 
  ylab("AUDPC pycnidia") +
  scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=20, hjust=1)) +
  geom_text(aes(x=Genotype, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0,500))
ggsave("~/Documents/Sobia/Watkins/Raw figures/pycnidia_by_isolate.pdf")


# Compare pycnidia AUDPC to Stigg and Longbow
# ===========================================================
list<-list()
for(i in AUDPC_pycnidia$Treatment){
  data<-AUDPC_pycnidia[(AUDPC_pycnidia$Treatment==i),]
  Chlor<-pycnidia[(pycnidia$Treatment==i),]
  pp_Chlor<-pairwise.wilcox.test(data$AUDPC,data$Variety,
                               p.adjust.method = "BH")
  mymat_Chlor<-tri.to.squ(pp_Chlor$p.value)
  myletters_Chlor<-as.data.frame(multcompLetters(mymat_Chlor,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_Chlor)<-"Letters"
  myletters_Chlor$Variety<-row.names(myletters_Chlor)
  py<-merge(Chlor, myletters_Chlor,by.x="Genotype", by.y="Variety")
  list[[length(list)+1]]<-py
}

all_pycnidia_var<-as.data.frame(do.call(rbind, list))
all_pycnidia_var$Factor<-paste(all_pycnidia_var$Treatment, all_pycnidia_var$Genotype)
all_pycnidia_var<-all_pycnidia_var[!(duplicated(all_pycnidia_var$Factor)),]
all_pycnidia_var$Treatment<-factor(all_pycnidia_var$Treatment, levels=c("Tween 20", "IPO323" ,  "IPO88004", "IPO89011","IPO90012" ,"IPO94269", "Isolate average"))

ggplot(all_pycnidia_var[!(all_pycnidia_var$Treatment=="Isolate average"),], aes(x=Genotype, y=AUDPC)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.9))+ 
  theme_bw() + theme(legend.position = "right") + 
  facet_wrap(~Treatment, ncol=2) +
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE),position=position_dodge(width=0.9), width=0.3) + 
  ylab("AUDPC pycnidia") +
  scale_fill_brewer(palette = "Set2")+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=20, hjust=1)) +
  geom_text(aes(x=Genotype, y=AUDPC+SE, label=Letters), position=position_dodge(width=0.9), vjust=-1, size=4)+
  coord_cartesian(ylim=c(0,500))
# ============================================================================


# ===========================================================


# Stigg and Longbow only Chlornidia
# ===========================================================
d

# ============================================================================