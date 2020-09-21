library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(multcompView)
library(readxl)
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

Part_1_Watkins <- read.csv(file="~/Documents/Sobia/Watkins/Data/Watkins_phenotyping_1.csv")
Part_2_Watkins <- read.csv(file="~/Documents/Sobia/Watkins/Data/Watkins_phenotyping_2.csv", na.strings = "-")
watkins_both<-read.csv(file="~/Documents/Sobia/Watkins/Data/Watkins_phenotyping_both.csv", na.strings = "-")
# ============================================================================
# necrosis analysis general
# Create a new column of "factor" which describes genotype and treatme

data<-watkins_both
data$Percentage_of_Longbow<-as.numeric(data$Percentage_of_Longbow)
data$Factor<-
  paste(data$Accession, data$Strain)
# convert to factor
data$Factor<-as.factor(data$Factor)

# aggregate data for plotting
data<-na.omit(data)
AUDPC<-
  aggregate(data$Percentage_of_Longbow, 
            by=list(data$Accession, data$Strain, data$Species), FUN=mean)
AUDPC$SD<- 
  aggregate(data$Percentage_of_Longbow, 
            by=list(data$Accession, data$Strain, data$Species), FUN=sd)$x
AUDPC$n<-
  aggregate(data$Percentage_of_Longbow, 
            by=list(data$Accession, data$Strain, data$Species), FUN=length)$x

AUDPC$SE<-AUDPC$SD/sqrt(AUDPC$n)
colnames(AUDPC)<-c("Accession", "Strain","Species" ,"AUDPC", "SD", "n", "SE")

# test for normality
ks.test(AUDPC_total$AUDPC, "pnorm", 
        alternative = "two.sided")
shapiro.test(AUDPC_total$AUDPC)

# do pairwise stats on each Variety
list<-list()
for(i in unique(watkins_both$Strain)){
  data<-watkins_both[(watkins_both$Strain==i),]
  data$Percentage_of_Longbow<-as.numeric(data$Percentage_of_Longbow)
  d2<-AUDPC[(AUDPC$Strain==i),]
  pp_necro<-pairwise.wilcox.test(data$Percentage_of_Longbow,data$Accession,
                                 p.adjust.method = "BH")
  mymat_necro<-tri.to.squ(pp_necro$p.value)
  myletters_necro<-as.data.frame(multcompLetters(mymat_necro,compare="<=",threshold=0.05,Letters=letters)$Letters)
  colnames(myletters_necro)<-"Letters"
  myletters_necro$Accession<-row.names(myletters_necro)
  d3<-merge(d2, myletters_necro, by="Accession")
  list[[length(list)+1]]<-d3
}

all<-as.data.frame(do.call(rbind, list))

# isolate average
AUDPC2<-
  aggregate(data$Percentage_of_Longbow, 
            by=list(data$Accession, data$Species), FUN=mean)
AUDPC2$SD<- 
  aggregate(data$Percentage_of_Longbow, 
            by=list(data$Accession, data$Species), FUN=sd)$x
AUDPC2$n<-
  aggregate(data$Percentage_of_Longbow, 
            by=list(data$Accession, data$Species), FUN=length)$x

AUDPC2$SE<-AUDPC2$SD/sqrt(AUDPC2$n)
colnames(AUDPC2)<-c("Accession","Species" ,"AUDPC", "SD", "n", "SE")



data<-watkins_both
data$Percentage_of_Longbow<-as.numeric(data$Percentage_of_Longbow)
d2<-AUDPC2
pp_necro<-pairwise.wilcox.test(data$Percentage_of_Longbow,data$Accession,
                               p.adjust.method = "BH")
mymat_necro<-tri.to.squ(pp_necro$p.value)
myletters_necro<-as.data.frame(multcompLetters(mymat_necro,compare="<=",threshold=0.05,Letters=letters)$Letters)
colnames(myletters_necro)<-"Letters"
myletters_necro$Accession<-row.names(myletters_necro)
d3<-merge(d2, myletters_necro, by="Accession")
isolate_average<-d3


order<-
 c( "WAT1190182",
"T1010003",
"T1010012",
"T2220009",
"T1010004",
"Stigg",
"WAT1190912",
"WAT1190110",
"WAT1190337",
"WAT1190402",
"T2220033",
"W10052",
"WBCDB0009",
"T2220053",
"T2220018",
"WAT1190158",
"WAT1190756",
"WAT1190482",
"WAT1190621",
"T1010011",
"WAT1190451",
"WAT1190105",
"Longbow",
"WAT1190450",
"WAT1190371",
"WAT1190363",
"T2220012")

all$Accession<-factor(all$Accession, levels=order)
isolate_average$Accession<-factor(isolate_average$Accession, levels=order)
write.csv(isolate_average, "~/Documents/Sobia/Watkins/Data/isolate.average_stats.csv")
ggplot(all, aes(x=Accession, y=AUDPC)) + 
  geom_bar(stat="identity", aes(fill=Species))+ 
  theme_bw() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE), width=0.3) + 
  ylab("AUDPC") +
  scale_fill_manual(values=c("#66C2A5" ,"#FC8D62", "#8DA0CB" ,"#E78AC3"), 
                    labels=c(expression(paste(italic("Aegilops tauschii"))), 
                             expression(paste(italic("Triticum aestivum"))), 
                             expression(paste(italic("Triticum macha"))), 
                             expression(paste(italic("Triticum urartu")))))+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1),
        axis.text.y = element_text(colour="black", size=15)) +
  facet_wrap(~Strain) +
  geom_hline(aes(yintercept=100), linetype="dashed", alpha=0.7)
ggsave("~/Documents/Sobia/Watkins/Raw figures/exp.1.pdf")


isolate_average<-read.csv("~/Documents/Sobia/Watkins/Data/isolate.average_stats.csv")
isolate_average$Accession<-factor(isolate_average$Accession, levels=order)

ggplot(isolate_average, aes(x=Accession, y=AUDPC)) + 
  geom_bar(stat="identity", aes(fill=Species))+ 
  theme_bw() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE), width=0.3) + 
  ylab("AUDPC") +
  scale_fill_manual(values=c("#66C2A5" ,"#FC8D62", "#8DA0CB" ,"#E78AC3"), 
                    labels=c(expression(paste(italic("Aegilops tauschii"))), 
                             expression(paste(italic("Triticum aestivum"))), 
                             expression(paste(italic("Triticum macha"))), 
                             expression(paste(italic("Triticum urartu")))))+
  theme(text = element_text(size=20, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=45, hjust=1),
        axis.text.y = element_text(colour="black", size=15)) +
  geom_hline(aes(yintercept=100), linetype="dashed", alpha=0.7)+  
  geom_text(aes(x=Accession, y=AUDPC+SE, label=Significant), position=position_dodge(width=0.9), vjust=-1, size=6)+
ggsave("~/Documents/Sobia/Watkins/Raw figures/isolate_average.pdf")

