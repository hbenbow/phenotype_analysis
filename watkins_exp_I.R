library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(multcompView)
library(readxl)
# Read data in and wrangle
# ============================================================================

watkins_both<-read.csv(file="~/Documents/Sobia/Watkins/Data/Sobia_original_watkins_spss.csv", na.strings = "-")
watkins_both$Percentage_of_Longbow<-as.numeric(watkins_both$Percentage_of_Longbow)
watkins_both$Accession<-as.factor(watkins_both$Accession)
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


list<-list()
data0<-watkins_both
for(i in unique(data0$Strain)){
  data<-data0[(data0$Strain==i),]
  d2<-AUDPC[(AUDPC$Strain==i),]
  pp_necro<-dunnTest(data$Percentage_of_Longbow,data$Accession,
                               method = "bonferroni")$res
  comparison<-as.data.frame(do.call(rbind, strsplit(pp_necro$Comparison, split=" - ")))
  pp_necro$C1<-comparison$V1
  pp_necro$C2<-comparison$V2
  sigs<-pp_necro[(pp_necro$P.adj<0.05),]
  sigs$Strain<-paste(i)
  list[[length(list)+1]]<-sigs
}
all1<-as.data.frame(do.call(rbind, list))
sig_Longbow<-all1[(all1$C1=="Longbow"),]
sig_Longbow$Significant<-"*"
sig_Longbow$Factor<-paste(sig_Longbow$C2, sig_Longbow$Strain)
sig_Longbow_varieties<-sig_Longbow[,c(6:9)]
AUDPC$Factor<-paste(AUDPC$Accession, AUDPC$Strain)
AUDPC_by_strain<-left_join(AUDPC, sig_Longbow_varieties, by="Factor")


# isolate average
data<-watkins_both
data<-na.omit(data)
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
pp_necro<-dunnTest(data$Percentage_of_Longbow,data$Accession,
                   method = "bonferroni")$res
comparison<-as.data.frame(do.call(rbind, strsplit(pp_necro$Comparison, split=" - ")))
pp_necro$C1<-comparison$V1
pp_necro$C2<-comparison$V2
sigs<-pp_necro[(pp_necro$P.adj<0.05),]

sig_Longbow2<-sigs[(sigs$C1=="Longbow"),]
sig_Longbow2$Significant<-"*"
AUDPC3<-merge(AUDPC2, sig_Longbow2, by.x = "Accession", by.y="C2", all.x=T)
chosen<-as.data.frame(c(
  "T1010004", 
  "T2220033",
  "WBCDB0009",
  "WAT1190182", 
  "WAT1190912", 
  "WAT1190337"
))
chosen$symbol<-"+"
colnames(chosen) <-c("Accession", "Chosen")

AUDPC3<-merge(AUDPC3, chosen, by= "Accession", all.x=T)


order<-
  c( "WAT1190182",
     "WBCDB0056",
     "T1010012",
     "T2220009",
     "T1010004",
     "T1010003",
     "WAT1190912",
     "WAT1190402",
     "T2220053",
     "WAT1190337",
     "WAT1190110",
     "T2220033",
     "Stigg",
     "WAT1190756",
     "WAT1190482",
     "WBCDB0009",
     "T2220018",
     "WAT1190158",
     "WAT1190621",
     "T1010011",
     "WAT1190450",
     "WAT1190105",
     "WAT1190451",
     "Longbow",
     "WAT1190371",
     "WAT1190363",
     "T2220012")

AUDPC3$Accession<-factor(AUDPC3$Accession, levels=order)
AUDPC_by_strain$Accession<-factor(AUDPC_by_strain$Accession, levels=order)




ggplot(AUDPC_by_strain, aes(x=Accession, y=AUDPC)) + 
  geom_bar(stat="identity", aes(fill=Species))+ 
  theme_bw() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE), width=0.3) + 
  ylab("AUDPC (% of cv. Longbow)") +
  scale_fill_manual(values=c("#66C2A5" ,"#FC8D62", "#8DA0CB" ,"#E78AC3","#A6D854" ),
                    labels=c(expression(paste(italic("Aegilops tauschii"))),
                             expression(paste(italic("Triticum aestivum"))),
                             expression(paste(italic("Triticum durum"))),                           
                             expression(paste(italic("Triticum macha"))),
                             expression(paste(italic("Triticum urartu")))))+
  theme(text = element_text(size=25, colour="black"), 
        axis.text.x = element_text(colour="black", size=13, angle=45, hjust=1),
        axis.text.y = element_text(colour="black", size=15)) +
  facet_wrap(~Strain.x) +
  geom_hline(aes(yintercept=100), linetype="dashed", alpha=0.7)+
  geom_text(aes(x=Accession, y=AUDPC/2, label=Significant), position=position_dodge(width=0.9), vjust=-1, size=8)
ggsave("~/Documents/Sobia/Watkins/Raw figures/exp.1.pdf")



ggplot(AUDPC3, aes(x=Accession, y=AUDPC)) + 
  geom_bar(stat="identity", aes(fill=Species))+ 
  theme_bw() + theme(legend.position = "right") + 
  geom_errorbar(aes(ymin=AUDPC - SE, ymax= AUDPC+SE), width=0.3) + 
  ylab("AUDPC (% of cv. Longbow)") +
  scale_fill_manual(values=c("#66C2A5" ,"#FC8D62", "#8DA0CB" ,"#E78AC3","#A6D854"),
                    labels=c(expression(paste(italic("Aegilops tauschii"))),
                             expression(paste(italic("Triticum aestivum"))),
                             expression(paste(italic("Triticum durum"))),                           
                             expression(paste(italic("Triticum macha"))),
                             expression(paste(italic("Triticum urartu"))))) +
  theme(text = element_text(size=25, colour="black"), 
        axis.text.x = element_text(colour="black", size=15, angle=45, hjust=1),
        axis.text.y = element_text(colour="black", size=15)) +
  geom_hline(aes(yintercept=100), linetype="dashed", alpha=0.7)+  
  geom_text(aes(x=Accession, y=AUDPC+SE, label=Significant), position=position_dodge(width=0.9), vjust=-1, size=8)+
  geom_text(aes(x=Accession, y=0, label=Chosen), position=position_dodge(width=0.9), vjust=-1, size=8)
  
  ggsave("~/Documents/Sobia/Watkins/Raw figures/isolate_average.pdf")

