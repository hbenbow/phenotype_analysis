library(jtools)
library(interactions)
library(ggplot2)
library('knitr')
library("ggpubr")

setwd("~/Documents/Sobia/Watkins/Data/")
AUDPC_Pycnidia <-
  read.csv("~/Documents/Sobia/Watkins/Data/AUDPC_Pycnidia.csv")
AUDPC_Chlorosis <-
  read.csv("~/Documents/Sobia/Watkins/Data/AUDPC_Chlorosis.csv")
AUDPC_Chlorosis <- na.omit(AUDPC_Chlorosis)
AUDPC_Pycnidia <- na.omit(AUDPC_Pycnidia)
AUDPC_Chlorosis$Factor <-
  paste(AUDPC_Chlorosis$Variety, AUDPC_Chlorosis$Treatment)
# convert to factor
AUDPC_Chlorosis$Factor <- as.factor(AUDPC_Chlorosis$Factor)
AUDPC_Pycnidia$Factor <-
  paste(AUDPC_Pycnidia$Variety, AUDPC_Pycnidia$Treatment)
# convert to factor
AUDPC_Pycnidia$Factor <- as.factor(AUDPC_Pycnidia$Factor)

# Sources of variance etc

AUDPC_Chlorosis$ID <-
  paste(AUDPC_Chlorosis$Factor,
        AUDPC_Chlorosis$Replication,
        AUDPC_Chlorosis$Plant.no)
AUDPC_Pycnidia$ID <-
  paste(AUDPC_Pycnidia$Factor,
        AUDPC_Pycnidia$Replication,
        AUDPC_Pycnidia$Plant.no)
chlor.and.pyc <- merge(AUDPC_Chlorosis, AUDPC_Pycnidia, by = "ID")
chlor.and.pyc <- chlor.and.pyc[c(1, 2, 3, 4, 5, 6, 8, 7, 14)]
colnames(chlor.and.pyc) <-
  c(
    "ID",
    "Variety",
    "Genotype",
    "Treatment",
    "Plant.no",
    "Replication",
    "Factor",
    "AUDPCnecrosis",
    "AUDPCpycnidia"
  )


ggplot(chlor.and.pyc, aes(x = AUDPCnecrosis, y = AUDPCpycnidia)) +
  geom_smooth(method="lm", aes(colour=Treatment)) +
  facet_wrap(~Variety)

ggline(chlor.and.pyc, x = "AUDPCnecrosis", y = "AUDPCpycnidia", color="Variety",
       add = c("mean_se", "jitter"),
       ylab = "Weight", xlab = "Treatment")

av_necrosis <- aov(AUDPCnecrosis ~   Replication + Plant.no + Treatment * Variety, data = chlor.and.pyc)
summary(av_necrosis)
av_pycnidia <- aov(AUDPCpycnidia ~   Replication + Plant.no + Treatment * Variety, data = chlor.and.pyc)
summary(av_pycnidia)
knit(summary(av_necrosis))
interact_plot(fiti, pred = AUDPCnecrosis, modx = Treatment)

i="Longbow"

for (i in unique(chlor.and.pyc$Variety)) {
  data <- chlor.and.pyc[(chlor.and.pyc$Variety == i), ]
  data <- data[!(data$Treatment == "Isolate average"), ]
  fiti <- lm(AUDPCpycnidia ~ AUDPCnecrosis *Treatment, data = data)
  assign(paste(i, "_lm", sep=""), summary(fiti))
  corrs<-list()
  r2l<-list()
  for(iso in unique(data$Treatment)){
    data2<-data[(data$Treatment==iso),]
    corr<-as.data.frame(cor(data2$AUDPCnecrosis, data2$AUDPCpycnidia, method="spearman"))
    corr$isolate<-paste(iso)
    corrs[[(length(corrs)+1)]]<-corr
    fi2 <- lm(AUDPCpycnidia ~ AUDPCnecrosis, data = data2)
    r2<-as.data.frame(summary(lm2)$adj.r.squared)
    r2$isolate<-paste(iso)
    r2l[[(length(r2l)+1)]]<-r2
    }
  correl<-as.data.frame(do.call(rbind, corrs))
  r2s<-as.data.frame(do.call(rbind, r2l))
  r2s$Variety<-paste(i)
  correl$Variety<-paste(i)
  assign(paste(i, "_corr", sep=""), correl)
  assign(paste(i, "_R2", sep=""), r2s)
  ggplot(data, aes(x = AUDPCnecrosis, y = AUDPCpycnidia)) +
    geom_smooth(method="lm", aes(colour=Treatment), se=F, fullrange=T, level=0.9 ) +
    theme_classic() +
    theme(text = element_text(size=20, colour="black"), 
          axis.text.x = element_text(colour="black", size=15),
          axis.text.y = element_text(colour="black", size=15))
  ggsave(paste("~/Documents/Sobia/Watkins/Raw figures/",i, ".pdf", sep=""))
}


