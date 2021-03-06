library(jtools)
library(interactions)
library(ggplot2)
library(ggpubr)


AUDPC_Pycnidia <-
  read.csv("~/Documents/Sobia/Watkins/Data/AUDPC_Pycnidia.csv")
AUDPC_Chlorosis <-
  read.csv("~/Documents/Sobia/Watkins/Data/AUDPC_necrosis.csv")
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
  geom_smooth(aes(colour = Variety))

ggline(chlor.and.pyc, x = "AUDPCnecrosis", y = "AUDPCpycnidia", color="Variety",
       add = c("mean_se", "jitter"),
       ylab = "Weight", xlab = "Treatment")

av_necrosis <- aov(AUDPCnecrosis ~   Replication + Plant.no + Treatment * Variety, data = chlor.and.pyc)
summary(av_necrosis)
av_pycnidia <- aov(AUDPCpycnidia ~   Replication + Plant.no + Treatment * Variety, data = chlor.and.pyc)
summary(av_pycnidia)
knit(summary(av_necrosis))
interact_plot(fiti, pred = AUDPCnecrosis, modx = Treatment)


for (i in chlor.and.pyc$Variety) {
  data <- chlor.and.pyc[(chlor.and.pyc$Variety == i), ]
  data <- data[!(data$Treatment == "Isolate average"), ]
  fiti <- lm(AUDPCpycnidia ~ AUDPCnecrosis *Treatment, data = data)
  assign(paste(i, "lm"), summary(fiti))
  interact_plot(fiti, pred = AUDPCnecrosis, modx = Treatment) + ggtitle(paste(i)) +
    theme(text = element_text(size = 20, colour = "black")) +
    theme_classic()
  ggsave(paste(i, ".pdf"))
}


list<-list()
for(i in unique(chlor.and.pyc$Variety)){
    data <- chlor.and.pyc[(chlor.and.pyc$Variety == i), ]
    for(t in unique(data$Treatment)){
    data2 <- data[(data$Treatment == t), ]
    cor<-as.data.frame(cor(data2$AUDPCnecrosis, data2$AUDPCpycnidia))
    cor$Treatment<-paste(t)
    cor$Genotype<-paste(i)
    list[[length(list)+1]]<-cor
  } 
}
all1<-as.data.frame(do.call(rbind, list))

