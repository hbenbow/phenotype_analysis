library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(multcompView)
library(readxl)
library(VennDiagram)


# TIlling mutations

mutations<-read.csv("~/Documents/Sobia/Watkins/Data/TILLING/selectedMutations.csv")
table<-spread(as.data.frame(table(mutations$gene, mutations$line)), key = "Var2", value="Freq")
table<-table[2:9693,]
hist(table$Cadenza0449)
hist(table$Cadenza0451)
commons<-table[(table$Cadenza0449 >=1),]
commons<-commons[(commons$Cadenza0451 >=1),]

interesting_consequences<-c(
"missense_variant",
"missense_variant&splice_region_variant",     
"stop_gained",
"splice_donor_variant&non_coding_transcript_variant",
"splice_region_variant&intron_variant",                                            
"splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant",
"splice_donor_variant",                                                       
"splice_region_variant&synonymous_variant",                                         
"splice_acceptor_variant",                                        
"start_lost",                                             
"splice_region_variant&5_prime_UTR_variant",                                           
"stop_gained&splice_region_variant",                                                
"splice_region_variant&intron_variant&non_coding_transcript_variant",               
"start_lost&splice_region_variant")

interesting_mutations<-subset(mutations, mutations$consequence %in% interesting_consequences)
interesting_table<-spread(as.data.frame(table(interesting_mutations$gene, interesting_mutations$line)), key = "Var2", value="Freq")
i_commons<-interesting_table[(interesting_table$Cadenza0449 >=1),]
i_commons<-i_commons[(i_commons$Cadenza0451 >=1),]
