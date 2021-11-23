#Author: Ziyan Lin
#Description: this is a script of genearting the subset (Severe and Mild/Moderate) of phyloseq object (Genus will full level names)
######################
#r/4.0.3
#library package
library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)
library(plotly)
library(dplyr)

options(scipen = 999) #That would make all the numbers to appear as decimals, convert scientific
######################
#input data
######################
load("my.physeq.Robj")

sam <- subset(physeq@sam_data,Visit=="V1")
all_patient <- sam$SubjectID

sub_sam <- subset(physeq@sam_data,Visit=="V1"& Perio == "Severe")
patient <- sub_sam$SubjectID
Sphyseq <- subset_samples(physeq,SubjectID %in% patient)

save(Sphyseq,file="my.Sphyseq.Robj")

patient2 <- all_patient[-which(all_patient %in% patient)]
Mphyseq <- subset_samples(physeq,SubjectID %in% patient2)
save(Mphyseq,file="my.Mphyseq.Robj")


