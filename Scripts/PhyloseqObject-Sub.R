#Author: Ziyan Lin
#Description: this is a script of genearting the subset of phyloseq object (Severe and Mild/Moderate separately)
######################
#r/4.0.3
#library package
library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)
library(dplyr)

######################
#input data
######################
load("my.physeq.Robj")

sam <- subset(physeq@sam_data,Visit=="V1")
all_patient <- sam$SubjectID

sub_sam <- subset(physeq@sam_data,Visit=="V1"& Perio == "Severe")
patient <- sub_sam$SubjectID
sub_physeq <- subset_samples(physeq,SubjectID %in% patient)

Sphyseq <- sub_physeq
save(Sphyseq,file="my.Sphyseq.Robj")

patient2 <- all_patient[-which(all_patient %in% patient)]
sub_physeq2 <- subset_samples(physeq,SubjectID %in% patient2)
Mphyseq <- sub_physeq2

save(Mphyseq,file="my.Mphyseq.Robj")
