#Author: Ziyan Lin
#Description: this is the script for plotting Figure 3A
library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)
library(plotly)
library("RColorBrewer")
library(dplyr)
library("PairedData")

source("Functions/Saliva-Barplots-Function.R")

load("my.Sphyseq.Robj")
physeq <- Sphyseq
#########################################
#Barplot
#########################################
taxa = "Phylum"
group = "Group"
wrap = "Visit"

#########################################
#Specofic Taxa Barplot
#########################################
#taxa_list
#[1] "Firmicutes"               "Proteobacteria"
#[3] "Spirochaetes"             "Others"
#[5] "Saccharibacteria_(TM7)"   "Absconditabacteria_(SR1)"
#[7] "Bacteroidetes"            "Actinobacteria"
#[9] "Fusobacteria"

taxa_list <- c("Firmicutes","Proteobacteria","Spirochaetes","Saccharibacteria_(TM7)","Absconditabacteria_(SR1)","Bacteroidetes","Actinobacteria","Fusobacteria","Others")

for (TAXA in taxa_list){
    TAXA = TAXA
p <- Saliva_Taxo_Barplot(physeq,TAXA=TAXA,taxa_level=taxa,min_pro=0.1,group= group,x_font=12,y_font= 10,y_title_font=18,title_font=16,wrap=wrap)

png(paste0("Barplot-",taxa,"-",TAXA,".png"),height=5,width=5,unit="in",res=380)
print(p)
dev.off()

pdf(paste0("Barplot-",taxa,"-",TAXA,".pdf"),height=5,width=5)
print(p)
dev.off()
}
