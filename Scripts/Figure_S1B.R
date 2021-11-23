#Author: Ziyan Lin
#Description: this is the script for plotting Figure S1B
######################
#library package
library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)
library(plotly)
library("RColorBrewer")
library(dplyr)
library(scatterplot3d)

options(scipen = 999) #That would make all the numbers to appear as decimals, convert scientific

source("Functions/Saliva-PCoA-Function.R")

###########################################################
#First Line 
load("my.Sphyseq.Robj")
physeq <- Sphyseq
############################################
#PCoA
############################################
group = "Visit"
d_method = "wunifrac"
add_pvalue = TRUE
px = 0.0025
py = 0.004

group_list <- c("CS","ES","NS")
myColor <- c("tomato","mediumpurple")

for (i in 1:length(group_list)){
    grp <- group_list[i]
    sub_physeq <- subset_samples(physeq,Group==grp) #select interested group
    
    p <- Saliva_Beta(sub_physeq,color= myColor,group = group,d_method = d_method,dot_size=3,x_font=13,y_font=13,x_title_font=15,y_title_font=15,add_pvalue=TRUE,px = px,py = py,return_dist = FALSE,return_coord = FALSE,return_pvalue = FALSE)

    png(paste0("PCoA-",d_method,"-",group,"-",grp,".png"),height=5,width=6,unit="in",res=680)
    print(p)
    dev.off()

    pdf(paste0("PCoA-",d_method,"-",group,"-",grp,".pdf"),height=5,width=6)
    print(p)
    dev.off()

    #output PCoA p-value
    pval <- Saliva_Beta(sub_physeq,group = group,d_method = d_method,return_pvalue = TRUE)
    write.csv(pval,paste0("PCoA-",d_method,"-",group,"-",grp,"-pval.csv"))
}

###########################################################
#Second Line
load("my.Mphyseq.Robj")
physeq <- Mphyseq
############################################
#PCoA
############################################
group = "Visit"
d_method = "wunifrac"
add_pvalue = TRUE
px = 0.0025
py = 0.004

group_list <- c("CS","ES","NS")
myColor <- c("tomato","mediumpurple")

for (i in 1:length(group_list)){
    grp <- group_list[i]
    sub_physeq <- subset_samples(physeq,Group==grp) #select interested group
    
    p <- Saliva_Beta(sub_physeq,color= myColor,group = group,d_method = d_method,dot_size=3,x_font=13,y_font=13,x_title_font=15,y_title_font=15,add_pvalue=TRUE,px = px,py = py,return_dist = FALSE,return_coord = FALSE,return_pvalue = FALSE)

    png(paste0("PCoA-",d_method,"-",group,"-",grp,".png"),height=5,width=6,unit="in",res=680)
    print(p)
    dev.off()

    pdf(paste0("PCoA-",d_method,"-",group,"-",grp,".pdf"),height=5,width=6)
    print(p)
    dev.off()

    #output PCoA p-value
    pval <- Saliva_Beta(sub_physeq,group = group,d_method = d_method,return_pvalue = TRUE)
    write.csv(pval,paste0("PCoA-",d_method,"-",group,"-",grp,"-pval.csv"))
}

