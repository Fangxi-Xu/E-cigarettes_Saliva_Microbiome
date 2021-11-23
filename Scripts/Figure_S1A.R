#Author: Ziyan Lin
#Description: this is the script for plotting Figure S1A
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

load("my.physeq.Robj")

patient_list_S <- readLines("V1-S.txt")

physeq@sam_data$PerioV1 <- "M"
physeq@sam_data$PerioV1[which(physeq@sam_data$SubjectID %in% patient_list_S)] <- "S"

group = "PerioV1"
d_method = "wunifrac"
add_pvalue = TRUE
px = 0.0025
py = 0.004

myColor <- c("darkolivegreen","indianred4")
group_list <- c("V1","V2")
group_list2 <- c("CS","ES","NS")

for (i in 1:length(group_list)){
    grp <- group_list[i]
    for (j in 1:length(group_list2)){
        grp2 <- group_list2[j]
        sub_physeq <- subset_samples(physeq,Visit==grp & Group==grp2) #select interested group

        p <- Saliva_Beta(sub_physeq,color= myColor,group = group,d_method = d_method,dot_size=3,x_font=13,y_font=13,x_title_font=15,y_title_font=15,add_pvalue=TRUE,px = px,py = py,return_dist = FALSE,return_coord = FALSE,return_pvalue = FALSE)

        png(paste0("PCoA-",d_method,"-",group,"-",grp,"-",grp2,".png"),height=5,width=6,unit="in",res=680)
        print(p)
        dev.off()

        pdf(paste0("PCoA-",d_method,"-",group,"-",grp,"-",grp2,".pdf"),height=5,width=6)
        print(p)
        dev.off()

        #output PCoA p-value
        pval <- Saliva_Beta(sub_physeq,group = group,d_method = d_method,return_pvalue = TRUE)
        write.csv(pval,paste0("PCoA-",d_method,"-",group,"-",grp,"-",grp2,"-pval.csv"))
    }
}

