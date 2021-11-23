#Author: Ziyan Lin
#Description: this is the script for plotting Figure 2A-B
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

load("my.Sphyseq.Robj")
physeq <- Sphyseq
############################################
#PCoA
############################################
group = "Group"
d_method = "wunifrac"
add_pvalue = TRUE
px = 0.0025
py = 0.004

group_list <- c("V1","V2")
myColor <- c("firebrick3","goldenrod2","dodgerblue2")
for (i in 1:length(group_list)){
    grp <- group_list[i]
    sub_physeq <- subset_samples(physeq,Visit==grp) #select interested group
    
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


groups_combo <- list(c("CS","ES"),c("CS","NS"),c("ES","NS"))
Color_combo <- list(c("firebrick3","goldenrod2"),c("firebrick3","dodgerblue2"),c("goldenrod2","dodgerblue2"))

#Each Visit
visit_list <- c("V1","V2")
for (myVisit in visit_list){
    
    my_physeq <- subset_samples(physeq,Visit==myVisit)
    
    groups_combo <- groups_combo
    for(i in 1:length(groups_combo)){
        Color_combo <- Color_combo
        myColor <- Color_combo[[i]]
        group_list <- groups_combo[[i]]
        plot_title <- paste0(group_list[1],"_vs_",group_list[2])
        
        sub_physeq <- subset_samples(my_physeq,Group==group_list[1]|Group==group_list[2])
        
        p <- Saliva_Beta(sub_physeq,color= myColor,group = group,d_method = d_method,dot_size=3,x_font=13,y_font=13,x_title_font=15,y_title_font=15,add_pvalue=TRUE,px = px,py = py,return_dist = FALSE,return_coord = FALSE,return_pvalue = FALSE)

        png(paste0("PCoA-",d_method,"-",plot_title,"-",myVisit,".png"),height=5,width=6,unit="in",res=680)
        print(p)
        dev.off()

        pdf(paste0("PCoA-",d_method,"-",plot_title,"-",myVisit,".pdf"),height=5,width=6)
        print(p)
        dev.off()

        #output PCoA p-value
        pval <- Saliva_Beta(sub_physeq,group = group,d_method = d_method,return_pvalue = TRUE)
        write.csv(pval,paste0("PCoA-",d_method,"-",plot_title,"-",myVisit,"-pval.csv"))
    }

}
