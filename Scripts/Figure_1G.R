#Author: Ziyan Lin
#Description: this is the script for plotting Figure 2C
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
library(picante)#calculate PD (alpha diversity)
library(scatterplot3d)
library(rstatix)
library(tidyr)

options(scipen = 999) #That would make all the numbers to appear as decimals, convert scientific

load("my.physeq.Robj")
source("Functions/Saliva-AlphaDiversity-Functions.R")

alpha.theme <- theme(axis.title.y=element_text(size=16, color = "black"),
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y=element_text(size=14, color="black"),
                     strip.text = element_text(size = 12,  color="black"),
                     plot.title = element_text(hjust = 0.5, size=12, face = "bold"),
                     legend.position="top",
                     legend.text=element_text(size=16, color="black"),
                     legend.title = element_blank(),
                     plot.background = element_rect(fill = "transparent",colour = NA), )
######################
#Alpha Diversity
######################
SampleID = "SampleID"
group = "Visit"
pvalue_grp = "Visit"
paired_grp= "SubjectID"
wrap = "Group"
measures=c("Shannon","Observed")
PD = TRUE
facet_order=c("Observed","Shannon","PD")
add_pvalue = TRUE
t_paired = TRUE

variable.labs <- c("Observed ASVs","Shannon","PD")
names(variable.labs) <- facet_order

group_list <- c("CS","ES","NS")
myColor <- c("gray","black")
myColor2 <- c("firebrick3","goldenrod2","dodgerblue2")

color_label <- "V1Perio"
#######################################
#add variable
sam <- subset(physeq@sam_data,Visit=="V1")
all_patient <- sam$SubjectID
sub_sam <- subset(physeq@sam_data,Visit=="V1"& Perio == "Severe")
patient <- sub_sam$SubjectID

physeq@sam_data$V1Perio <- "M"
physeq@sam_data$V1Perio[which(physeq@sam_data$SubjectID %in% patient)] <- "S"
#######################################

for (i in 1:length(group_list)){
    
    grp <- group_list[i]
    fill_col <- myColor2[i]
    sub_physeq <- subset_samples(physeq,Group==grp) #select interested group

    #Get alpha diversity box plot
    p <- Saliva_Alpha_Sub(sub_physeq,color= myColor,SampleID = SampleID,group= group,wrap = wrap,measures=measures,PD = PD,facet_order=facet_order,x_font=12,y_font= 10,y_title_font=18,return_data=FALSE,pvalue_grp = pvalue_grp, pvalue_data = FALSE, paired_grp=paired_grp,t_paired = t_paired,add_pvalue=add_pvalue,p_style = "p.signif",box_width=0.6,fill_col=fill_col,color_label=color_label)
    
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + alpha.theme
    
    p <- p + facet_wrap(~variable,,scales = "free",labeller = labeller(variable = variable.labs))


    png(paste0("Alp-",group,"-",grp,".png"),height=4,width=5,unit="in",res=680)
    print(p)
    dev.off()

    pdf(paste0("Alp-",group,"-",grp,".pdf"),height=4,width=5)
    print(p)
    dev.off()

    #output alpha diversity info
    alp_data <- Saliva_Alpha(sub_physeq,return_data=TRUE,SampleID = SampleID,group= group,wrap = wrap,measures=measures,PD = PD,facet_order=facet_order,pvalue_grp = pvalue_grp, pvalue_data = FALSE, paired_grp=paired_grp,t_paired = t_paired)

    write.csv(alp_data,paste0("Alp-",group,"-",grp,"-data.csv"))

    #output alpha p-value
    pval_tbl <- Saliva_Alpha(sub_physeq,pvalue_data = TRUE,return_data=FALSE,SampleID = SampleID,group= group,wrap = wrap,measures=measures,PD = PD,facet_order=facet_order,pvalue_grp = pvalue_grp, paired_grp=paired_grp,t_paired = t_paired)

    write.csv(pval_tbl,paste0("Alp-",group,"-",grp,"-pval.csv"))
}


