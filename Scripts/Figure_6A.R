#Author: Ziyan Lin
#Description: this is the script for plotting Figure 6A
######################
#library package
library(phyloseq)
library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)
library(plotly)
library("RColorBrewer")
library(dplyr)
library(rstatix)
library(tidyr)
library(reshape)

options(scipen = 999) #That would make all the numbers to appear as decimals, convert scientific

######################
#input data
V1_msd <- read.csv("data/V1_101_clean.csv",row.names=1)
V1_msd$Visit <- "V1"
V2_msd <- read.csv("data/V2_101_clean.csv",row.names=1)
V2_msd$Visit <- "V2"

patient_list_S <- readLines("data/V1-S.txt")
patient_list_M <- readLines("data/V1-M.txt")

df <- rbind(V1_msd,V2_msd)
df$Perio <- "M"
df$Perio[which(df$SubjectID %in% patient_list_S)] <- "S"

tbl <- melt(df)
tbl <- tbl[order(tbl$variable,tbl$SubjectID),] #order for paired wilcoxon
tbl$variable <- gsub("\\.","-",tbl$variable)
######################
#Visit compare
SampleID = "SampleID"
group = "Visit"
pvalue_grp = "Visit"
paired_grp= "SubjectID"
wrap = "Group"
add_pvalue = TRUE
wilcox_paired = TRUE
box_width=0.8
x_font=12
y_font= 10
y_title_font=18
p_style = "p.signif"

MSDDodge <- function(tbl,variable="IL-6",color=c("black","gray"),color2 = c("firebrick3","goldenrod2","dodgerblue2"),group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE, paired_grp="SubjectID",add_pvalue = FALSE,dot_size=3, p_style ="p.signif",...) {
    
    tbl = tbl
    variable <- variable
    myColor = color
    myColor2 = color2
    group = group
    wrap = wrap
    x_font=x_font
    y_font=y_font
    title_font =title_font
    wilcox_paired = wilcox_paired
    add_pvalue = add_pvalue
    p_style = p_style
    paired_grp=paired_grp
    dot_size = dot_size
    
    subset_df <- tbl[tbl$variable==variable,]
    
    subset_df$Label <- as.vector(as.matrix(subset_df)[,group])
    subset_df$Wrap <- as.vector(as.matrix(subset_df)[,wrap])
    subset_df$Paired <- as.vector(as.matrix(subset_df)[,paired_grp])
    subset_df$Color_label <- as.vector(as.matrix(subset_df)[,"Perio"])

    subset_df <- subset_df[order(subset_df$variable,subset_df$Paired),] #order for paired wilcoxon

    #genere box plot
    set.seed(123)
    p<-NULL
    p <- ggplot(subset_df,aes(x=Label, y=value))+ #,color=Wrap,
        geom_boxplot(aes(color=Wrap,fill=Label),position=position_dodge(width=0.9),outlier.shape = NA) +
        geom_point(aes(color=Color_label, group = Wrap),position = position_jitterdodge(dodge.width=0.9),size=dot_size) +
        theme_classic() +
        ggtitle(variable) +
        xlab("")+
        ylab("Concentration (pg/mL)")+
        theme(plot.title = element_text(face = "bold.italic",hjust = 0.5,size=title_font),
            axis.text.x=element_text(angle=0,size=x_font,color="black"),
            axis.text.y=element_text(size=y_font,color="black"))

        p <- p + scale_color_manual(values=myColor) +
                scale_fill_manual(values=myColor2,
                                    name=group)
                            
        if (add_pvalue == TRUE){
            f <- as.formula(paste0("value ~ ",wrap))
            stat.test <- subset_df %>%
                group_by(Label) %>%
                wilcox_test(f,paired = wilcox_paired) %>%
                add_significance("p")

                # Add p-values onto the box plots
            stat.test <- stat.test %>%
                add_xy_position(x = "Label", dodge = 0.8)
                
          # stat.test$y.position[2:3] <- as.numeric(c("20","20"))
           # stat.test$y.position <- as.numeric(c("3","3","3"))

        p <- p + stat_pvalue_manual(stat.test,label = p_style, hide.ns = TRUE)#tip.length = 0.002/0.005
            return(p)
        
        } else {
            return(p)
            
        }
        
}

myColor <- c("gray","black","black","black")
myColor2 <- c("firebrick3","goldenrod2","dodgerblue2")

p1 <- MSDDodge(tbl= tbl,variable = "IFN-γ",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)

p2  <- MSDDodge(tbl= tbl,variable = "IL-10",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)
#p2 <- p2 +  scale_y_continuous(limit= c(0,3.2))

p3 <- MSDDodge(tbl= tbl,variable = "IL-12p70",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)

p4 <- MSDDodge(tbl= tbl,variable = "IL-13",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)

p5 <- MSDDodge(tbl= tbl,variable = "IL-1β",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)

p6 <- MSDDodge(tbl= tbl,variable = "IL-2",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)

p7 <- MSDDodge(tbl= tbl,variable = "IL-4",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)

p8 <- MSDDodge(tbl= tbl,variable = "IL-6",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)
#p8 <- p8 +  scale_y_continuous(limit= c(0,20))

p9 <- MSDDodge(tbl= tbl,variable = "IL-8",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)

p10 <- MSDDodge(tbl= tbl,variable ="TNF-α",color=myColor,color2=myColor2,group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,wilcox_paired = TRUE,add_pvalue = add_pvalue,p_style="p.signif",paired_grp=paired_grp,dot_size=0.5)


library(gridExtra)
library(ggpubr)
library(Cairo)

png(paste0("MSD-All.png"),height=7,width=10,unit="in",res=680)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,ncol=5, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

CairoPDF(paste0("MSD-All.pdf"),height=8,width=12)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,ncol=5, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()



