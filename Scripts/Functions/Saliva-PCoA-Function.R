Saliva_Beta <- function(phyloseqObject,color= rainbow(50),group= "Group",d_method = "bray",dot_size =3,x_font=13,y_font=13,x_title_font=15,y_title_font=15,add_pvalue=TRUE,px = 0.45,py = 0.45,return_dist = FALSE,return_coord = FALSE,return_pvalue = FALSE){
    
    #load parameter
    physeq = phyloseqObject
    myColor = color
    group = group
    d_method = d_method #"wunifrac"/"unifrac"/"bray"
    add_pvalue = add_pvalue
    px = px
    py = py
    dot_size = dot_size
    x_font = x_font
    y_font = y_font
    x_title_font = x_title_font
    y_title_font = y_title_font
    return_dist = return_dist
    return_coord = return_coord
    return_pvalue = return_pvalue

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
    
    options(scipen = 999)
    
    #Get normalized object
    physeq_norm  = transform_sample_counts(physeq, function(x) x / sum(x))
    
    #calculate distance
    iDist <- distance(physeq_norm, method=d_method)
    
    #get coordinate
    iPCoA <- ordinate(physeq_norm, "PCoA", distance=iDist)
    
    sample_data(physeq_norm)$Label <- as.vector(as.matrix(sample_data(physeq_norm))[,group])
    
    #adonis(p-value)
    set.seed(123)
    stats<- adonis(formula = iDist ~ sample_data(physeq_norm)$Label)
    p_val<-stats$aov.tab$"Pr(>F)"[1]
    
    #generate PCoA plot
    p <- NULL
    p <- plot_ordination(physeq_norm, iPCoA,color=group) #shape=
    
    if (return_dist == TRUE){
        return(as.matrix(iDist))
    } else if(return_coord == TRUE){
        return(as.data.frame(p$data))
    } else if (return_pvalue == TRUE){
        return(as.matrix(stats$aov.tab))
    } else{
        pc1 <- gsub("Axis.1  ","PCo1",p$labels$x)
        pc2 <- gsub("Axis.2  ","PCo2",p$labels$y)

        df <- p$data
        p <- ggplot(df,aes(x=Axis.1,y=Axis.2,color=Label,fill=Label)) +
                 geom_point(size=dot_size) +
                 ggtitle(paste0("PCoA (", d_method,") --- ",group)) +
                 xlab(paste0(pc1)) +
                 ylab(paste0(pc2)) +
                 stat_ellipse(geom = "polygon",type = "norm",level = 0.95,size=1.2,alpha=0.1) + theme_bw()+ #stat_ellipse(type = "t")
                 theme_classic()+
                 theme(axis.text.x=element_text(angle=0,size=x_font,color="black"),
                       axis.text.y=element_text(size=y_font,color="black"),
                       axis.title.x=element_text(size=x_title_font),
                       axis.title.y=element_text(size=y_title_font))

        p <- p + scale_color_manual(values = myColor,
                                    name = group) +
                 scale_fill_manual(values = alpha(myColor),
                                    name  = group)
                                    
        if(add_pvalue == TRUE){
            p <- p + annotate("text", x = px, y = py, label = paste0("p = ",p_val),size=5,fontface = 'italic')
            return(p)
            
        } else {
            return(p)
        }
                                
    }
}
