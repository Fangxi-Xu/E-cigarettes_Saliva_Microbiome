Saliva_Alpha <- function(phyloseqObject,color= rainbow(50),SampleID = "SampleID",group= "Group",wrap = "Visit",measures=c("Shannon","Observed","Chao1"),PD = TRUE,facet_order=c("Observed","Shannon","Chao1","PD"),x_font=12,y_font= 10,y_title_font=18,return_data=FALSE,pvalue_grp = "Group", pvalue_data = FALSE, paired_grp="SubjectID",t_paired = FALSE,add_pvalue = FALSE,p_style = "p.signif",box_width=0.8){

    #load parameter
    physeq = phyloseqObject
    myColor = color
    measures = measures
    SampleID = SampleID
    group = group
    x_font = x_font
    y_font = y_font
    y_title_font = y_title_font
    return_data=return_data
    pvalue_grp = pvalue_grp
    pvalue_data = pvalue_data
    t_paired = t_paired
    PD = PD
    facet_order = facet_order
    paired_grp = paired_grp
    wrap = wrap
    add_pvalue = add_pvalue
    p_style = p_style
    box_width = box_width
    
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
    
    options(scipen = 999)
    
    #Calculate Alpha Diversity
    phy_alp  <- prune_taxa(taxa_sums(physeq) > 0, physeq) #prune OTUs that are not present in any of the samples
    
    p <- NULL
    p <- plot_richness(phy_alp,measures=measures)
    variables <- c(SampleID,group,wrap,paired_grp,"variable","value","se")
    df1 <- subset(p$data,select=variables)
    
    if (PD == TRUE){
        #Calculate PD
        phy_alp_otu <- as.data.frame(phy_alp@otu_table)
        phy_alp_tree <- phy_alp@phy_tree# it is a rooted tree
        df.pd.temp <- pd(t(phy_alp_otu), phy_alp_tree,include.root=T) # t(ou_table) transposes the table for use in picante
    
        #format PD result
        grps <- c(SampleID,group,wrap,paired_grp)
        df.pd <- subset(phy_alp@sam_data,select=grps)
        df.pd$variable <- "PD"
        df.pd$value <- df.pd.temp$PD
        df.pd$se <- "NA"

    #combine all method
    alp_data <- rbind(df1,df.pd)
    } else { alp_data <- df1 }
    
    alp_data$variable <- factor(alp_data$variable,levels=facet_order)#customize facet order
    
    if (return_data == TRUE){
        return(as.data.frame(alp_data))
    } else if (pvalue_data == TRUE) {
        alp_data$Paired <- as.vector(as.matrix(alp_data)[,paired_grp])
        alp_data <- alp_data[order(alp_data$variable,alp_data$Paired),] #order for paired wilcoxon
        f <- as.formula(paste0("value ~ ",group))
        pval_tbl <- compare_means(formula=f,data=alp_data,method="t.test",group.by="variable",paired = t_paired,p.adjust.method = "bonferroni")
        return(as.data.frame(pval_tbl))
    } else {
        
        alp_data$Paired <- as.vector(as.matrix(alp_data)[,paired_grp])
        alp_data <- alp_data[order(alp_data$variable,alp_data$Paired),] #order for paired wilcoxon

        alp_data$Label <- as.vector(as.matrix(alp_data)[,group])
        alp_data$Wrap <- as.vector(as.matrix(alp_data)[,wrap])

        
    }
    
    
        if(t_paired == TRUE){
            set.seed(123)
            p<-NULL
            p <- ggplot(alp_data,aes(x=Label, y=value))+
            facet_wrap(~variable,scales = "free")+
            geom_boxplot(width=box_width,outlier.shape = NA) +
            geom_point(aes(color=Label,group=Paired), position = position_dodge(0.2)) +
            geom_line(aes(group=Paired),color="gray",linetype = "dashed", position = position_dodge(0.2)) +
            theme_classic() +
            xlab("")+
            ylab("Alpha Diversity Measure")+
            theme(axis.text.x=element_text(angle=0,size=x_font,color="black"),
                  axis.text.y=element_text(size=y_font,color="black"))+
            theme(strip.background = element_rect(fill=NULL,linetype="blank"), #facet background
                  strip.text.x = element_text(size=10,face = "bold.italic")) #facet font size
            p <- p + scale_color_manual(values=myColor,
                                        name=group)
        }
        else {
            #genere box plot
            set.seed(123)
            p<-NULL
            p <- ggplot(alp_data,aes(x=Label, y=value))+
            facet_wrap(~variable,scales = "free")+
            geom_boxplot(width=box_width,outlier.shape = NA) +
            geom_jitter(aes(color=Label),width=.15) +
            theme_classic() +
            xlab("")+
            ylab("Alpha Diversity Measure")+
            theme(axis.text.x=element_text(angle=0,size=x_font,color="black"),
                axis.text.y=element_text(size=y_font,color="black"))+
            theme(strip.background = element_rect(fill=NULL,linetype="blank"), #facet background
                strip.text.x = element_text(size=10,face = "bold.italic")) #facet font size
        
            p <- p + scale_color_manual(values=myColor,
                                  name=group)
            }

        #prepare for p-value (t.test)
        lev <- levels(as.factor(alp_data$Label)) # get the variables
        L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])# make a pairwise list that we want to compare.

        if (add_pvalue == TRUE){
        #add p value(significant level)
        p <- p + stat_compare_means(
                    comparisons = L.pairs,
                    hide.ns = T,
                    label = p_style, #"p.signif",#"p.format"
                    method = "t.test",
                    paired = t_paired)
        return(p)
            
        } else {return(p)}
        
}


Saliva_Alpha_Sub <- function(phyloseqObject,color= rainbow(50),SampleID = "SampleID",group= "Group",wrap = "Visit",measures=c("Shannon","Observed","Chao1"),PD = TRUE,facet_order=c("Observed","Shannon","Chao1","PD"),x_font=12,y_font= 10,y_title_font=18,return_data=FALSE,pvalue_grp = "Group", pvalue_data = FALSE, paired_grp="SubjectID",t_paired = FALSE,add_pvalue = FALSE,p_style = "p.signif",box_width=0.8,fill_col="white",color_label){

    #load parameter
    physeq = phyloseqObject
    myColor = color
    measures = measures
    SampleID = SampleID
    group = group
    x_font = x_font
    y_font = y_font
    y_title_font = y_title_font
    return_data=return_data
    pvalue_grp = pvalue_grp
    pvalue_data = pvalue_data
    t_paired = t_paired
    PD = PD
    facet_order = facet_order
    paired_grp = paired_grp
    wrap = wrap
    add_pvalue = add_pvalue
    p_style = p_style
    box_width = box_width
    fill_col = fill_col
    color_label = color_label
    
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
    
    options(scipen = 999)
    
    #Calculate Alpha Diversity
    phy_alp  <- prune_taxa(taxa_sums(physeq) > 0, physeq) #prune OTUs that are not present in any of the samples
    
    p <- NULL
    p <- plot_richness(phy_alp,measures=measures)
    variables <- c(SampleID,group,wrap,paired_grp,color_label,"variable","value","se")
    df1 <- subset(p$data,select=variables)
    
    if (PD == TRUE){
        #Calculate PD
        phy_alp_otu <- as.data.frame(phy_alp@otu_table)
        phy_alp_tree <- phy_alp@phy_tree# it is a rooted tree
        df.pd.temp <- pd(t(phy_alp_otu), phy_alp_tree,include.root=T) # t(ou_table) transposes the table for use in picante
    
        #format PD result
        grps <- c(SampleID,group,wrap,paired_grp,color_label)
        df.pd <- subset(phy_alp@sam_data,select=grps)
        df.pd$variable <- "PD"
        df.pd$value <- df.pd.temp$PD
        df.pd$se <- "NA"

    #combine all method
    alp_data <- rbind(df1,df.pd)
    } else { alp_data <- df1 }
    
    alp_data$variable <- factor(alp_data$variable,levels=facet_order)#customize facet order
    
    if (return_data == TRUE){
        return(as.data.frame(alp_data))
    } else if (pvalue_data == TRUE) {
        alp_data$Paired <- as.vector(as.matrix(alp_data)[,paired_grp])
        alp_data <- alp_data[order(alp_data$variable,alp_data$Paired),] #order for paired wilcoxon
        f <- as.formula(paste0("value ~ ",group))
        pval_tbl <- compare_means(formula=f,data=alp_data,method="t.test",group.by="variable",paired = t_paired,p.adjust.method = "bonferroni")
        return(as.data.frame(pval_tbl))
    } else {
        
        alp_data$Paired <- as.vector(as.matrix(alp_data)[,paired_grp])
        alp_data <- alp_data[order(alp_data$variable,alp_data$Paired),] #order for paired wilcoxon

        alp_data$Label <- as.vector(as.matrix(alp_data)[,group])
        alp_data$Wrap <- as.vector(as.matrix(alp_data)[,wrap])
        alp_data$ColorLabel <- as.vector(as.matrix(alp_data)[,color_label])

    }
    
    
        if(t_paired == TRUE){
            set.seed(123)
            p<-NULL
            p <- ggplot(alp_data,aes(x=Label, y=value,fill=Wrap))+
            facet_wrap(~variable,scales = "free")+
            geom_boxplot(width=box_width,outlier.shape = NA) +
            geom_point(aes(color=ColorLabel,group=Paired), position = position_dodge(0.2)) +
            geom_line(aes(group=Paired),color="gray",linetype = "dashed", position = position_dodge(0.2)) +
            theme_classic() +
            xlab("")+
            ylab("Alpha Diversity Measure")+
            theme(axis.text.x=element_text(angle=0,size=x_font,color="black"),
                  axis.text.y=element_text(size=y_font,color="black"))+
            theme(strip.background = element_rect(fill=NULL,linetype="blank"), #facet background
                  strip.text.x = element_text(size=10,face = "bold.italic")) #facet font size
            p <- p + scale_color_manual(values=myColor,
                                        name=color_label) +
                    scale_fill_manual(values=fill_col,
                                     name=wrap)
                                        
        }
        else {
            #genere box plot
            set.seed(123)
            p<-NULL
            p <- ggplot(alp_data,aes(x=Label, y=value,fill=Wrap))+
            facet_wrap(~variable,scales = "free")+
            geom_boxplot(width=box_width,outlier.shape = NA) +
            geom_jitter(aes(color=ColorLabel),width=.15) +
            theme_classic() +
            xlab("")+
            ylab("Alpha Diversity Measure")+
            theme(axis.text.x=element_text(angle=0,size=x_font,color="black"),
                axis.text.y=element_text(size=y_font,color="black"))+
            theme(strip.background = element_rect(fill=NULL,linetype="blank"), #facet background
                strip.text.x = element_text(size=10,face = "bold.italic")) #facet font size
        
            p <- p + scale_color_manual(values=myColor,
                                  name=color_label) +
                    scale_fill_manual(values=fill_col,
                                      name=wrap)
                
            }

        #prepare for p-value (t.test)
        lev <- levels(as.factor(alp_data$Label)) # get the variables
        L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])# make a pairwise list that we want to compare.

        if (add_pvalue == TRUE){
        #add p value(significant level)
        p <- p + stat_compare_means(
                    comparisons = L.pairs,
                    hide.ns = T,
                    label = p_style, #"p.signif",#"p.format"
                    method = "t.test",
                    paired = t_paired)
        return(p)
            
        } else {return(p)}
        
}



AlphaDodge <- function(alp_data,variable="Observed",color=c("black","gray"),color2 = c("firebrick3","dodgerblue2","goldenrod2"),group="Group",wrap="Visit",x_font=12,y_font= 10,title_font=12,t_paired = FALSE, paired_grp="SubjectID",add_pvalue = FALSE,p_style ="p.signif",...) {
    
    alp_data = alp_data
    variable <- variable
    myColor = color
    myColor2 = color2
    group = group
    wrap = wrap
    x_font=x_font
    y_font=y_font
    title_font =title_font
    t_paired = t_paired
    add_pvalue = add_pvalue
    p_style = p_style
    paired_grp=paired_grp
    
    subset_df <- alp_data[alp_data$variable==variable,]
    
    subset_df$Label <- as.vector(as.matrix(subset_df)[,group])
    subset_df$Wrap <- as.vector(as.matrix(subset_df)[,wrap])
    subset_df$Paired <- as.vector(as.matrix(subset_df)[,paired_grp])
    
    subset_df <- subset_df[order(subset_df$variable,subset_df$Paired),] #order for paired wilcoxon

    #genere box plot
    set.seed(123)
    p<-NULL
    p <- ggplot(subset_df,aes(x=Label, y=value,color=Wrap,fill=Label))+
        geom_boxplot(position=position_dodge(width=0.9),outlier.shape = NA) +
        geom_point(position = position_jitterdodge(dodge.width=0.9)) +
        theme_classic() +
        ggtitle(variable) +
        xlab("")+
        ylab("Alpha Diversity Measure")+
        theme(plot.title = element_text(face = "bold.italic",hjust = 0.5,size=title_font),
            axis.text.x=element_text(angle=0,size=x_font,color="black"),
            axis.text.y=element_text(size=y_font,color="black"))

        p <- p + scale_color_manual(values=myColor,
                            name=wrap) +
                scale_fill_manual(values=myColor2,
                                    name=group)
                            
        if (add_pvalue == TRUE){
            f <- as.formula(paste0("value ~ ",wrap))
            stat.test <- subset_df %>%
                group_by(Label) %>%
                t_test(f,paired = t_paired) %>%
                add_significance("p")

                # Add p-values onto the box plots
            stat.test <- stat.test %>%
                add_xy_position(x = "Label", dodge = 0.8)

            p <- p + stat_pvalue_manual(stat.test,label = p_style, hide.ns = TRUE)#tip.length = 0
            return(p)
        
        } else {
            return(p)
            
        }
        
}

