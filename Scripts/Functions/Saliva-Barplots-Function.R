Saliva_Barplot <- function(phyloseqObject,taxa_level="Phylum",color= rainbow(50),min_pro=0.1,group= "Group",x_font=12,y_font= 10,y_title_font=18,return_data=FALSE){
    #load parameter
    physeq = phyloseqObject
    taxa_level = taxa_level
    myColor = color
    min_pro = min_pro
    group = group
    x_font = x_font
    y_font = y_font
    y_title_font = y_title_font
    return_data=return_data
    
    rm_name = paste0(tolower(strsplit(taxa_level,"")[[1]][1]),"__")
    
    #library package
    library(phyloseq)
    library(biomformat)
    library(ggplot2)
    library(vegan)
    library(ggsignif)
    library(ggpubr)
    #library(ggrepel)
    #library(cowplot)
    library(plotly)
    library("RColorBrewer")
    library(dplyr)
    library(picante)#calculate PD (alpha diversity)
    library(scatterplot3d)
    library("PairedData")

    options(scipen = 999)
    
    #Get Percentage
    physeq_norm  = transform_sample_counts(physeq, function(x) x / sum(x))
    
    phy_taxo <- tax_glom(physeq_norm,taxa_level) #add up the otu proportion based on selected taxo
    
    df <- psmelt(phy_taxo)
    df$Abundance <- df$Abundance*100 #change to percentage(%)

    my_group <- c(taxa_level,group)
    taxo_filter <- df %>%
               group_by_at(my_group) %>%  # the grouping variable
               summarise(mean_Abun = mean(Abundance),
               n_Abun = n()) %>%
               filter(mean_Abun < min_pro) %>% #obtain phlyums Others
               summarise(num = sum(n_Abun), .groups = 'drop') %>%
               filter(num >= nrow(sample_data(physeq_norm))) #remove phylum min_pro (default=0.1%) in all conditions
               
    df[,taxa_level][which(df[,taxa_level] %in% as.matrix(taxo_filter[,taxa_level]))] <- "Others" #assign Others Taxo to others

    df.stats <- df %>%
                group_by_at(taxa_level) %>% summarise(mean_Abun = mean(Abundance))

    taxo_order_tbl <- df.stats[order(df.stats$mean_Abun,decreasing=T),]
    taxo_order <- as.vector(as.matrix(taxo_order_tbl)[,taxa_level])
    taxo_name <- c(gsub(rm_name,"",taxo_order[-length(taxo_order)]),"Others") #formatted name for Taxo
    names(taxo_name) <- taxo_order
    
    #calculate standard error
    df_summary <- df %>%
                  group_by_at(my_group) %>%  # the grouping variable
                  summarise(mean_Abun = mean(Abundance),# calculates the mean of each group
                  sd_Abun = sd(Abundance), # calculates the standard deviation of each group
                  n_Abun = n(),  # calculates the sample size per group
                  SE_Abun = sd(Abundance)/sqrt(n()))# calculates the standard error of each group
 
    #output data in csv format
    if (return_data== TRUE){
    return(as.data.frame(df_summary))} else{
 
 
    df_summary$TAXO <- as.vector(as.matrix(df_summary)[,taxa_level])
    df_summary$TAXO <- factor(df_summary$TAXO,levels=taxo_order,labels=taxo_name) # set facet order
    
    df_summary$Label <- as.vector(as.matrix(df_summary)[,group])

    #generate barplot
    p<-NULL
    p <- ggplot(df_summary, aes(fill=TAXO, x=Label, y=mean_Abun)) +
                  facet_wrap(~TAXO,scales = "free") +
    geom_bar(aes(color=TAXO,fill=TAXO),stat="identity",color="black",width=0.55) +
                  ylab("Relative Abundance(%)") +
                  xlab("") +
                  geom_errorbar(aes(ymin=mean_Abun,ymax=mean_Abun+SE_Abun),width=.2) +#one side
                  theme_classic()+
                  theme(axis.text.x=element_text(size=x_font,color="black"),
                        axis.text.y=element_text(size=y_font,color="black"))+
                  theme(axis.title.y=element_text(size=y_title_font,color="black"))+
                  theme(strip.background = element_rect(fill=NULL,linetype="blank"), #facet background
                        strip.text.x = element_text(size=10,face = "bold.italic")) #facet font size
    # theme_classic()

    p <- p + scale_fill_manual(values=myColor,
                               breaks=taxo_name,
                               name = taxa_level)
                               
    return(p)
    }

    
}


Saliva_Barplot_Pvalue <- function(phyloseqObject,taxa_level="Phylum",min_pro=0.1,group= "Group",pvalue_grp = "Group", pairwise = FALSE,paired_grp="SubjectID", wilcox_paired = FALSE){
    #load parameter
    physeq = phyloseqObject
    taxa_level = taxa_level
    min_pro = min_pro
    group = group
    pvalue_grp = pvalue_grp
    pairwise = pairwise
    paired_grp = paired_grp
    wilcox_paired = wilcox_paired
    
    rm_name = paste0(tolower(strsplit(taxa_level,"")[[1]][1]),"__")
    
    #library package
    library(phyloseq)
    library(vegan)
    library("RColorBrewer")
    library(dplyr)
    library("PairedData")

    options(scipen = 999)
    
    #Get Percentage
    physeq_norm  = transform_sample_counts(physeq, function(x) x / sum(x))
    
    phy_taxo <- tax_glom(physeq_norm,taxa_level) #add up the otu proportion based on selected taxo
    
    df <- psmelt(phy_taxo)
    df$Abundance <- df$Abundance*100 #change to percentage(%)

    my_group <- c(taxa_level,group)
    taxo_filter <- df %>%
               group_by_at(my_group) %>%  # the grouping variable
               summarise(mean_Abun = mean(Abundance),
               n_Abun = n()) %>%
               filter(mean_Abun < min_pro) %>% #obtain phlyums Others
               summarise(num = sum(n_Abun), .groups = 'drop') %>%
               filter(num >= nrow(sample_data(physeq_norm))) #remove phylum min_pro (default=0.1%) in all conditions
               
    df[,taxa_level][which(df[,taxa_level] %in% as.matrix(taxo_filter[,taxa_level]))] <- "Others" #assign Others Taxo to others
    
    #order table for paired wilcoxon test
    df$Paired <- as.vector(as.matrix(df)[,paired_grp])
    df$Label <- as.vector(as.matrix(df)[,group])

    df <- df[order(df$Label),]
    df <- df[order(df$Paired,df$OTU),]

    f <- as.formula(paste0("Abundance ~ ",pvalue_grp))
    
    if (pairwise == TRUE){ #calculate p-value using wilcox (pairwise)
        pval_tbl_pairwise <- compare_means(formula=f,data=df,method="wilcox.test",group.by=taxa_level,paired=wilcox_paired)
        #output p-value(pairwise)
        return(as.data.frame(pval_tbl_pairwise))
        } else { #calculate p-value using Kruskal-Wallis
            
        pval_tbl <- compare_means(formula=f,data=df,method="kruskal.test",group.by=taxa_level,paired=wilcox_paired)
        #output p-value
        return(as.data.frame(pval_tbl))
    }

}


Saliva_Stacked_Barplot <- function(phyloseqObject,taxa_level="Phylum",color= rainbow(50),min_pro=0.1,group= "Group",wrap="Visit",x_font=12,y_font= 10,y_title_font=18,facet_font=15){
    #load parameter
    physeq = phyloseqObject
    taxa_level = taxa_level
    myColor = color
    min_pro = min_pro
    group = group
    x_font = x_font
    y_font = y_font
    y_title_font = y_title_font
    wrap = wrap
    facet_font = facet_font
    rm_name = paste0(tolower(strsplit(taxa_level,"")[[1]][1]),"__")
    
    #library package
    library(phyloseq)
    library(biomformat)
    library(ggplot2)
    library(vegan)
    library(ggsignif)
    library(ggpubr)
    #library(ggrepel)
    #library(cowplot)
    library(plotly)
    library("RColorBrewer")
    library(dplyr)
    library(picante)#calculate PD (alpha diversity)
    library(scatterplot3d)
    library("PairedData")

    options(scipen = 999)
    
    #Get Percentage
    physeq_norm  = transform_sample_counts(physeq, function(x) x / sum(x))
    
    phy_taxo <- tax_glom(physeq_norm,taxa_level) #add up the otu proportion based on selected taxo
    
    df <- psmelt(phy_taxo)
    df$Abundance <- df$Abundance*100 #change to percentage(%)

    my_group <- c(taxa_level,group)
    taxo_filter <- df %>%
               group_by_at(my_group) %>%  # the grouping variable
               summarise(mean_Abun = mean(Abundance),
               n_Abun = n()) %>%
               filter(mean_Abun < min_pro) %>% #obtain phlyums Others
               summarise(num = sum(n_Abun), .groups = 'drop') %>%
               filter(num >= nrow(sample_data(physeq_norm))) #remove phylum min_pro (default=0.1%) in all conditions
               
    df[,taxa_level][which(df[,taxa_level] %in% as.matrix(taxo_filter[,taxa_level]))] <- "Others" #assign Others Taxo to others
    
    df.stats <- df %>%
                group_by_at(taxa_level) %>% summarise(mean_Abun = mean(Abundance))

    taxo_order_tbl <- df.stats[order(df.stats$mean_Abun,decreasing=T),]
    taxo_order <- as.vector(as.matrix(taxo_order_tbl)[,taxa_level])
    taxo_name <- c(gsub(rm_name,"",taxo_order[-length(taxo_order)]),"Others") #formatted name for Taxo
    names(taxo_name) <- taxo_order
    
    my_group <- c(taxa_level,group,wrap)
    df_summary <- df %>%
                  group_by_at(my_group)  %>% # the grouping variable
                  summarise(mean_Abun = mean(Abundance))

    df_summary$TAXO <- as.vector(as.matrix(df_summary)[,taxa_level])
    df_summary$TAXO <- factor(df_summary$TAXO,levels=taxo_order,labels=taxo_name) # set facet order
    
    df_summary$Label <- as.vector(as.matrix(df_summary)[,group])
    
    
    #Get plot
    p <- NULL
    p <- ggplot(df_summary, aes(fill=TAXO, x=Label, y=mean_Abun)) +
                  facet_wrap(as.formula(paste0("~",wrap)),scales = "free") +
    geom_bar(aes(color=TAXO,fill=TAXO),stat="identity",color="black",width=0.65) +
                  ylab("Relative Abundance(%)") +
                  xlab("") +
                  theme_classic()+
                  theme(axis.text.x=element_text(angle=0,size=x_font,color="black"),
                        axis.text.y=element_text(color="black",size=y_font))+
                  theme(axis.title.y=element_text(size=y_title_font,color="black"))+
                  theme(strip.background = element_rect(fill=NULL,linetype="blank"), #facet background
                        strip.text.x = element_text(size=facet_font,face="bold")) #facet font size

    p <- p + scale_fill_manual(values=myColor,
                               breaks=taxo_name,
                               name = taxa_level)
                               
    return(p)
    
}


Saliva_Taxo_Barplot <- function(phyloseqObject,taxa_level="Phylum",TAXA,min_pro=0.1,group= "Group",wrap="Visit",x_font=12,y_font= 10,y_title_font=18,title_font=16){
    #load parameter
    physeq = phyloseqObject
    taxa_level = taxa_level
    TAXA = TAXA
    min_pro = min_pro
    group = group
    x_font = x_font
    y_font = y_font
    y_title_font = y_title_font
    title_font = title_font
    wrap = wrap
    rm_name = paste0(tolower(strsplit(taxa_level,"")[[1]][1]),"__")
    
    #library package
    library(phyloseq)
    library(biomformat)
    library(ggplot2)
    library(vegan)
    library(ggsignif)
    library(ggpubr)
    #library(ggrepel)
    #library(cowplot)
    library(plotly)
    library("RColorBrewer")
    library(dplyr)
    library(picante)#calculate PD (alpha diversity)
    library(scatterplot3d)
    library("PairedData")

    options(scipen = 999)
    
    #Get Percentage
    physeq_norm  = transform_sample_counts(physeq, function(x) x / sum(x))
    
    phy_taxo <- tax_glom(physeq_norm,taxa_level) #add up the otu proportion based on selected taxo
    
    df <- psmelt(phy_taxo)
    df$Abundance <- df$Abundance*100 #change to percentage(%)

    my_group <- c(taxa_level,group)
    taxo_filter <- df %>%
               group_by_at(my_group) %>%  # the grouping variable
               summarise(mean_Abun = mean(Abundance),
               n_Abun = n()) %>%
               filter(mean_Abun < min_pro) %>% #obtain phlyums Others
               summarise(num = sum(n_Abun), .groups = 'drop') %>%
               filter(num >= nrow(sample_data(physeq_norm))) #remove phylum min_pro (default=0.1%) in all conditions
               
    df[,taxa_level][which(df[,taxa_level] %in% as.matrix(taxo_filter[,taxa_level]))] <- "Others" #assign Others Taxo to others
    
    df.stats <- df %>%
                group_by_at(taxa_level) %>% summarise(mean_Abun = mean(Abundance))

    taxo_order_tbl <- df.stats[order(df.stats$mean_Abun,decreasing=T),]
    taxo_order <- as.vector(as.matrix(taxo_order_tbl)[,taxa_level])
    taxo_name <- c(gsub(rm_name,"",taxo_order[-length(taxo_order)]),"Others") #formatted name for Taxo
    names(taxo_name) <- taxo_order
    
    my_group <- c(taxa_level,group,wrap)
    df_summary <- df %>%
                  group_by_at(my_group)  %>% # the grouping variable
                  summarise(mean_Abun = mean(Abundance),
                  sd_Abun = sd(Abundance), # calculates the standard deviation of each group
                  n_Abun = n(),  # calculates the sample size per group
                  SE_Abun = sd(Abundance)/sqrt(n()))# calculates the standard error of each group)

    df_summary$TAXO <- as.vector(as.matrix(df_summary)[,taxa_level])
    df_summary$TAXO <- factor(df_summary$TAXO,levels=taxo_order,labels=taxo_name) # set facet order
    
    df_summary$Label <- as.vector(as.matrix(df_summary)[,group])
    df_summary$Wrap <- as.vector(as.matrix(df_summary)[,wrap])

    df_taxo <- df_summary %>% filter(TAXO == TAXA)

    #Get plot
    p <- NULL
    p <- ggplot(df_taxo, aes(fill=Wrap, x=Label, y=mean_Abun)) +
    geom_bar(stat="identity",color="black",width=0.65, position=position_dodge()) +
                  ylab("Relative Abundance(%)") +
                  xlab("") +
                  ggtitle(TAXA) + 
                  geom_errorbar(aes(ymin=mean_Abun,ymax=mean_Abun+SE_Abun),width=.2,position=position_dodge(.65)) +#one side
                  theme_classic()+
                  theme(axis.text.x=element_text(size=x_font,color="black"),
                        axis.text.y=element_text(color="black",size=y_font))+
                  theme(axis.title.y=element_text(size=y_title_font,color="black"),
                        plot.title = element_text(face = "bold.italic",size=title_font))

    p <- p + scale_fill_manual(values=c('#999999','#E69F00'),
                               name=paste(wrap))
    p <- p + theme(plot.title = element_text(hjust = 0.5)) #center title

                               
    return(p)
    
}
