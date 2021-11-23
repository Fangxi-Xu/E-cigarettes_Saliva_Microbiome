Saliva_Heatmap_tbl <- function(phyloseqObject,n_top = 20,group= "Group",taxa_level = "Genus"){

    #load parameter
    physeq = phyloseqObject
    n_top = n_top
    group = group
    taxa_level = taxa_level
    
    #library package
    library(phyloseq)
    library(pheatmap)
    library(dplyr)
    
    options(scipen = 999)
    
    physeq_norm  = transform_sample_counts(physeq, function(x) x / sum(x))

    phy_taxo <- tax_glom(physeq_norm,taxa_level) #add up the otu proportion based on selected taxa (group otu)
    top_taxa <- prune_taxa(names(sort(taxa_sums(phy_taxo),TRUE)[1:n_top]), phy_taxo)
    
    df <- psmelt(top_taxa)
    
    df$myName <- as.vector(as.matrix(df)[,taxa_level])
    df$Label <- as.vector(as.matrix(df)[,group])
    
    df_summary <- df %>%
                  group_by(myName,Label) %>%
                  summarise(mean_Abun = mean(Abundance))
                
    ######################
    #formatting the plot data
    tbl <- matrix(ncol=length(unique(df_summary$Label)),nrow=length(unique(df_summary$myName)))
    df_summary <- df_summary[order(df_summary$Label,df_summary$myName),]
    
    colnames(tbl) <- unique(df_summary$Label)
    rownames(tbl) <- unique(df_summary$myName)
    
    tbl <- tbl[order(rownames(tbl)),order(colnames(tbl))]

    tbl[rownames(tbl),colnames(tbl)] <- df_summary$mean_Abun #assign the Abundance value
    
    return(as.matrix(tbl))

}
