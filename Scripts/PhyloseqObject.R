#Author: Ziyan Lin
#Description: this is a script of genearting the phyloseq object
######################
#r/4.0.3
#library package
library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(ggsignif)
library(ggpubr)
library(dplyr)
######################
#input data
######################
taxo <- read.csv("data_qiime2/taxo.csv",header=T,row.names=1)
#rarefied feature table
biom_tbl <- read_biom("data_qiime2/feature-table.biom")
feature_tbl <- as.data.frame(as.matrix(biom_data(biom_tbl)))
#metadata
meta <- import_qiime_sample_data("data_qiime2/metadata.tsv")

OTU <- otu_table(feature_tbl,taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(taxo))
samples <- sample_data(meta)
phylodata <- phyloseq(OTU,TAX,samples)

#import tree
rooted_tree <- read_tree("data_qiime2/rooted_tree.nwk")
physeq <- merge_phyloseq(phylodata,rooted_tree)

save(physeq,file="my.physeq.Robj")
