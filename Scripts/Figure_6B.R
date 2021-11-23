#Author: Ziyan Lin
#Description: this is the script for plotting Figure 6B
######################
#library package
library(phyloseq)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(Hmisc)
library(reshape2)
library(corrplot)
library(RColorBrewer)

options(scipen = 999) #That would make all the numbers to appear as decimals, convert scientific

load("my.physeq.Robj")
sam_data <- physeq@sam_data[order(physeq@sam_data$Visit,physeq@sam_data$SubjectID),]
######################
#input data
V1_msd <- read.csv("data/V1_101_clean.csv",row.names=1)
V1_msd <- V1_msd[,-ncol(V1_msd)]
V1_msd <- V1_msd[order(V1_msd$SubjectID),]
rownames(V1_msd) <- sam_data$SampleID[which(sam_data$Visit == "V1")]

V2_msd <- read.csv("data/V2_101_clean.csv",row.names=1)
V2_msd <- V2_msd[,-ncol(V2_msd)]
V2_msd <- V2_msd[order(V2_msd$SubjectID),]
rownames(V2_msd) <- sam_data$SampleID[which(sam_data$Visit == "V2")]

MSD <- rbind(V1_msd[,-1],V2_msd[,-1])
######################
n_top = 40

physeq_norm  = transform_sample_counts(physeq, function(x) x / sum(x))
phy_taxo <- tax_glom(physeq_norm,"Genus")
top_taxa <- prune_taxa(names(sort(taxa_sums(phy_taxo),TRUE)[1:n_top]), phy_taxo)
df <- psmelt(top_taxa)

df_summary <- df %>%
              group_by(Genus,SampleID) %>%
              summarise(mean_Abun = mean(Abundance))
df_summary$Genus <- gsub(".*\\|g__","",df_summary$Genus)

tbl <- dcast(df_summary,SampleID ~ Genus)
tbl <- tbl[order(tbl$SampleID),]
rownames(tbl) <- tbl$SampleID
tbl <- tbl[,-1]

MSD <- MSD[order(rownames(MSD)),]

mat <- cbind(MSD,tbl)
res <- rcorr(as.matrix(mat),type="spearman")
cor <- res$r
cor <- cor[colnames(tbl),colnames(MSD)]
colnames(cor) <- gsub("\\.","-",colnames(cor))

cell_colors = rev(brewer.pal(11,"Spectral"))

show_rowname = TRUE #TRUE/FALSE
fontsize_row = 12 #row font
fontsize_col = 12 #colname font
cluster_rows = TRUE #TRUE/FALSE
cluster_cols = FALSE #TRUE/FALSE
scale = "none" #"none","row","column"

png(paste0("Cor-",n_top,".png"),height=6,width=6,unit="in",res=680)
pheatmap(cor, color = cell_colors, border_color = NA, scale = scale, cluster_rows =cluster_rows, cluster_cols = cluster_cols, main = "Heatmap", fontsize_row = fontsize_row, fontsize_col = fontsize_col, show_rownames = show_rowname)
dev.off()


library(Cairo)
CairoPDF(paste0("Cor-",n_top,".pdf"),height=6,width=6)
pheatmap(cor, color = cell_colors, border_color = NA, scale = scale, cluster_rows =cluster_rows, cluster_cols = cluster_cols, main = "Heatmap", fontsize_row = fontsize_row, fontsize_col = fontsize_col, show_rownames = show_rowname)
dev.off()
