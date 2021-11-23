#Author: Ziyan Lin
#Description: this is the script for plotting Figure 4A

library(phyloseq)
library(pheatmap)
library(dplyr)

source("../Functions/Saliva-Heatmap-Functions.R")
load("my.Sphyseq.Robj")

physeq <- Sphyseq
######################
#Heatmap (by Condition)
######################
physeq@sam_data$Condition <- paste0(physeq@sam_data$Group,"|",physeq@sam_data$Visit)
n_top = 40
group = "Condition"

#Get Heatmap table
tbl <- Saliva_Heatmap_tbl(physeq, n_top = n_top,group= group)

#formatted genus names
rownames(tbl) <- gsub(".*\\|g__","",rownames(tbl))

#get heatmap
#set parameter
cell_colors = colorRampPalette(c("#043177", "#244B88", "#FAFAFA", "#C62E2E", "#BF0F0F"))(50)
show_rowname = TRUE #TRUE/FALSE
fontsize_row = 8 #row font
fontsize_col = 12 #colname font
cluster_rows = TRUE #TRUE/FALSE
cluster_cols = FALSE #TRUE/FALSE
scale = "row" #"none","row","column"

png(paste0("Heatmap-top",n_top,"-",group,".png"),width =6, height = 6,unit="in",res=680)
pheatmap(tbl, color = cell_colors, border_color = NA, scale = scale, cluster_rows =cluster_rows, cluster_cols = cluster_cols, main = "Heatmap", fontsize_row = fontsize_row, fontsize_col = fontsize_col, show_rownames = show_rowname)
dev.off()

pheatmap(tbl, color = cell_colors, border_color = NA, scale = scale, cluster_rows =cluster_rows, cluster_cols = cluster_cols, main = "Heatmap", fontsize_row = fontsize_row, fontsize_col = fontsize_col, show_rownames = show_rowname,width =6,height = 6,filename=paste0("Heatmap-top",n_top,"-",group,".pdf"))
