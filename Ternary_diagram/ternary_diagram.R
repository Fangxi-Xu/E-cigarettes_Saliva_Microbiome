setwd("/Volumes/Samsung_T5/IECOH_SALv1v2/R-Fangxi")
# create plot folder
path <- "./plot"
dir.create(path)
#################################################################################
#use phylo_object.script to import data
#################################################################################
library(microbiomeutilities)
library(microbiome)
library(knitr)
library(tibble)
library(dplyr)
#101 subjects all together
phylo_relative = transform_sample_counts(phylo_clean, function(x) {(x/sum(x))*100} )
phylo_relative
phylo_relative_v1 <-subset_samples(phylo_relative, Visit == "V1")
phylo_relative_v2 <-subset_samples(phylo_relative, Visit == "V2")
#################################################################################
phylo.SR1.v1 = subset_taxa(phylo_relative_v1, Phylum == "p__Absconditabacteria_(SR1)")
phylo.SR1.v2 = subset_taxa(phylo_relative_v2, Phylum == "p__Absconditabacteria_(SR1)")

tern_df_SR1_v1 <- prep_ternary(phylo.SR1.v1, group="Group", level= "Species")
head(tern_df_SR1_v1)
tern_df_SR1_v2 <- prep_ternary(phylo.SR1.v2, group="Group", level= "Species")
head(tern_df_SR1_v2)
#################################################################################
phylo.TM7.v1 = subset_taxa(phylo_relative_v1, Phylum == "p__Saccharibacteria_(TM7)")
phylo.TM7.v2 = subset_taxa(phylo_relative_v2, Phylum == "p__Saccharibacteria_(TM7)")

tern_df_TM7_v1 <- prep_ternary(phylo.TM7.v1, group="Group", level= "Species")
view(tern_df_TM7_v1)
tern_df_TM7_v1 <- tern_df_TM7_v1[tern_df_TM7_v1$Phylum != "Unknown", ]       
tern_df_TM7_v2 <- prep_ternary(phylo.TM7.v2, group="Group", level= "Species")
view(tern_df_TM7_v2)
tern_df_TM7_v2 <- tern_df_TM7_v2[tern_df_TM7_v2$Phylum != "Unknown", ]  
#################################################################################
phylo.bacteroidetes.v1 = subset_taxa(phylo_relative_v1, Phylum == "p__Bacteroidetes")
phylo.bacteroidetes.v2 = subset_taxa(phylo_relative_v2, Phylum == "p__Bacteroidetes")

tern_df_bacteroidetes_v1 <- prep_ternary(phylo.bacteroidetes.v1, group="Group", level= "Species")
tern_df_bacteroidetes_v1 <- tern_df_bacteroidetes_v1[tern_df_bacteroidetes_v1$Phylum != "Unknown", ]     

tern_df_bacteroidetes_v2 <- prep_ternary(phylo.bacteroidetes.v2, group="Group", level= "Species")
tern_df_bacteroidetes_v2 <- tern_df_bacteroidetes_v2[tern_df_bacteroidetes_v2$Phylum != "Unknown", ]
#################################################################################
phylo.fuso.v1 = subset_taxa(phylo_relative_v1, Phylum == "p__Fusobacteria")
phylo.fuso.v2 = subset_taxa(phylo_relative_v2, Phylum == "p__Fusobacteria")

tern_df_fuso_v1 <- prep_ternary(phylo.fuso.v1, group="Group", level= "Species")
tern_df_fuso_v1 <- tern_df_fuso_v1 [tern_df_fuso_v1 $Phylum != "Unknown", ]     

tern_df_fuso_v2 <- prep_ternary(phylo.fuso.v2, group="Group", level= "Species")
tern_df_fuso_v2 <- tern_df_fuso_v2[tern_df_fuso_v2$Phylum != "Unknown", ]
#################################################################################
phylo.spirochaetes.v1 = subset_taxa(phylo_relative_v1, Phylum == "p__Spirochaetes")
phylo.spirochaetes.v2 = subset_taxa(phylo_relative_v2, Phylum == "p__Spirochaetes")

tern_df_spirochaetes_v1 <- prep_ternary(phylo.spirochaetes.v1, group="Group", level= "Species")
tern_df_spirochaetes_v1 <- tern_df_spirochaetes_v1 [tern_df_spirochaetes_v1 $Phylum != "Unknown", ]     

tern_df_spirochaetes_v2 <- prep_ternary(phylo.spirochaetes.v2, group="Group", level= "Species")
tern_df_spirochaetes_v2 <- tern_df_spirochaetes_v2[tern_df_spirochaetes_v2$Phylum != "Unknown", ]

#################################################################################
phylo.firmicutes.v1 = subset_taxa(phylo_relative_v1, Phylum == "p__Firmicutes")
phylo.firmicutes.v2 = subset_taxa(phylo_relative_v2, Phylum == "p__Firmicutes")

tern_df_firmicutes_v1 <- prep_ternary(phylo.firmicutes.v1, group="Group", level= "Species")
tern_df_firmicutes_v1 <- tern_df_firmicutes_v1 [tern_df_firmicutes_v1 $Phylum != "Unknown", ]     

tern_df_firmicutes_v2 <- prep_ternary(phylo.firmicutes.v2, group="Group", level= "Species")
tern_df_firmicutes_v2 <- tern_df_firmicutes_v2[tern_df_firmicutes_v2$Phylum != "Unknown", ]

#################################################################################
phylo.actino.v1 = subset_taxa(phylo_relative_v1, Phylum == "p__Actinobacteria")
phylo.actino.v2 = subset_taxa(phylo_relative_v2, Phylum == "p__Actinobacteria")

tern_df_actino_v1 <- prep_ternary(phylo.actino.v1, group="Group", level= "Species")
tern_df_actino_v1 <- tern_df_actino_v1 [tern_df_actino_v1 $Phylum != "Unknown", ]     

tern_df_actino_v2 <- prep_ternary(phylo.actino.v2, group="Group", level= "Species")
tern_df_actino_v2 <- tern_df_actino_v2[tern_df_actino_v2$Phylum != "Unknown", ]
#################################################################################
phylo.proteo.v1 = subset_taxa(phylo_relative_v1, Phylum == "p__Proteobacteria")
phylo.proteo.v2 = subset_taxa(phylo_relative_v2, Phylum == "p__Proteobacteria")

tern_df_proteo_v1 <- prep_ternary(phylo.proteo.v1, group="Group", level= "Species")
tern_df_proteo_v1 <- tern_df_proteo_v1 [tern_df_proteo_v1 $Phylum != "Unknown", ]     

tern_df_proteo_v2 <- prep_ternary(phylo.proteo.v2, group="Group", level= "Species")
tern_df_proteo_v2 <- tern_df_proteo_v2[tern_df_proteo_v2$Phylum != "Unknown", ]
#################################################################################
#tern_df$size <- (apply(tern_df[2:4], 1, mean))
pdf("tern_df_SR1_v1.pdf", width = 4, height = 4)
ggtern(tern_df_SR1_v1, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum),
             size=5,
             alpha=0.7, 
             show.legend=T) + scale_size(range=c(0, 6)) +  geom_mask() + 
  geom_text(aes(label = Species), vjust=1,size=2) +
  scale_color_manual(values = c("#FA6B09FF"))+ tern.theme
dev.off()

pdf("tern_df_SR1_v2.pdf", width = 4, height = 4)
ggtern(tern_df_SR1_v2, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum),
             size=5,
             alpha=0.7, 
             show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#FA6B09FF"))+ tern.theme
dev.off()
#################################################################################
pdf("tern_df_TM7_v1.pdf", width = 4, height = 4)
ggtern(tern_df_TM7_v1, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum),
             size=5,
             alpha=0.7, 
             show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size=2) +
  scale_color_manual(values = c("#4B0082"))+tern.theme
dev.off()

pdf("tern_df_TM7_v2.pdf", width = 4, height = 4)
ggtern(tern_df_TM7_v2, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum),
             size=5,
             alpha=0.7, 
             show.legend=T) +
  scale_size(range=c(0, 6)) +  geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size=2) +
  scale_color_manual(values = c("#4B0082"))+tern.theme
dev.off()
#################################################################################
pdf("tern_df_bacteroidetes_v1.pdf", width = 4, height = 4)
ggtern(tern_df_bacteroidetes_v1, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Genus), vjust=1, size = 2) +
  scale_color_manual(values = c("#FFD700"))+tern.theme
dev.off()

pdf("tern_df_bacteroidetes_v2.pdf", width = 4, height = 4)
ggtern(tern_df_bacteroidetes_v2, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = unique), vjust=1, size = 2) +
  scale_color_manual(values = c("#FFD700"))+tern.theme
dev.off()
#################################################################################
pdf("tern_df_fuso_v1.pdf", width = 4, height = 4)
ggtern(tern_df_fuso_v1, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#FF0000"))+tern.theme
dev.off()

pdf("tern_df_fuso_v2.pdf", width = 4, height = 4)
ggtern(tern_df_fuso_v2, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#FF0000"))+tern.theme
dev.off()
#################################################################################
pdf("tern_df_spirochaetes_v1.pdf", width = 4, height = 4)
ggtern(tern_df_spirochaetes_v1, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#FFFF00"))+tern.theme
dev.off()

pdf("tern_df_spirochaetes_v2.pdf", width = 4, height = 4)
ggtern(tern_df_spirochaetes_v2, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#FFFF00"))+tern.theme
dev.off()
############################################################################################

pdf("tern_df_firmicutes_v1.pdf", width = 4, height = 4)
ggtern(tern_df_firmicutes_v1, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#008000"))+tern.theme
dev.off()

pdf("tern_df_firmicutes_v2.pdf", width = 4, height = 4)
ggtern(tern_df_firmicutes_v2, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#008000"))+tern.theme
dev.off()
############################################################################################
pdf("tern_df_actino_v1.pdf", width = 4, height = 4)
ggtern(tern_df_actino_v1, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#FF1493"))+tern.theme
dev.off()

pdf("tern_df_actino_v2.pdf", width = 4, height = 4)
ggtern(tern_df_actino_v2, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#FF1493"))+tern.theme
dev.off()
############################################################################################
pdf("tern_df_proteo_v1.pdf", width = 4, height = 4)
ggtern(tern_df_proteo_v1, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#1E90FF"))+tern.theme
dev.off()

pdf("tern_df_proteo_v2.pdf", width = 4, height = 4, bg = "transparent")
ggtern(tern_df_proteo_v2, mapping= aes(x=CS, y=ES, z=NS)) + 
  geom_point(aes(color= Phylum), size=5,alpha=0.7,show.legend=T) +
  scale_size(range=c(0, 6)) + geom_mask() + 
  geom_text(aes(label = Species), vjust=1, size = 2) +
  scale_color_manual(values = c("#1E90FF"))+tern.theme
dev.off()