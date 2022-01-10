setwd("/Volumes/Samsung_T5/IECOH_SALv1v2/R-Fangxi")
# create plot folder
path <- "./plot"
dir.create(path)
#Alpha Diversity
# Setup environment
library(microbiome) # data analysis and visualization
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(microbiomeutilities) # some utility tools 
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(DT) # interactive tables in html and markdown
library(data.table) # alternative to data.frame
library(dplyr)# data handling 
library(rstatix) #stats test
library(ggpubr)
#######################################################################################
#V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1V1
#######################################################################################
#plot all 101 patients V1 alpha
all_v1 <- alpha_index(phylo_v1_clean)
all_v1_levene <- alpha_levene(all_v1)#accept
all_v1_normality <-alpha_normality(all_v1)#accept
all_v1_anova <- anova_stats(all_v1)#observed and pd reject, shannon accept
all_v1_tukey <- tukey_stats(all_v1)#pd reject: adjusted p = 0.0267
all_v1_tukey_sig <- all_v1_tukey%>%filter(p.adj.signif != "ns") #filter non-sig pairs
pdf(file.path(path, "101_V1_alpha.pdf"),width=6,height=4)
ggboxplot(all_v1, x = "Group", y = "value", fill = "Group", color = "Group", order = c("CS","ES","NS"),
          legend= "right",  add = "jitter", facet.by = "variable", scales = "free")+
  stat_pvalue_manual(all_v1_tukey_sig, label = "p.adj.signif", y.position =c(23.5) , bracket.size=0.6, size=6 )+
  scale_fill_manual(values=c("white", "white", "white"))+
  scale_color_manual(values = c("CS"="firebrick3", "ES"="goldenrod2", "NS"= "dodgerblue2"))+
  ylab("Alpha Diversity Measure")+
  alpha.theme
dev.off()
write.csv(all_v1_anova, file = "/Volumes/Samsung_T5/IECOH_SALv1v2/R-Fangxi/plot/alpha_div/alpha_stats/all_v1_anova.csv")
write.csv(all_v1_tukey, file = "/Volumes/Samsung_T5/IECOH_SALv1v2/R-Fangxi/plot/alpha_div/alpha_stats/all_v1_tukey.csv")
#plot severe V1 alpha
se_v1 <- alpha_index(phylo_se_v1_clean)
se_v1_levene <- alpha_levene(se_v1)#accept
se_v1_normality <-alpha_normality(se_v1)#accept
se_v1_anova <- anova_stats(se_v1)#reject
se_v1_tukey <- tukey_stats(se_v1)#ES NS reject 
se_v1_tukey_sig <- se_v1_tukey %>%filter(p.adj.signif != "ns") #filter non-sig pairs
pdf(file.path(path, "se_V1_alpha.pdf"),width=6,height=4)
ggboxplot(se_v1, x = "Group", y = "value", fill = "Group", color = "Group", order = c("CS","ES","NS"),
          legend= "right",  add = "jitter", facet.by = "variable", scales = "free")+
  stat_pvalue_manual(se_v1_tukey_sig, label = "p.adj.signif", y.position =c(380,4.7, 24) , bracket.size=0.6, size=6 )+
  scale_fill_manual(values=c("white", "white", "white"))+
  scale_color_manual(values = c("CS"="firebrick3", "ES"="goldenrod2", "NS"= "dodgerblue2"))+
  ylab("Alpha Diversity Measure")+
  alpha.theme
dev.off()
write.csv(se_v1_anova, file = "/Volumes/Samsung_T5/IECOH_SALv1v2/R-Fangxi/plot/alpha_div/alpha_stats/se_v1_anova.csv")
write.csv(se_v1_tukey, file = "/Volumes/Samsung_T5/IECOH_SALv1v2/R-Fangxi/plot/alpha_div/alpha_stats/se_v1_tukey.csv")
#plot mild/moderate(mm) V1 alpha
mm_v1 <- alpha_index(phylo_mm_v1_clean)
mm_v1_levene <- alpha_levene(mm_v1)#accept
mm_v1_normality <-alpha_normality(mm_v1)#accept
mm_v1_anova <- anova_stats(mm_v1)#accept
#no post hoc conducted because no differnece identified by anova
pdf(file.path(path, "mm_V1_alpha.pdf"),width=6,height=4)
ggboxplot(mm_v1, x = "Group", y = "value", fill = "Group", color = "Group", order = c("CS","ES","NS"),
          legend= "right",  add = "jitter", facet.by = "variable", scales = "free")+
  scale_fill_manual(values=c("white", "white", "white"))+
  scale_color_manual(values = c("CS"="firebrick3", "ES"="goldenrod2", "NS"= "dodgerblue2"))+
  ylab("Alpha Diversity Measure")+
  alpha.theme
dev.off()
write.csv(mm_v1_anova, file = "/Volumes/Samsung_T5/IECOH_SALv1v2/R-Fangxi/plot/alpha_div/alpha_stats/mm_v1_anova.csv")

#######################################################################################
#V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2V2
#######################################################################################
#plot all 101 patients V2 alpha
all_v2 <- alpha_index(phylo_v2_clean)
all_v2_levene <- alpha_levene(all_v2)#accept
all_v2_normality <-alpha_normality(all_v2)#pd reject
all_v2_anova <- anova_stats(all_v2)#accept
all_v2_kw <-kw_stats(all_v2)#accept
pdf(file.path(path, "101_V2_alpha.pdf"),width=6,height=4)
ggboxplot(all_v2, x = "Group", y = "value", fill = "Group", color = "Group", order = c("CS","ES","NS"),
          legend= "right",  add = "jitter", facet.by = "variable", scales = "free")+
  scale_fill_manual(values=c("white", "white", "white"))+
  scale_color_manual(values = c("CS"="firebrick3", "ES"="goldenrod2", "NS"= "dodgerblue2"))+
  ylab("Alpha Diversity Measure")+
  alpha.theme
dev.off()
#plot severe V2 alpha
se_v2 <- alpha_index(phylo_se_v2_clean)
se_v2_levene <- alpha_levene(se_v2)#accept
se_v2_normality <-alpha_normality(se_v2)#pd reject
se_v2_anova <- anova_stats(se_v2)#accept
se_v2_kw <-kw_stats(se_v2)#accept
pdf(file.path(path, "se_V2_alpha.pdf"),width=6,height=4)
ggboxplot(se_v2, x = "Group", y = "value", fill = "Group", color = "Group", order = c("CS","ES","NS"),
          legend= "right",  add = "jitter", facet.by = "variable", scales = "free")+
  scale_fill_manual(values=c("white", "white", "white"))+
  scale_color_manual(values = c("CS"="firebrick3", "ES"="goldenrod2", "NS"= "dodgerblue2"))+
  ylab("Alpha Diversity Measure")+
  alpha.theme
dev.off()

#plot mild/moderate V2 alpha
mm_v2 <- alpha_index(phylo_mm_v2_clean)
mm_v2_levene <- alpha_levene(mm_v2)#accept
mm_v2_normality <-alpha_normality(mm_v2)#pd reject
mm_v2_anova <- anova_stats(mm_v2)#accept
mm_v2_kw <-kw_stats(mm_v2)#accept
pdf(file.path(path, "mm_V2_alpha.pdf"),width=6,height=4)
ggboxplot(phylo_df_melt, x = "Group", y = "value", fill = "Group", color = "Group", order = c("CS","ES","NS"),
          legend= "right",  add = "jitter", facet.by = "variable", scales = "free")+
  scale_fill_manual(values=c("white", "white", "white"))+
  scale_color_manual(values = c("CS"="firebrick3", "ES"="goldenrod2", "NS"= "dodgerblue2"))+
  ylab("Alpha Diversity Measure")+
  alpha.theme
dev.off()
####################################################################################################
################################################################################################
#paired sample compare visits
#Ziyan did this section
#paired sample t test
stat.test.tt <- phylo_df_melt %>%
  group_by(variable) %>%
  t_test(value ~ Visit,p.adjust.method = "bonferroni", paired = TRUE ) %>%
  add_significance("p") %>%
  adjust_pvalue() %>%
  add_significance() %>%
  #filter(p.signif != "ns") #filter non-sig pairs
  
  stat.test.tt
#paired sample wilcoxon signed rank test
stat.test.wil <- phylo_df_melt %>%
  group_by(variable) %>%
  wilcox_test(value ~ Visit,p.adjust.method = "bonferroni", paired = TRUE ) %>%
  add_significance("p") %>%
  adjust_pvalue() %>%
  add_significance() %>%
  filter(p.signif != "ns") #filter non-sig pairs
stat.test.wil

#paired sample visit compare
#Ziyan did it - please add script