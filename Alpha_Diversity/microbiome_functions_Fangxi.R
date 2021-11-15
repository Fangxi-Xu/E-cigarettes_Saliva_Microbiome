#Author: Fangxi Xu
#calculate alpha diversity index: Observed; Shannon; PD
alpha_index <- function(file) {
phylo_object <- file #input data has to be a phyloseq object
phylo.div <- alpha(phylo_object, index = "all") #calculate for all index
phylo.meta <- meta(phylo_object) #add metadata
phylo.meta$SampleID <- rownames(phylo.meta) #Add the rownames as a new colum for easy integration later
phylo.div$SampleID <- rownames(phylo.div) # Add the rownames to table
phylo.df <- merge(phylo.div,phylo.meta, by = "SampleID")
library("picante") #package for calculate PD
phy_alp_otu <- as.data.frame(phylo_object@otu_table)
phy_alp_tree <- phylo_object@phy_tree
df.pd.temp <- pd(t(phy_alp_otu), phy_alp_tree,include.root=T)
df.pd.temp$SampleID <- rownames(df.pd.temp)
phylo.df <- merge(phylo.df,df.pd.temp, by = "SampleID") ##merge these two data frames into one
phylo.df2 <-phylo.df[, c("Group", "Visit","observed", "diversity_shannon", "PD")]
colnames(phylo.df2) <- c("Group", "Visit","Observed ASVs", "Shannon", "PD")
phylo_df_melt <- reshape2::melt(phylo.df2)
}

#test
v1_101_alpha <- alpha_index(phylo_v1_clean)

#function of Levene's test for homogeneity of variance across groups
alpha_levene <- function(file) {
  alpha_index <- file # input data is output of alpha_index function
  alpha_observed <- subset(alpha_index, variable == "Observed ASVs")
  alpha_shannon <- subset(alpha_index, variable == "Shannon")
  alpha_pd <- subset(alpha_index, variable == "PD")
  library("car")
  leveneTest_observed <- alpha_observed%>% levene_test(value ~ Group)
  leveneTest_shannon <- alpha_shannon%>% levene_test(value ~ Group)
  leveneTest_pd <- alpha_pd%>% levene_test(value ~ Group)
  levene_results <- cbind(leveneTest_observed$p, leveneTest_shannon$p,leveneTest_pd$p)
  colnames(levene_results) <- c("Observed ASVs", "Shannon", "PD")
  levene_results <- data.frame(levene_results)
}

#test
v1_101_alpha_levene <- alpha_levene(v1_101_alpha)

#function of Shapiro-Wilk test of normality
alpha_normality <- function(file) {
  alpha_index <- file # input data is output of alpha_index function
  alpha_observed <- subset(alpha_index, variable == "Observed ASVs")
  alpha_shannon <- subset(alpha_index, variable == "Shannon")
  alpha_pd <- subset(alpha_index, variable == "PD")
  normality_observed <- shapiro_test(alpha_observed$value)
  normality_shannon <- shapiro_test(alpha_shannon$value)
  normality_pd <- shapiro_test(alpha_pd$value)
  normality_results <- rbind(normality_observed, normality_shannon, normality_pd)
}

#test
v1_101_alpha_normality <- alpha_normality(v1_101_alpha)

#Next step for analysis: Based on the results of the homogeneity of variance test and normality test, choose a proper test to use
#ANOVA assumes that the data is normally distributed and homogeneity of variance. Therefore, if use ANOVA and tukey post hoc if accept null hypo
#if reject null hypo, use non-parametric version of ANOVA - Dunn's post hoc
#same assumption apply for t-test and wilcoxon test if only have 2 groups

#anova function
anova_stats <- function(file) {
  alpha_index <- file # input data is output of alpha_index function
  stat.test.anova <- alpha_index %>%
    group_by(variable) %>%
     anova_test(value ~ Group) %>%
      add_significance("p")
  stat.test.anova$method <- "ANOVA"
  stat.test.anova <-data.frame(stat.test.anova)
}

#test
v1_101_anova <- anova_stats(v1_101_alpha)


#Kruskal–Wallis one-way analysis of variance: rank sum based
kw_stats <- function(file) {
  alpha_index <- file # input data is output of alpha_index function
  stat.test.kw <- alpha_index %>%
    group_by(variable) %>%
    kruskal_test(value ~ Group) %>%
    add_significance("p")
  stat.test.kw$method <- "Kruskal–Wallis"
  stat.test.kw <- data.frame(stat.test.kw)
}

#test
v1_101_kw <- kw_stats(v1_101_alpha)

##correct for multiple comparisons: tukey_hsd following anova
tukey_stats <-function(file) {
  alpha_index <- file # input data is output of alpha_index function
  stat.test.tukey <-  alpha_index %>%
    group_by(variable) %>%
    tukey_hsd(value ~ Group) 
    stat.test.tukey$method <- "Tukey_HSD"
    stat.test.tukey <- data.frame(stat.test.tukey)
}

#test
v1_101_tukey <- tukey_stats(v1_101_alpha)

#if only have 2 comparative groups to analyze, use t test or wilcox based on normality
tat.test.tt <- phylo_df_melt %>%
  group_by(variable) %>%
  t_test(value ~ Visit,p.adjust.method = "bonferroni", paired = TRUE ) %>%
  add_significance("p") %>%
  adjust_pvalue() %>%
  add_significance() %>%
  #filter(p.signif != "ns") #filter non-sig pairs
  
  stat.test.tt
#plot alpha diversity box-plot using output of alpha_index function!!!(filter ns stats before plotting)




