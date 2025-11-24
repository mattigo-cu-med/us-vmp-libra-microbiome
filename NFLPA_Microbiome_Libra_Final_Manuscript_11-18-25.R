#NFLPA VMP metagenomic analysis
setwd("/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/NFLPA_Microbiome/Analysis")
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ANCOMBC)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(vegan)
library(knitr)
library(forcats)
library(tidyverse)
library(DT)
library(outliers)
library(beeswarm)
library(WGCNA)
library(flextable)
library(magrittr)
library(dplyr)
library(tidyheatmaps)
library(devtools)
library(ggvegan)
theme_set(theme_bw())
#####Demographics#####
#summarySE funciton 
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#Data sets
#nfl <- read.csv("Metadata/NFL-VMP_timeMatched-demos-alpha.csv")
#libra <- read.csv("LIBRA/LIBRA_variables_noNAs_012325.csv", stringsAsFactors = T)
#a.dat <- read.table("Metagenomic_data/alpha diversity/alphadiversity_NFL.txt", sep = "\t", header = T, 
#                   stringsAsFactors = T)
#libra.alpha <- merge(x = libra, y = a.dat, by.x = "SampleID", by.y = "Sample_ID")

######################################################
#after all these data were combined which was done in the **metadata-wrangling.R** 
# file, this was written out to check for consistency and make the Libra variable
#write.csv(libra.alpha.all, "LIBRA/LIBRA_variables_noNAs_021125-alpha.csv", row.names = F)
libra.alpha.all <- read.csv("LIBRA_variables_noNAs_021225-alpha.csv", stringsAsFactors = T)
libra.alpha.all$LIBRA.score.quartile <- as.factor(libra.alpha.all$LIBRA.score.quartile)
#Examining the distribution of the libra variable
#this has already been added to the metadata and will be commented out 
lib.var <- libra.alpha.all$LIBRA.score
lib.var <- data.frame(lib.var)
myCol <- lapply(lib.var, function(x) cut(x, breaks = quantile(x, na.rm = T), labels = FALSE))
myCol
myCol$lib.var <- as.factor(myCol$lib.var)
#the equation could not comute 2 values that should be in Q1 FYI
#adding in quartile data and writing out for future use
#write.csv(libra.alpha.all, "LIBRA/LIBRA_variables_noNAs_021125-alpha1.csv", row.names = T)
#libra.alpha.all$LIBRA.score.quartile <- myCol$lib.var
#table(myCol)
########################################################
#Making some metadata summary plots
library(table1)
#p-value function
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    ano <- aov(y ~ g)
    p <- summary(ano)[[1]][[5]][1]
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=2, eps=0.001)))
  
}
p <- aov(Age ~ as.factor(LIBRA.score.quartile) , data = libra.alpha.all)
summary(p)
#Demographics
libra_table_2 = table1(~ Age + Race + Ethnicity + BMI + Homeless.ever + Homeless.Times + Multiple.Mild.TBIs, data = libra.alpha.all)
libra_table_3 = table1(~ Age + Race + Ethnicity + BMI + Homeless.ever + Homeless.Times + Multiple.Mild.TBIs | LIBRA.score.quartile, 
                       data = libra.alpha.all, overall = T,extra.col=list(`P-value`=pvalue))

#Libra vaiables
libra_table_2a = table1(~ LIBRA.score +  Obese +  Lifetime_smoking + Diabetes + Hypertension + 
                          Hypercholesterolemia + Coronary.heart.disease + 
                          Chronic_kidney_disease +  Current_alcohol + Depressed + Healthy_diet + Physical.activity.H.L, data = libra.alpha.all)


t1flex(libra_table_2) %>% 
  save_as_docx(path="libra_table_2.docx")

t1flex(libra_table_2a) %>% 
  save_as_docx(path="libra_table_2a.docx")


t1flex(libra_table_3) %>% 
  save_as_docx(path="libra_table_3.docx")

#### 

install.packages("ggbeeswarm")
library(ggbeeswarm)
p_swarm <- ggplot(data = libra.alpha.all, aes(y = LIBRA.score, x="",color = LIBRA.score.quartile)) +
  geom_beeswarm(cex = 3)+
  scale_y_continuous(breaks = seq(0, 10, by=2))+
  scale_color_manual(values = colvec)+
  guides(color = guide_legend(reverse = TRUE))+
  labs(y="LIBRA Score", x="LIBRA Score Distribution", color="LIBRA Score Quartile")+
  theme_bw()

ggsave("Figures/S1_LIBRA_swarm_figure.jpeg", p_swarm, w=6, h=6, units="in")


####Figure 1 Alpha and Beta
##############Alpha#####
#Alpha plots
#Note color Is differen
#Facet titles are different
#Botom Title different
#Add Letter code
libra.alpha.all$LIBRA.score.quartile <- factor(libra.alpha.all$LIBRA.score.quartile)
melt.libra.alpha <- melt(libra.alpha.all, id.vars = c("SampleID", "LIBRA.score", "LIBRA.score.quartile"), 
                         variable.name = "alpha", measure.vars = c("observed_otus", "shannon"))
new_labels <- c("observed_otus" = "Observed OTUs", "shannon" = "Shannon Index")
colvec <- c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A")

p1b <- ggplot(data = melt.libra.alpha %>% dplyr::filter(alpha=="shannon"), aes(x = LIBRA.score.quartile, y = value)) +
  geom_boxplot(aes(fill = LIBRA.score.quartile)) + 
  geom_jitter(aes(color = LIBRA.score.quartile)) + 
  scale_color_manual(values=(colvec)) + 
  scale_fill_manual(values=(colvec))+
  theme_bw() + #facet_wrap(~alpha, scales = "free", labeller = labeller(alpha = new_labels)) + 
  #coord_cartesian(clip="off") +
  #annotate("text", x= 0, y=98, label = "bold(A)" )
  labs(tag = "B", x="LIBRA Score Quartile", y="Shannon Index") +
  stat_compare_means(method = "kruskal.test", label.y.npc="top", label.x = 2.75, size = 3) + # Add Kruskal-Wallis p-value
  theme(plot.tag = element_text(size = 12, face = "bold"),
        legend.position = "none") 

?labs

ggsave(paste0("Figures/Figure_1b_Alpha_Boxplots", Sys.Date(),".jpeg"),p1a, height = 10, width = 14, units="in" )
#ggsave("LIBRA/Figures/alpha/LIBRAscore-quaritles.pdf")
#Change Location
#Change Titles
p1a <- ggplot(data = melt.libra.alpha %>% dplyr::filter(alpha=="shannon"), aes(x = LIBRA.score, y = value)) +
  geom_point() + 
  geom_smooth(se = F, method = "lm") +
  stat_cor(label.y.npc="top", label.x = 5.25, size = 3)+
  theme_bw() + #facet_wrap(~alpha, scales = "free", labeller = labeller(alpha = new_labels)) +
  labs(tag = "A", x="LIBRA Score", y="Shannon Index") +
  theme(plot.tag = element_text(size = 12, face = "bold"))


ggsave(paste0("Figures/Figure_1a_Alpha_Points", Sys.Date(),".jpeg"),p1a, height = 10, width = 14, units="in" )

#####Beta######
x <- read.table("bray_curtis_dm_NFL.txt", check.names = F)
#next line will determine which metadata is being piped in so make sure to check before proceeding 
meta.in <- read.csv("LIBRA_variables_CCA_noNAs_021225.csv", stringsAsFactors = T)
#matching up metadata with the distance matrix
ind <- match(x = rownames(x), table = meta.in$SampleID)
ind1 <- !is.na(ind)
x1 <- x[ind1,ind1]
ind2 <- match(x = rownames(x1), table = meta.in$SampleID)
meta.in1 <- meta.in[ind2,]
#this next line should be all true - if not, then something is wrong
table(rownames(x1) == meta.in1$SampleID) 
rownames(x1)
colnames(x1)
meta.in1$SampleID
#read in taxa
taxa <- read.table(file = "taxatable-species-short-absolute_NFL.tsv", 
                   sep = "\t", row.names = 1, header = T, check.names = F)

Taxa <- t(taxa)
head(Taxa)

ind1 <- match(x = colnames(x1), table = rownames(Taxa))
Taxa1 <- Taxa[ind1,]
head(Taxa1)

#**Traditional** CCA start
pca_w=capscale(x1~1, data=meta.in1[,-1])
plot(pca_w)
biplot(pca_w, type=c("points") )
fit=envfit(pca_w, meta.in1[,-1], na.rm = T)
#p.max in the code below will determine which variables end up on the graph - adjust for more or less
plot(fit, col="blue", p.max=0.5)
fit
#p.max in the code below will determine which taxa end up on the graph - adjust for more or less
fit=envfit(pca_w, Taxa1, na.rm = T)
plot(fit, col="red",p.max=0.92,)
fit

#**Special plot** CCA start
#this code is for coloring dots on biplot, if not wanting that then next 3 lines can be removed
meta.in$LIBRA.score.quartile <- as.factor(meta.in$LIBRA.score.quartile)
with(meta.in, levels(LIBRA.score.quartile))
scl <- 2
#colvec <- c("red", "black", "lightgreen", "royalblue")
#colvec <- c("royalblue", "lightgreen", "black", "red")

#CCA start
pca_w=capscale(x1~1, data=meta.in1[,-1])
pca_test = data.frame(pca_w$Ybar)
test_filter_pca = pca_test %>% filter(Dim1 == min(Dim1))
test_filter_pca2 = pca_test %>% slice_max(Dim1, n = 3)

plot(pca_w, scaling = 2, type = "n")
biplot(pca_w, type = c("points"),scaling = scl)
with(meta.in1, points(pca_w, col = colvec[LIBRA.score.quartile], scaling = scl, pch = 21, bg = colvec[LIBRA.score.quartile]))
head(with(meta.in1, colvec[LIBRA.score.quartile]))
meta.in1$LIBRA.score.quartile <- as.factor(meta.in1$LIBRA.score.quartile)

with(meta.in1, legend("topright", legend = levels(LIBRA.score.quartile), bty = "n", col = colvec, pch = 21, pt.bg = colvec))
fit=envfit(pca_w, meta.in1[,-1], na.rm = T)

?capscale
colvec <- c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A")

ggplot_beta = fortify(pca_w)
ggplot_beta = ggplot_beta %>% rename(SampleID = label)
beta_data_set = merge(ggplot_beta, meta.in, by="SampleID")
p1c = ggplot(beta_data_set, aes(x=MDS1, y=MDS2, color = factor(LIBRA.score.quartile, levels = c("4", "3", "2", "1")))) +
  geom_point()+
  theme_bw()+
  coord_cartesian()+
  scale_color_manual(values=rev(colvec))+
  stat_ellipse(type = "norm", linetype = 2) +
  labs(color = "Libra Quartile", tag = "C") +
  theme(plot.tag = element_text(size = 12, face = "bold")) +
  annotate("text",label=paste0("PERMANOVA = ",round(result$`Pr(>F)`[1],2)), x=1.55,y=1.9, size =3.5)


ggsave(paste0("Figures/Figure_1c_Beta", Sys.Date(),".jpeg"),p1c, height = 10, width = 14, units="in" )

gg_top = ggarrange(p1a, p1b, ncol=2)
p1 = ggarrange(gg_top, p1c, nrow=2)

ggsave(paste0("Figures/Figure_1_Alpha_Beta_Diversity", Sys.Date(),".jpeg"),p1, height = 8, width = 7, units="in" )

#####Figure 2 Heatmap and DEA#####
#MaAsLin2 for significant species
#MaAsLin for taxa 
library(maaslin3)
BiocManager::install("Maaslin2")
library(Maaslin2)

meta <- read.csv("LIBRA_variables_noNAs_021225-alpha.csv", header = T, check.names = F, 
                 row.names = 1, stringsAsFactors = T)
table(meta$Diabetes)
table(meta$LIBRA.score.quartile, useNA = "ifany")
table(meta$Age.baseline, useNA = "ifany")
table(meta$Race, useNA = "ifany")
############
#note that this code is repeated for each taxonomic level so that the input_data 
#and the output destination were changed to match

input_data <- read.table("taxatable-species-short-absolute_NFL.tsv", 
                         header = T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)
table(meta$LIBRA.score, useNA = "ifany")
names(meta)

gen.fit_LIBRA <- Maaslin2(input_data = input_data, input_metadata = meta, min_prevalence = .5, 
                          min_abundance = 0.01, normalization = "TSS", 
                          output = "MaAsLin2_Out_final", 
                          fixed_effects = c("LIBRA.score","Age","BMI"),
                          plot_scatter = F,
                          plot_heatmap = F,
                          correction="BH",
                          analysis_method = "LM",
                          transform = "LOG",
                          max_significance = 0.3)


#####Heatmap for relative abundances
libra_mas = data.frame(gen.fit_LIBRA$results)
mas.sig = libra_mas %>% dplyr::filter(pval < 0.05 & value == "LIBRA.score")
mas.sig_qval = libra_mas %>% dplyr::filter(qval < 0.25 & value == "LIBRA.score")

#####Code run for less stringent output
dendo = as.dendrogram((hclust(dist(t(relative_taxa_adbund), method = "euclidian"), method="complete")))
?hclust
plot(dendo)
row.ord <- order.dendrogram(dendo)
reordered_data <- t(relative_taxa_adbund)[row.ord, ]
updated_relative_taxa = as.data.frame(t(reordered_data))

relative_taxa_adbund_filtered = updated_relative_taxa %>% dplyr::filter(row.names(.) %in% mas.sig$feature)

relative_taxa_adbund_filtered_update <- relative_taxa_adbund_filtered %>%
  rownames_to_column(var = "Species")

relative_taxa_adbund_filtered_update2 <- data.frame(t(relative_taxa_adbund_filtered_update[-1]))
colnames(relative_taxa_adbund_filtered_update2) <- relative_taxa_adbund_filtered_update[, 1]

relative_taxa_adbund_filtered_update3 <- relative_taxa_adbund_filtered_update2 %>%
  rownames_to_column(var = "SampleID")


#reordered_df <- relative_taxa_adbund_filtered_update3[match(row.ord, relative_taxa_adbund_filtered_update3$SampleID), ]

relative_taxa_adbund_filtered_update3_melt = melt(relative_taxa_adbund_filtered_update3, id.var="SampleID", variable.name="Species", value.name = "Abundance") 


all_taxa_meta_data = plyr::join(relative_taxa_adbund_filtered_update3_melt, meta.in, by="SampleID") %>% 
  mutate(Species = str_remove(Species, "s__")) %>% 
  mutate(Species = str_replace(Species, "_", " "))

all_taxa_meta_data2 = all_taxa_meta_data %>% mutate(Abundance_log = log10(Abundance + 1E-10))
all_taxa_meta_data2$LIBRA.score.quartile = factor(all_taxa_meta_data2$LIBRA.score.quartile, levels = c("4", "3", "2", "1"))
all_taxa_meta_data3 = all_taxa_meta_data2 %>% arrange(Species, (LIBRA.score)) %>%
  mutate(Libra.score.half = case_when(LIBRA.score.quartile == "1" | LIBRA.score.quartile == "2" ~ "Lower",
                                      T~"Upper"))


all_taxa_meta_data_test = all_taxa_meta_data2 %>% 
  arrange(factor(Species, levels =  mas_sig_list$feature), LIBRA.score)

colvec <- c("#4F94CD", "#FFC125","#CD3333", "#4A4A4A")


cols = list( `LIBRA Score Quartile`  = c( "4" = "#4A4A4A", "3" = "#CD3333","2" = "#FFC125","1" = "#4F94CD" ))


#####Final heatmap
all_taxa_meta_data_test = all_taxa_meta_data_test %>%      
  rename(
    `LIBRA Score` = LIBRA.score,
    `LIBRA Score Quartile` = LIBRA.score.quartile)

library(grid)
#setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
tdy_hm_rel2=tidyheatmap(df = all_taxa_meta_data_test,
                        rows = Species,
                        columns = SampleID,
                        values = Abundance_log,
                        #cluster_cols = T,
                        #clustering_distance_cols = "euclidean",
                        #clustering_method = "complete",
                        annotation_col = c(`LIBRA Score Quartile`,  `LIBRA Score`) ,
                        annotation_colors = cols,
                        colors = c("#ffffff","#145afc"),
                        scale= "none",
                        display_numbers = F,
                        gaps_col = `LIBRA Score Quartile`,
                        show_colnames = F)
#setHook("grid.newpage", NULL, "replace")
#setHook("grid.newpage", NULL, "replace")
#grid.text("xlabel example", y=-0.07, gp=gpar(fontsize=16))
dev.off()
#my_gtable = tdy_hm_rel2$gtable
tdy_hm_rel2$gtable$grobs[[2]]$gp = gpar(fontsize = 8, fontface="italic")# assuming that the xlabels are in the third grob



install.packages("ggplotify")
library(ggplotify)
test_hm = as.ggplot(tdy_hm_rel2) 
test_hm = test_hm + labs(tag = "B") +
  annotate("text", x = 0.325, y = 0.00001, label = "Samples") +  
  coord_cartesian(ylim = c(0, 1), clip = "off")+
  theme(plot.margin = margin(t = 0, r = -2.4, b = 0.3, l = 0, unit = "cm"),
        plot.tag = element_text(size = 12, face = "bold"))

ggsave(paste0("Figures/Figure_3_Libra_Quartile_Heatmap_",Sys.Date(),".jpeg"), test_hm, jpeg, height = 8, width = 16, units = "in" )



#####figure 2b
######LFC of significant species####3
gen.fit_LIBRA
mas.sig2 = mas.sig %>% mutate(Up_Down = case_when(coef > 0 ~ "Up",
                                                  T ~ "Down")) %>%
  arrange((qval))%>% 
  mutate(feature = str_remove(feature, "s__")) %>% 
  mutate(feature = str_replace(feature, "_", " "))

mas.sig2$Up_Down = factor(mas.sig2$Up_Down, levels = c("Up", "Down")) 
#Creates figure for horizontal LFC
mas_sig_list = mas.sig2 %>% arrange(desc(coef))

p3 = ggplot(mas.sig2, aes(x=reorder(feature,coef), y=coef, fill=Up_Down)) +
  geom_bar(stat='identity') +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))+
  scale_x_discrete(position = "top")+
  #theme(axis.text.y = element_text(face = "italic")) + # Assign specific colors to groups
  coord_flip() +
  labs(#title="Species LFC by LIBRA Score", 
    x="Species", 
    y="Log Fold Change") +
  scale_fill_manual(values = c("Up" = "darkred", "Down" = "skyblue"))
?theme


ggsave(paste0("Figures/Figure_2_LFC_Significant_Species_",Sys.Date(),".jpeg"), p3, jpeg, height = 8, width = 12, units = "in" )

tdy_hm_rel2

p3 = p3 + labs(tag = "A") +
  theme(plot.tag = element_text(size = 12, face = "bold"))
p2 = ggarrange(p3, test_hm, ncol=1)
ggsave(paste0("Figures/Figure_2_Significant_Species_LFC_HM_",Sys.Date(),".jpeg"), p2, jpeg, height = 10, width = 12, units = "in" )

mas.sig2_outtable = mas.sig2 %>% dplyr::select(feature, coef, stderr, pval, qval) %>% arrange(desc(coef))

sup_table_2 = flextable(mas.sig2_outtable)
sup_table_2 = set_caption(sup_table_2, caption = "Supplementary Table 2. Species associated with LIBRA Scores") %>% autofit()
save_as_docx(sup_table_2, path = "Figures/Supplementary_Table_2_MaAsLin_Results.docx")




#####Pathway Analysis
#Figure 3
####Volcano
#DeSeq for genes and pathways
#Differential Abundance with DESeq2 - done at genus level
BiocManager::install("DESeq2")
detach("package:phyloseq", unload = TRUE)
library(DESeq2)
library(pheatmap)
BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)
BiocManager::install("DOSE")

library(DOSE)
BiocManager::install("pathview")

library(pathview)
BiocManager::install("clusterProfiler")

library(clusterProfiler)
BiocManager::install("AnnotationHub")

library(AnnotationHub)
BiocManager::install("ensembldb")


library(ensembldb)
library(tidyverse)
BiocManager::install("apeglm")
library(apeglm)
dat <- read.csv("uniref-pathway-filtered-absolute_NFL1.csv", 
                header = TRUE, row.names = 1, check.names = F)


#dat[-1] <-sweep(dat[-1], 2, colSums(dat[,-1]), `/`) * 10000
#Round data so that all the data are integers
dat <- dat %>% 
  mutate_if(is.numeric, round)
keep_genes <- rowSums( dat > 100 ) >= 42
dat <- dat[keep_genes,]

dat <- as.matrix(dat)
meta <- read.csv("LIBRA_variables_noNAs_021225-alpha.csv", header = T, 
                 row.names = 1, stringsAsFactors = T, check.names = F)

#align data
ind <- match(x = colnames(dat), table = rownames(meta))
meta1 <- meta[ind,] 
meta1[1:10, 1:10]
#Check that data are aligned
all(colnames(dat) == rownames(meta_transformed))
#make deseq object
meta1
library(caret)
meta1_scale_factors = subset(meta1, select = c(Age, BMI, LIBRA.score))
meta1_scaled  = scale(meta1_scale_factors, scale = T, center = T)


dds <- DESeqDataSetFromMatrix(countData = dat, colData = meta1_scaled, design = ~ Age + BMI + LIBRA.score)
?DESeqDataSetFromMatrix
#generate normalized counts 
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
?counts
?estimateSizeFactors()
#write out table for use in the future
#write.csv(normalized_counts, file="LIBRA/Metagenomic-results/Pathways/multiple-mild/uniref-pathway-filtered-absolute_NFL-normalized_counts.csv",
#          quote=F)
#QC for DE analysis using DESeq2
#If input file is very large > 800 MB then, use the below command
rld <- varianceStabilizingTransformation(dds, blind=T)
#rld <- rlog(dds, blind=TRUE)
#PCA 
plotPCA(rld, intgroup="LIBRA.score")
#Hierarchical Clustering
### Extract the rlog matrix from the object
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor)
#differential expression analysis
dds <- DESeq(dds)
plotDispEsts(dds)
res <- results(dds)
dds_shrink = lfcShrink(dds = dds, res = res, coef="LIBRA.score")
?lfcShrink()
summary(res)
summary(dds_shrink)
summary(dds)

plotMA(res, ylim=c(-2,2))
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(dds_shrink))
all_genes <- as.character(rownames(res))

# Extract significant results

signif_res <- dds_shrink[dds_shrink$padj < 0.3 & dds_shrink$pvalue < 0.05 & !is.na(dds_shrink$padj), ]
signif_res <- res[res$padj < 0.25 & !is.na(res$padj), ]

plotMA(signif_res, ylim=c(-4,4))
#This will be empty if no genes were significant
signif_genes <- as.character(rownames(signif_res))
signif_genes[1:10,1:10]

signif_genes <- as.character(rownames(signif_res))
all_genes[1:10,1:10]


#Below is old code to export sig-taxa
sigtab = cbind(as(res, "data.frame"), all_genes)

sigtab = cbind(as(dds_shrink, "data.frame"), all_genes)


#scale_fill_discrete <- function(palname = "Set1", ...) {scale_fill_brewer(palette = palname, ...)}
#x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
#x = sort(x, TRUE)
#sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
#ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +ggtitle("Food Desert")+theme_bw()
#write.csv(sigtab, "LIBRA/Metagenomic-results/Genes/libra.score/LIBRA-score-sig-pathways1.csv", 
#         row.names = F) 

volcplot <- function(data, padj_threshold = 0.05, fc = 1.0354, plot_title = 'Volcano Plot', plot_subtitle = NULL, genelist_vector = genes$all_genes, genelist_filter = FALSE) {
  
  # Set the fold-change thresholds
  neg_log2fc <- -log2(fc)
  pos_log2fc <- log2(fc)
  
  # Make a dataset for plotting, add the status as a new column
  plot_ready_data <- data %>%
    mutate_at('padj', ~replace(.x, is.na(.x), 1)) %>%
    mutate_at('log2FoldChange', ~replace(.x, is.na(.x), 0)) %>%
    mutate(
      log2fc_threshold = ifelse(log2FoldChange >= pos_log2fc & pvalue <= padj_threshold, 'up',
                                ifelse(log2FoldChange <= neg_log2fc & pvalue <= padj_threshold, 'down', 'ns')
      )
    ) %>%
    mutate(all_genes = replace_na(all_genes, 'none'))
  
  if (genelist_filter) {
    plot_ready_data <- plot_ready_data %>% dplyr::filter(all_genes %in% genelist_vector)
  }
  
  if(!is.null(genelist_vector)) {
    plot_ready_data <- plot_ready_data %>% mutate(all_genes = ifelse(all_genes %in% genelist_vector & pvalue < padj_threshold & log2fc_threshold != 'ns', all_genes, ''))
  }
  
  # Get the number of up, down, and unchanged genes
  up_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'up') %>% nrow()
  down_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'down') %>% nrow()
  unchanged_genes <- plot_ready_data %>% dplyr::filter(log2fc_threshold == 'ns') %>% nrow()
  
  # Make the labels for the legend
  legend_labels <- c(
    str_c('Up: ', up_genes),
    str_c('NS: ', unchanged_genes),
    str_c('Down: ', down_genes)
  )
  
  # Set the x axis limits, rounded to the next even number
  x_axis_limits <- DescTools::RoundTo(
    max(abs(plot_ready_data$log2FoldChange)),
    0.5,
    ceiling
  )
  
  # Set the plot colors
  plot_colors <- c(
    'up' = 'darkred',
    'ns' = 'gray',
    'down' = 'dodgerblue1'
  )
  
  
  # Make the plot, these options are a reasonable strting point
  plot <- ggplot(plot_ready_data) +
    geom_point(
      alpha = 0.25,
      size = 1.5
    ) +
    aes(
      x = log2FoldChange,
      y = -log10(pvalue),
      color = log2fc_threshold,
      label = all_genes
    ) +
    geom_vline(
      xintercept = c(neg_log2fc, pos_log2fc),
      linetype = 'dashed'
    ) +
    geom_hline(
      yintercept = -log10(padj_threshold),
      linetype = 'dashed'
    ) +
    scale_x_continuous(
      'log2(FC)',
      limits = c(-x_axis_limits, x_axis_limits)
    ) +
    scale_color_manual(
      values = plot_colors#,
      #labels = legend_labels
    ) +
    #labs(
    #  color = str_c(fc, '-fold, padj ≤', padj_threshold),
    #  title = plot_title,
    #  subtitle = plot_subtitle
    #) +
    theme_bw(base_size = 24) +
    theme(legend.position = "none",
          aspect.ratio = 1,
          axis.text = element_text(color = 'black'),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),  # Reduces dead area around legend
          legend.spacing.x = unit(0.2, 'cm'),
          legend.title = element_blank()
    )
  
  # Add gene labels if needed
  if (!is.null(genelist_vector)) {
    plot <- plot +
      geom_label_repel(
        size = 2.5,
        force = 1,
        max.overlaps = 100000,
        nudge_x = -.05,
        nudge_y = -.20,
        segment.color = 'black',
        min.segment.length = 0,
        show.legend = FALSE
      )
  }
  plot
}
deseq_results <- as.data.frame(mas.sig2) %>%
  rownamnes_to_column(var = 'ensembl_id') %>%
  left_join({ensembl_id_signif_genes})

library(ggrepel)
volc = volcplot(sigtab2)
volc2 = volcplot(sigtab2)

sigtab2 = na.omit(sigtab)
rownames(sigtab2) <- NULL

genes = sigtab2 %>% dplyr::filter(log2FoldChange > log2(1.0354) | log2FoldChange< -log2(1.0354) ) %>%
  arrange(desc(abs(log2FoldChange))) %>% slice_head(n=6)

out_table = sigtab2 %>% dplyr::filter(padj < 0.3 & pvalue < 0.05) %>%
  arrange(desc(log2FoldChange)) 

ggsave(paste0("Figures/Figure_3_Pathway_Volcano_Plot_",Sys.Date(),".jpeg"), volc, jpeg, height = 8, width = 12, units = "in" )


out_table_pathways = out_table %>% dplyr::select(all_genes, log2FoldChange, lfcSE, pvalue, padj) %>% arrange(desc(log2FoldChange))

sup_table_3 = flextable(out_table_pathways)
sup_table_3 = set_caption(sup_table_3, caption = "Supplementary Table 3. Metabolic Pathways associated with LIBRA Scores") %>% autofit()
save_as_docx(sup_table_3, path = "Figures/Supplementary_Table_3_DESeq_Results.docx")




######Extra for MaAsLin
mal_fun= function(effect,ref){
  fit = Maaslin2(input_data = input_data, input_metadata = meta, min_prevalence = .5, 
                 min_abundance = .01, normalization = "TSS", 
                 output = paste0("MaAsLin_Out_species-level_",effect,"_",Sys.Date()), 
                 fixed_effects = effect,
                 reference = ref,
                 plot_scatter = F,
                 plot_heatmap = F,                          
                 correction="BH",
                 analysis_method = "LM",
                 transform = "LOG",
                 max_significance = 0.3)#,
  #reference = "Current_alcohol,No")}
  return(fit)
}

meta$Physical.inactivity_lib.weight = as.factor(meta$Physical.inactivity_lib.weight)
meta$Healthy_diet
depressed_mal = mal_fun("Depressed", "Depressed,No")

obese_mal = mal_fun("Obese", "Obese,No")

alc_mal = mal_fun("Current_alcohol","Current_alcohol,No")

chd_mal = mal_fun("Coronary.heart.disease", "Coronary.heart.disease,No")

cdd_mal = mal_fun("Chronic_kidney_disease", "Chronic_kidney_disease,No")

diabetes_mal = mal_fun("Diabetes","Diabetes,No")

Hypercholesterolemia_mal = mal_fun("Hypercholesterolemia", "Hypercholesterolemia,No")

Lifetime_smoking_mal = mal_fun("Lifetime_smoking", "Lifetime_smoking,No")

Hypertension_mal = mal_fun("Hypertension","Hypertension,No")

Physical.inactivity_lib.weight_mal = mal_fun("Physical.inactivity_lib.weight", "Physical.inactivity_lib.weight,0")

Healthy_diet_mal = mal_fun("Healthy_diet", "Healthy_diet,Yes")








###Extra
########Figure 2#######
####relative abundance of significant species######
relative_taxa_adbund <- read.table(file = "taxatable-species-short-relative_NFL.tsv", 
                                   sep = "\t", row.names = 1, header = T, check.names = F)

mas.sig <- read.table("significant_results.tsv", 
                      check.names = F, sep = "\t", header = T)

mas.sig_filter = mas.sig %>% filter(qval < 0.1)


relative_taxa_adbund_filtered = relative_taxa_adbund %>% filter(row.names(.) %in% mas.sig$feature)

relative_taxa_adbund_filtered_update <- relative_taxa_adbund_filtered %>%
  rownames_to_column(var = "Species")

relative_taxa_adbund_filtered_update2 <- data.frame(t(relative_taxa_adbund_filtered_update[-1]))
colnames(relative_taxa_adbund_filtered_update2) <- relative_taxa_adbund_filtered_update[, 1]

relative_taxa_adbund_filtered_update3 <- relative_taxa_adbund_filtered_update2 %>%
  rownames_to_column(var = "SampleID")

relative_taxa_adbund_filtered_update3_melt = melt(relative_taxa_adbund_filtered_update3, id.var="SampleID", variable.name="Species", value.name = "Abundance") 


all_taxa_meta_data = merge(relative_taxa_adbund_filtered_update3_melt, meta.in, by="SampleID")

abundance_plot = ggplot(all_taxa_meta_data, aes(y=Abundance, x=LIBRA.score)) + 
  geom_point() +
  facet_wrap(~Species, scales="free") +
  geom_smooth(method="lm") +
  theme(strip.text.x = element_text(face = "italic")) # For x-axis facet titles)
ggsave(paste0("Figures/Figure_2_Abundance_Change_Libra_Significant_Species_",Sys.Date(),".jpeg"), abundance_plot, jpeg, height = 20, width = 25, units = "in" )



