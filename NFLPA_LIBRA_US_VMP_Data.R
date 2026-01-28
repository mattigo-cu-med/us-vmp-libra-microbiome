#US-VMP LIBRA RData
setwd("/Volumes/BrainHealth/B4 Health Lab/Matt_Folder/NFLPA_Microbiome/Analysis")

libra.alpha.all <- read.csv("LIBRA_variables_noNAs_021225-alpha.csv", stringsAsFactors = T)

x_beta <- read.table("bray_curtis_dm_NFL.txt", check.names = F)

meta.in <- read.csv("LIBRA_variables_CCA_noNAs_021225.csv", stringsAsFactors = T)


taxa_species <- read.table(file = "taxatable-species-short-absolute_NFL.tsv", 
                   sep = "\t", row.names = 1, header = T, check.names = F)


meta <- read.csv("LIBRA_variables_noNAs_021225-alpha.csv", header = T, check.names = F, 
                 row.names = 1, stringsAsFactors = T)


input_data <- read.table("taxatable-species-short-absolute_NFL.tsv", 
                         header = T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)


dat <- read.csv("uniref-pathway-filtered-absolute_NFL1.csv", 
                header = TRUE, row.names = 1, check.names = F)


meta <- read.csv("LIBRA_variables_noNAs_021225-alpha.csv", header = T, 
                 row.names = 1, stringsAsFactors = T, check.names = F)


input_data_phylum <- read.table("/Volumes/B4_Backup/Stamper_Mac_VA_files_060325/NFL-grant/Metagenomic_data/filtered taxonomy tables/taxatables-aggregated-absolute/taxatable-phylum-short-absolute_NFL.tsv", 
                                header = T, sep = "\t", row.names = 1, stringsAsFactors = F, check.names = F)


#These need to get renamed
x_pathway <- read.table('/Volumes/B4_Backup/Stamper_Mac_VA_files_060325/NFL-grant/Metagenomic_data/beta diversity/bray_curtis_pathway_dm_NFL.txt', check.names = F)
pathway_abundance <- read.csv(file = '/Volumes/B4_Backup/Stamper_Mac_VA_files_060325/NFL-grant/Metagenomic_data/filtered functional tables/uniref-pathway-filtered-absolute_NFL1.csv', 
                 row.names = 1, header = T, check.names = F)


save(file="NFLPA_LIBRA_US_VMP_Microbiome.RData")