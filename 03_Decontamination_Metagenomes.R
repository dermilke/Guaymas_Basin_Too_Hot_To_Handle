#### This script identifies potential contaminants in metagenome data ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)
library(furrr)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

# Read processed gene-library in datalist format 
datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Full_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv")) %>%
  filter_taxa_datalist((!is.na(KEGG)) | (!is.na(CAZy_ID)) | (Kingdom != "Unclassified"))

# Identify potential contaminants from control-samples. Genes that are slightly enriched
# in control-samples are treated as potential contaminants. Genes are slightly enriched
# in control samples if their abundance is higher than their overall 5% quantile.
datalist_contaminant <- datalist %>%
  filter_station_datalist(Type == "Control") %>%
  singleton_filter(min_station = 1, min_count = 1)

remove <- datalist_contaminant$Count_Data$Gene_ID

datalist_reduced <- datalist %>%
  filter_taxa_datalist(Gene_ID %in% remove)

plan(multisession, workers = 6)

remove_tested <- furrr::future_map_lgl(remove, function(remove_gene) {
  abundance_tmp <- datalist_reduced$Count_Data %>%
    filter(Gene_ID == remove_gene) %>%
    #select_if(is.numeric) %>%
    select(-c(1:12)) %>%
    unlist()
  
  quantile_tmp <- quantile(abundance_tmp[abundance_tmp > 0], 0.05)
  
  test <- datalist_contaminant$Count_Data %>%
    filter(Gene_ID == remove_gene) %>%
    #select_if(is.numeric) %>%
    select(-c(1:12)) %>%
    unlist() > quantile_tmp
  
  return(ifelse(sum(test) == 0, FALSE, TRUE))
}, .progress = T)

plan(sequential)

remove_taxa_class <- read_csv("data/Contaminant_Taxa.csv")$Class %>% unique()

# Also, remove genes as contaminants that were classified as likely contaminant
# taxa. We selected taxonomic groups as potential contaminants that are either
# known contaminant taxa (e.g. Propionibacteriales) or taxa whose typical lifestyle
# are not compatible with hydrothermal, deep sediment  habitat (see manuscript).
removed_total <- datalist %>%
  filter_taxa_datalist(!(Gene_ID %in% remove[remove_tested])) %>%
  filter_taxa_datalist(!(Class %in% remove_taxa_class)) %>% 
  filter_taxa_datalist(!(Phylum %in% c("Bacillota", "Cyanobacteriota", "Bacteroidota", "Chlamydiota") | Order == "Propionibacteriales")) %>%
  filter_taxa_datalist(!(Class %in% read_csv("data/Contaminant_Taxa.csv")$Class)) %>%
  .$Count_Data %>%
  select_if(is.numeric) %>% colSums()

# Visualize proportion of contaminants per sample
enframe(removed_total/(datalist$Count_Data %>% select_if(is.numeric) %>% colSums())) %>%
  dplyr::rename("Sample_ID" = "name") %>%
  left_join(datalist$Meta_Data, by = "Sample_ID") %>%
  filter(Type == "Metagenome") %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Temperature)$Sample_ID))) %>%
  ggplot(aes(x = Sample_ID, y = value*100)) +
    geom_bar(stat = "identity") +
    facet_grid(~Core, scales = "free_x", space = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = .5, vjust = 1)) +
    labs(x = "", y = "Non-contaminated genes (%)") +
    ylim(c(0, 100))

ggsave("fig/Contaminant_removal.png", width = 6, height = 5, dpi = 300)

# Save decontaminated count-table
datalist %>%
  filter_taxa_datalist(!(Gene_ID %in% remove[remove_tested])) %>%
  filter_taxa_datalist(!(Class %in% remove_taxa_class)) %>% 
  filter_taxa_datalist(!(Phylum %in% c("Bacillota", "Cyanobacteriota", "Bacteroidota", "Chlamydiota") | Order == "Propionibacteriales")) %>%
  filter_taxa_datalist(!(Class %in% read_csv("data/Contaminant_Taxa.csv")$Class)) %>%
  .$Count_Data %>%
  write_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv")

