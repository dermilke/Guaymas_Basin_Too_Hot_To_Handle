#### This script prepares output of bioinformatic pipeline for further analyses ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

# Read raw taxonomic classification and merge with count table from read recruitment
# Normalization of counts was already done as part of bioinformatic pipeline (see manuscript methods)
taxonomy_raw <- read_tsv("data/Gene_Library/raw/Complete_Classification_Export.tsv")

counts_raw <- read_tsv("data/Gene_Library/raw/Count_Table_normalized.tsv")

count_data_raw <- right_join(taxonomy_raw, counts_raw, by = "Gene_ID")

# Prepare associated sample meta-data according to count-table sample order
meta_data_raw <- read.csv2("data/meta_data_complete_MG.csv") %>%
  as_tibble() %>%
  mutate_at(2:15, as.numeric)

sample_accession <- read_csv2("data/NCBI_SRA_Accession_Table.csv") %>%
  select(Run, `Sample Name`)

sample_accession$`Sample Name`[match(colnames(count_data_raw), sample_accession$Run)][-c(1:14)]

count_data_processed <- count_data_raw %>%
  select(-evalue, -score) %>%
  magrittr::set_colnames(c("Gene_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
                           "UniRef_ID", "CAZy_ID", "CAZy_Family", "KEGG", 
                           sample_accession$`Sample Name`[match(colnames(count_data_raw), sample_accession$Run)][-c(1:14)])) %>%
  mutate_at(2:8, function(x) ifelse(is.na(x), "Unclassified", x))

meta_data_processed <- meta_data_raw %>%
  dplyr::slice(match(names(count_data_processed), Sample_ID)) %>%
  mutate(Type = ifelse(str_detect(string = Sample_ID, pattern = "_MT$"), "Metatranscriptome", 
                       ifelse(str_detect(string = Sample_ID, pattern = "^U15"), "Metagenome", "Control"))) %>%
  mutate(Core = str_split_fixed(string = Sample_ID, pattern = "_", n = 2)[,1])

# Save processed data
write_tsv(count_data_processed, "data/Gene_Library/Processed/Full_Gene_Count.tsv")
write_tsv(meta_data_processed, "data/Gene_Library/Meta_Data/Meta_Data.tsv")
