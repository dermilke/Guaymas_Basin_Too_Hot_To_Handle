#### This script analyses output of bioinformatic pipeline. It analyses gene-origins ####
#### by assigning each gene to the sample it was found in. Gene cluster containing several ####
#### genes are associated to consensus sequence - which is the longest sequence here ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/unfold_kegg.R")
source("R/Create_Datatable_mod.R")
source("R/datatable_color.R")

# First load the complete gene table into R - this takes some time and will
# likely flood your memory as the table has several gigabyte. Hence,
# we only take the gene-origin sample_ID for each gene_ID and save it in an extra
# file that is much smaller. Then remove the complete_gene_table from memory.
data <- read_tsv("data/Gene_Library/complete_gene_table_filtered.tsv")

meta_data <- read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv") %>%
  left_join(., read_csv2("data/NCBI_SRA_Accession_Table.csv") %>%
              select(Run, 'Sample Name') %>%
              dplyr::rename("Sample_ID" = "Sample Name"), by = "Sample_ID") %>%
  filter(!is.na(Run)) %>%
  filter(Type == "Metagenome")

data %>%
  select(Gene_ID, Run) %>%
  write_tsv("data/Gene_Origin.tsv")

rm("data")

# The resulting Gene_Origin table is much smaller and better to handle
gene_origin_table <- read_tsv("data/Gene_Origin.tsv")

# Load metagenome datalist
datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv")) %>%
  filter_taxa_datalist(Kingdom %in% c("Archaea", "Bacteria")) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_station_datalist(Type == "Metagenome")

# Assign an origin-sample to each gene_ID
datalist$Count_Data <- datalist$Count_Data %>%
  left_join(., gene_origin_table, by = "Gene_ID") %>%
  left_join(., select(meta_data, Run, Sample_ID), by = "Run")  %>%
  dplyr::rename("Origin" = "Sample_ID") %>%
  mutate_if(is.character, function(x) str_replace_all(x, pattern = "^Unknown .*", replacement = "Unknown"))

# Free memory of large table
rm(gene_origin_table)

# Transform count-table into flat datatable
datatable_origin <- datalist %>%
  create_datatable(grpBy = Origin, upper_grp = Origin, otherThreshold = 0, addOthers = F, addColorScheme = F)

# Define temperature-intervals (according to cluster-distribution)
breaks <- c(0, 20, 45, 70)

datatable_origin_mod <- datatable_origin %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., select(meta_data, Sample_ID, Type, Temperature) %>%
              dplyr::rename("Group" = "Sample_ID", "Origin_Type" = "Type", "Temp_Origin" = "Temperature"), by = "Group") %>%
  mutate(Origin_Temp_Grp = cut(Temp_Origin, breaks = breaks, include.lowest = TRUE, labels = c("0-20", "20-45", ">45")))

# Visualize gene origins according to temperature-intervals in a stacked barplot
datatable_origin_mod %>%
  filter(!is.na(Origin_Temp_Grp)) %>%
  filter(!is.na(Temperature)) %>%
  group_by(Sample_ID) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Temperature)$Sample_ID),
                             labels = distinct(arrange(select(., Sample_ID, Temperature), Temperature))$Temperature)) %>%
  ggplot(aes(x = Sample_ID, y = Abundance*100, fill = Origin_Temp_Grp)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cbbPalette <- c("grey", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                                             "#CC79A7", "darkred", "forestgreen")[c(1, 7, 4)]) +
  facet_grid(.~Core, scales = "free_x", space = "free") +
  labs(x = "Temperature (°C)", fill = "Gene Origin:", y = "Rel. Abundance (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text = element_text(size = 5))

ggsave("fig/Gene_Origins.png", width = 7.5, height = 5.5, dpi = 300)
