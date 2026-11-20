#### This script visualizes gene-sequence properties from Prodigal ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                "#CC79A7", "darkred", "forestgreen")

# Import gene-table from bioinformatic pipeline containing output of Prodigal
gene_table <- read_tsv("data/complete_gene_table_filtered.tsv") %>%
  mutate(Gene_length = Right_coord - Left_coord) %>%
  mutate(Contig_Name = str_replace_all(string = Header, pattern = " # .*$", replacement = "")) %>%
  left_join(., select(read_csv2("data/NCBI_SRA_Accession_Table.csv"), Run, `Sample Name`),  by = "Run") %>%
  dplyr::rename("Sample_ID" = "Sample Name") %>%
  left_join(., datalist$Meta_Data, by = "Sample_ID") %>%
  mutate(Temperature = ordered(Temperature)) %>%
  filter(!is.na(Temperature)) %>%
  mutate(Core = str_replace_all(Core, pattern = "^U", replacement = "") %>%
           str_replace_all(pattern = "B$", replacement = ""))

# Keep only decontaminated genes in each sample
gene_table_decontaminated <- purrr::map(datalist$Meta_Data$Sample_ID, function(sample_tmp) {
  
  genes_decontaminated <- datalist$Count_Data %>%
    select(Gene_ID, !!sample_tmp) %>%
    magrittr::set_colnames(c("Gene_ID", "Abundance")) %>%
    filter(Abundance > 0) %>%
    .$Gene_ID
  
  gene_table_tmp <- gene_table %>%
    filter(Sample_ID == sample_tmp) %>%
    filter(Gene_ID %in% genes_decontaminated)
  
  return(gene_table_tmp)
  
}) %>%
  bind_rows()

# Visualize start-codon distribution
gene_table_decontaminated %>%
  group_by(Sample_ID) %>%
  summarize(TTG = sum(Start_type == "TTG"),
            ATG = sum(Start_type == "ATG"),
            GTG = sum(Start_type == "GTG"),
            Edge = sum(Start_type == "Edge")) %>%
  rowwise() %>%
  mutate(TTG_prop = TTG/sum(c(TTG, ATG, GTG))) %>%
  mutate(ATG_prop = ATG/sum(c(TTG, ATG, GTG))) %>%
  mutate(GTG_prop = GTG/sum(c(TTG, ATG, GTG))) %>%
  select(-c(2:5)) %>%
  reshape2::melt() %>%
  magrittr::set_colnames(c("Sample_ID", "Start_Type", "Prop")) %>%
  left_join(., datalist$Meta_Data, by = "Sample_ID") %>%
  filter(!is.na(Temperature)) %>%
  mutate(Start_Type = ordered(Start_Type, levels = c("ATG_prop", "TTG_prop", "GTG_prop"), labels = c("ATG", "TTG", "GTG"))) %>%
  
  ggplot(aes(x = Temperature, y = Prop, col = Core)) +
  geom_point() +
  geom_smooth(se = F, span = .9) +
  scale_color_manual(values = cbbPalette[-c(2, 6)]) +
  facet_grid(Start_Type~., scales = "free") +
  labs(x = "Temperature (°C)", y = "Proportion of all start codons", col = "Core") +
  theme_bw() +
  theme(strip.text = element_text(size = 5),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = .5))

ggsave("fig/Start_Codons_Genes.png", width = 4, height = 5, dpi = 300)

# Visualize GC-content per sample
gene_table_decontaminated %>%
  ggplot(aes(x = Temperature, y = GC_cont)) +
  geom_boxplot() +
  facet_grid(.~Core, scales = "free_x", space = "free") +
  labs(x = "Temperature (°C)", y = "GC content of genes") +
  theme_bw() +
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 9))

ggsave("fig/GC_Content_Genes.png", width = 5, height = 4, dpi = 300)

# Visualize gene-length per sample
gene_table_decontaminated %>%
  ggplot(aes(x = Temperature, y = Gene_length)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_grid(.~Core, scales = "free_x", space = "free") +
  labs(x = "Temperature (°C)", y = "Gene length (bp)") +
  theme_bw() +
  theme(strip.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 9))

ggsave("fig/Length_Genes.png", width = 5, height = 4, dpi = 300)