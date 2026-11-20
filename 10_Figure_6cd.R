#### This script visualizes composition of lipid degradation genes and ####
#### Lipid composition in DOM data based on compound classification against LipidMap DB ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Create_Datatable_mod.R")

# Load DOM and gene datalists
datalist_DOM <- list(Count_Data = read_tsv("data/DOM/Final_Count_Data.tsv"),
                     Meta_Data = read_tsv("data/DOM/Final_Meta_Data.tsv"))

datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv"))

# Select all genes associated to KEGG module M00086 and M00087 (modules for lipid degradation)
# that were not associated to multiple KEGG modules.
lipid_degradation_modules <- read_tsv("data/KEGG_DB/Mapping_Module_KO.tsv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  filter(Module %in%  c("M00086", "M00087")) %>%
  left_join(., select(read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv"), KO, KO_D) %>% distinct())

# Transform selected functional genes into flat datatable and standardize
# functional abundances by the energy-metabolism category (see script for Figure 5). 
# Resulting ratios are compositionally more stable than pure rel. abundances.
datatable_lipids <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(KEGG %in% lipid_degradation_modules$KO) %>%
  create_datatable(grpBy = KEGG, upper_grp = KEGG, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., dplyr::rename(lipid_degradation_modules, "Group" = "KO")) %>%
  left_join(., read_tsv("data/Energy_Metabolism_per_Sample.tsv")) %>%
  mutate(Abundance = Abundance / Abundance_Energy)

# Visualize lipid degradation genes in stacked barplot
datatable_lipids %>%
  filter(Core %in% c("U1545B", "U1547B")) %>%
  mutate(Core = str_replace_all(Core, pattern = "B$", replacement = "")) %>%
  filter(!is.na(Temperature)) %>%
  mutate(KO_D = str_replace_all(KO_D, pattern = ";.*$", replacement = "")) %>%
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = KO_D)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = read_csv("../ExCom/Data/Colors/Province_Colour.csv")$Colour[-2]) +
  theme_bw() +
  facet_grid(.~Core, scale = "free_x", space = "free_x") +
  labs(x = "Temperature (°C)", y = "Rel. Abundance (%)", title = "Lipid Degradation", fill = "") +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = .5),
        axis.title = element_text(size = 12)) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

ggsave("fig/Lipid_Degradation.png", width = 6, height = 3, dpi = 300)

# DOM compound classification based on LIPID-MAPS. Use only lipid sumformulas that
# were present only once in the database.
lipid_DB <- read_tsv("data/lipid_formulas_LIPID_MAPS.tsv") %>%
  distinct() %>%
  filter(Formula %in% {table(Formula) %>% enframe() %>% filter(value < 2) %>% .$name})

# Classify DOM compounds by sumformulas
datalist_DOM$Count_Data <- datalist_DOM$Count_Data %>%
  left_join(., lipid_DB %>% dplyr::rename("OTU_ID" = "Formula")) %>%
  filter(!is.na(Lipid_Category))

# Transform intensity matrix to flat datatable format
datatable_DOM_lipids <- list(Count_Data = datalist_DOM$Count_Data,
                             Meta_Data = datalist_DOM$Meta_Data) %>%
  create_datatable(grpBy = Lipid_Category, upper_grp = Lipid_Category, otherThreshold = 0, addColorScheme = F, addOthers = F) %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  filter(Core %in% c("U1545", "U1547"))

# Visualize composition of lipids in stacked barplot
datatable_DOM_lipids %>%
  filter(!(Group %in% c("Others", "NA"))) %>%
  
  ggplot(aes(x = as.factor(Temperature_C), y = Abundance*100, fill = Group)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(.~Core, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c("#0077BB", "#33BBEE", "#009988", "#EE7733", "#CC3311", "#EE3377", "#AA4499", "#FFC20A")) +
    theme_bw() +
    theme(strip.text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = .5),
          axis.title = element_text(size = 12)) +
    labs(x = "Temperature (°C)", y = "Rel. Intensity (%)", fill = "Lipid Category")

ggsave("fig/Lipid_DOM_Compounds.png", width = 8, height = 4.5, dpi = 300)
