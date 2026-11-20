#### This script finds DOM compounds matching selected aromatics degradation pathways ####
#### Pathways, including compounds and KOs, were manually selected based on KEGG description ####
#### and detailed pathways are shown in manuscript figures ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Create_Datatable_mod.R")

# Import datalists for DOM and functional genes
datalist_DOM <- list(Count_Data = read_tsv("data/DOM/Final_Count_Data.tsv"),
                     Meta_Data = read_tsv("data/DOM/Final_Meta_Data.tsv"))

datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv")) %>%
  filter_station_datalist(Type == "Metagenome") %>%
  filter_taxa_datalist(Kingdom %in% c("Bacteria", "Archaea")) %>%
  mutate_count_datalist(function(x) x/sum(x)) 

# Load manually curated table with aromatic degradation pathways, including KEGG IDs,
# Modules and Compounds (substrates and products)
aromatics_degradation_pathways <- read_csv("data/Selected_Aromatic_Pathways.csv")

# Filter DOM data for compounds by matching sumformulas and transform table into flat datatable format
datatable_DOM_aromatics <- datalist_DOM %>%
  filter_taxa_datalist(OTU_ID %in% unique(c(aromatics_degradation_pathways$Substrate, aromatics_degradation_pathways$Product))) %>%
  create_datatable(grpBy = OTU_ID, upper_grp = OTU_ID, otherThreshold = 0, addColorScheme = F, addOthers = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = ".* - ", replacement = "")) %>%
  left_join(., compounds_aromatics %>% dplyr::rename("Group" = "Formula"))

# Iterate over each module and plot intensity of compound matches in a scatter plot
purrr::map(unique(aromatics_degradation_pathways$Module), function(module_tmp) {
  
  compounds_tmp <- unique(c(filter(aromatics_degradation_pathways, Module == !!module_tmp)$Product, 
                            filter(aromatics_degradation_pathways, Module == !!module_tmp)$Substrate))
  
  datatable_DOM_aromatics %>%
    filter(Group %in% compounds_tmp) %>%
    filter(Core %in% c("U1545", "U1547")) %>%
    
    ggplot(aes(x = (Temperature_C), y = Abundance*100)) +
      geom_smooth(se = F, col = "darkgrey") +
      geom_point(aes(fill = Core), pch = 21, size = 2.5) +
      scale_fill_manual(values = read_csv("../ExCom/Data/Colors/Province_Colour.csv")$Colour[-c(1)]) +
      facet_grid(.~Group, scales = "free_x", space = "free_x") +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            strip.text = element_text(size = 12, angle = 0)) +
      labs(x = "Temperature (°C)", y = "Rel. Intensity (%)")
  
  ggsave(paste0("fig/Aromatics_Degradation_Compounds_Module_", module_tmp, ".png"), 
         width = 5, height = 5, dpi = 300)
  
})

# Filter functional-genes datalist for selected KEGG-Orthologues
datatable_aromatics <- datalist %>%
  filter_taxa_datalist(KEGG %in% aromatics_degradation_pathways$KO) %>%
  create_datatable(grpBy = KEGG, upper_grp = KEGG, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., aromatics_degradation_pathways %>% dplyr::rename("Group" = "KO"))

# Iterate over all aromatics degradation modules and visualize rel. abundance of KOs associated to
# those modules in a stacked barplot
purrr::map(unique(aromatics_degradation_pathways$Module), function(module_tmp) {
  datatable_aromatics %>%
    filter(Module == "M00534") %>%
    filter(Core %in% c("U1545B", "U1547B")) %>%
    ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = Group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = read_csv("../ExCom/Data/Colors/Province_Colour.csv")$Colour[-9]) +
    facet_grid(.~Core, scales = "free_x", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    labs(x = "Temperature (°C)", y = "Abundance (%)", fill = "KO")
  
  ggsave(paste0("fig/Toluene_aerobic_Degradation_Genes_Module_", module_tmp, ".png"), width = 6, height = 3.5, dpi = 300)
})
