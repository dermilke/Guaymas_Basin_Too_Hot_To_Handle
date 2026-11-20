#### This script visualizes composition of aromatic degradation genes and ####
#### composition of DOM compounds that are associated to aromatic degradation modules within KEGG ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)
library(KEGGREST)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Create_Datatable_mod.R")
source("R/datatable_color.R")
source("R/datatable_color.R")
source("R/Colors_Script.R")

find_matching_compounds <- function(Modules) {
  
  potential_compounds <- read_tsv("data/KEGG_DB/Mapping_Compound_Module.tsv") %>%
    filter(Module %in% Modules)
  
  compounds <- map(seq(1, nrow(potential_compounds)), function(ind) {
    tmp <- keggGet(potential_compounds$Compound[ind])
    Sys.sleep(0.1)
    
    return(tibble(Module = potential_compounds$Module[ind],
                  Formula = tmp[[1]]$FORMULA))
    
  }) %>%
    bind_rows() %>%
    left_join(., select(read_tsv("data/KEGG_DB/KEGG_Module_Hierarchy.tsv"), Level_D, Module) %>% distinct()) %>%
    filter(!is.na(Level_D)) %>%
    mutate(Level_D = str_replace_all(string = Level_D, pattern = ", .*", replacement = ""))
    
  return(compounds)
  
}

datalist_DOM <- list(Count_Data = read_tsv("data/DOM/Final_Count_Data.tsv"),
                     Meta_Data = read_tsv("data/DOM/Final_Meta_Data.tsv"))

datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv"))

# Select all genes associated to KEGG modules from map01220 (degradation of aromatic compounds)
# that were not associated to multiple KEGG modules.
aromatics_degradation_modules <- read_tsv("data/KEGG_DB/Mapping_Module_KO.tsv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  filter(Module %in%  c("M00538",       "M00537",       "M00419",       "M00547",       "M00548",       "M00551",       "M00568",       "M00569",       "M00539",       "M00543",       "M00544",
                        "M00418",       "M00541",       "M00540",       "M00534",       "M00624",       "M00623",       "M00636",       "M00545")) %>%
  left_join(., select(read_tsv("data/KEGG_DB/KEGG_Module_Hierarchy.tsv"), Level_D, Module) %>% distinct()) %>%
  mutate(Level_D = str_replace_all(string = Level_D, pattern = ", .*", replacement = ""))

# Find matching DOM compounds based on sumformulas. Use KEGGREST to scan KEGG database.
compounds_aromatics <- find_matching_compounds(aromatics_degradation_modules$Module) %>%
  filter(Formula %in% {table(.$Formula) %>% enframe() %>% filter(value < 2) %>% .$name})

# Filter out DOM compound matches with KEGG database and transform table into flat datatable format
datatable_DOM_aromatics <- datalist_DOM %>%
  filter_taxa_datalist(OTU_ID %in% compounds_aromatics$Formula) %>%
  create_datatable(grpBy = OTU_ID, upper_grp = OTU_ID, otherThreshold = 0, addColorScheme = F, addOthers = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = ".* - ", replacement = "")) %>%
  left_join(., compounds_aromatics %>% dplyr::rename("Group" = "Formula"))

# Define overarching colorscheme for both, DOM and gene data
aromatics_colors <- tibble(Group = sort(unique(aromatics_degradation_modules$Level_D)),
                           Color = read_csv("../ExCom/Data/Colors/Province_Colour.csv")$Colour[1:16]) %>%
  filter(Group != "Benzoyl-CoA degradation")

# Add vertical lines to DOM plot do mark temperature thresholds between clusters
vlines <- tibble(x = c(7.5, 17.5, 25.5, 2.5, 7.5, 11.5),
                 Core = c("U1545", "U1545", "U1545", "U1547", "U1547", "U1547"))

# Visualize composition of DOM compounds associated to aromatics degradation as stacked barplot
datatable_DOM_aromatics %>%
  filter(Core %in% c("U1545", "U1547")) %>%
  ggplot(aes(x = as.factor(Temperature_C), y = Abundance*100, fill = Level_D)) +
  geom_vline(aes(xintercept = x), data = vlines, col = "black", linetype = "dashed") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values =  dplyr::slice(aromatics_colors, match(sort(unique(datatable_DOM_aromatics$Level_D)), Group))$Color) +
  facet_grid(.~Core, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = .5),
        axis.title = element_text(size = 12)) +
  labs(x = "Temperature (°C)", y = "Rel. Intensity (%)", fill = "")

ggsave("fig/DOM_Aromatics_Degradation.png", width = 8, height = 4.5, dpi = 300)

# Transform selected functional genes into flat datatable and standardize
# functional abundances by the energy-metabolism category (see script for Figure 5). 
# Resulting ratios are compositionally more stable than pure rel. abundances.
datatable_aromatics <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(KEGG %in% aromatics_degradation_modules$KO) %>%
  create_datatable(grpBy = KEGG, upper_grp = KEGG, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., aromatics_degradation_modules %>% dplyr::rename("Group" = "KO")) %>%
  left_join(., read_tsv("data/Energy_Metabolism_per_Sample.tsv")) %>%
  mutate(Abundance = Abundance / Abundance_Energy)

# Visualize aromatics degradation genes in stacked barplot
datatable_aromatics %>%
  filter(Core %in% c("U1545B", "U1547B")) %>%
  filter(!(Level_D %in% c("Benzoyl-CoA degradation", "Homoprotocatechuate degradation", "Methane oxidation",
                          "Trans-cinnamate degradation"))) %>%
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = Level_D)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = read_csv("../ExCom/Data/Colors/Province_Colour.csv")$Colour[c(-16, -17)]) +
  theme_bw() +
  facet_grid(.~Core, scale = "free", space = "free_x") +
  labs(x = "Temperature (°C)", y = "Rel. Abundance (%)", fill = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text.x = element_text(size = 7, angle = 0),
        strip.text.y = element_text(size = 7, angle = 0)) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

ggsave("fig/Aromatic_Degradation_Genes.png", width = 7, height = 4, dpi = 300)
