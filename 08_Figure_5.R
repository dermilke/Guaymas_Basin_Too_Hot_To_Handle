#### This script infers abundances of selected functional genes and visualizes them ####
## The following functions are analyzed here:
## 1) General KEGG metabolic categories
## 2) Wood-Ljungdahl Pathway Genes
## 3) Alkane Degradation Genes
## 4) CAZyme Genes
## 5) Methanogenesis Genes
## 6) Sulfate Metabolism Genes
## 7) Hydrogen Metabolism Genes
##
## For each group we selected only those genes that are known core-genes of the specific
## pathway and are not associated to other metabolic roles (only associated to one KEGG module).
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)
library(furrr)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Create_Datatable_mod.R")
source("R/datatable_color.R")

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Import decontaminated datalist and select only functional genes with taxonomic annotation on 
# Kingdom level. Standardize each sample into rel. abundances.
datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv")) %>%
  filter_taxa_datalist(Kingdom %in% c("Archaea", "Bacteria")) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_station_datalist(Type == "Metagenome") 

## 1) General Metabolic Categories
# General categories were defined as second hierarchical level of KEGG orthology
# and only functional genes associated to a single KEGG module were considered here
general_modules <- read_tsv("data/KEGG_DB/Mapping_Module_KO.tsv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  left_join(., select(read_tsv("data/KEGG_DB/KEGG_Module_Hierarchy.tsv"), Level_B, Module) %>% distinct() %>%
              filter(!(Level_B %in% c("Module set", "Gene set", "Xenobiotics biodegradation"))))

# Transform selected functional genes into flat datatable
datatable_general <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(KEGG %in% general_modules$KO) %>%
  create_datatable(grpBy = KEGG, upper_grp = KEGG, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., dplyr::rename(general_modules, "Group" = "KO"))

# Visualize general KEGG categories
# Slightly change the temperature of one sample so that it will be rounded up
# for visualization purposes
p1 <- datatable_general %>%
  filter(!is.na(Level_B)) %>%
  filter(!is.na(Temperature)) %>%
  filter(Core %in% c("U1545B", "U1547B")) %>%
  mutate(Temperature = ifelse(Temperature == 17.5, 17.4, Temperature)) %>%
  mutate(Temperature = round(Temperature)) %>%
  
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = Level_B)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = read_csv("../ExCom/Data/Colors/Province_Colour.csv")$Colour) +
    theme_bw() +
    facet_grid(.~Core, scale = "free", space = "free_x") +
    labs(x = "Temperature (°C)", y = "Rel. Abundance (%)", fill = "", title = "General Categories") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          strip.text.x = element_text(size = 7, angle = 0),
          strip.text.y = element_text(size = 7, angle = 0)) +
    guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

## 2) Wood-Ljungdahl-Pathway
# Select all genes associated to KEGG module M00377 (module for Wood-Ljungdahl-Pathway)
# that were not associated to multiple KEGG modules.
wood_modules <- read_tsv("data/KEGG_DB/KEGG_Module_KO_list.csv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  filter(Module == "M00377")  %>%
  left_join(., select(read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv"), KO, KO_D) %>% distinct())

# Transform selected functional genes into flat datatable and standardize
# functional abundances by the energy-metabolism category from 1) above. 
# Resulting ratios are compositionally more stable than pure rel. abundances.
datatable_wood <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(KEGG %in% wood_modules$KO) %>%
  create_datatable(grpBy = KEGG, upper_grp = KEGG, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., wood_modules %>% dplyr::rename("Group" = "KO")) %>%
  left_join(., read_tsv("data/Energy_Metabolism_per_Sample.tsv")) %>%
  mutate(Abundance = Abundance / Abundance_Energy)

# Visualize Wood-Ljungdahl genes in stacked barplot
p2 <- datatable_wood %>%
  filter(!is.na(Temperature)) %>%
  filter(Core %in% c("U1545B", "U1547B")) %>%
  mutate(KO_D = str_replace_all(KO_D, pattern = " \\/ .*$", replacement = "")) %>%
  mutate(Temperature = ifelse(Temperature == 17.5, 17.4, Temperature)) %>%
  mutate(Temperature = round(Temperature)) %>%
  mutate(KO_D = str_replace_all(KO_D, pattern = "; .*$", replacement = "")) %>%
  
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = KO_D)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(cbbPalette, "darkred",  "forestgreen", "steelblue", "salmon")) +
  theme_bw() +
  facet_grid(.~Core, scale = "free", space = "free_x") +
  labs(x = "Temperature (°C)", y = "Energy Metabolism Genes (%)",title = "Wood Ljungdahl", fill = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text.x = element_text(size = 7, angle = 0),
        strip.text.y = element_text(size = 7, angle = 0)) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

## 2) Alkane-degradation genes
# Select known core-genes in alkane degradation pathways: 
# NAO, dhaA, ssuD, alkB1_2, ladA, pyruvate formate lyase
alkane_degradation_modules <- read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv") %>%
  filter(str_detect(string = .$KO_D, pattern = "NAO; nitroalkane oxidase") |
           str_detect(string = .$KO_D, pattern = "dhaA") |
           str_detect(string = .$KO_D, pattern = "ssuD") |
           str_detect(string = .$KO_D, pattern = "alkB1_2, alkM; alkane 1-monooxygenase") |
           str_detect(string = .$KO_D, pattern = "ladA") |
           str_detect(string = .$KO_D, pattern = "pyruvate formate lyase activating enzyme"))

# Transform selected functional genes into flat datatable and standardize
# functional abundances by the energy-metabolism category from 1) above. 
datatable_alkane <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(KEGG %in% alkane_degradation_modules$KO) %>%
  create_datatable(grpBy = KEGG, upper_grp = KEGG, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., alkane_degradation_modules %>% dplyr::rename("Group" = "KO")) %>%
  left_join(., read_tsv("data/Energy_Metabolism_per_Sample.tsv")) %>%
  mutate(Abundance = Abundance / Abundance_Energy)

# Visualize Alkane degradation genes in stacked barplot
p3 <- datatable_alkane %>%
  filter(!is.na(Temperature)) %>%
  filter(Core %in% c("U1545B", "U1547B")) %>%
  mutate(Temperature = ifelse(Temperature == 17.5, 17.4, Temperature)) %>%
  mutate(Temperature = round(Temperature)) %>%
  mutate(KO_D = str_replace_all(KO_D, pattern = "; .*$", replacement = "")) %>%
  mutate(KO_D = ifelse(KO_D == "pflA, pflC, pflE", "assD2", KO_D)) %>%
  
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = KO_D)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = read_csv("../ExCom/Data/Colors/Province_Colour.csv")$Colour[c(2, 1, 4, 6, 7, 5, 3)]) +
    theme_bw() +
    facet_grid(.~Core, scale = "free", space = "free_x") +
    labs(x = "Temperature (°C)", y = "Energy Metabolism Genes (%)",title = "Alkane Degradation", fill = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          strip.text.x = element_text(size = 7, angle = 0),
          strip.text.y = element_text(size = 7, angle = 0)) +
    guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

## 4) CAZyme Genes
# Collapse original datalist into CAZy-datalist based on CAZy-Families and
# transform it into flat datatable and standardize gene abundance by dividing
# it against total abundances of energy-metabolism genes
datalist_CAZY <- datalist
datalist_CAZY$Count_Data <- datalist_CAZY$Count_Data %>%
  mutate(CAZy_Family = str_replace_all(CAZy_Family, pattern = "^.*.fasta\\|", replacement = ""))

datatable_CAZy <- datalist_CAZY %>%
  filter_taxa_datalist(!is.na(CAZy_Family)) %>%
  create_datatable(grpBy = CAZy_Family, upper_grp = CAZy_Family, otherThreshold = 0, addOthers = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., read_tsv("data/Energy_Metabolism_per_Sample.tsv")) %>%
  mutate(Abundance = Abundance / Abundance_Energy)

# Visualize as stacked barplot
p4 <- datatable_CAZy %>%
  filter(!is.na(Temperature)) %>%
  filter(Core %in% c("U1545B", "U1547B")) %>%
  mutate(Temperature = ifelse(Temperature == 17.5, 17.4, Temperature)) %>%
  mutate(Temperature = round(Temperature)) %>%
  filter(Group %in% c("AA", "CBM", "CE", "GH", "GZ", "PL")) %>%
  
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = Group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("steelblue", cbbPalette[4], "orange", cbbPalette[9], "#486663", "#27dc9a")) +
    theme_bw() +
    facet_grid(.~Core, scale = "free", space = "free_x") +
    labs(x = "Temperature (°C)", y = "Energy Metabolism Genes (%)",title = "CAZyme Families", fill = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          strip.text.x = element_text(size = 5, angle = 0),
          strip.text.y = element_text(size = 7, angle = 0)) +
    guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

## 5) Methanogenesis genes
# Select core-genes in methanogenesis mcrA,B,C,D,G based on KEGG-IDs
methanogenesis_table <-  read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv") %>%
  filter(KO %in% c("K00399", "K00400", "K00401", "K00402", "K03421", "K03422")) %>%
  mutate(KO_D = str_replace_all(KO_D, pattern = ";.*", replacement = "")) %>%
  mutate(KO_D = ifelse(KO_D == "K00400", "mcrA2", KO_D))

# Transform selected functional genes into flat datatable and standardize
# functional abundances by the energy-metabolism category from 1) above. 
# Resulting ratios are compositionally more stable than pure rel. abundances.
datatable_methanogen <- datalist %>%
  filter_taxa_datalist(KEGG %in% methanogenesis_table$KO) %>%
  create_datatable(grpBy = KEGG, upper_grp = KEGG, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., dplyr::rename(methanogenesis_table, "Group" = "KO"))  %>%
  left_join(., read_tsv("data/Energy_Metabolism_per_Sample.tsv")) %>%
  mutate(Abundance = Abundance / Abundance_Energy)

# Visualize as stacked barplot
p5 <- datatable_methanogen %>%
  filter(!is.na(Temperature)) %>%
  filter(Core %in% c("U1545B", "U1547B")) %>%
  mutate(Temperature = ifelse(Temperature == 17.5, 17.4, Temperature)) %>%
  mutate(Temperature = round(Temperature)) %>%
  
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = KO_D)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("steelblue", "forestgreen", "darkred", "salmon", "red", "orange")) +
  theme_bw() +
  facet_grid(.~Core, scale = "free", space = "free_x") +
  labs(x = "Temperature (°C)", y = "Energy Metabolism Genes (%)",title = "Methanogenesis", fill = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text.x = element_text(size = 5, angle = 0),
        strip.text.y = element_text(size = 7, angle = 0)) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

## 6) Sulfate metabolism genes
# Select the three core genes aprA, dsrA and dsrB for sulfate metabolism.
sulfate_reduction_modules <- read_tsv("data/KEGG_DB/Mapping_Module_KO.tsv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  left_join(., read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv")) %>%
  mutate(KO_D = str_replace_all(KO_D, pattern = ";.*$", replacement = "")) %>%
  filter(KO_D %in% c("aprA", "dsrA", "dsrB", "qmoB")) 

# Transform selected functional genes into flat datatable and standardize
# functional abundances by the energy-metabolism category from 1) above. 
# Resulting ratios are compositionally more stable than pure rel. abundances.
datatable_sulfate <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(KEGG %in% sulfate_reduction_modules$KO) %>%
  create_datatable(grpBy = KEGG, upper_grp = KEGG, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., dplyr::rename(sulfate_reduction_modules, "Group" = "KO")) %>%
  left_join(., read_tsv("data/Energy_Metabolism_per_Sample.tsv")) %>%
  mutate(Abundance = Abundance / Abundance_Energy)

# Visualize in stacked barplot
p6 <- datatable_sulfate %>%
  filter(Core %in% c("U1545B", "U1547B")) %>%
  filter(!is.na(Temperature)) %>%
  left_join(., read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv") %>% dplyr::rename("Group" = "KO"))  %>%
  mutate(Temperature = ifelse(Temperature == 17.5, 17.4, Temperature)) %>%
  mutate(Temperature = round(Temperature)) %>%
  
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = KO_D)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = cbbPalette[c(5, 4, 6)]) +
    theme_bw() +
    facet_grid(.~Core, scale = "free", space = "free_x") +
    labs(x = "Temperature (°C)", y = "Energy Metabolism Genes (%)", fill = "", title = "Sulfate Metabolism") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          strip.text.x = element_text(size = 7, angle = 0),
          strip.text.y = element_text(size = 7, angle = 0)) +
    guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

## 7) Hydrogen metabolism genes
# Select the different KOs for iron-dependent hydrogenase genes
hydrogen_modules <- read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  filter(KO %in% c("K00437", "K18008", "K17997", "K17998", "K17999")) %>%
  #mutate(KO_D = ifelse(KO == "K00437", paste0("[NiFe] hydrogenase ", KO_D), KO_D)) %>%
  #mutate(KO_D = ifelse(KO == "K18008", paste0("[NiFe] hydrogenase ", KO_D), KO_D)) %>%
  #mutate(KO_D = ifelse(KO == "K17997", paste0("[FeFe] hydrogenase ", KO_D), KO_D)) %>%
  #mutate(KO_D = ifelse(KO == "K17998", paste0("[FeFe] hydrogenase ", KO_D), KO_D)) %>%
  #mutate(KO_D = ifelse(KO == "K17999", paste0("[FeFe] hydrogenase ", KO_D), KO_D)) %>%
  mutate(KO_D = str_replace_all(KO_D, pattern = ";.*", replacement = ""))

# Transform selected functional genes into flat datatable and standardize
# functional abundances by the energy-metabolism category from 1) above. 
# Resulting ratios are compositionally more stable than pure rel. abundances.
datatable_hydrogen <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(KEGG %in% hydrogen_modules$KO) %>%
  create_datatable(grpBy = KEGG, upper_grp = KEGG, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others") %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*", replacement = "")) %>%
  left_join(., hydrogen_modules %>% dplyr::rename("Group" = "KO")) %>%
  left_join(., read_tsv("data/Energy_Metabolism_per_Sample.tsv")) %>%
  mutate(Abundance = Abundance / Abundance_Energy)

# Visualize in stacked barplot
p7 <- datatable_hydrogen %>%
  filter(!is.na(Temperature)) %>%
  filter(Core %in% c("U1545B", "U1547B")) %>%
  mutate(Temperature = ifelse(Temperature == 17.5, 17.4, Temperature)) %>%
  mutate(Temperature = round(Temperature)) %>%
  
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = KO_D)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("steelblue", "forestgreen", "darkred", cbbPalette)) +
  theme_bw() +
  facet_grid(.~Core, scale = "free", space = "free_x") +
  labs(x = "Temperature (°C)", y = "Energy Metabolism Genes (%)",title = "Hydrogen Metabolism", fill = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text.x = element_text(size = 7, angle = 0),
        strip.text.y = element_text(size = 7, angle = 0)) +
  guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))

cowplot::plot_grid(p1 + theme(legend.position = "none"),
                   cowplot::plot_grid(p2 + theme(legend.position = "none"), p3  + theme(legend.position = "none"), p4  + theme(legend.position = "none"),
                                      p5  + theme(legend.position = "none"), p6  + theme(legend.position = "none"), p7  + theme(legend.position = "none"), ncol = 3),
                   rel_widths = c(0.3, 0.7))

ggsave("fig/Function_Composite.png", width = 10, height = 6, dpi = 300)

ggsave(plot = cowplot::get_legend(p1), filename = "fig/Legend_p1.png")
ggsave(plot = cowplot::get_legend(p2), filename = "fig/Legend_p2.png")
ggsave(plot = cowplot::get_legend(p3), filename = "fig/Legend_p3.png")
ggsave(plot = cowplot::get_legend(p4), filename = "fig/Legend_p4.png")
ggsave(plot = cowplot::get_legend(p5), filename = "fig/Legend_p5.png")
ggsave(plot = cowplot::get_legend(p6), filename = "fig/Legend_p6.png")
ggsave(plot = cowplot::get_legend(p7), filename = "fig/Legend_p7.png")