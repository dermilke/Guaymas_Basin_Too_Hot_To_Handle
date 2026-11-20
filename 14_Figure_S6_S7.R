#### This script visualizes taxonomic composition of selected functional genes ####
#### of Alkane-degradation and Aromatics-degradation pathways ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Create_Datatable_mod.R")

# Read gene-datalist
datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv")) %>%
  filter_taxa_datalist(Kingdom %in% c("Archaea", "Bacteria")) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_station_datalist(Type == "Metagenome")

# Select core-genes of the alkane-degradation pathway (see Figure 5)
alkane_KOs <- tibble(KO = c("K00496", "K20938", "K19823", "K04091", "K01563",
                            "K04069"),
                     Alkane_Type = c("alkB\nalkane-1\nmonooxygenase",
                              "ladA\nlong-chain alkane\nmonooxygenase",
                              "NAO\nnitroalkane\noxidase",
                              "ssuD\nalkanesulfonate\nmonooxygenase", 
                              "dhaA\nhaloalkane\ndehalogenase",
                              "assD2\nAlkylsuccinate\nsynthase (II) "))

# Calculate for each selected core-gene the taxonomic composition along the
# temperature gradient
datatable_alkanes <- purrr::map(alkane_KOs$KO, function(KO_tmp) {
  
  datatable_tmp <- datalist %>%
    filter_taxa_datalist(KEGG == KO_tmp) %>%
    create_datatable(grpBy = Class, upper_grp = Phylum, otherThreshold = 0.03, addColorScheme = T) 
  
  datatable_tmp$table <- datatable_tmp$table %>%
    mutate(KO_Type = KO_tmp) %>%
    mutate(Group = as.character(Group))
  
  return(datatable_tmp)
  
})

# Create a color-vector for the taxonomic composition
colors_new <- purrr::map(datatable_alkanes, "table") %>%
  bind_rows() %>%
  mutate(Group = str_replace_all(Group, pattern = " - Unknown .*", replacement = " - Unknown")) %>%
  mutate(Group = str_replace_all(Group, pattern = "^Unknown .* - ", replacement = "Unknown - ")) %>%
  add_colorscheme(upper_grp = "Phylum") %>%
  mutate(Group = ifelse(Phylum == "Other Phylum", "Other Phylum", paste0(Phylum, " - ", lower_grp)))

# Combine the taxonomic composition for all alkane-degradation genes along the
# temperature gradient and adapt names for plotting
table_alkanes <- purrr::map(datatable_alkanes, "table") %>%
  bind_rows() %>%
  left_join(., dplyr::rename(alkane_KOs, "KO_Type" = "KO"), by = "KO_Type") %>%
  filter(!is.na(Temperature)) %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Temperature)$Sample_ID))) %>%
  mutate(Group = str_replace_all(Group, pattern = " - Unknown .*", replacement = " - Unknown")) %>%
  mutate(Group = str_replace_all(Group, pattern = "^Unknown .* - ", replacement = "Unknown - ")) %>%
  mutate(Group = ordered(Group, levels = colors_new$Group)) %>%
  mutate(Alkane_Type = ordered(Alkane_Type, levels = group_by(., Alkane_Type) %>% summarize(Abundance = mean(Abundance)) %>%
                                 arrange(desc(Abundance)) %>% .$Alkane_Type))

# Visualize alkane-degradation taxonomy with stacked barplots
table_alkanes %>%
    ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = Group)) +
    geom_bar(stat = "identity", position = "stack", col = "black", linewidth = .1) +
    scale_fill_manual(values = colors_new$Color) +
    facet_grid(Alkane_Type~Core, scales = "free", space = "free_x") +
    theme_bw() +
    guides(fill = guide_legend(ncol = 1, 
                               keywidth = 0.8, keyheight = 0.7)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          legend.text = element_text(size = 8),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(angle = 0)) +
    labs(x = "Temperature (°C)", y = "Abundance (%)", title = "Alkane degradation genes",
         fill = "Phylum - Class")
  
ggsave("fig/Alkane_Genes_Taxonomy.png", width = 12, height = 8, dpi = 300)

# Select core-genes of the aromatics-degradation for taxonomic analyses
aromatic_KOs <- tibble(KO = c("K03381", "K14581", "K14578", "K14579", "K14580",
                              "K11946", "K14585", "K18257", "K07540",
                              "K00446", "K07104", "K16270"),
                       Enzyme_Type = c("catA\ncatechol\ndioxygenase", 
                                       "nahAa\nnaphthalene\ndioxygenase",
                                       "nahAb\nnaphthalene\ndioxygenase",
                                       "nahAc\nnaphthalene\ndioxygenase",
                                       "nahAd\nnaphthalene\ndioxygenase",
                                       "phdG\nhydratase\naldolase",
                                       "nahE\nhydroxybenzylidenepyruvate\nhydratase-aldolase",
                                       "phdE\ndihydrophenanthrene\ndehydrogenase",
                                       "benzylsuccinate\nsynthase",
                                       "dmpB\ncatechol\ndioxygenase",
                                       "catE\ncatechol\ndioxygenase",
                                       "todE\nmethylcatechol\ndioxygenase"))

# Calculate for each selected core-gene the taxonomic composition along the
# temperature gradient
datatable_aromatics <- purrr::map(aromatic_KOs$KO, function(KO_tmp) {
  
  datalist_tmp <- datalist %>%
    filter_taxa_datalist(KEGG == KO_tmp)
  
  datalist_tmp$Count_Data <- datalist_tmp$Count_Data %>%
    mutate_if(is.character, function(x) str_replace_all(x, pattern = "^Unknown .*", replacement = "Unknown")) %>%
    mutate_if(is.character, function(x) str_replace_all(x, pattern = "^candidate", replacement = "Candidate"))
  
  if (nrow(datalist_tmp$Count_Data) > 0) {
    
    datatable_tmp <- create_datatable(datalist_tmp, grpBy = Class, upper_grp = Phylum, otherThreshold = 0.05, addColorScheme = T) 
    
    datatable_tmp$table <- datatable_tmp$table %>%
      mutate(KO_Type = KO_tmp) %>%
      mutate(Group = as.character(Group))
    
  } else {
    datatable_tmp = NULL
  }   
  
  return(datatable_tmp)
  
})

# Create a color-vector for the taxonomic composition
colors_new <- purrr::map(datatable_aromatics, "table") %>%
  bind_rows() %>%
  add_colorscheme(upper_grp = "Phylum") %>%
  mutate(Group = ifelse(Phylum == "Other Phylum", "Other Phylum", paste0(Phylum, " - ", lower_grp))) 

# Combine the taxonomic composition for all aromatics-degradation genes along the
# temperature gradient and adapt names for plotting
table_aromatics <- purrr::map(datatable_aromatics, "table") %>%
  bind_rows() %>%
  left_join(., dplyr::rename(aromatic_KOs, "KO_Type" = "KO"), by = "KO_Type") %>%
  filter(!is.na(Temperature)) %>%
  mutate(Sample_ID = ordered(Sample_ID, levels = unique(arrange(., Temperature)$Sample_ID))) %>%
  mutate(Group = ordered(Group, levels = colors_new$Group)) %>%
  mutate(Enzyme_Type = ordered(Enzyme_Type, levels = {group_by(., Enzyme_Type) %>% summarize(Abundance = mean(Abundance)) %>%
                                 arrange(desc(Abundance)) %>% .$Enzyme_Type}[c(1, 8, 5, 9, 4, 7, 2, 6, 3, 10)]))

# Visualize aromatics-degradation taxonomy with stacked barplot
table_aromatics %>%
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", col = "black", linewidth = .1) +
  scale_fill_manual(values = colors_new$Color) +
  facet_grid(Enzyme_Type~Core, scales = "free", space = "free_x") +
  theme_bw() +
  guides(fill = guide_legend(ncol = 1, 
                             keywidth = 0.8, keyheight = 0.7)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.text = element_text(size = 8),
        strip.text.x = element_text(size = 5),
        strip.text.y = element_text(angle = 0)) +
  labs(x = "Temperature (°C)", y = "Abundance (%)", title = "Aromatics degradation genes",
       fill = "Phylum - Class")

ggsave("fig/Aromatics_Genes_Taxonomy.png", width = 12, height = 10, dpi = 300)