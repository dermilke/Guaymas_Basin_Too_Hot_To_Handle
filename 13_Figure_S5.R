#### This script visualizes taxonomic composition of all functional gene categories ####
#### used in Figure 5 and 6, summarized by temperature-intervals matching cluster distribution ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Create_Datatable_mod.R")
source("R/datatable_color.R")
source("R/create_others.R")
source("R/datatable_color.R")

# Import gene datalist
datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv")) %>%
  filter_taxa_datalist(Kingdom %in% c("Archaea", "Bacteria")) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_station_datalist(Type == "Metagenome") %>%
  mutate_meta_datalist(Temp_Grp = cut(Temperature, c(0, 20, 45, 70)))

# Adapt taxonomic names: Collapse all "Unknown [...]" taxa into "Unknown" groups 
datalist$Count_Data <- datalist$Count_Data %>%
  mutate_if(is.character, function(x) str_replace_all(x, pattern = "Unknown .*$", replacement = "Unknown"))

# Select core functional genes based on the same gene selection as done
# in Figure 5 for all different functional categories:
general_modules <- read_tsv("data/KEGG_DB/Mapping_Module_KO.tsv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  left_join(., select(read_tsv("data/KEGG_DB/KEGG_Module_Hierarchy.tsv"), Level_B, Module) %>% distinct()) %>%
  mutate(Type = "General")

sulfate_reduction_modules <- read_tsv("data/KEGG_DB/Mapping_Module_KO.tsv") %>%
  filter(Module %in% c("M00176", "M00596", "M00986", "M00595", "M00985", "M00990",
                       "M00991", "M00992", "M00993")) %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  left_join(., select(read_tsv("data/KEGG_DB/KEGG_Module_Hierarchy.tsv"), Level_D, Module) %>% distinct()) %>%
  mutate(Level_D = str_replace_all(string = Level_D, pattern = ", .*", replacement = "")) %>%
  mutate(Type = "Sulfate")

aromatics_degradation_modules <- read_tsv("data/KEGG_DB/Mapping_Module_KO.tsv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  filter(Module %in%  c("M00538",       "M00537",       "M00419",       "M00547",       "M00548",       "M00551",       "M00568",       "M00569",       "M00539",       "M00543",       "M00544",
                        "M00418",       "M00541",       "M00540",       "M00534",       "M00624",       "M00623",       "M00636",       "M00545")) %>%
  left_join(., select(read_tsv("data/KEGG_DB/KEGG_Module_Hierarchy.tsv"), Level_D, Module) %>% distinct()) %>%
  mutate(Level_D = str_replace_all(string = Level_D, pattern = ", .*", replacement = "")) %>%
  mutate(Type = "Aromatics")

alkane_degradation_modules <- read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv") %>%
  filter(str_detect(string = .$KO_D, pattern = "NAO; nitroalkane oxidase") |
           str_detect(string = .$KO_D, pattern = "dhaA") |
           str_detect(string = .$KO_D, pattern = "ssuD") |
           str_detect(string = .$KO_D, pattern = "alkB1_2, alkM; alkane 1-monooxygenase") |
           str_detect(string = .$KO_D, pattern = "ladA") |
           str_detect(string = .$KO_D, pattern = "pyruvate formate lyase activating enzyme")) %>%
  mutate(Type = "Alkane")

hydrogen_modules <- read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  filter(KO %in% c("K00437", "K18008", "K17997", "K17998", "K17999")) %>%
  mutate(KO_D = str_replace_all(KO_D, pattern = ";.*", replacement = "")) %>%
  mutate(Type = "Hydrogen")

lipid_degradation_modules <- read_tsv("data/KEGG_DB/Mapping_Module_KO.tsv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  filter(Module %in%  c("M00086", "M00087")) %>%
  left_join(., select(read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv"), KO, KO_D) %>% distinct()) %>%
  mutate(Type = "Lipids")

methanogenesis_table <-  read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv") %>%
  filter(KO %in% c("K00399", "K00400", "K00401", "K00402", "K03421", "K03422")) %>%
  mutate(KO_D = str_replace_all(KO_D, pattern = ";.*", replacement = "")) %>%
  mutate(KO_D = ifelse(KO_D == "K00400", "mcrA2", KO_D)) %>%
  mutate(Type = "Methanogenesis")

wood_modules <- read_tsv("data/KEGG_DB/KEGG_Module_KO_list.csv") %>%
  filter(KO %in% {enframe(table(.$KO)) %>% filter(value < 2) %>% .$name}) %>%
  filter(Module == "M00377")  %>%
  left_join(., select(read_tsv("data/KEGG_DB/Mapping_KO_Description.tsv"), KO, KO_D) %>% distinct()) %>%
  mutate(Type = "Wood")

gene_selection_list <- list(general_modules, sulfate_reduction_modules, aromatics_degradation_modules,
                            alkane_degradation_modules, hydrogen_modules, methanogenesis_table,
                            wood_modules, lipid_degradation_modules)

# Transform the datalist-table of all functional categories into flat datatable format
datatable_list <- purrr::map(gene_selection_list, function(gene_table) {
  
  datatable <- datalist %>%
    filter_taxa_datalist(KEGG %in% gene_table$KO) %>%
    mutate_count_datalist(function(x) x/sum(x)) %>%
    mutate_meta_datalist(Temp_Grp = as.numeric(Temp_Grp)) %>%
    summarize_by_param(Temp_Grp) %>%
    filter_station_datalist(!is.na(Temp_Grp)) %>%
    mutate_meta_datalist(Type = gene_table$Type[1]) %>%
    create_datatable(grpBy = Class, upper_grp = Phylum, otherThreshold = 0.01, addOthers = T, addColorScheme = F) %>%
    mutate(Group = as.character(Group))
  
  return(datatable)
  
})

# And extra for CAZymes since they are encoded in their own column (CAZY_Family)
datalist_CAZY <- datalist
datalist_CAZY$Count_Data <- datalist_CAZY$Count_Data %>%
  mutate(CAZy_Family = str_replace_all(CAZy_Family, pattern = "^.*.fasta\\|", replacement = ""))

# Also transform CAZyme datalist into flat datatable
datatable_cazy <- datalist_CAZY %>%
  filter_taxa_datalist(CAZy_Family %in% c("AA", "CBM", "CE", "GH", "PL")) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  mutate_meta_datalist(Temp_Grp = as.numeric(Temp_Grp)) %>%
  summarize_by_param(Temp_Grp) %>%
  filter_station_datalist(!is.na(Temp_Grp)) %>%
  create_datatable(grpBy = Class, upper_grp = Phylum, otherThreshold = 0.01, addOthers = T, addColorScheme = F) %>%
  mutate(Group = as.character(Group)) %>%
  mutate(Type = "CAZyme")

# combine all datatables for all functional categories and group low-abundance
# taxonomies into "others"
new_datatable <- bind_rows(datatable_list)  %>%
  bind_rows(datatable_cazy) %>%
  create_others(otherThreshold = 0.1, Type, Temp_Grp) %>%
  filter(!is.na(Temp_Grp)) %>%
  mutate(Temp_Grp = ifelse(Temp_Grp == 1, "0-20°C",
                           ifelse(Temp_Grp == 2, "20-45°C", ">45°C"))) %>%
  mutate(Temp_Grp = ordered(Temp_Grp, levels = c("0-20°C", "20-45°C", ">45°C"))) %>%
  mutate(Group = ifelse(Group == "Others", "Other Phylum", Group)) %>%
  mutate(Group = str_replace_all(Group, pattern = "candidate", replacement = "Candidate")) %>%
  mutate(Group = ordered(Group, levels = c(sort(unique(Group[Group != "Other Phylum"])), "Other Phylum")))

# Create new color-bar for the complete taxonomy
new_colors <- add_colorscheme(new_datatable, upper_grp = "Phylum") %>%
  mutate(Group = paste0(Phylum, " - ", lower_grp)) %>%
  mutate(Group = ifelse(Phylum == "Other Phylum", "Other Phylum", Group)) %>%
  mutate(Group = ifelse(Phylum == "Others", "Other Phylum", Group))

# Define a function for plotting stacked barplots with ggplot
plot_function <- function(datatable, title, colors) {
  p1 <- datatable %>%
    ggplot(aes(x = as.factor(Temp_Grp), y = Abundance*100, fill = Group)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors) +
      theme_bw() +
      labs(x = "Temperature (°C)", y = "Rel. Abundance (%)", fill = "Phylum - Class", title = title) +
      theme(strip.text.x = element_text(size = 7, angle = 0),
            strip.text.y = element_text(size = 7, angle = 0)) +
      guides(fill = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.7))
  return(p1)
}

# Plot all the different functional categories with the previously defined
# plotting-function
p1 <- new_datatable %>%
  filter(Type == "General") %>%
  plot_function(title = "General Categories", filter(new_colors, Group %in% .$Group)$Color)

p2 <- new_datatable %>%
  filter(Type == "Wood") %>%
  plot_function(title = "Wood Ljungdahl", filter(new_colors, Group %in% .$Group)$Color)

p3 <- new_datatable %>%
  filter(Type == "Alkane") %>%
  plot_function(title = "Alkane Degradation", filter(new_colors, Group %in% .$Group)$Color)

p4 <- new_datatable %>%
  filter(Type == "CAZyme") %>%
  plot_function(title = "CAZyme Families", filter(new_colors, Group %in% .$Group)$Color)

p5 <- new_datatable %>%
  filter(Type == "Methanogenesis") %>%
  plot_function(title = "Methanogenesis", filter(new_colors, Group %in% .$Group)$Color)

p6 <- new_datatable %>%
  filter(Type == "Sulfate") %>%
  plot_function(title = "Sulfate Metabolism", filter(new_colors, Group %in% .$Group)$Color)

p7 <- new_datatable %>%
  filter(Type == "Hydrogen") %>%
  plot_function(title = "Hydrogen Metabolism", filter(new_colors, Group %in% .$Group)$Color)

p8 <- new_datatable %>%
  filter(Type == "Lipids") %>%
  plot_function(title = "Lipid Degradation", filter(new_colors, Group %in% .$Group)$Color)

p9 <- new_datatable %>%
  filter(Type == "Aromatics") %>%
  plot_function(title = "Aromatics Degradation", filter(new_colors, Group %in% .$Group)$Color)

# Combine all individual plots into a single plot
cowplot::plot_grid(p1 + theme(legend.position = "none"),
                   cowplot::plot_grid(p2 + theme(legend.position = "none"), p3  + theme(legend.position = "none"), p4  + theme(legend.position = "none"), p5  + theme(legend.position = "none"),
                                      p6  + theme(legend.position = "none"), p7  + theme(legend.position = "none"), p8  + theme(legend.position = "none"), p9  + theme(legend.position = "none"),  ncol = 4),
                   rel_widths = c(0.25, 0.75)
)

ggsave("fig/Taxonomy_Composite.png", width = 13, height = 8, dpi = 300)

# And save color-legend extra
new_datatable %>% plot_function(title = "Legend", new_colors$Color) %>%
  cowplot::get_legend() %>%
  ggsave(plot = ., filename = "fig/Taxonomy_Legend.png", width = 5, height = 12, dpi = 300)

