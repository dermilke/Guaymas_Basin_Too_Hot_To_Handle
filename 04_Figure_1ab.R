#### This script visualizes taxonomic composition of functional genes ####
#### and cell-abundances vs. temperatures per drilling core ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Create_Datatable_mod.R")

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                "#CC79A7", "darkred", "forestgreen")

# Read metagenomes in datalist format, discard control- and transcriptome-samples
# and transform counts into rel. abundance data (TSS).
datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv")) %>%
  filter_station_datalist(Type == "Metagenome") %>%
  filter_taxa_datalist(Kingdom %in% c("Bacteria", "Archaea")) %>%
  mutate_count_datalist(function(x) x/sum(x)) 

# Collapse gene data into taxonomic abundances
datalist$Count_Data <- datalist$Count_Data %>%
  mutate_if(is.character, function(x) str_replace_all(x, pattern = "Unknown .*$", replacement = "Unknown")) %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>%
  summarize_if(is.numeric, sum) %>%
  ungroup()
  
# Transform count-matrix into flat datatable
datatable <- datalist %>%
  create_datatable(grpBy = Class, upper_grp = Phylum, otherThreshold = 0.005, addColorScheme = T, addOthers = T)

# Visualize taxonomic composition
datatable$table %>%
  mutate(Group = as.character(Group)) %>%
  mutate(Group = ifelse(Group == "Unknown - Unknown", "Unclassified", Group)) %>%
  mutate(Group = ifelse(Group == "candidate division Zixibacteria - Unknown", "Candidate division Zixibacteria - Unknown", Group)) %>% 
  mutate(Group = ordered(Group, levels = c(sort(unique(Group)[!(unique(Group) %in% c("Other Phylum", "Unclassified"))]), "Other Phylum", "Unclassified"))) %>%
  ggplot(aes(x = as.factor(Temperature), y = Abundance*100, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", col = "black", lwd = 0.1) +
  facet_grid(~Core, scales = "free_x", space = "free") +
  scale_fill_manual(values = c(datatable$color[1:(length(datatable$color)-2)], "steelblue", "darkgrey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 5)) +
  labs(x = "Temperature (°C)", y = "Rel. abundance (%)", fill = "Phylum - Class") +
  guides(fill = guide_legend(ncol = 1,
                             keywidth = 0.8, keyheight = 0.7))

ggsave("fig/01_Taxonomic_Composition_Class.png", width = 10.3, height = 7.5, dpi = 300)

### Infer and visualize cell-counts per drilling core ###
# Problem: We only have cell-counts and depth but lack temperature.
# Because temperature increases linearly with depth, we can infer it from depth.

# First infer linear models of Depth - Temperature
temperature_data <- read_csv("data/Depth_Temperature.csv") %>%
  mutate(Site = str_replace_all(Site, pattern = "^385-", replacement = "") %>%
           str_replace_all(pattern = "-", replacement = "")) %>%
  mutate(Site = str_replace_all(Site, pattern = "[A,B]$", replacement = "A,B")) %>%
  filter(Depth <= 350)

lm_list <- purrr::map(unique(temperature_data$Site), function(site_tmp) {
  temperature_data %>%
    filter(Site == site_tmp) %>%
    with(., lm(Temperature ~ Depth, data = .))
}) %>%
  magrittr::set_names(unique(temperature_data$Site))

# Then predict temperature based on depth of cell-counts samples
cell_abundance <- read_csv("data/Cell_Abundance.csv") %>%
  mutate(Site = str_replace_all(Sample_ID, pattern = "_.*$", replacement = "")) %>%
  mutate(Site = str_replace_all(Site, pattern = "[A,B]$", replacement = "A,B"))

cell_abundance_plotting_data <- purrr::map(unique(cell_abundance$Site), function(site_tmp) {
  cell_abundance %>%
    filter(Site == site_tmp) %>%
    mutate(Temperature = predict(lm_list[[site_tmp]], newdata = data.frame(Depth = Depth)))
}) %>%
  bind_rows() %>%
  mutate(Site = str_replace_all(Site, pattern = "A,B", replacement = ""))

# Visualize cell-counts per temperature
axis_label <- expression(10^4, 10^6, 10^8)
axis_title <- expression(Cell~abundance~(cells~cm^{-3}))

cell_abundance_plotting_data %>%
  ggplot(aes(x = Temperature, y = Cell_Abundance, col = Site)) +
  geom_point(size = 2) +
  scale_y_log10(breaks = c(10^4, 10^6, 10^8),
                labels = axis_label) +
  geom_smooth(method = "lm", se = F, lty = "dashed") +
  scale_color_manual(values = cbbPalette[-c(2, 6)]) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = axis_title,
       col = "Core")

ggsave("fig/Cell_Abundance.png", width = 4, height = 5, dpi = 300)

