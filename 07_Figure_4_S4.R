#### This script Imports FastSpar Co-occurrences of DOM compounds from the complete dataset ####
#### And analyses clusters within this co-occurrence network ####
## Relies on output of script 01_FastSpar_DOM.R
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)
library(igraph)
library(cluster)
library(gghalves)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Cluster_Pipeline.R")
source("R/Import_SparCC_Network.R")

add_compound_info <- function(datalist, atom, .position = 2) {
  
  datalist$Count_Data <- datalist$Count_Data %>%
    mutate(!!atom := {str_extract(string = OTU_ID, pattern =  paste0(atom, "[0-9]*")) %>%
        str_replace_all(string = ., pattern = paste0("^", atom, "$"), replacement = "1") %>%
        str_replace_all(string = ., pattern = atom, replacement = "") %>%
        unlist() %>% as.numeric() %>% ifelse(is.na(.), 0, .)}, .before = .position) 
  
  return(datalist)
  
}

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                "#CC79A7", "darkred", "forestgreen")

# Import DOM data
DOM_datalist <- list(Count_Data = read_tsv("data/DOM/Final_Count_Data.tsv"),
                     Meta_Data = read_tsv("data/DOM/Final_Meta_Data.tsv"))

# Transform intensity-matrix into flat datatable
DOM_datatable <- DOM_datalist %>%
  create_datatable(grpBy = DOM_Classification, upper_grp = DOM_Classification, otherThreshold = 0, addColorScheme = F, addOthers = F)

# Infer Co-occurrence Clusters based on consensus-clustering pipeline from 
# Milke et al. (2026). For that, we first import the correlation-matrix and
# filter out correlation with p-value < 0.05. Next, we transform the matrix
# into flat table format and filter out co-occurrences with correlation values
# below min_r.
network_data_complete <- import_sparcc_network("output/FastSpar/fastspar_cor.tsv",
                                               "output/FastSpar/fastspar_pvalues.tsv",
                                               min_r = 0.5, min_p = 0.05) %>%
  filter(Type == "Positive")

# Then we use this flat co-occurrence table as input for the network-clustering
# pipeline. We determined the optimal cluster-resolution interval using the
# find_parameters = T mode, which yielded res_param = c(0.6, 1.3). 
# Cluster_num is overwritten by the optimal cluster-number that is automatically
# determined using the "consensus_network" approach. Disconnected networks with
# fewer than 5 nodes are discarded and clusters with less than 10 compounds as well.
cluster_data <- network_data_complete %>%
  complete_cluster_pipeline(., find_parameters = F,
                            vertex_num_drop = 5, res_param = c(0.6, 1.3), 
                            cluster_method = "Consensus_network", cluster_num = 6,
                            min_cluster_size = 10)

# Save compound-cluster-table
write_csv(cluster_data$Cluster_Tab, "output/FastSpar/Network_Cluster_DOM.csv")

# Since cluster-label are arbitrary, we change their names so that they 
# better fit their order along the temperature gradient (cluster 1 = low temperature,
# cluster 2 = intermediate temperature, cluster 3 = high temperature)
Cluster_Tab <- read_csv("output/FastSpar/Network_Cluster_DOM.csv") %>%
  mutate(Cluster = ifelse(Cluster == 1, 4, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 3, 1, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 4, 3, Cluster))

# Visualize the final co-occurrence network
cluster_colors <- cbbPalette[c( 9, 1, 5, 6, 2, 3)]
node_size <- ifelse((log(degree(cluster_data$Network))) < 3 , 1.5, (log(degree(cluster_data$Network))*0.8))
node_size[node_size > 2] <- 2
node_color <- cluster_colors[c(3, 2, 1)][vertex_attr(cluster_data$Network, "Cluster")]
layout_tmp <- layout_with_fr(cluster_data$Network, grid = "nogrid")

png(filename = "fig/03_DOM_Network.png", width = 10, height = 10, res = 600, units = "in")
par(mfrow = c(1,1))
plot(cluster_data$Network, vertex.label = NA, vertex.color = node_color, 
     vertex.frame.color = NA, vertex.size = 1,
     layout = layout_tmp, edge.width = 0.3)
dev.off()

## Visualize cluster composition in each sample as barplot
# Collapse DOM_datalist compound intensity-matrix into cluster-intensities
datalist_DOM_cluster <- DOM_datalist

datalist_DOM_cluster$Count_Data <- right_join(dplyr::rename(Cluster_Tab, "OTU_ID" = "ASV"), DOM_datalist$Count_Data, by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  mutate(Cluster = paste0("Cluster ", Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum)

# Transform cluster-intensity-matrix into flat datatable format
datatable_DOM_cluster <- datalist_DOM_cluster %>%
  create_datatable(grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0, addOthers = F)

# Visualize cluster composition in barplots
datatable_DOM_cluster %>%
  filter(Group != "Others") %>%
  group_by(Temperature_C, Group, Core) %>%
  summarize_if(is.numeric, mean) %>%
  mutate(Temperature_C = round(Temperature_C, digits = 1)) %>%
  ggplot(aes(x = as.factor(Temperature_C), y = Abundance * 100, fill = Group)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(.~Core, scales = "free_x", space = "free") +
    scale_fill_manual(values = cbbPalette[c( 9, 1, 5, 6, 2, 3)]) +
    scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 3)]) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    labs(y = "Proportion (%)", x = "Temperature (°C)", fill = "")

ggsave("fig/03_DOM_CoOccurrence_Cluster_Composition2.png", width = 10, height = 4, dpi = 300)

## Characterize DOM co-occurrence cluster based on Van-Krevelen diagrams

# Compute atomic composition of each DOM compound
DOM_datalist_atomic <- DOM_datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  add_compound_info(atom = "S") %>%
  add_compound_info(atom = "C") %>%
  add_compound_info(atom = "H") %>%
  add_compound_info(atom = "O") %>%
  add_compound_info(atom = "N") %>%
  add_compound_info(atom = "P") 

# Visualize number of atoms for each compound in H/C - O/C space and analyse 
# compounds according to cluster-association
DOM_datalist_atomic$Count_Data %>%
  left_join(dplyr::rename(Cluster_Tab, "OTU_ID" = "ASV"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  mutate(Cluster = paste0("Cluster ", Cluster)) %>%
  select(Cluster, P, N, O, H, C, S, DOM_Classification) %>%
  mutate(DOM_Classification = str_replace_all(DOM_Classification, pattern = "_", replacement = " ")) %>%
  mutate(H_C = H/C, O_C = O/C) %>%
  filter(DOM_Classification != "Unclassified") %>%
  ggplot(aes(x = O_C, y = H_C, col = DOM_Classification)) +
    geom_point(size = 0.3) +
    scale_color_manual(values = c(cbbPalette, "salmon")[-7],
                       guide = guide_legend(override.aes = list(size = 1.5))) +
    facet_grid(.~Cluster) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    labs(x = "O/C", y = "H/C", col = "DOM Classification")

ggsave("fig/03_DOM_CoOccurrence_Cluster_VanKrevelen.png", width = 8, height = 3, dpi = 300)

## Analyse sample-wise ratios of cluster intensities. This metric is robust to compositional
## data (log-ratios). We correlate this ratio to temperature to see how temperature 
## affects DOM cluster composition.

# Calculate sample-wise log-ratios of cluster_2 : cluster_1
ratio_data <- datatable_DOM_cluster %>%
  group_by(Sample_ID) %>%
  mutate(Ratio = log(Abundance[Group == "Cluster 2 - Cluster 2"] / Abundance[Group == "Cluster 1 - Cluster 1"])) %>%
  ungroup()

# Estimate linear-model performance
summary(lm(ratio_data$Ratio ~ ratio_data$Temperature_C))
  
# Visualize linear dependency
ggplot(ratio_data, aes(y = Ratio, x = Temperature_C)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, col = "darkred") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    labs(x = "Temperature (°C)", y = "Intensity ratio\nlog( DOM Cluster 2 : DOM Cluster 1 )")

ggsave("fig/03_DOM_CoOccurrence_Cluster_Ratio.png", width = 3.7, height = 3.7, dpi = 300)  


## Compute Elemental composition of each cluster within each sample. Use intensity-weighted
## average number of atoms per sample for that. Finally, visualize the resulting
## average number of atoms per sample and cluster in a boxplot.
elemental_cluster_composition <- DOM_datalist_atomic$Count_Data %>%
  left_join(dplyr::rename(Cluster_Tab, "OTU_ID" = "ASV"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  mutate(Cluster = paste0("Cluster ", Cluster)) %>%
  pivot_longer(starts_with("MGC"),
               names_to = "Sample", values_to = "Intensity") %>%
  filter(!is.na(Intensity), Intensity > 0) %>%
  group_by(Sample, Cluster) %>%
  summarise(
    `C Atoms` = weighted.mean(C, Intensity),
    `H Atoms` = weighted.mean(H, Intensity),
    `O Atoms` = weighted.mean(O, Intensity),
    `P Atoms` = weighted.mean(P, Intensity),
    `N Atoms` = weighted.mean(N, Intensity),
    `S Atoms` = weighted.mean(S, Intensity),
    .groups = "drop"
  ) %>%
  mutate(`H/C` = `H Atoms` / `C Atoms`) %>%
  mutate(`O/C` = `O Atoms` / `C Atoms`) %>%
  select(-Sample) %>%
  reshape2::melt() %>%
  magrittr::set_colnames(c("Cluster", "Atom", "Abundance")) %>%
  mutate(Cluster = ordered(Cluster, levels = c("Cluster 1", "Cluster 2", "Cluster 3"))) 

elemental_cluster_composition %>%
  ggplot(aes(x = Cluster, y = Abundance)) +
    geom_boxplot(width = 0.25, 
                 position = position_nudge(x = -0.15)) +
    geom_jitter(aes(x = as.numeric(Cluster) + 0.15, colour = Cluster),
                width = 0.07, height = 0, alpha = 0.6) +
    scale_color_manual(values = cbbPalette[c( 9, 1, 5, 6, 2, 3)]) +
    facet_wrap(~Atom, scales = "free_y", nrow = 2) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = "", y = "Average number of atoms\nweighted by compound intensity")

ggsave("fig/03_DOM_CoOccurrence_Cluster_Elements.png", width = 9, height = 5, dpi = 300)

elemental_cluster_composition %>%
  filter(Atom == "S Atoms") %>%
  ggplot(aes(x = Cluster, y = Abundance)) +
  geom_boxplot(width = 0.25, 
               position = position_nudge(x = -0.15)) +
  geom_jitter(aes(x = as.numeric(Cluster) + 0.15, colour = Cluster),
              width = 0.07, height = 0, alpha = 0.6) +
  scale_color_manual(values = cbbPalette[c( 9, 1, 5, 6, 2, 3)]) +
  facet_wrap(~Atom, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(x = "", y = "Average number of atoms\nweighted by compound intensity")

ggsave("fig/03_DOM_CoOccurrence_Cluster_Elements_S_atoms.png", width = 3.5, height = 3.6, dpi = 300)
