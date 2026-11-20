#### This script independently clusters DOM compounds and functional genes ####
#### based on individual temperature patterns and compares resulting cluster-distribution ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)
library(furrr)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/unfold_kegg.R")
source("R/Create_Datatable_mod.R")
source("R/loess_clustering.R")

cbbPalette <- c("grey", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                "#CC79A7", "darkred", "forestgreen")

# Import gene datalist and remove functions with unknown metabolic roles 
# i.e. second hierarchical level of KEGG orthology classified as "unclassified"
datalist_genes <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                       Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv")) %>%
  filter_taxa_datalist(Kingdom %in% c("Bacteria", "Archaea")) %>%
  unfold_kegg("../MG_Pipeline/ressources/KO_Orthology.tsv") %>%
  filter_taxa_datalist(!str_detect(string = KO_B, pattern = "^Unclassified")) %>%
  filter_station_datalist(Type == "Metagenome") %>%
  mutate_count_datalist(function(x) x/sum(x))

# Collapse functional genes based on their functional description 
# (fourth level of KEGG orthology)
datalist_genes$Count_Data <- datalist_genes$Count_Data %>%
  select(-c(2:12, 17)) %>%
  group_by(KO_A, KO_B, KO_C, KO_D) %>%
  summarize_if(is.numeric, sum) %>%
  ungroup()

# Transform gene-count matrix into flat datatable
datatable_genes <- datalist_genes %>%
  create_datatable(grpBy = KO_D, upper_grp = KO_D, otherThreshold = 0, addColorScheme = F, addOthers = F) %>%
  filter(Group != "Others")

# Import DOM datalist
datalist_DOM <- list(Count_Data = read_tsv("data/DOM/Final_Count_Data.tsv"),
                     Meta_Data = read_tsv("data/DOM/Final_Meta_Data.tsv")) %>%
  filter_station_datalist(Temperature_C <= max(datalist_genes$Meta_Data$Temperature))

# And transform intensity-matrix into flat datatable
datatable_DOM <- datalist_DOM %>%
  create_datatable(grpBy = OTU_ID, upper_grp = OTU_ID, otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  filter(Group != "Others")

#### Clustering of genes and compounds ####
# Adjust maximum memory size per furrr-multisession 
options(future.globals.maxSize = 1.5 * 1024^3)  

# Compute easy loess-models to smooth data for clustering, visualize them
# and save cluster-table. The script can be run with sil_analysis = T to 
# identify optimal number of clusters (then k_cluster is overridden by optimal number)
cluster_tab_temp_genes <- loess_cluster(datatable_genes, x_parameter = Temperature, k_clusters = 3, loess_span = 0.75,
                                        sil_analysis = T, threads = 8)

ggsave("fig/Cluster_Genes.png", width = 7, height = 4, dpi = 300)

write_csv(cluster_tab_temp_genes$Cluster_Table, "output/Gene_Cluster_hclust.csv")

# And the same clustering for DOM data
cluster_tab_temp_DOM <- loess_cluster(datatable_DOM, x_parameter = Temperature_C, k_clusters = 4, loess_span = 0.75,
                                      sil_analysis = T, threads = 8)

ggsave("fig/Cluster_DOM.png", width = 7, height = 4, dpi = 300)

write_csv(cluster_tab_temp_DOM$Cluster_Table, "output/DOM_Cluster_hclust.csv")

# Combine visualization of DOM and functional gene clusters based on
# clusters identified above.
temperature_vec_combined <- seq(max(range(datatable_DOM$Temperature_C)[1], range(datatable_genes$Temperature)[1]), 
                                min(range(datatable_DOM$Temperature_C)[2], range(datatable_genes$Temperature)[2]),
                                length.out = 100)

# Compute loess-prediction of abundances along fixed temperature grid for
# genes and DOM compounds
loess_smooth_genes_comb <- furrr::future_map(unique(datatable_genes$Group), function(KO_D_tmp) {
  tmp_data <- datatable_genes %>%
    filter(Group == KO_D_tmp) %>%
    mutate(Abundance = (Abundance-mean(Abundance))/sd(Abundance))
  
  tmp_model <- loess(Abundance ~ Temperature, data = tmp_data)
  
  tmp_prediction <- predict(tmp_model, data.frame(Temperature = temperature_vec_combined))
  
  return(tibble(KO_D = KO_D_tmp, Abundance_pred = tmp_prediction, Temperature = temperature_vec_combined))
}) %>%
  bind_rows()

loess_smooth_DOM_comb <- furrr::future_map(unique(datatable_DOM$Group), function(Compound_tmp) {
  tmp_data <- datatable_DOM %>%
    filter(Group == Compound_tmp) %>%
    mutate(Abundance = (Abundance-mean(Abundance))/sd(Abundance))
  
  tmp_model <- loess(Abundance ~ Temperature_C, data = tmp_data)
  
  tmp_prediction <- predict(tmp_model, data.frame(Temperature_C = temperature_vec_combined))
  
  return(tibble(Compound = Compound_tmp, Abundance_pred = tmp_prediction, Temperature = temperature_vec_combined))
}) %>%
  bind_rows()

# Combine cluster-associations of genes and DOM compounds into one table
cluster_table_comb <- bind_rows(
  cluster_tab_temp_genes$Cluster_Table %>%
    mutate(Case = "Genes"),
  cluster_tab_temp_DOM$Cluster_Table %>%
    mutate(Case = "DOM")
)

# Because cluster labels are arbitrary, we harmonize both clustering approaches by
# re-labeling cluster-names to make them comparable between DOM and functional genes.
# (so that DOM cluster 1 is similar to gene cluster 1)
cluster_table_comb <- cluster_table_comb %>%
  mutate(Cluster = ifelse(Cluster == 1 & Case == "DOM", 5, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 3 & Case == "DOM", 1, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 4 & Case == "DOM", 3, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 2 & Case == "DOM", 4, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 5 & Case == "DOM", 2, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 1 & Case == "Genes", 4, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 2 & Case == "Genes", 1, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 3 & Case == "Genes", 2, Cluster)) %>%
  mutate(Cluster = ifelse(Cluster == 4 & Case == "Genes", 3, Cluster)) %>%
  mutate(Cluster = paste0("Cluster: ", Cluster))

# Combine loess-predictions with cluster assignments
loess_smooth_combined_comb <- bind_rows(
  loess_smooth_genes_comb %>%
    dplyr::rename("Group" = "KO_D") %>%
    mutate(Group = as.character(Group)) %>%
    mutate(Case = "Genes"),
  loess_smooth_DOM_comb %>%
    dplyr::rename("Group" = "Compound") %>%
    mutate(Group = as.character(Group)) %>%
    mutate(Case = "DOM")
) %>%
  left_join(., cluster_table_comb, by = c("Group", "Case")) %>%
  mutate(Case = ordered(Case, levels = c("DOM", "Genes")))

# Compute for each cluster a mean temperature-abundance-pattern for visualization
mean_smooth_combined_comb <- loess_smooth_combined_comb %>%
  group_by(Temperature, Cluster, Case) %>%
  summarize(Abundance_pred_mean = mean(Abundance_pred)) %>%
  mutate(Case = ordered(Case, levels = c("DOM", "Genes")))

# Visualize everything together in one plot
ggplot() +
  geom_line(aes(x = Temperature, y = Abundance_pred, group = Group), col = "darkgrey", alpha = .05,
            data = loess_smooth_combined_comb) +
  geom_line(aes(x = Temperature, y = Abundance_pred_mean, group = 1), col = "darkred", lwd = 1.5,
            data = mean_smooth_combined_comb) +
  facet_grid(Case~Cluster, scales = "free_y") +
  labs(x = "Temperature (°C)", y = "Abundance (z-scored)") +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave("fig/Gene_DOM_Cluster_Comparison.png", width = 7, height = 6, dpi = 300)

#### COmposition per sample ####
# This part of the script was for internal control and not used for the publication.
# We check the abundance of each inferred cluster in the total sample.

datalist_genes_cluster <- list(Count_Data = datalist_genes$Count_Data %>%
    left_join(., dplyr::rename(cluster_table_comb, "KO_D" = "Group") %>% select(-Case) %>%
                mutate(KO_D = str_replace_all(KO_D, pattern = " - .*$", replacement = ""))) %>%
    group_by(Cluster) %>%
    summarize_if(is.numeric, sum),
  Meta_Data = datalist_genes$Meta_Data)

datalist_DOM_cluster <- list(Count_Data = datalist_DOM$Count_Data %>%
                               left_join(., dplyr::rename(cluster_table_comb, "OTU_ID" = "Group") %>% select(-Case) %>%
                                           mutate(OTU_ID = str_replace_all(OTU_ID, pattern = " - .*$", replacement = ""))) %>%
                               group_by(Cluster) %>%
                               summarize_if(is.numeric, sum),
                             Meta_Data = datalist_DOM$Meta_Data)

temp_interval <- seq(0, 65, 5)

datatable_gene_cluster <- datalist_genes_cluster %>%
  create_datatable(grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0, addOthers = F) %>%
  filter(Type == "Metagenome") %>%
  group_by(Temperature, Core, Group) %>%
  summarize(Abundance = mean(Abundance)) %>%
  filter(Group != "Others") %>%
  mutate(Case = "Genes") %>%
  mutate(Group = as.character(Group))

datatable_DOM_cluster <- datalist_DOM_cluster %>%
  create_datatable(grpBy = Cluster, upper_grp = Cluster, otherThreshold = 0, addOthers = F) %>%
  dplyr::rename("Temperature" = "Temperature_C") %>%
  group_by(Temperature, Core, Group) %>%
  summarize(Abundance = mean(Abundance)) %>%
  filter(Group != "Others") %>%
  mutate(Case = "DOM") %>%
  mutate(Group = as.character(Group))

bind_rows(datatable_gene_cluster, datatable_DOM_cluster) %>%
  mutate(Temperature_int = temp_interval[findInterval(Temperature, temp_interval)]) %>%
  group_by(Case, Temperature_int, Group) %>%
  summarize(Abundance = mean(Abundance)) %>%
  mutate(Group = str_replace_all(Group, pattern = " - .*$", replacement = "")) %>%
  mutate(Case = ordered(Case, levels = c("DOM", "Genes"))) %>%
  ggplot(aes(x = as.factor(Temperature_int), y = Abundance*100, fill = Group)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = cbbPalette[-1]) +
    facet_grid(Case~.) +
    labs(x = "Temperature group (°C)", y = "Abundance (%)", fill = "") +
    theme_bw() +
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave("fig/Cluster_Composition_Comparison.png", width = 6, height = 6, dpi = 300)
