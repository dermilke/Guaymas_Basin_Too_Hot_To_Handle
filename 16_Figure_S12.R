#### This script visualizes taxonomic composition of functional genes ####
#### and cell-abundances vs. temperatures per drilling core ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)
library(furrr)
library(vegan)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/Create_Datatable_mod.R")
source("R/datatable_color.R")
source("R/create_others.R")

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
                "#CC79A7", "darkred", "forestgreen")

add_compound_info <- function(datalist, atom, .position = 2) {
  
  datalist$Count_Data <- datalist$Count_Data %>%
    mutate(!!atom := {str_extract(string = OTU_ID, pattern =  paste0(atom, "[0-9]*")) %>%
        str_replace_all(string = ., pattern = paste0("^", atom, "$"), replacement = "1") %>%
        str_replace_all(string = ., pattern = atom, replacement = "") %>%
        unlist() %>% as.numeric() %>% ifelse(is.na(.), 0, .)}, .before = .position) 
  
  return(datalist)
  
}

# Read datalists for DOM and genes and add meta-data column for Temperature-Intervals 
# (according to cluster distribution)
datalist <- list(Count_Data = read_tsv("data/Gene_Library/Processed/Decontaminated_Gene_Count.tsv"),
                 Meta_Data = read_tsv("data/Gene_Library/Meta_Data/Meta_Data.tsv")) %>%
  filter_taxa_datalist(Kingdom %in% c("Archaea", "Bacteria")) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_station_datalist(Type == "Metagenome") %>%
  mutate_meta_datalist(Temp_Grp = ifelse(Temperature < 20, "0-20",
                           ifelse(Temperature > 20 & Temperature <= 45, "20-45",
                                  ifelse(Temperature > 45, ">45", "Control")))) %>%
  mutate_meta_datalist(Temp_Grp = ordered(Temp_Grp, levels = c("0-20", "20-45", ">45"),
                            labels = c("0-20°C", "20-45°C", ">45°C"))) 

datalist_DOM <- list(Count_Data = read_tsv("data/DOM/Final_Count_Data.tsv"),
                     Meta_Data = read_tsv("data/DOM/Final_Meta_Data.tsv")) %>%
  mutate_meta_datalist(Temp_Grp = ifelse(Temperature_C < 20, "0-20",
                                         ifelse(Temperature_C > 20 & Temperature_C <= 45, "20-45",
                                                ifelse(Temperature_C > 45, ">45", "Control")))) %>%
  mutate_meta_datalist(Temp_Grp = ordered(Temp_Grp, levels = c("0-20", "20-45", ">45"),
                                          labels = c("0-20°C", "20-45°C", ">45°C"))) 

# Collapse gene datalist by KEGG function
datalist_KEGG <- datalist %>%
  summarize_by_taxa(tax_lvl = KEGG)

# Collapse gene datalist by CAZy function
datalist_CAZY <- datalist %>%
  summarize_by_taxa(tax_lvl = CAZy_ID)

# Extract count table and transpose
KEGG_cca_count_data <- datalist_KEGG %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  t()

CAZY_cca_count_data <- datalist_CAZY %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  t()

# Extract Meta-Data and select only potential explanatory parameter
# Parameter pre-selection was done based on which parameter were measured from the same
# sampling material as the metagenomes were taken from
KEGG_cca_meta_data <- datalist_KEGG %>%
  .$Meta_Data %>%
  select(-c("Sample_ID", "Type", "Core", "Temp_Grp", "Depth"))

CAZY_cca_meta_data <- datalist_CAZY %>%
  .$Meta_Data %>%
  select(-c("Sample_ID", "Type", "Core", "Temp_Grp", "Depth"))

# Run CCA for both, minimal model (~1) and maximal model (~.)
KEGG_model_min <- vegan::cca(KEGG_cca_count_data ~ 1, KEGG_cca_meta_data)
KEGG_model_max <- vegan::cca(KEGG_cca_count_data ~ ., KEGG_cca_meta_data)

CAZY_model_min <- vegan::cca(CAZY_cca_count_data ~ 1, CAZY_cca_meta_data)
CAZY_model_max <- vegan::cca(CAZY_cca_count_data ~ ., CAZY_cca_meta_data)

# Run forward-selection ordistep, start from minimal model, scope is maximal model
KEGG_model_final <- vegan::ordistep(KEGG_model_min, scope = formula(KEGG_model_max), direction = c("both"))
CAZY_model_final <- vegan::ordistep(CAZY_model_min, scope = formula(CAZY_model_max), direction = c("both"))

cca_CAZY <- vegan::cca(datalist_CAZY %>%
                         filter_station_datalist(Type == "Metagenome") %>%
                         .$Count_Data %>%
                         select_if(is.numeric) %>%
                         t(),
                       datalist_CAZY %>%
                         filter_station_datalist(Type == "Metagenome") %>%
                         .$Meta_Data %>%
                         select(Alkalinity, SO4, PO4, H2S, NH4, CH4, CO, DOC, TOC, TN))

# Extract CCA scores from final model and merge with meta-data table
cca_tab_sites_KEGG <- KEGG_model_final$CA$u[,1:2] %>%
  as_tibble(rownames = "Sample_ID") %>%
  left_join(datalist_KEGG$Meta_Data)

cca_tab_sites_CAZY <- CAZY_model_final$CA$u[,1:2] %>%
  as_tibble(rownames = "Sample_ID") %>%
  left_join(datalist_CAZY$Meta_Data)

# Extract DOM data for CCA
DOM_for_cca <- datalist_DOM %>%
  filter_station_datalist(!is.na(DOC_mM)) %>%
  filter_station_datalist(!is.na(BC))

DOM_cca_count_data <- DOM_for_cca %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  t()

# And pre-select environmental parameter that were measured at the same location.
# We include here parameter that were derived from the DOM data itself. 
# We also highlight that in the manuscript.
DOM_cca_meta_data <- DOM_for_cca %>%
  .$Meta_Data %>%
  select("Temperature_C", "DOC_mM", "HC", "C", "N", "S", "P", "ideg", "CRAM", "NOSC", 
         "Unsaturated_wm", "Saturated_wm", "BC", "Alkalinity", "SO42", "NH4")

# Run CCA for both, minimal model (~1) and maximal model (~.)
DOM_model_min <- vegan::cca(DOM_cca_count_data ~ 1, DOM_cca_meta_data)
DOM_model_max <- vegan::cca(DOM_cca_count_data ~ ., DOM_cca_meta_data)

# Run forward-selection ordistep, start from minimal model, scope is maximal model
DOM_model_final <- vegan::ordistep(DOM_model_min, scope = formula(DOM_model_max), direction = c("both"))

# Compute significance of individual terms and axes
anova.cca(KEGG_model_final, by = "term")
anova.cca(KEGG_model_final, by = "axis")

anova.cca(CAZY_model_final, by = "term")
anova.cca(CAZY_model_final, by = "axis")

anova.cca(DOM_model_final, by = "term")
anova.cca(DOM_model_final, by = "axis")

# Define plot-size
plot_width <- 7
plot_height <- 5

# Plot CCA biplots for each of KEGG, CAZy, and DOM CCAs individually and save as PNG files
png(filename = "fig/CCA_KEGG.png", width = plot_width, height = plot_height, res = 300, units = "in")

plot(KEGG_model_final, sit.par = list(type = "p", cex = 2, pch = 21, bg = cbbPalette[cca_tab_sites_KEGG$Temp_Grp]),
     spe.par = list(cex = .6, col = adjustcolor("darkgrey", alpha.f = 0.4), pch = 10),
     bip.par = list(type = "text", cex = .7, col = "black", font = 2),
     xlab = paste0("CCA1 (", round(KEGG_model_final$CCA$eig[1] / KEGG_model_final$tot.chi * 100, 1), "%)"),
     ylab = paste0("CCA2 (", round(KEGG_model_final$CCA$eig[2] / KEGG_model_final$tot.chi * 100, 1), "%)"))

title(sub = paste0("Constrained = ",
                   round(KEGG_model_final$CCA$tot.chi / KEGG_model_final$tot.chi * 100, 1), "%;  ",
                   "Unconstrained = ",
                   round(KEGG_model_final$CA$tot.chi / KEGG_model_final$tot.chi * 100, 1), "%"))

legend("bottomleft",
       legend = levels(cca_tab_sites_KEGG$Temp_Grp),
       pt.bg = cbbPalette[seq_along(levels(cca_tab_sites_KEGG$Temp_Grp))],
       pch = 21,
       pt.cex = 2,
       bty = "n")

dev.off()

png(filename = "fig/CCA_CAZY.png", width = plot_width, height = plot_height, res = 300, units = "in")

plot(CAZY_model_final, sit.par = list(type = "p", cex = 2, pch = 21, bg = cbbPalette[cca_tab_sites_CAZY$Temp_Grp]),
     spe.par = list(cex = .6, col = adjustcolor("darkgrey", alpha.f = 0.4), pch = 10),
     bip.par = list(type = "text", cex = .7, col = "black", font = 2),
     xlab = paste0("CCA1 (", round(CAZY_model_final$CCA$eig[1] / CAZY_model_final$tot.chi * 100, 1), "%)"),
     ylab = paste0("CCA2 (", round(CAZY_model_final$CCA$eig[2] / CAZY_model_final$tot.chi * 100, 1), "%)"))


title(sub = paste0("Constrained = ",
                   round(CAZY_model_final$CCA$tot.chi / CAZY_model_final$tot.chi * 100, 1), "%;  ",
                   "Unconstrained = ",
                   round(CAZY_model_final$CA$tot.chi / CAZY_model_final$tot.chi * 100, 1), "%"))

legend("bottomleft",
       legend = levels(cca_tab_sites_CAZY$Temp_Grp),
       pt.bg = cbbPalette[seq_along(levels(cca_tab_sites_CAZY$Temp_Grp))],
       pch = 21,
       pt.cex = 2,
       bty = "n")

dev.off()

png(filename = "fig/CCA_DOM.png", width = plot_width, height = plot_height, res = 300, units = "in")

plot(DOM_model_final, sit.par = list(type = "p", cex = 2, pch = 21, bg = cbbPalette[DOM_for_cca$Meta_Data$Temp_Grp]),
     spe.par = list(cex = .6, col = adjustcolor("darkgrey", alpha.f = 0.4), pch = 10),
     bip.par = list(type = "text", cex = .7, col = "black", font = 2),
     xlab = paste0("CCA1 (", round(DOM_model_final$CCA$eig[1] / DOM_model_final$tot.chi * 100, 1), "%)"),
     ylab = paste0("CCA2 (", round(DOM_model_final$CCA$eig[2] / DOM_model_final$tot.chi * 100, 1), "%)"))

title(sub = paste0("Constrained = ",
                   round(DOM_model_final$CCA$tot.chi / DOM_model_final$tot.chi * 100, 1), "%;  ",
                   "Unconstrained = ",
                   round(DOM_model_final$CA$tot.chi / DOM_model_final$tot.chi * 100, 1), "%"))

legend("bottomleft",
       legend = levels(DOM_for_cca$Meta_Data$Temp_Grp),
       pt.bg = cbbPalette[seq_along(levels(DOM_for_cca$Meta_Data$Temp_Grp))],
       pch = 21,
       pt.cex = 2,
       bty = "n")

dev.off()