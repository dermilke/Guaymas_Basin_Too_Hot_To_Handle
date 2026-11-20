#### This script visualizes selected environmental variables ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)
library(ggtext)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

# This function calculates the number of atoms per sum-formula and stores them in 
# individual columns of the count-matrix
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

# Import DOM data to infer sample-wise DOM metrics
datalist_DOM <- list(Count_Data = read_tsv("data/DOM/Final_Count_Data.tsv"),
                     Meta_Data = read_tsv("data/DOM/Final_Meta_Data.tsv"))

# Compute atom numbers for each compound
datalist_DOM_tmp <- datalist_DOM %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  add_compound_info(atom = "S") %>%
  add_compound_info(atom = "C") %>%
  add_compound_info(atom = "H") %>%
  add_compound_info(atom = "O") %>%
  add_compound_info(atom = "N") %>%
  add_compound_info(atom = "P") 

# Calculate compound-intensity weighted number of atoms for each sample
# for S, C, H, O, N, compound-length, O/C, H/C
datatable_DOM <- datalist_DOM$Meta_Data %>%
  mutate(S_avg = {datalist_DOM_tmp$Count_Data %>%
      select(-c(1:8)) %>%
      mutate_if(is.numeric, function(x) x*datalist_DOM_tmp$Count_Data$S) %>%
      colSums()}) %>%
  mutate(C_avg = {datalist_DOM_tmp$Count_Data %>%
      select(-c(1:8)) %>%
      mutate_if(is.numeric, function(x) x*datalist_DOM_tmp$Count_Data$C) %>%
      colSums()}) %>%
  mutate(H_avg = {datalist_DOM_tmp$Count_Data %>%
      select(-c(1:8)) %>%
      mutate_if(is.numeric, function(x) x*datalist_DOM_tmp$Count_Data$H) %>%
      colSums()}) %>%
  mutate(O_avg = {datalist_DOM_tmp$Count_Data %>%
      select(-c(1:8)) %>%
      mutate_if(is.numeric, function(x) x*datalist_DOM_tmp$Count_Data$O) %>%
      colSums()}) %>%
  mutate(N_avg = {datalist_DOM_tmp$Count_Data %>%
      select(-c(1:8)) %>%
      mutate_if(is.numeric, function(x) x*datalist_DOM_tmp$Count_Data$N) %>%
      colSums()}) %>%
  mutate(length_avg = {datalist_DOM_tmp$Count_Data %>%
      select(-c(1:8)) %>%
      mutate_if(is.numeric, function(x) x*(datalist_DOM_tmp$Count_Data$S +
                                             datalist_DOM_tmp$Count_Data$C +
                                             datalist_DOM_tmp$Count_Data$H +
                                             datalist_DOM_tmp$Count_Data$O +
                                             datalist_DOM_tmp$Count_Data$N)) %>%
      colSums()}) %>%
  mutate(o_c_avg = {datalist_DOM_tmp$Count_Data %>%
      select(-c(1:8)) %>%
      mutate_if(is.numeric, function(x) x*datalist_DOM_tmp$Count_Data$O / datalist_DOM_tmp$Count_Data$C) %>%
      colSums()})  %>%
  mutate(h_c_avg = {datalist_DOM_tmp$Count_Data %>%
      select(-c(1:8)) %>%
      mutate_if(is.numeric, function(x) x*datalist_DOM_tmp$Count_Data$H / datalist_DOM_tmp$Count_Data$C) %>%
      colSums()})

### Abiotic parameter ###
# Select abiotic parameter for visualization
abiotic_parameter <- c("Sample_ID", "Core", "Temperature_C", "BC", "TDN_mM", "TN_weight",
                       "Alkalinity", "TS_weight", "SO42", "H2S", "NH4", "Fe")

# Visualize those abiotic parameter and format parameter names together with measurement units
datatable_DOM %>%
  select(abiotic_parameter) %>%
  reshape2::melt(id.vars = c("Sample_ID", "Core", "Temperature_C")) %>%
  magrittr::set_colnames(c("Sample_ID", "Core", "Temperature_C", "Parameter", "Value")) %>%
  as_tibble() %>%
  filter(Core %in% c("U1545", "U1547")) %>%
  mutate(Parameter = str_replace_all(Parameter, pattern = "BC", replacement = "Dissolved black<br>carbon (µM)") %>%
           str_replace_all(pattern = "TDN_mM", "Total dissolved<br>nitrogen (mM)") %>%
           str_replace_all(pattern = "Fe", "Fe (µM)") %>%
           str_replace_all(pattern = "H2S", "H<sub>2</sub>S (µM)") %>%
           str_replace_all(pattern = "TN_weight", "Total<br>nitrogen (wt %)") %>%
           str_replace_all(pattern = "TS_weight", "Total<br>sulfur (wt %)") %>%
           str_replace_all(pattern = "SO42", "SO<sub>4</sub><sup>2−</sup> (mM)") %>%
           str_replace_all(pattern = "NH4", "NH<sub>4</sub> (mM)")) %>%
  
  ggplot(aes(x = Temperature_C, y = Value, col = Parameter)) +
    geom_point() +
    geom_smooth(se = F) +
    geom_vline(xintercept = c(20, 45, 60), linetype = 2) +
    facet_grid(Parameter~Core, scales = "free_y") +
    scale_color_manual(values = cbbPalette) +
    labs(x = "Temperature (°C)", y = "") +
    theme_bw() +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 12, angle = 0),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text.y = element_markdown(size = 10, angle = 0))

ggsave("fig/Abiotic_Parameter.png", width = 5, height = 8, dpi = 300)

### Biotic parameter ###

# Select the subset of relevant parameter for visualization, format their names
# and visualize with ggplot
datatable_DOM %>%
  select(Sample_ID, Temperature_C, Core, C, H, O, HC, OC, N, P, S, CRAM, 
         Polycyclic_Aromatics_all, Highly_unsaturated_wm, Unsaturated_wm,
         Saturated_wm, mz) %>%
  reshape2::melt(id.vars = c("Sample_ID", "Core", "Temperature_C")) %>%
  magrittr::set_colnames(c("Sample_ID", "Core", "Temperature_C", "Parameter", "Value")) %>%
  as_tibble() %>%
  filter(Core %in% c("U1545", "U1547")) %>%
  mutate(Parameter = str_replace_all(Parameter, pattern = "Polycyclic_Aromatics_all", replacement = "Polycyclic\nAromatics") %>%
           str_replace_all(pattern = "Highly_unsaturated_wm", "Highly\nUnsaturated") %>%
           str_replace_all(pattern = "Unsaturated_wm", "Unsaturated") %>%
           str_replace_all(pattern = "Saturated_wm", "Saturated") %>%
           str_replace_all(pattern = "mz", "m/z (Da)") %>%
           str_replace_all(pattern = "HC", "H/C") %>%
           str_replace_all(pattern = "DOC_mM", "DOC (mM)") %>%
           str_replace_all(pattern = "^OC$", "O/C")) %>%
  mutate(Parameter = ordered(Parameter, levels = c("m/z (Da)", "H/C", "O/C", "C", "H", "O", 
                                                   "N", "S", "P", "CRAM", "Polycyclic\nAromatics",
                                                   "Highly\nUnsaturated", "Unsaturated", "Saturated"))) %>%
  filter(!(Parameter %in% c("DOC (mM)", "H/C", "ideg", "NOSC"))) %>%
  
  ggplot(aes(x = Temperature_C, y = Value, col = Parameter)) +
    geom_vline(xintercept = c(20, 45, 60), linetype = 2) +
    geom_point() +
    geom_smooth(se = F) +
    facet_grid(Parameter~Core, scales = "free_y") +
    scale_color_manual(values = c(cbbPalette, "#66cdaa", "#ff6666", "#ffff66", "#bada55", "#0a75ad", "#008000")) +
    labs(x = "Temperature (°C)", y = "") +
    theme_bw() +
    theme(legend.position = "none",
          strip.text.x = element_text(size = 12, angle = 0),
          strip.text.y = element_text(size = 12, angle = 0),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))

ggsave("fig/Biotic_Parameter_supplement.png", width = 6, height = 11, dpi = 300)

#### Biotic parameter - Subset for Figure 2 ####
# Parameter subset for manuscript main figure
biotic_parameter_red <- c("Sample_ID", "Temperature_C", "Core", "DOC_mM", "HC", "NOSC", "ideg")

# Names-formating and visualization
datatable_DOM %>%
  select(biotic_parameter_red) %>%
  reshape2::melt(id.vars = c("Sample_ID", "Core", "Temperature_C")) %>%
  magrittr::set_colnames(c("Sample_ID", "Core", "Temperature_C", "Parameter", "Value")) %>%
  as_tibble() %>%
  filter(Core %in% c("U1545", "U1547")) %>%
  mutate(Parameter = str_replace_all(Parameter, pattern = "Polycyclic_Aromatics_all", replacement = "Polycyclic\nAromatics") %>%
           str_replace_all(pattern = "Highly_unsaturated_wm", "Highly\nUnsaturated") %>%
           str_replace_all(pattern = "Unsaturated_wm", "Unsaturated") %>%
           str_replace_all(pattern = "Saturated_wm", "Saturated") %>%
           str_replace_all(pattern = "mz", "m/z (Da)") %>%
           str_replace_all(pattern = "HC", "H/C") %>%
           str_replace_all(pattern = "DOC_mM", "DOC (mM)") %>%
           str_replace_all(pattern = "^OC$", "O/C")) %>%
  mutate(Parameter = ifelse(Parameter == "ideg", "I<sub>DEG</sub>", Parameter)) %>%
  mutate(Parameter = ordered(Parameter, levels = c("DOC (mM)", "m/z (Da)", "H/C", "O/C", "N", "S", "NOSC", "I<sub>DEG</sub>"))) %>%
  filter((Parameter %in% c("DOC (mM)", "H/C","NOSC", "I<sub>DEG</sub>"))) %>%
  
  ggplot(aes(x = Temperature_C, y = Value, col = Parameter)) +
    geom_vline(xintercept = c(20, 45, 60), linetype = 2) +
    geom_point() +
    geom_smooth(se = F) +
    facet_grid(Parameter~Core, scales = "free_y") +
    scale_color_manual(values = c(cbbPalette[-c(4:6)], "#66cdaa", "#ff6666", "#ffff66", "#bada55", "#0a75ad", "#008000")) +
    labs(x = "Temperature (°C)", y = "") +
    theme_bw() +
    theme(legend.position = "none",
          strip.text.y = element_markdown(size = 9))

ggsave("fig/Biotic_Parameter_reduced2.png", width = 5.5, height = 5, dpi = 300)
