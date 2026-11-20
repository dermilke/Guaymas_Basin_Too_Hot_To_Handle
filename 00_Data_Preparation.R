#### This script prepares raw DOM data for further analyses ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

DOM_taxonomy_maker <- function(taxonomy) {
  
  taxonomy %>%
    mutate(Type = ifelse(Aromatic == 1, "Aromatic", "Unclassified")) %>%
    mutate(Type = ifelse(Aromatic_O_rich == 1, "Aromatic_O_rich", Type)) %>%
    mutate(Type = ifelse(Aromatic_O_poor == 1, "Aromatic_O_poor", Type)) %>%
    mutate(Type = ifelse(Highly_unsaturated == 1, "Highly_unsaturated", Type)) %>%
    mutate(Type = ifelse(Unsaturated == 1, "Unsaturated", Type)) %>%
    mutate(Type = ifelse(Unsaturated_O_rich == 1, "Unsaturated_O_rich", Type)) %>%
    mutate(Type = ifelse(Unsaturated_O_poor == 1, "Unsaturated_O_poor", Type)) %>%
    mutate(Type = ifelse(Unsaturated_with_N == 1, "Unsaturated_with_N", Type)) %>%
    mutate(Type = ifelse(Highly_unsaturated_O_rich == 1 & Unsaturated_with_N == 1, "Highly_unsaturated_O_rich_with_N", Type)) %>%
    mutate(Type = ifelse(Highly_unsaturated_O_rich == 1 & Unsaturated_with_N == 0, "Highly_unsaturated_O_rich", Type)) %>%
    mutate(Type = ifelse(Highly_unsaturated_O_poor == 1 & Unsaturated_with_N == 1, "Highly_unsaturated_O_poor_with_N", Type)) %>%
    mutate(Type = ifelse(Highly_unsaturated_O_poor == 1 & Unsaturated_with_N == 0, "Highly_unsaturated_O_poor", Type)) %>%
    mutate(Type = ifelse(Saturated == 1, "Saturated", Type)) %>%
    mutate(Type = ifelse(Highly_unsaturated_O_poor == 1, "Highly_unsaturated_O_poor", Type)) %>%
    mutate(Type = ifelse(Saturated_O_rich == 1, "Saturated_O_rich", Type)) %>%
    mutate(Type = ifelse(Saturated_O_poor == 1, "Saturated_O_poor", Type)) %>%
    dplyr::rename("DOM_Classification" = "Type") %>%
    return()
}

Meta_Data <- read_csv("data/DOM/Meta_Data_raw_More.csv")

DOM_data <- read_csv("data/DOM/DOM_Data_raw.csv") %>%
  mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>% 
  select_if(names(.) %in% c("formula_isotopefree", Meta_Data$Sample_ID)) %>%
  dplyr::rename("OTU_ID" = "formula_isotopefree") %>%
  group_by(OTU_ID) %>%
  summarize_if(is.numeric, mean) %>%
  mutate_if(is.numeric, function(x) x/sum(x))

DOM_datalist <- list(Count_Data = DOM_data,
                     Meta_Data = Meta_Data %>%
                       dplyr::slice(match(names(DOM_data), Sample_ID)))

DOM_datalist %>%
  mutate_count_datalist(function(x) round(x*1000000)) %>%
  singleton_filter(min_count = 10, min_station = 15) %>%
  .$Count_Data %>%
  dplyr::rename("#OTU ID" = "OTU_ID") %>%
  write_tsv("output/FastSpar/FastSpar_Dataprep_DOM.tsv")

taxonomy <- read_csv("data/DOM/DOM_Data_raw.csv") %>%
  mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>% 
  dplyr::rename("OTU_ID" = "formula_isotopefree") %>%
  select(1:27) %>%
  group_by(OTU_ID) %>%
  summarize_if(is.numeric, mean) %>%
  ungroup() %>%
  DOM_taxonomy_maker() %>%
  select(-c("Aromatic", "Aromatic_O_rich", "Aromatic_O_poor", "Highly_unsaturated", 
            "Highly_unsaturated_O_rich", "Highly_unsaturated_O_poor", "Unsaturated",
            "Unsaturated_O_rich", "Unsaturated_O_poor", "Unsaturated_with_N", "Saturated",
            "Saturated_O_rich", "Saturated_O_poor")) %>%
  select(-mz, -mean_signal_to_MDL, -DBE)

DOM_datalist <- list(Count_Data = bind_cols(select(taxonomy, OTU_ID, DOM_Classification), select(DOM_data, -OTU_ID)),
                     Meta_Data = DOM_datalist$Meta_Data)

write_tsv(DOM_datalist$Count_Data, "data/DOM/Final_Count_Data.tsv")
write_tsv(DOM_datalist$Meta_Data, "data/DOM/Final_Meta_Data.tsv")
