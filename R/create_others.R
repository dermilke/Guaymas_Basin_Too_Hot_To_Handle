create_others <- function(datatable, otherThreshold = 0.0005, ...) {
  
  additional_Vars = enquos(...)
  
  datatable_tmp <- datatable %>%
    mutate(upper_groups = str_split_fixed(Group, pattern = " - ", n = 2)[,1]) %>%
    mutate(lower_groups = str_split_fixed(Group, pattern = " - ", n = 2)[,2]) 
  
  upper_groups_keep <- datatable_tmp %>%
    group_by(upper_groups) %>%
    summarize(Abundance = sum(Abundance)) %>%
    filter(Abundance > otherThreshold) %>%
    .$upper_groups
  
  groups_keep <- datatable_tmp %>%
    filter(upper_groups %in% upper_groups_keep) %>%
    group_by(Group) %>%
    summarize(Abundance = sum(Abundance)) %>%
    filter(Abundance >= otherThreshold/10) %>%
    .$Group
  
  max_abundance_table <- datatable_tmp %>%
    group_by(Sample_ID, !!!additional_Vars) %>%
    summarize(Total_Abundance = sum(Abundance)) %>%
    ungroup()
  
  filtered_abundance_table <- datatable_tmp %>%
    filter(upper_groups %in% upper_groups_keep) %>% 
    group_by(Sample_ID, !!!additional_Vars) %>% 
    summarize(Filtered_Abundance = sum(Abundance)) %>%
    ungroup()
  
  abundance_others_upper <- left_join(max_abundance_table, filtered_abundance_table) %>% 
    mutate(Filtered_Abundance = ifelse(is.na(Filtered_Abundance), 0, Filtered_Abundance)) %>%
    mutate(Abundance = Total_Abundance - Filtered_Abundance) %>%
    mutate(Group = "Others") %>%
    select(-Total_Abundance, -Filtered_Abundance)
  
  max_abundance_table_lower <- datatable_tmp %>%
    filter(upper_groups %in% upper_groups_keep) %>% 
    group_by(Sample_ID, Group, !!!additional_Vars) %>%
    summarize(Total_Abundance = sum(Abundance)) %>%
    ungroup()
  
  filtered_abundance_table_lower <- datatable_tmp %>%
    filter(upper_groups %in% upper_groups_keep) %>% 
    filter(Group %in% groups_keep) %>% 
    group_by(Sample_ID, Group, !!!additional_Vars) %>% 
    summarize(Filtered_Abundance = sum(Abundance)) %>%
    ungroup()
  
  abundance_others_lower <- left_join(max_abundance_table_lower, filtered_abundance_table_lower) %>% 
    mutate(Filtered_Abundance = ifelse(is.na(Filtered_Abundance), 0, Filtered_Abundance)) %>%
    mutate(Abundance = Total_Abundance - Filtered_Abundance) %>% 
    mutate(upper_group = str_split_fixed(Group, pattern = " - ", n = 2)[,1]) %>%
    group_by(Sample_ID, !!!additional_Vars, upper_group) %>%
    summarize(Abundance = sum(Abundance)) %>%
    mutate(Group = paste0(upper_group, " - Others")) %>%
    select(-upper_group)

  bind_rows(filter(datatable, Group %in% groups_keep) %>% select(Sample_ID, !!!additional_Vars, Abundance, Group),
            abundance_others_upper, abundance_others_lower) %>%
    filter(Abundance > 0) %>%
    left_join(., select(datatable, -Abundance, -Group) %>% distinct()) %>%
    return()

}
