add_colorscheme = function(datatable, upper_grp) {

  colorTaxonomy <- read_csv("~/PhD/Projects/ExCom/Data/Colors/Taxonomy_Colour.csv")

  taxa_table <- str_split_fixed(datatable$Group, " - ", 2) %>%
    unique() %>%
    as_tibble() %>%
    magrittr::set_colnames(c(upper_grp, "Lower_Grp")) %>%
    arrange(!!sym(upper_grp))
  
  main_color_table <- taxa_table[,1] %>%
    table() %>%
    as_tibble() %>%
    dplyr::rename("Group" = upper_grp) %>%
    mutate(Group = ordered(Group, levels = c(unique(Group[Group != paste0("Other ", upper_grp)]), 
                                             paste0("Other ", upper_grp)))) %>%
    arrange(Group) %>%
    mutate(main_color = colorTaxonomy$Colour[match(Group, colorTaxonomy$Group)] %>%
             ifelse(is.na(.), "grey", .)) %>%
    mutate(rowNumber = 1:n())
    
  derivative_color <- NULL
      
  for (i in seq_along(row_number(main_color_table))) {
        
    derivative_color <- c(derivative_color,
                          colorspace::lighten(main_color_table$main_color[i], 
                                              amount = {if (main_color_table$n[i] == 1) 0
                                                    else if (main_color_table$n[i] > 6) (seq_len(main_color_table$n[i])/main_color_table$n[i])*1.6-0.8
                                                    else (seq_len(main_color_table$n[i])/main_color_table$n[i])-0.5}))
        
  }
  
  tibble(!!upper_grp := rep(main_color_table$Group, main_color_table$n),
         Color = derivative_color) %>%
    with(., {if (paste0("Other ", upper_grp) %in% taxa_table$Phylum) {
      mutate(., lower_grp = taxa_table %>% filter(!!sym(upper_grp) != paste0("Other ", upper_grp)) %>% 
             arrange(!!sym(upper_grp)) %>% add_case(!!sym(upper_grp) := paste0("Other ", upper_grp)) %>% .$Lower_Grp) %>%
        slice(-which(!!sym(upper_grp) == paste0("Other ", upper_grp))) %>%
        add_case(!!sym(upper_grp) := paste0("Other ", upper_grp), Color = "#BEBEBE") 
      } else mutate(., lower_grp = arrange(taxa_table, !!sym(upper_grp))$Lower_Grp)}) %>%
    return()

}
