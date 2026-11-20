create_datatable <- function(datalist, grpBy = Family, upper_grp = Class, otherThreshold = 0.01, selectGrps = NULL, 
         addColorScheme = F, addOthers = T) {
  
  make_proportion <- function(x) {
    result <- x/sum(x, na.rm = T)
    return(ifelse(is.finite(result), result, 0))
  }
  
  grpBy <- enquo(grpBy)
  upper_grp <- enquo(upper_grp)
  
  datalist$Count_Data <- mutate(datalist$Count_Data, !!quo_name(grpBy) := paste(!!upper_grp, !!grpBy, sep = " - "))
  
  if (is.null(selectGrps)) {
    
    if (addOthers) {
      
      tmp_Upper <- datalist$Count_Data %>%
        group_by(!!upper_grp) %>%
        select_if(is.numeric) %>%
        summarize_all(sum, na.rm = T)
      
      tmp_Lower <- datalist$Count_Data %>%
        group_by(!!grpBy) %>%
        select_if(is.numeric) %>%
        summarize_all(sum, na.rm = T) 
      
      if (nrow(tmp_Lower) > 1) 
        tmp_Lower <- filter(tmp_Lower, apply(select(tmp_Lower, -1), 2, make_proportion) %>% apply(., 1, mean, na.rm = T) > otherThreshold) 
      
      tmp_both <- tmp_Lower %>%
        mutate(!!quo_name(upper_grp) := pull(left_join(tmp_Lower, distinct(select(datalist$Count_Data, !!upper_grp, !!grpBy))), !!upper_grp), .before = 1) 
      
      tmp_Upper_red <- tmp_both %>%
        group_by(!!upper_grp) %>%
        select_if(is.numeric) %>%
        summarize_all(sum, na.rm = T)
      
      tmp <- (select_if(filter(tmp_Upper, !!upper_grp %in% pull(tmp_Upper_red, !!upper_grp)), is.numeric) - 
                select_if(tmp_Upper_red, is.numeric)) %>%
        as_tibble() %>%
        mutate(!!quo_name(grpBy) := paste0(pull(tmp_Upper_red, !!upper_grp), " - Others"), .before = 1) %>%
        rbind(., tmp_Lower) %>%
        mutate_at(1, ~replace(., is.na(.), "")) %>%
        mutate_at(2:ncol(.), as.numeric) %>%
        select_if(~sum(is.na(.)) == 0) %>%
        add_case({colSums(select_if(datalist$Count_Data, is.numeric)) -
            summarise(., across(where(is.numeric), sum))} %>%
              mutate(!!quo_name(grpBy) := paste0("Other ", as_label(upper_grp)))) %>%
        filter(rowSums(select_if(., is.numeric)) > 10^(-16))
      
    } else {
      
      tmp <- datalist$Count_Data %>%
        group_by(!!grpBy) %>%
        select_if(is.numeric) %>%
        summarize_all(sum, na.rm = T) %>%
        filter(with(., {if (nrow(.) > 1) {
          select(., -1) %>%
            apply(2, make_proportion) %>% 
            apply(1, mean, na.rm = T)} else rowMeans(.[,-1])}) > otherThreshold) %>%
        rbind(., c("Others", apply(select_if(datalist$Count_Data, is.numeric), 2, sum) - apply(select(., -1), 2, sum))) %>%
        mutate_at(1, ~replace(., is.na(.), "")) %>%
        mutate_at(2:ncol(.), as.numeric) %>%
        select_if(~sum(is.na(.)) == 0)
      
    }
    
    
  } else {
    
    tmp <- datalist$Count_Data %>%
      group_by(!! grpBy) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T) %>%
      filter(!! grpBy %in% selectGrps) %>%
      rbind(., c(paste0("Others (<", otherThreshold * 100,"%)"), apply(select_if(datalist$Count_Data, is.numeric), 2, sum) - apply(select(., -1), 2, sum))) %>%
      mutate_at(1, ~replace(., is.na(.), "")) %>%
      mutate_at(2:ncol(.), as.numeric) %>%
      select_if(~sum(is.na(.)) == 0)
    
  }
  
  if (addOthers) {
    
    group_label <- NA
    
    for (i in 1:length(pull(tmp_Upper_red, !!upper_grp))) {
      
      group_label <- c(group_label, pull(filter(tmp_both, !!upper_grp == pull(tmp_Upper_red, !!upper_grp)[i]), !!grpBy), 
                       paste0(pull(tmp_Upper_red, !!upper_grp)[i], " - Others"))
      
    }
    
    group_label <- c(group_label[-1], paste0("Other ", as_label(upper_grp)))
    
    datatable <- tmp %>%
      reshape2::melt(.) %>%
      dplyr::rename("Sample_ID" = "variable", "Abundance" = "value") %>%
      left_join(., distinct(datalist$Meta_Data), by = "Sample_ID") %>%
      as_tibble() %>%
      dplyr::rename("Group" = as_label(grpBy)) %>%
      mutate(Group = ordered(Group, levels = group_label)) %>%
      arrange(Group)
    
  } else {
    
    tmp_richness <- datalist$Count_Data %>%
      mutate_if(is.numeric, function(x) as.numeric(x > 0)) %>%
      group_by(!!grpBy) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T) %>%
      filter(!!grpBy %in% deframe(select(tmp, as_label(grpBy)))) %>%
      rbind(., c("Others", apply(select_if(datalist$Count_Data, is.numeric) > 0, 2, sum) - apply(select(., -1), 2, sum))) %>%
      mutate_at(1, ~replace(., is.na(.), "")) %>%
      mutate_at(2:ncol(.), as.numeric) %>%
      select_if(~sum(is.na(.)) == 0) 
    
    datatable <- tmp %>%
      reshape2::melt(.) %>%
      dplyr::rename("Sample_ID" = "variable", "Abundance" = "value") %>%
      left_join(., {tmp_richness %>% reshape2::melt(.) %>% dplyr::rename("Sample_ID" = "variable", "Richness" = "value")}, by = c("Sample_ID", as_label(grpBy))) %>%
      left_join(., distinct(datalist$Meta_Data), by = "Sample_ID") %>%
      as_tibble() %>%
      dplyr::rename("Group" = as_label(grpBy)) %>%
      mutate(Group = ordered(Group, levels = c(unique(Group)[-which(Group == "Others")], "Others"))) %>%
      arrange(Group)
    
  }
  
  if (addColorScheme) {
    
    #colorTaxonomy <- read_csv("https://raw.githubusercontent.com/dermilke/ExCom/master/Data/Colors/Taxonomy_Colour.csv")
    colorTaxonomy <- read_csv("~/PhD/Projects/ExCom/Data/Colors/Taxonomy_Colour.csv")
    
    derivative_color = str_split_fixed(datatable$Group, " - ", 2) %>%
      unique() %>%
      .[,1] %>%
      table() %>%
      as_tibble() %>%
      dplyr::rename("Group" = ".") %>%
      mutate(Group = ordered(Group, levels = c(unique(Group[Group != paste0("Other ", as_label(upper_grp))]), 
                                               paste0("Other ", as_label(upper_grp))))) %>%
      arrange(Group) %>%
      mutate(main_color = colorTaxonomy$Colour[match(Group, colorTaxonomy$Group)] %>%
               ifelse(is.na(.), "grey", .)) %>%
      mutate(rowNumber = 1:n()) %>%
      with(., {
        derivative_color <- NULL
        for (i in rowNumber) {
          
          derivative_color <- c(derivative_color,
                                colorspace::lighten(main_color[i], 
                                                    amount = {if (n[i] == 1) 0
                                                      else if (n[i] > 6) (seq_len(n[i])/n[i])*1.6-0.8
                                                      else (seq_len(n[i])/n[i])-0.5}
                                )
          )
          
        }
        derivative_color
      })
    
    result <- list(table = datatable, color = derivative_color)
    
  } else {
    
    result <- datatable
    
  }
  
  return(result)
  
}
