import_sparcc_network <- function(cor_file, pval_file, min_r, min_p) {
  
  cor <- read.delim(cor_file, header = T,
                    check.names = F, sep = "\t") %>%
    magrittr::set_rownames(.[,1]) %>%
    .[,-1] %>%
    as.matrix()
  
  pval <- read.delim(pval_file, header = T, check.names = F, sep = "\t") %>%
    magrittr::set_rownames(.[,1]) %>%
    .[,-1] %>%
    as.matrix() %>%
    reshape2::melt(as.is = T) %>%
    dplyr::rename("From" = "Var1", "To" = "Var2", "pvalue" = "value")
  
  edges <- cor %>%
    reshape2::melt(as.is = T) %>%
    dplyr::rename("From" = "Var1", "To" = "Var2", "weight" = "value") %>%
    as_tibble() %>%
    left_join(., pval, by = c("From", "To")) %>%
    filter(abs(weight) >= min_r & pvalue <= min_p) %>%
    filter(weight != 0) %>%
    mutate(Type = ifelse(weight < 0, "Negative", "Positive"))
  
  return(edges)
  
}

create_network <- function(datalist_input, cor_mat, r_threshold = 0.5, type = "Positive",
                           plot_location = NULL) {
  
  Max_Count <- datalist_input %>%
    mutate_count_datalist(function(x) x/sum(x)) %>%
    .$Count_Data %>%
    select_if(is.numeric) %>%
    mutate(Max = apply(., 1, which.max)) %>%
    mutate(OTU_ID = datalist_input$Count_Data$OTU_ID) %>%
    select(OTU_ID, Max) %>%
    cbind(., datalist_input$Meta_Data[.$Max,]) %>%
    as_tibble()
  
  network_data <- map(cor_mat, function(x) {
    import_sparcc_network(cor_file = x, 
                          pval_file = str_replace_all(x, pattern = "cor", replacement = "pvalues"),
                          min_r = r_threshold, min_p = 0.05)
  }) %>%
    bind_rows() %>%
    filter(Type == type) %>%
    group_by(From, To) %>%
    summarize(weight = max(weight)) 
  
  network <- network_data %>%
    filter(weight != 1) %>%
    #select(-weight) %>%
    igraph::graph_from_data_frame(d = ., directed = F,
                                  vertices = slice(Max_Count, match(unique(c(pull(., From), 
                                                                             pull(., To))),
                                                                    OTU_ID)) %>%
                                    mutate_if(is.factor, as.character) %>%
                                    select_if(function(x) is.character(x) | is.numeric(x))) %>%
    igraph::set_vertex_attr(graph = ., name = "label", value = NA) %>%
    igraph::simplify(edge.attr.comb = "first")
  
  deg_tmp <- igraph::degree(network, mode = "all")
  
  colors_cluster <- read_csv("data/Colors_Cluster.csv")$Colour
  
  colors_cluster_ramp <- colorRampPalette(colors = colors_cluster)
  
  cluster_obj <- cluster_fast_greedy(network)
  
  cluster_mem <- tibble(OTU_ID = cluster_obj$names, Cluster = cluster_obj$membership) %>%
    group_by(Cluster) %>%
    mutate(N = n()) %>%
    ungroup() %>%
    filter(N >= 6) %>%
    mutate(Cluster = ordered(Cluster, levels = unique(Cluster)) %>% as.numeric()) %>%
    mutate(Color = with(., if (length(unique(.$Cluster)) > 16) colors_cluster_ramp(length(unique(.$Cluster)))[Cluster] 
                        else colors_cluster[Cluster])) 
  
  vertex_colors <- tibble(OTU_ID = names(V(network))) %>%
    left_join(., cluster_mem, by = "OTU_ID") %>%
    .$Color
  
  network <- igraph::set.vertex.attribute(network, name = "Plot_Size", value =  ifelse((log(deg_tmp)) < 3 , 2.5, (log(deg_tmp)*1.8)))
  
  network <- igraph::delete_vertices(network, V(network)[!(V(network)$name %in% cluster_mem$OTU_ID)])
  network <- igraph::set.vertex.attribute(network, name = "Cluster", value =  cluster_mem$Cluster)
  network <- igraph::set.vertex.attribute(network, name = "Color", value =  cluster_mem$Color)
  
  layout_network <- igraph::layout_nicely(network)
  
  cluster_table <- tibble(OTU_ID = igraph::get.vertex.attribute(network, name = "name"),
                          Cluster = igraph::get.vertex.attribute(network, name = "Cluster"),
                          Color = igraph::get.vertex.attribute(network, name = "Color")) %>%
    filter(!is.na(Color)) %>%
    mutate(Cluster = ordered(Cluster, levels = unique(Cluster), labels = seq(1, length(unique(Cluster))))) 
  
  return(list(network_data = network_data, network = network, 
              cluster_obj = cluster_obj, cluster_table = cluster_table))
  
}
