#### This script contains the function to compute the optimal number of loess-cluster ####
## It is flexible towards different environmental gradients and can be run in 2
## modes: with or without shilouette width computation. Typically, you need to run
## this analysis with sil_analysis = T to identify the optimal number of clusters 
## (highest shilouette-width). Once you know it, you can omit this step and 
## pre-determine the cluster-number with k_clusters.
## Input parameter:
# datatable - Input datatable containing observational units in Group column
#             and Abundances of these units in the Abundance column and the
#             the environmental gradient along which the analysis is executed
# x_parameter - Name of the environmental parameter used to compute loess-cluster
# k_clusters - pre-determined number of clusters (use only once you know the optimal number)
# sil_analysis - if TRUE: computes shilouette-width to determine optimal number of clusters
# loess_span - the span used for loess-model computation. Higher values leads to stronger smoothing
# threads - Number of CPUs for multisessions
#
# Author: Felix Milke
# Date: 19.11.2026

loess_cluster <- function(datatable, x_parameter, k_clusters, sil_analysis = F, loess_span = 0.5, threads = 4) {
  
  require(rlang)
  
  x_parameter <- enquo(x_parameter)
  
  plan(multisession, workers = threads)
  
  # Define environmental space
  param_vec <- seq(range(pull(datatable, !!x_parameter))[1], range(pull(datatable, !!x_parameter))[2], length.out = 100)
  
  # Compute loess-models along environmental axis
  loess_smooth <- furrr::future_map(unique(datatable$Group), function(Group_tmp) {
    tmp_data <- datatable %>%
      filter(Group == Group_tmp) %>%
      mutate(Abundance = (Abundance-mean(Abundance))/sd(Abundance))
    
    tmp_model <- loess(as.formula(paste("Abundance ~", as_name(x_parameter))), data = tmp_data,
                       span = loess_span)
    
    tmp_prediction <- predict(tmp_model, tibble(!!as_name(x_parameter) := param_vec))
    
    return(tibble(Group = Group_tmp, Abundance_pred = tmp_prediction, !!as_name(x_parameter) := param_vec))
  }) %>%
    bind_rows()
  
  # Compute Distance-matrix of predicted abundance-matrix for clustering
  bc_mat <- loess_smooth %>%
    reshape2::dcast(formula = as.formula(paste0("Group ~ ", as_name(x_parameter))), value.var = "Abundance_pred") %>% 
    magrittr::set_rownames(.[,1]) %>%
    select(-1) %>%
    vegan::vegdist(method = "euclidean")
  
  # If sil_analysis = T: Compute shilouette-width for 2-10 clusters
  if (sil_analysis) {
    sil_width <- purrr::map_dbl(seq(2, 10), function(cluster_num) {
      hclust(bc_mat) %>%
        cutree(k = cluster_num) %>%
        cluster::silhouette(., bc_mat) %>%
        .[,3] %>% mean()
    })
    plot(seq(2, 10), sil_width)
    k_clusters <- as.numeric(readline(prompt = "Please enter number of clusters: "))
  }
  
  okay = FALSE
  
  while (!okay) {
    
    # Cluster dissimilarity-matrix based on hierarchical clustering and optimal number of clusters
    cluster_nums <- hclust(bc_mat, method = "complete") %>%
      cutree(k = k_clusters)
    
    cluster_table <- tibble(Group = names(cluster_nums),
                            Cluster = cluster_nums)
    
    # Gather data for plotting
    loess_smooth_mod <- loess_smooth %>%
      left_join(., cluster_table, by = "Group")
    
    mean_smooth <- loess_smooth_mod %>%
      group_by(!!x_parameter, Cluster) %>%
      summarize(Abundance_pred_mean = mean(Abundance_pred))
    
    p1 <- ggplot() +
      geom_line(aes(x = !!x_parameter, y = Abundance_pred, group = Group), col = "darkgrey", alpha = .05,
                data = loess_smooth_mod) +
      geom_line(aes(x = !!x_parameter, y = Abundance_pred_mean, group = 1), col = "darkred", lwd = 1.5,
                data = mean_smooth) +
      facet_grid(~Cluster, scales = "free_y") +
      labs(x = as_name(x_parameter), y = "Abundance (z-scored)")
    
    print(p1)
    
    k_clusters <- as.numeric(readline(prompt = "Enter new number of clusters (or leave empty if okay): "))
    
    okay <- ifelse(is.na(k_clusters), TRUE, FALSE)
    
  }
  
  return(list(Cluster_Table = cluster_table, p1))
  
}