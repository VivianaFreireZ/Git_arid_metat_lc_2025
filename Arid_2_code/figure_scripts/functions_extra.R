# -------------------------------------------------------------------------
plot_dwt_clusters <- function(dwt_res, nrow = NULL, ncol = NULL){

  clusts <- tibble(FeatureID = names(dwt_res@cluster),
                   cluster = paste0('Cluster ', dwt_res@cluster))
  
  feature_dwt_clusts <- as.data.frame(dwt_res@datalist) %>% 
    rowid_to_column(var = 'time') %>% 
    mutate(time = case_when(time == 1 ~ 'D0',
                            time == 2 ~ 'D3',
                            time == 3 ~ 'D7',
                            time == 4 ~ 'D14')) %>% 
    mutate(time = factor(time, levels = c('D0', 'D3', 'D7', 'D14'))) %>% 
    pivot_longer(!time, names_to = 'FeatureID', values_to = 'value') %>% 
    left_join(clusts, by = 'FeatureID')
  
  rename_vec <- set_names(1:length(unique(dwt_res@cluster)),
                          paste0('Cluster ', unique(dwt_res@cluster)))
  
  dwt_centroids <- as.data.frame(dwt_res@centroids) %>% 
    rename(rename_vec) %>% 
    rowid_to_column(var = 'time') %>% 
    pivot_longer(!time, names_to = 'cluster', values_to = 'value')
  
  feature_dwt_plot <- ggplot(feature_dwt_clusts) +
    geom_line(aes(x = time,
                  y = value,
                  group = FeatureID),
              size = 0.5) +
    geom_line(data = dwt_centroids,
              aes(y = value,
                  x = time,
                  group = 1),
              linetype = 'dashed',
              color = 'red',
              size = 1) +
    facet_wrap(~cluster, nrow = nrow, ncol = ncol) +
    theme_bw() 
  
  return(feature_dwt_plot)
}