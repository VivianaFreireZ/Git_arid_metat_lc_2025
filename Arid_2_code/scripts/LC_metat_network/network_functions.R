# Get kegg network plot of specified compounds
create_kegg_network <- function(compound_list,
                                ko_list,
                                ko_abundances = NULL,
                                cid_abundances = NULL,
                                pathways,
                                prevalence_per_reaction = .5,
                                expand_search = 1,
                                generate_plot = TRUE){
  
  if(!is.null(ko_abundances)){
    stopifnot(
      'ko_abundances needs columns ko_id and getmm'='ko_id' %in% colnames(ko_abundances)
    )
    stopifnot(
      'all KO ids in ko_abundances must be in ko_list'=all(ko_abundances$ko_id %in% ko_list)
    )
  }
  
  # Loading kegg tables
  suppressMessages({
    ko_rx <- read_csv("scripts/LC_metat_network/output/KO_to_reaction_2025.csv") %>% 
      mutate(ko_id = str_remove(ko_id, 'ko:'),
             reaction = str_remove(reaction, 'rn:'))
    
    rx_compound <- read_csv("scripts/LC_metat_network/output/Reaction_to_compound_2025.csv") %>% 
      mutate(compound_id = str_remove(compound_id, 'cpd:'),
             reaction = str_remove(reaction, 'rn:'))
    
    path_ko <- read_csv('scripts/LC_metat_network/output/Pathway_to_ko_2025.csv') %>% 
      mutate(pathway = str_remove(pathway, 'path:'),
             ko_id = str_remove(ko_id, 'ko:'))
    
    path_rx <- read_csv("scripts/LC_metat_network/output/Pathway_to_reaction_2025.csv") %>% 
      mutate(pathway = str_remove(pathway, 'path:'),
             reaction = str_remove(reaction, 'rn:'))
    
    path_compound <- read_csv("scripts/LC_metat_network/output/Pathway_to_compound_2025.csv") %>% 
      mutate(pathway = str_remove(pathway, 'path:'),
             compound_id = str_remove(compound_id, 'cpd:'))
  }) 
  
  compounds_to_omit <- c("C00828",  # Menaquinone
                         "C00534",  # Pyridoxamine
                         "C00006",  # NADP+
                         "C00003",  # NAD+
                         "C00002",  # ATP
                         "C00314",  # Pyridoxine
                         "C00864",  # Pantothenate
                         "C00504",  # Folate
                         "C00032",  # Heme
                         "C05443",  # Vitamin D3
                         "C00253",  # Nicotinate
                         "C00250",  # Pyridoxal
                         "C11378",  # Ubiquinone-10
                         "C05777",  # Coenzyme F430
                         "C00072",  # Ascorbate
                         "C00378",  # Thiamine
                         "C00101",  # Tetrahydrofolate
                         "C00029",  # UDP-glucose
                         "C00068",  # Thiamin diphosphate
                         "C00061",  # FMN
                         "C00063",  # CTP
                         "C05776",  # Vitamin B12
                         "C00113",  # PQQ
                         "C18237",  # Molybdoenzyme molybdenum cofactor
                         "C00051",  # Glutathione
                         "C00010",  # CoA
                         "C00016",  # FAD
                         "C00018",  # Pyridoxal phosphate
                         "C00019",  # S-Adenosyl-L-methionine
                         "C00153",  # Nicotinamide
                         "C04628",  # Coenzyme B
                         "C00862",  # Methanofuran
                         "C15672",  # Heme O
                         "C15670",  # Heme A
                         "C02059",  # Phylloquinone
                         "C03576",  # Coenzyme M
                         "C05441",  # Vitamin D2
                         "C00272",  # Tetrahydrobiopterin
                         "C02477",  # alpha-Tocopherol
                         "C00473",  # Retinol
                         "C00120",  # Biotin
                         "C00725",  # Lipoate
                         "C00053",  # 3'-Phosphoadenylyl sulfate
                         "C00194",  # Cobamide coenzyme
                         "C00255",  # Riboflavin
                         'C00001',  # H2O
                         'C00008',  # ADP
                         'C00013',  # Diphosphate
                         'C00004',  # NADH
                         'C00005',  # NADPH
                         'C00080',  # H+
                         'C00009',  # Orthophosphate
                         'C00008',  # ADP
                         'C00004',  # NADH
                         'C00020',  # AMP
                         'C00007',  # Oxygen
                         'C00015')
  # Getting colnames
  if(!is.null(ko_abundances)){
    column_names <- ko_abundances %>% 
      ungroup() %>% 
      select(-ko_id) %>% 
      colnames()
  } else {
    column_names <- 'no_ko_abundances'
  }
  
  # Get reactions of selected pathways
  sel_rx_path <-  path_rx %>% 
    filter(pathway %in% pathways)
  
  # Get compounds of selected pathways
  sel_cpd_path <- path_compound %>% 
    filter(pathway %in% pathways)
  
  # Get kos of selected pathways
  sel_ko_path <- path_ko %>% 
    filter(pathway %in% pathways)
  
  # Get reactions that match my CPDs and selected pathways
  my_metabolites_rx_from_cpd <- rx_compound %>% 
    filter(compound_id %in% compound_list,
           compound_id %in% sel_cpd_path$compound_id,
           reaction %in% sel_rx_path$reaction)
  
  # Get reactions that match KO ids and selected pathways
  my_ko_rx <- ko_rx %>% 
    group_by(reaction) %>% 
    mutate(kos_per_rx = n()) %>% 
    filter(ko_id %in% ko_list,
           ko_id %in% sel_ko_path$ko_id,
           reaction %in% sel_rx_path$reaction) %>% 
    mutate(kos_in_list = n()) %>% 
    filter(kos_in_list / kos_per_rx >= prevalence_per_reaction)
  
  # Getting reaction expression from kos
  if(!is.null(ko_abundances)){
    my_ko_rx_values <- ko_rx %>% 
      group_by(reaction) %>% 
      mutate(kos_per_rx = n()) %>% 
      inner_join(ko_abundances, by = 'ko_id') %>% 
      filter(ko_id %in% sel_ko_path$ko_id,
             reaction %in% sel_rx_path$reaction) %>%
      pivot_longer(all_of(column_names), names_to = 'id', values_to = 'value') %>% 
      filter(value > 0) %>% 
      group_by(reaction, id) %>% 
      mutate(kos_in_ids = n()) %>% 
      filter(kos_in_ids / kos_per_rx >= prevalence_per_reaction) %>% 
      summarise(getmm_rx = median(value)) %>% 
      pivot_wider(names_from = id, values_from = 'getmm_rx')
    
    min_getmm <- my_ko_rx_values %>% 
      ungroup %>% 
      select(-reaction) %>% 
      min(na.rm = TRUE)
    
    max_getmm <- my_ko_rx_values %>% 
      ungroup %>% 
      select(-reaction) %>% 
      max(na.rm = TRUE) 
  } else {
    min_getmm <- 0.001
    max_getmm <- 1000
  }
  
  
  # Search additional reactions and metabolites based on available kos
  my_metabolites_rx_from_ko <- rx_compound %>% 
    filter(compound_id %in% sel_cpd_path$compound_id,
           !(compound_id %in% compounds_to_omit),
           reaction %in% my_ko_rx$reaction)
  
  
  # Create semi-final table with reactions from compounds and KOs
  
  my_metabolites_rx <- rbind(my_metabolites_rx_from_cpd,
                             my_metabolites_rx_from_ko) %>% 
    filter(reaction %in% sel_rx_path$reaction)
  
  i <- 1
  while(i <= expand_search){
    extra_metabolites <- rx_compound %>% 
      filter(reaction %in% my_metabolites_rx$reaction,
             !(compound_id %in% compounds_to_omit)) %>% 
      pull(compound_id) %>% unique
    
    my_metabolites_rx <- rx_compound %>% 
      filter(compound_id %in% extra_metabolites,
             compound_id %in% sel_cpd_path$compound_id,
             reaction %in% sel_rx_path$reaction)
    i <- i + 1
  }
  
  # Creating a table of edges of compounds and reactions
  # and adding edges weights
  if(!is.null(ko_abundances)){
    edge_table <- my_metabolites_rx %>% 
      left_join(my_ko_rx_values, by = 'reaction') %>% 
      select(to = reaction, from = compound_id, all_of(column_names)) %>% 
      distinct()
  } else {
    edge_table <- my_metabolites_rx %>% 
      select(to = reaction, from = compound_id) %>% 
      mutate(no_ko_abundances = NA)
  }
  
  # Creating node table
  
  
  if(is.null(cid_abundances)){
    node_table <- tibble(name = unique(c(edge_table$to, 
                                         edge_table$from))) %>% 
      mutate(node_type = case_when(str_detect(name, 'C') & 
                                     name %in% compound_list ~ 'Detected metabolite',
                                   str_detect(name, 'C') & 
                                     !(name %in% compound_list) ~ 'Undetected metabolite',
                                   str_detect(name, 'R') & 
                                     name %in% my_ko_rx$reaction ~ 'Reaction present',
                                   str_detect(name, 'R') & 
                                     !(name %in% my_ko_rx$reaction) ~ 'Reaction absent'),
             no_cid_abundances_node = NA)
    
    min_abun <- 3
    max_abun <- 3
    
  } else {
    cid_abundances_temp <- cid_abundances %>% 
      rename_with(~paste0(.x, '_node'))
    
    node_table <- tibble(name = unique(c(edge_table$to, 
                                         edge_table$from))) %>% 
      mutate(node_type = case_when(str_detect(name, 'C') & 
                                     name %in% compound_list ~ 'Detected metabolite',
                                   str_detect(name, 'C') & 
                                     !(name %in% compound_list) ~ 'Undetected metabolite',
                                   str_detect(name, 'R') & 
                                     name %in% my_ko_rx$reaction ~ 'Reaction present',
                                   str_detect(name, 'R') & 
                                     !(name %in% my_ko_rx$reaction) ~ 'Reaction absent')) %>% 
      left_join(cid_abundances_temp, by = c('name' = 'kegg_cid_node'))
    
    min_abun <- cid_abundances_temp %>% 
      select(-kegg_cid_node) %>% 
      min()
    max_abun <- cid_abundances_temp %>% 
      select(-kegg_cid_node) %>% 
      max()
  }
  
  # Creating igraph object
  net_igraph <- graph_from_data_frame(edge_table, directed = FALSE, 
                                      vertices = node_table)
  
  #lay <- graphlayouts::layout_with_sparse_stress(net_igraph, pivots = 100)
  
  net_df <- ggnetwork(net_igraph)
  
  suppressMessages({
    compound_names <- read_csv('scripts/LC_metat_network/input/KEGG_compound_db.csv') %>% 
      select(KEGG_id, KEGG_name) %>% 
      mutate(KEGG_id = str_remove(KEGG_id, 'cpd:'),
             KEGG_name = str_remove(KEGG_name, ';.*'))
  })
  
  path_names <- tibble(pathway_id = pathways) %>% 
    mutate(path_name = map_chr(pathway_id, function(id){
      keggGet(id)[[1]]$NAME[1]
    }))
  
  all_paths <- rbind(sel_cpd_path %>% rename(name = compound_id),
                     sel_rx_path %>% rename(name = reaction))
  
  temp_df <- net_df %>% 
    left_join(compound_names, by = c('name' = 'KEGG_id')) %>%
    mutate(KEGG_name = map2_chr(name, KEGG_name, function(id, kid){
      
      if(is.na(kid) & str_detect(id, 'C')){
        return(keggGet(id)[[1]]$NAME[1])
      } else {
        return(kid)
      }
    })) %>% 
    left_join(all_paths, by = 'name') %>%
    left_join(path_names, by = c('pathway' = 'pathway_id')) %>% 
    mutate(node_type = factor(node_type,
                              levels = c('Detected metabolite',
                                         'Undetected metabolite',
                                         'Reaction present',
                                         'Reaction absent')))
  
  if(!generate_plot){
    return(temp_df)
    
  } else {
    
    net_plots <- purrr::map(column_names, function(id){
      
      id <- rlang::as_string(rlang::ensym(id))
      
      if(is.null(cid_abundances)){
        plot_df <- temp_df %>% 
          select(x, y, name, node_type, xend, yend, KEGG_name, 
                 pathway, path_name, any_of(id), no_cid_abundances_node) %>% 
          rename(getmm_rx = all_of(id),
                 abun_node = no_cid_abundances_node)
      } else {
        plot_df <- temp_df %>% 
          select(x, y, name, node_type, xend, yend, KEGG_name, 
                 pathway, path_name, any_of(c(id, paste0(id, '_node')))) %>% 
          rename(getmm_rx = all_of(id),
                 abun_node = all_of(paste0(id, "_node")))
      }
      
      plot <- plot_df %>% 
        ggplot(aes(x = x, y = y, xend = xend, yend = yend))+
        geom_edges(data = . %>% filter(is.na(getmm_rx)),
                   color = 'gray20',
                   linewidth = .5) +
        geom_edges(data = . %>% filter(!is.na(getmm_rx)),
                   aes(color = getmm_rx,
                       linewidth = getmm_rx)) +
        scale_color_viridis_c(option = 'A', direction = -1,
                              limits = c(min_getmm, max_getmm)) +
        labs(color = 'geTMM') +
        ggnewscale::new_scale_color() +
        geom_nodes(data = . %>% filter(str_detect(node_type, 'metabolite')),
                   aes(shape = node_type,
                       color = node_type),
                   size = 3,
                   fill = 'white') +
        geom_nodes(data = . %>% filter(str_detect(node_type, 'Reaction')),
                   aes(shape = node_type,
                       color = node_type),
                   size = 2,
                   fill = 'white') +
        geom_nodes(data = . %>% filter(!is.na(abun_node)),
                   aes(shape = node_type,
                       color = node_type,
                       size = abun_node),
                   fill = 'white',
                   show.legend = FALSE) +
        geom_nodetext(data = . %>% filter(node_type == 'Undetected metabolite'),
                      aes(label = KEGG_name),
                      size = 1.5) +
        geom_nodetext_repel(data = . %>% filter(node_type == 'Detected metabolite'),
                            aes(label = KEGG_name),
                            fontface = 'bold',
                            size = 2.5,
                            max.overlaps = 100) +
        geom_nodetext(data = . %>% filter(node_type == 'Reaction present'),
                      aes(label = name),
                      size = 1.5) +
        scale_shape_manual(values = c('Detected metabolite' = 19, 
                                      'Undetected metabolite' = 21,
                                      'Reaction present' = 17,
                                      'Reaction absent' = 24)) +
        scale_size_continuous(range = c(8, 20), limits = c(min_abun, max_abun)) +
        # scale_size_manual(values = c('Detected metabolite' = 7, 
        #                              'Undetected metabolite' = 5,
        #                              'Reaction present' = 3,
        #                              'Reaction absent' = 2)) +
        scale_linewidth_continuous(range = c(1, 4)) +
        scale_color_manual(values = c('blue4', 'blue4', '#a3b19f', '#a3b19f')) +
        #ggpubr::get_palette('Set1', nrow(path_names))
        labs(color = 'Node type',
             title = id,
             shape = 'Node type') +
        guides(linewidth = 'none') +
        theme_blank() +
        theme(legend.position = 'bottom',
              plot.title = element_text(face = 'bold', hjust = 0.5, size = 8),
              legend.title.position = 'top')
      
      return(plot)
    })
    
    names(net_plots) <- column_names
    
    return(net_plots)
    
  }
}

