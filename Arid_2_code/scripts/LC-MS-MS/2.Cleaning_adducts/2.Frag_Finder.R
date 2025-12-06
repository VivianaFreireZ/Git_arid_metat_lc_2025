library(tidyverse)
library(psych)

## Author: John Bouranis - Tfaily lab

## Function

find_ppm <- function(mo, me){
  ppm <- (mo - me)/me * 10^6
  return(ppm)
}

read_mgf <- function(path){
  #Read in the mgf file
  l <- readLines(path)
  #Find the start of each feature
  ti <- which(str_detect(l, 'FEATURE_ID'))
  #Find the end of each feature
  tt <- which(str_detect(l, 'END IONS'))
  ck <- data.frame()
  #Pull out the necessary information
  for(i in 1:length(ti)){
    tmp <- tibble(Compound_ID = str_split(l[ti[i]], '=')[[1]][2],
                  ions = list(as.numeric(sapply(str_split(l[(ti[i]+8):(tt[i]-1)], ' '), function(x) x[1]))))
    ck <- rbind(ck, tmp)
  }
  #Check to see if "FT_" is at the start of feature names
  if(str_detect(ck$Compound_ID[1], 'FT_', negate = TRUE)){
    #If not, add it on
  frags <- ck %>%
    modify_at('Compound_ID', ~paste0('FT_', .x))
  } else {
    frags <- ck
  }
    
  return(frags)
}

find_ISFrags <- function(link_data, intensity_data, mgf_file,
                         time_tol = 0.005, ppm_tol = 5, cor_cutoff = 0.9){
  
  #Check data integrity:
  if(sum(colnames(link_data) %in% c('Compound_ID', 'm/z', 'RT [min]')) != 3){
    stop('Check Link Data, it must contain the following columns named exactly: Compound_ID, RT [min], m/z')
  }
  
  #Format the fragment table to be correct
  fragtab <- link_data %>%
    dplyr::select(Compound_ID, `m/z`, `RT [min]`) %>%
    dplyr::rename('mz' = `m/z`, 'rt' = `RT [min]`) 
  
  #Calculate the difference between each rt and all others:
  #Output the results into a list
  olist <- list()
  for(i in 1:nrow(fragtab)){
    ovec <- fragtab$rt - fragtab$rt[i]
    names(ovec) <- fragtab$Compound_ID
    olist[[i]] <- ovec
  }
  
  #Now trim each vector to only places where the rt dif <= 1 and pull the names
  olist_clean <- lapply(olist, function(x) names(x[abs(x) <= time_tol]))
  
  #Remove features that are a length of 1 - These have no matches
  olist_trimmed <- olist_clean[sapply(olist_clean, length) != 1]
  #And each set will show up twice so trim it down
  olist_final <- olist_trimmed[!duplicated(olist_trimmed)]
  
  #Now turn it into a dataframe
  #Making a unique group for each set
  dfo <- tibble()
  for(i in 1:length(olist_final)){
    df <- tibble(group = paste0('G', i),
                 matches = list(olist_final[[i]]))
    dfo <- rbind(dfo, df)
  }
  
  #Now make a table with the groups
  group_tab <- dfo %>%
    unnest('matches') %>%
    dplyr::rename('Compound_ID' = matches)
  
  #Make a table of intensities
  int_table <- intensity_data %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('Compound_ID')
  
  #Make a table with a correlation matrix of each group
  big_group <- group_tab %>%
    #Join on teh intensity data
    left_join(int_table) %>%
    #Make long
    pivot_longer(cols = c(-group, -Compound_ID), names_to = 'sample') %>%
    #Group
    group_by(group) %>%
    nest() %>%
    #Make a correlation matrix for each set of candidates
    mutate(cor_tab = purrr::map(data, function(x){
      ct <- x %>%
        pivot_wider(names_from = Compound_ID) %>%
        column_to_rownames('sample') %>%
        cor()
      return(ct)
    })) 
  
  #Filter the candidates based on correlation info:
  IS_cand <- big_group %>%
    mutate(cands = purrr::map_chr(cor_tab, function(x){
      x %>%
        as.data.frame() %>%
        rownames_to_column('var1') %>%
        #Make the correlation table long for ease
        pivot_longer(starts_with('FT_'), names_to = 'var2') %>%
        #Remove the diag
        filter(var1 != var2) %>%
        #Grab only ones above the cutoff
        filter(value >= cor_cutoff) %>%
        #Pull out var1 since it will include both sides
        pull(var1) %>%
        #Remove duplicates
        unique() %>%
        paste0(collapse = '; ')
    })) %>%
    ungroup() %>%
    #Remove candidates where nothing passed
    filter(cands != '') %>%
    #Clear out the duplicates
    filter(!duplicated(cands))
  
  #Now using the grouping information prep for MS2 data
  grp <- IS_cand %>%
    #Clean up the fragment data
    mutate(split = str_split(cands, pattern = '; ')) %>%
    dplyr::select(group, split) %>%
    unnest('split') %>%
    #Rename for merging
    dplyr::rename('Compound_ID' = split) %>%
    #Merge in the fragment data
    left_join(fragtab) %>%
    group_by(group) %>%
    #Keep the highest mz as the parent
    mutate(type = ifelse(max(mz) == mz, 'parent', 'frag'))
  
  frags <- read_mgf(mgf_file)
  
  #Merge in the MS2 data
  grp_frag <- grp %>%
    left_join(frags) %>%
    #Grab the necessary information from the fragment tab
    reframe(parent_mz = mz[type == 'parent'],
            parent_ft = Compound_ID[type == 'parent'],
            frag_mz = list(mz[type == 'frag']),
            frag_fts = list(Compound_ID[type == 'frag']),
            ions_parent = ions[type == 'parent']) %>%
    #Calculate PPMs of each fragment candidate against the MS2 data of parent fragment
    mutate(ppm_vals = map2(frag_mz, ions_parent, function(x,y){
      purrr::map(x, function(a){
        find_ppm(a,y)
      })
    })) %>%
    mutate(match = purrr::map2(ppm_vals, frag_fts, function(x,y){
      lgls <- sapply(x, function(x) any(abs(x) <= ppm_tol))
      data.frame(frag_cand = y,
                 match = lgls)
    })) %>%
    dplyr::select(parent_ft, match) %>%
    unnest('match') %>%
    filter(match)
  
  final_output <- grp_frag %>%
    dplyr::select(-match) %>%
    group_by(parent_ft) %>%
    summarise(fragments = paste0(frag_cand, collapse = '; '))

  return(final_output)
}

#link_data contains columns with m/z, RT [min.], and Compound_ID
#Intensity data in a numeric matrix of RAW (or adjusted for things like sample weight) of feature intensities
#Each column is a feature and each row is a sample with sample names as rownames ONLY
#mgf_file is the path of the mgf file you gave to SIRIUS
#time_tol is the tolerance for finding initial candidates - Defaults to 0.005
#cor_cutoff is the cutoff for correlations between two features to be considered a fragment - Default is 0.9
#ppm_tol is the ppm tolerance between candidate fragments (which pass the correlation test) and 
#the parent (heaviest) ion MS2 fragments - Default is 5


## Loading RP and HILIC RAW matrix

rp_matrix <- read_csv("1.Pre-processing/output/RP_matrix_clean_2025.csv") %>% 
  select(!c(FeatureID_mgf, Name, Formula, `Calc. MW`, `m/z`, `RT [min]`, 
            `Annotation MW`,`mzVault Best Match`, `mzCloud Best Match`, 
            `KEGG ID`,`Annot. Source: Predicted Compositions`, 
            `Annot. Source: mzCloud Search`,
            `Annot. Source: mzVault Search`, 
            `# mzCloud Results`, `# mzVault Results`))

hilic_matrix <- read_csv("1.Pre-processing/output/HILIC_matrix_clean_2025.csv") %>% 
  select(!c(FeatureID_mgf, Name, Formula, `Calc. MW`, `m/z`, `RT [min]`, 
            `Annotation MW`,`mzVault Best Match`, `mzCloud Best Match`, 
            `KEGG ID`,`Annot. Source: Predicted Compositions`, 
            `Annot. Source: mzCloud Search`,
            `Annot. Source: mzVault Search`, 
            `# mzCloud Results`, `# mzVault Results`))

## Creating link data

## RP

link_rp <- rp_matrix %>% 
  select(FeatureID) %>% 
  separate_wider_delim(FeatureID, delim = "_", names = c("FT", "RT [min]", "m/z"), cols_remove = FALSE) %>% 
  select(!FT) %>% 
  rename(Compound_ID = FeatureID) %>% 
  mutate(`m/z` = as.numeric(`m/z`),
         `RT [min]` = as.numeric(`RT [min]`))

## HILIC

link_hilic <- hilic_matrix %>% 
  select(FeatureID) %>% 
  separate_wider_delim(FeatureID, delim = "_", names = c("FT", "RT [min]", "m/z"), cols_remove = FALSE) %>% 
  select(!FT) %>% 
  rename(Compound_ID = FeatureID) %>% 
  mutate(`m/z` = as.numeric(`m/z`),
         `RT [min]` = as.numeric(`RT [min]`))

## RP 
rp_ready <- rp_matrix %>% 
  column_to_rownames(var = "FeatureID")

rp_ready <- as.data.frame(t(rp_ready))  

## HILIC

hilic_ready <- hilic_matrix %>% 
  column_to_rownames(var = "FeatureID")

hilic_ready <- as.data.frame(t(hilic_ready)) 

## Testing function HILIC

hilic_test <- find_ISFrags(link_data = link_hilic, intensity_data = hilic_ready, 
                        mgf_file = '2.Cleaning_adducts/input/HILIC_mgf_modified_matching_featureID_adduct.mgf',
                        time_tol = 0.015, cor_cutoff = 0.99)

write_csv(hilic_test, "2.Cleaning_adducts/output/Hilic_fragments_2025.csv")


## Testing function RP

rp_test <- find_ISFrags(link_data = link_rp, intensity_data = rp_ready, 
                           mgf_file = '2.Cleaning_adducts/input/RP_mgf_modified_matching_featureID_adduct.mgf',
                           time_tol = 0.015, cor_cutoff = 0.99)

write_csv(rp_test, "2.Cleaning_adducts/output/RP_fragments_2025.csv")
