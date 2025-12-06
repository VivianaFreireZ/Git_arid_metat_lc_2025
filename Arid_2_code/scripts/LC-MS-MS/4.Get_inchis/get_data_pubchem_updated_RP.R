# ***********************************************
# ** Script for querying PubChem using its API **
# ***********************************************

library(httr)
library(tidyverse)
library(classyfireR)

################
## FUNCTIONS ##
###############

#### Safe GET ####
# Function to avoid errors while accessing APIs

safe_GET <- safely(httr::GET)

#### Retrieve data from PubChem ####

# Function to get data from PubChem based on the compound name
# Inputs:
#   name: Vector with compound names to search
#   search_by: which type of identifier to use for searching
#   wait: Time in seconds to wait after finishing a search

retrieve_data_pubchem <- function(query, search_by = 'name', wait = NA){
  
  if(!(search_by %in% c('name', 'inchikey', 'inchi', 'smiles', 'cid'))){
    stop('search_by option should be one of: name, inchikey, inchi, smiles or cid')
  }
  
  new_name <- query %>% 
    str_replace_all(' ', '%20')
  
  url <- paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/',
                search_by,
                '/',
                new_name,
                '/property/Title,ExactMass,MolecularFormula,InChi,InChiKey,CanonicalSMILES/CSV')
  
  r <- safe_GET(url)
  
  if(is.null(r$error)){
    content <- try(httr::content(r$result, 'parsed', type = 'text/csv', encoding = "UTF-8"))
    
    if(ncol(content) < 3){
      content <- tibble("CID" = numeric(0),
                        "Title" = character(0),
                        "search" = character(0),
                        "ExactMass" = numeric(0),    
                        "MolecularFormula" = character(0),
                        "InChI" = character(0),
                        "InChIKey" = character(0),
                        "SMILES" = character(0)
      )
    }
    
  } else {
    content <- tibble("CID" = numeric(0),
                      "Title" = character(0),
                      "search" = character(0),
                      "ExactMass" = numeric(0),    
                      "MolecularFormula" = character(0),
                      "InChI" = character(0),
                      "InChIKey" = character(0),
                      "SMILES" = character(0)
    )
  }
  
  content <- content %>% 
    mutate(Query = query) %>% 
    select(Query, CID, Title_PubChem = Title, everything())
  
  if(!is.na(wait)){
    Sys.sleep(wait)
  }
  
  return(content)
  
}


#### Search Formulas ####

start_search_formulas_pubchem <- function(formula, wait = NA){
  
  #### Start PubChem search (using molecular formulas) ####
  
  # This function will start the molecular formula search and give you a listkey to access the search later
  # Inputs:
  #   formula: Molecular formula to search
  
  url <- paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/',
                formula,
                '/property/MolecularFormula/JSON')
  
  r <- safe_GET(url)
  
  if(is.null(r$error)){
    content <- try(httr::content(r$result, 'parsed', type = 'application/json', encoding = "UTF-8"))
    
    listkey <- content$Waiting$ListKey
    
  } else {
    listkey <- NA
    print('Search did not start')
  }
  
  if(!is.na(wait)){
    Sys.sleep(wait)
  }
  
  listkey_df <- tibble(Query_formula = formula,
                       listkey = listkey)
  
}


retrieve_search_formulas_pubchem <- function(key_df, wait = NA){
  
  #### Retrieve data from PubChem (using molecular formula) ####
  
  # Function to get data from PubChem based on the compound name
  # Inputs:
  #   key: Dataframe with key produced by previous search
  
  formula <- key_df$Query_formula
  key <- key_df$listkey
  
  url <- paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/',
                key,
                '/property/Title,ExactMass,MolecularFormula,InChi,InChiKey,CanonicalSMILES/CSV')
  
  r <- safe_GET(url)
  
  if(is.null(r$error)){
    content <- try(httr::content(r$result, 'parsed', type = 'text/csv', encoding = "UTF-8"))
    
    if(str_detect(content[1,1], 'Code: PUGREST.BadRequest')){
      content <- tibble("CID" = NA,
                        "Title" = as.character(content[2,1]),
                        "search" = NA,
                        "ExactMass" = NA,    
                        "MolecularFormula" = NA,
                        "InChI" = NA,
                        "InChIKey" = NA,
                        "SMILES" = NA
      )
    }
    
  } else {
    content <- tibble("CID" = numeric(0),
                      "Title" = character(0),
                      "search" = character(0),
                      "ExactMass" = numeric(0),    
                      "MolecularFormula" = character(0),
                      "InChI" = character(0),
                      "InChIKey" = character(0),
                      "SMILES" = character(0)
    )
  }
  
  content <- content %>% 
    mutate(Query_formula = formula) %>% 
    select(Query_formula, CID, Title_PubChem = Title, everything())
  
  if(!is.na(wait)){
    Sys.sleep(wait)
  }
  
  return(content)
}


############# Running in my data

## From name

rp_ready <- read_csv("scripts/LC-MS-MS/3.Cleanning_Annotation/output/RP_final_annotation_checked.csv")

rp_names <- rp_ready %>% 
  filter(Name != 'no name') %>% 
  pull(Name) %>% 
  unique()

res_list <- map(rp_names, ~ retrieve_data_pubchem(.x, search_by = 'name', wait = 1))

rp_inchi <- reduce(res_list, rbind)

## Retrieving Class using ClassyfireR

classification_list <- map(rp_inchi$InChIKey, get_classification)
names(classification_list) <- rp_inchi$Query

classyfire_df <- imap(classification_list, function(cf_obj, name){
  
  temp <- tibble(Level = c('kingdom', 'superclass', 'class', 'subclass', 'level 5'))
  
  if(!is.null(cf_obj)){
    df <- classification(cf_obj) %>% 
      right_join(temp, by = 'Level') %>% 
      mutate(Query = name) %>% 
      select(-CHEMONT) %>% 
      pivot_wider(names_from = Level, values_from = Classification) %>% 
      select(Query, kingdom, superclass, class, subclass, `level 5`)
  } else {
    df <- tibble(Query = name,
                 kingdom = NA,
                 superclass = NA,
                 class = NA,
                 subclass = NA,
                 `level 5` = NA)
  }
  return(df)
})

classyfire_df_ready <- reduce(classyfire_df, rbind)

## Joining to my data 

rp_ready_v1 <- rp_ready %>% 
  filter(Name != 'no name') %>% 
  left_join(classyfire_df_ready, by = c('Name' = 'Query')) %>% 
  left_join(rp_inchi, by = c('Name' = 'Query'))


write_csv(rp_ready_v1, 'scripts/LC-MS-MS/4.Get_inchis/output/RP_classyfire_from_name.csv')

## From Sirius

sirius_class <- read_tsv("scripts/LC-MS-MS/Sirius/RP/canopus_compound_summary.tsv") %>% 
  select(!id) %>% 
  distinct()

rp_sirius_class <- rp_ready %>% 
  filter(use_Sirius_class == 'yes') %>% 
  left_join(sirius_class, by = c('FeatureID' = 'featureId'))


write_csv(rp_sirius_class, 'scripts/LC-MS-MS/4.Get_inchis/output/RP_classyfire_canopus_v1.csv')


## Need to add from GUI

gui_features <- c('FT_1.262_380.20290', 'FT_1.264_374.12946', 'FT_1.360_350.19235',
                  'FT_1.416_275.12378', 'FT_1.428_448.16617', 'FT_11.232_326.12348',
                  'FT_8.272_360.20173', 'FT_9.687_356.09768')

rp_sirius_class <- rp_ready %>% 
  filter(FeatureID %in% gui_features) %>% 
  left_join(sirius_class, by = c('FeatureID' = 'featureId'))


write_csv(rp_sirius_class, 'scripts/LC-MS-MS/4.Get_inchis/output/RP_classyfire_canopus_v2.csv')




