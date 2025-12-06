# ***********************************************
# ** Script for querying PubChem using its API **
# ***********************************************

library(httr)
library(tidyverse)

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


########################
#### How to use it ####
#######################

# DATA RETRIEVAL BY NAMES, INCHI OR SMILES

# 1. Enter the identifier of the compound you want to search in the "retrieve_data_pubchem" function

## For example to search by name

res <- retrieve_data_pubchem('Glucose', search_by = 'name')

## For example to search by inchikey

res2 <- retrieve_data_pubchem('LCTONWCANYUPML-UHFFFAOYSA-M', search_by = 'inchikey')

# 2. To search multiple names at the same time

## 2.1 Create a vector of identifiers to search

compound_names <- c('Glucose', 'Fructose', 'Ethanol')

## 2.2 Alternatively if you have a dataframe, extract the column with the names

df <- tibble(Name = c('Glucose', 'Fructose', 'Ethanol'))
compound_names <- df$Name

## 2.3 Use the map function make multiple searches. 

## Pubchem does not like if a lot of searches are done at the same time.
## So if you are searching many compounds at the time, add some waiting time between searches

res_list <- map(compound_names, ~ retrieve_data_pubchem(.x, search_by = 'name', wait = 1))

## 2.4 Concatenate the results in a single dataframe

final_df <- do.call(rbind, res_list)


# ----------------------------------
# DATA RETRIEVAL BASED ON FORMULA

# 1. Enter the formula to search in the "start_search_formulas_pubchem" function

## For example to search by name

key <- start_search_formulas_pubchem('C88H3O93')

# 2. Wait about 2 to 5 minutes and then enter key into he "retrieve_search_formulas_pubchem" function

res <- retrieve_search_formulas_pubchem(key)

# 2. To search multiple names at the same time

## 2.1 Create a vector of identifiers to search

compound_mf <- c('C6H12O6', 'CH4', 'C27H30O16')

## 2.2 Alternatively if you have a dataframe, extract the column with the formulas

df <- tibble(Formula = c('C6H12O6', 'CH4', 'C27H30O16'))
compound_mf <- df$Formula

## 2.3 Use the map function make multiple searches. 

## Pubchem does not like if a lot of searches are done at the same time.
## So if you are searching many compounds at the time, add some waiting time between searches

key_list <- map(compound_mf, ~ start_search_formulas_pubchem(.x, wait = 1))
res_list <- map(key_list, ~ retrieve_search_formulas_pubchem(.x, wait = 1))

## 2.4 Concatenate the results in a single dataframe

final_df <- do.call(rbind, res_list)



############# Checking if predicted formulas are real

## HILIC

hilic_formula <- read_csv("scripts/LC-MS-MS/3.Cleanning_Annotation/output/HILIC_final_annotation.csv") %>% 
  filter(Name == 'no name' & Formula != 'no formula') %>% 
  mutate(Formula = str_remove_all(Formula, ' ')) %>% 
  select(Formula) %>%
  distinct() %>%
  pull()

## Formula seach in pubchem

## Breaking into steps

# 1 -50

key_list <- map(hilic_formula[1:50], ~ start_search_formulas_pubchem(.x, wait = 1))

res_list_1_50 <- map(key_list, ~ retrieve_search_formulas_pubchem(.x, wait = 1))


# 51 - 100

key_list_2 <- map(hilic_formula[51:100], ~ start_search_formulas_pubchem(.x, wait = 1))


res_list_51_100 <- map(key_list_2, ~ retrieve_search_formulas_pubchem(.x, wait = 1))


# 101 - 250

key_list_3 <- map(hilic_formula[101:250], ~ start_search_formulas_pubchem(.x, wait = 1))


res_list_101_250 <- map(key_list_3, ~ retrieve_search_formulas_pubchem(.x, wait = 1))


# 251 - 400

key_list_4 <- map(hilic_formula[251:400], ~ start_search_formulas_pubchem(.x, wait = 1))


res_list_251_400 <- map(key_list_4, ~ retrieve_search_formulas_pubchem(.x, wait = 1))


# 401 - 549

key_list_5 <- map(hilic_formula[401:549], ~ start_search_formulas_pubchem(.x, wait = 1))


res_list_401_549 <- map(key_list_5, ~ retrieve_search_formulas_pubchem(.x, wait = 1))


## Final table

join_table <- c(res_list_1_50, res_list_51_100,res_list_101_250, 
                res_list_251_400, res_list_401_549)

res_final <- reduce(join_table, plyr::rbind.fill)


## ID unreal formulas

delete_formula <- res_final %>% 
  filter(is.na(CID),
         Query_formula != 'C23H24O9') %>% 
  distinct()


## Filtered annot. table

## Ambiguous formula assigned as level 4

hilic_formula_check <- read_csv("scripts/LC-MS-MS/3.Cleanning_Annotation/output/HILIC_final_annotation.csv") %>% 
  mutate(Formula = str_remove_all(Formula, ' '),
         use_Sirius_class = ifelse(use_Sirius_class %in% c('yes', 'YES'), 'yes', 'no'),
         Annotation_level = ifelse(Formula %in% delete_formula$Query_formula & use_Sirius_class == 'no', 
                                   'L4', Annotation_level),
         Formula = ifelse(Formula %in% delete_formula$Query_formula & use_Sirius_class == 'no', 
                                   NA, Formula))

write_csv(hilic_formula_check, "scripts/LC-MS-MS/3.Cleanning_Annotation/output/HILIC_final_annotation_checked.csv")

################################################################################

## RP

rp_formula <- read_csv("scripts/LC-MS-MS/3.Cleanning_Annotation/output/RP_final_annotation.csv") %>% 
  filter(Name == 'no name' & Formula != 'no formula' & Formula != "halogenated,unknown"
         & Formula != 'na') %>% 
  mutate(Formula = str_remove_all(Formula, ' ')) %>% 
  select(Formula) %>%
  distinct() %>%
  pull()


# 1 -50

key_list_rp_1 <- map(rp_formula[1:50], ~ start_search_formulas_pubchem(.x, wait = 1))

res_list_rp_1 <- map(key_list_rp_1, ~ retrieve_search_formulas_pubchem(.x, wait = 1))



# 51 - 100

key_list_rp_2 <- map(rp_formula[51:100], ~ start_search_formulas_pubchem(.x, wait = 1))

res_list_rp_2 <- map(key_list_rp_2, ~ retrieve_search_formulas_pubchem(.x, wait = 1))


# 101 - 150

key_list_rp_3 <- map(rp_formula[101:150], ~ start_search_formulas_pubchem(.x, wait = 1))

res_list_rp_3 <- map(key_list_rp_3, ~ retrieve_search_formulas_pubchem(.x, wait = 1))


# 151 - 177

key_list_rp_4 <- map(rp_formula[151:177], ~ start_search_formulas_pubchem(.x, wait = 1))

res_list_rp_4 <- map(key_list_rp_4, ~ retrieve_search_formulas_pubchem(.x, wait = 1))

## Final table

join_table_rp <- c(res_list_rp_1, res_list_rp_2, res_list_rp_3, res_list_rp_4)

res_final_rp <- reduce(join_table_rp, plyr::rbind.fill)


## ID unreal formulas

delete_formula <- res_final_rp %>% 
  filter(is.na(CID),
         Query_formula != 'halogenated,unknown' & Query_formula != 'na') %>% 
  distinct()


## Filtered annot. table

## Ambiguous formula assigned as level 4

rp_formula_check <- read_csv("scripts/LC-MS-MS/3.Cleanning_Annotation/output/RP_final_annotation.csv") %>% 
  mutate(Formula = str_remove_all(Formula, ' '),
         use_Sirius_class = ifelse(use_Sirius_class %in% c('yes', 'YES'), 'yes', 'no'),
         Annotation_level = ifelse(Formula %in% delete_formula$Query_formula & use_Sirius_class == 'no', 
                                   'L4', Annotation_level),
         Formula = ifelse(Formula %in% delete_formula$Query_formula & use_Sirius_class == 'no',
                          NA, Formula))

write_csv(rp_formula_check, "scripts/LC-MS-MS/3.Cleanning_Annotation/output/RP_final_annotation_checked.csv")

