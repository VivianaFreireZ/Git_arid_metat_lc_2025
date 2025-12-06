
## Topic: Retrieve KEGG Modules definitions


library(KEGGREST)
library(jsonlite)
library(tidyverse)


## Use keggLink to match KOs with reactions and compounds

## Check KEGG databases 

listDatabases()

## From KOs to Reaction ID

ko_rex <- keggLink('ko', 'reaction')

## From Reactions to Compounds

rex_compound <- keggLink('reaction', 'compound')

## From Pathways to Compounds

path_compound <- keggLink('pathway', 'compound')

## From pathways to kos

path_ko <- keggLink('pathway', 'ko')

## Creating table with KO, RX, Compound

table_rex <- tibble(ko_id = ko_rex,
                    reaction = names(ko_rex))

write_csv(table_rex, "scripts/LC_metat_network/output/KO_to_reaction_2025.csv")


table_compound <- tibble(reaction = rex_compound,
                         compound_id = names(rex_compound))

write_csv(table_compound, "scripts/LC_metat_network/output/Reaction_to_compound_2025.csv")

table_compound_path <- tibble(pathway = path_compound,
                              compound_id = names(path_compound))

write_csv(table_compound_path, "scripts/LC_metat_network/output/Pathway_to_compound_2025.csv")


table_path_ko <- tibble(pathway = path_ko,
                        ko_id = names(path_ko)) %>% 
  filter(str_detect(pathway, 'map'))

write_csv(table_path_ko, "scripts/LC_metat_network/output/Pathway_to_ko_2025.csv")


## Retrieving compounds names 

compound_list <- keggList('compound')


compounds_names <- tibble(compound_name = compound_list,
                          compound_id = names(compound_list))

## Retrieving KO pathways


pathways <- keggLink('ko', 'pathway')

table_path <- tibble(ko_id = pathways,
                     pathway_id = names(pathways))


## Retrieving Reaction to pathway


pathways_2 <- keggLink('pathway', 'reaction')

table_path2 <- tibble(pathway = pathways_2,
                      reaction = names(pathways_2)) %>% 
  filter(str_detect(pathway, 'map'))

write_csv(table_path2, "scripts/LC_metat_network/output/Pathway_to_reaction_2025.csv")

