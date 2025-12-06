suppressPackageStartupMessages({
  library(tidyverse)
  library(KEGGREST)
}) 

#### List of compounds in the database ####

compound_list <- keggList('compound')

# Loop to get info about compounds (takes long (2-3 h) do not run unless you need it)

compound_df <- tibble(KEGG_id = NA, KEGG_name = NA, KEGG_formula = NA, KEGG_pathway_id = NA,
                      KEGG_pathway = NA, KEGG_module = NA, KEGG_brite = NA, KEGG_enzyme = NA, KEGG_reaction = NA, .rows = 0)


for(i in 1:length(compound_list)){
  id <- names(compound_list[i])
  cpd_info <- keggGet(id)
  for(j in length(cpd_info)){
    temp <- tibble(KEGG_id = id,
                   KEGG_name = paste(cpd_info[[j]]$NAME, collapse = ''),
                   KEGG_pathway_id = ifelse(!is.null(names(cpd_info[[j]]$PATHWAY)), paste(names(cpd_info[[j]]$PATHWAY), collapse = ';'), NA),
                   KEGG_pathway = ifelse(!is.null(cpd_info[[j]]$PATHWAY), paste(cpd_info[[j]]$PATHWAY, collapse = ';'), NA),
                   KEGG_module = ifelse(!is.null(cpd_info[[j]]$MODULE), paste(cpd_info[[j]]$MODULE, collapse = ';'), NA),
                   KEGG_brite = ifelse(!is.null(cpd_info[[j]]$BRITE), paste(cpd_info[[j]]$BRITE, collapse = ';'), NA),
                   KEGG_enzyme = ifelse(!is.null(cpd_info[[j]]$ENZYME), paste(cpd_info[[j]]$ENZYME, collapse = ';'), NA),
                   KEGG_reaction = ifelse(!is.null(cpd_info[[j]]$REACTION), paste(cpd_info[[j]]$REACTION, collapse = ';'), NA))
  }
  compound_df <- rbind(compound_df, temp)
  print(paste(Sys.time(), 'Compound No:', i))
}

filename <- file.path('KEGG_compound_db.csv')
write_csv(compound_df, filename)


