
# Script to calculate differently expressed genes using DESeq2
# Author: Viviana Freire Zapata

library(DESeq2)
library(ggpubr)
library(ggrepel)
library(clusterProfiler)
library(readxl)
library(patchwork)
library(rlist)
library(igraph)
library(ggh4x)
library(ggraph)
library(clusterProfiler)
library(tidyverse)


# Loading metadata

metadata <- read_csv('input/metadata_omics.csv') 

sra_table <- read_csv("input/SraRunTable.csv") %>% 
  select(Run, `Sample Name`) %>% 
  mutate(treatment_code = str_replace_all(`Sample Name`, '-', '_'))


# Loading metaT raw conts

rnaseq_counts <- read_tsv('input/rnaseq_counts_final_2025.txt', skip = 1) %>% 
  filter(!str_detect(Geneid, 'RNA'))

# Fixing names

raw_matrix <- rnaseq_counts %>% 
  select(!c(Chr, Start, End, Strand, Length)) %>% 
  column_to_rownames(var = "Geneid") %>% 
  rename_with(.fn = function(x) str_remove(x, '/xdisk/tfaily/vfreirezapata/metat_2025_final/bowtie_results/')) %>% 
  rename_with(.fn = function(x) str_remove(x, '_mapped_sorted.bam'))

# Summing counts per sample

matrix_sum <- raw_matrix %>% 
  rownames_to_column(var = 'Geneid') %>% 
  pivot_longer(!Geneid, names_to = 'Run', values_to = 'counts') %>% 
  left_join(sra_table, by = 'Run') %>% 
  group_by(treatment_code, Geneid) %>% 
  summarise(counts_sum = sum(counts))

## Reformat -- creating DeSeq object

countData <-  matrix_sum %>% 
  pivot_wider(names_from = treatment_code, values_from = counts_sum) %>% 
  column_to_rownames(var = 'Geneid')


## Creating list of MAGs

taxonomy <- read_csv("input/mag_taxonomy_denovo.csv")

## NOTE: A total of 3 MAGs have only 1 and 6 genes. Delete these MAGs

delete_mags <- c('bcn_bin_145_1', 'bco_bin_69_1', 'bco_bin_177_1')

MAG_list <- taxonomy %>% 
  dplyr::select(bin) %>%
  filter(!bin %in% delete_mags) %>% 
  unique() %>% 
  pull()

names(MAG_list) <- MAG_list

## Create paired vectors of comparisons ti vs tf

t_initial <- c('May', 'July_1', 'July_2', 'July_3', 'August', 'October')
t_final <- c('July_1', 'July_2', 'July_3', 'August', 'October', 'May_22')

deresults_titf <- map2(t_initial, t_final, function(ti, tf){
  
  all_samples <- metadata %>% 
    filter(month == ti | month == tf) %>% 
    pull(treatment_code)
  
  count_data_df <- map(MAG_list, function(MAG){
    
    count_data_ft <- countData %>% 
      rownames_to_column(var = "geneID") %>% 
      mutate(MAG_id = str_remove(geneID, "_k.*")) %>% 
      filter(MAG_id == MAG) %>% 
      dplyr::select(geneID, MAG_id, all_of(c(all_samples)))
    
    count_data_final <- count_data_ft %>%
      dplyr::select(-MAG_id) %>% 
      column_to_rownames(var = "geneID")
    
    return(count_data_final)
    
  })
  
  ############# DeSeq2 per MAG in time #######################
  
  
  colData <- metadata %>% 
    filter(month == ti| month == tf) %>% 
    column_to_rownames(var = 'treatment_code')
  
  DE_cluster_dds <- map(count_data_df, function(df){
    
    deseq_obj <- DESeqDataSetFromMatrix(countData = df,
                                        colData = colData,
                                        design = ~ month)
    deseq_obj$month <- relevel(deseq_obj$month, ref = ti)
    
    keep <- rowSums(counts(deseq_obj) >= 10) >= 4
    
    dds <- try(DESeq(deseq_obj[keep,]))
    
    if(is(dds, 'try-error')){
      dds <- 'All genes have at least one zero. No DE analysis possible'
    }
    
    return(dds)
    
  })
  
  ## Note: If MAGs have genes matrix with at least 1 zero, DESeq2 has an error
  ## 9 MAGs had this issue and the DE analyses wasn't performed. This MAGs were filtered out
  
  idx <- map(DE_cluster_dds, function(x){
    !is.character(x)
  }) %>% 
    unlist()
  
  ## Continue analysis with 273 MAGs
  
  DE_cluster_dds_filt <- DE_cluster_dds[idx]
  
})



## Saving DE results per MAG, to avoid running again

write_rds(deresults_titf, 'output/DE_month_10_2025.rds')


## LFC shrinkage

comparisons_names <- c('month_July_1_vs_May',
                       'month_July_2_vs_July_1',
                       'month_July_3_vs_July_2',
                       'month_August_vs_July_3',
                       'month_October_vs_August',
                       'month_May_22_vs_October')

names(comparisons_names) <- c('from_May_to_July_1',
                              'from_July_1_to_July_2',
                              'from_July_2_to_July_3',
                              'from_July_3_to_August',
                              'from_August_to_October',
                              'from_October_to_May_22')

names(deresults_titf) <- comparisons_names



# The res object has the dataframes for each month: res$July

lfc_results <- imap(deresults_titf, function(de_obj, comparison){
  
  lfc_mag <- map(de_obj, function(dds){
    
    # Recommended to shrink logfoldchange
    # Each month can be visualized as lfc_res$July
    
    res_shrink <- lfcShrink(dds,
                            coef = comparison,
                            type = 'apeglm')
    
    df <- as.data.frame(res_shrink) %>% 
      drop_na(padj) %>% 
      mutate(comparison = comparison) %>% 
      rownames_to_column(var = 'FeatureID')
    
    return(df)
    
  })
  
  res <- purrr::reduce(lfc_mag, rbind)
  
  return(res)
  
})

lfc_final <- purrr::reduce(lfc_results, rbind)


write_csv(lfc_final, "output/DE_result_month_slice_10_2025.csv")
