# Author: VFZ
# Topic: MAG Gene count normalization - TMM


library(edgeR)
library(tidyverse)



## Loading gene counts

mag_genes_counts <- read_tsv("input/metag_genecounts_2025.txt", skip = 1)

## Run table

sra_table <- read_csv("input/SraRun_metag.csv") %>% 
  dplyr::select(Run, `Sample Name`)

## annotation

annot <- read_csv("input/annot_with_cydbs_2.23.25.csv") %>% 
  rename(geneid = `...1`)

## Fixing names

raw_matrix <- mag_genes_counts %>% 
  dplyr::select(!c(Chr, Start, End, Strand, Length)) %>% 
  column_to_rownames(var = "Geneid") %>% 
  rename_with(.fn = function(x) str_remove(x, '/xdisk/tfaily/vfreirezapata/kraken_test/gene_count_map/')) %>% 
  rename_with(.fn = function(x) str_remove(x, '_mapped_sorted.bam'))

## Delete deep sequenced runs

delete_sr <- c('SRR31679670', 'SRR31679671', 'SRR31679672', 'SRR31679673', 
               'SRR31679674', 'SRR31679675')

# Summing counts per sample

matrix_sum <- raw_matrix %>% 
  rownames_to_column(var = 'Geneid') %>% 
  pivot_longer(!Geneid, names_to = 'Run', values_to = 'counts') %>% 
  left_join(sra_table, by = 'Run') %>% 
  filter(!Run %in% delete_sr) %>% 
  group_by(`Sample Name`, Geneid) %>% 
  summarise(counts_sum = sum(counts))

# Summing genes with same KO 

ko_matrix <- matrix_sum %>% 
  left_join(annot, by = c('Geneid' = 'geneid')) %>% 
  filter(!is.na(ko_id)) %>% 
  group_by(fasta, ko_id, `Sample Name`) %>% 
  summarise(sum_counts = sum(counts_sum)) 

## Reformat


countData_temp <- ko_matrix %>%
  pivot_wider(names_from = 'Sample Name', values_from = 'sum_counts', values_fill = 0) %>%
  ungroup() %>%
  unite(fasta_ko, fasta, ko_id, sep = '__') %>%
  column_to_rownames(var = 'fasta_ko')


## Removing genes with zero counts 

raw_matrix_ft <- countData_temp[rowSums(countData_temp >= 5) >= 3,]

# TMM normalization

## Create DGEList and filter

dge <- DGEList(counts=raw_matrix_ft)


## Compute TMM factors
dge <- calcNormFactors(dge, method="TMM") 


## Save TMM normalized counts

final_mx <- as.data.frame(dge$counts)%>% 
  rownames_to_column(var = 'mag_geneid')

write_csv(final_mx, "output/TMM_gene_counts_mags_2025.csv")
