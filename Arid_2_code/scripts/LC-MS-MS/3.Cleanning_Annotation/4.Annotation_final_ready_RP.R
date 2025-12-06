
library(tidyverse)

# Loading cleaned annotation files


L1 <- read_csv('scripts/LC-MS-MS/3.Cleanning_Annotation/input/manually_curated_tables/RP_level_1_mzvault_ready.csv')


L2 <- read_csv('scripts/LC-MS-MS/3.Cleanning_Annotation/input/manually_curated_tables/RP_level_2_flat_mzCloud_2024_ready.csv')


L2_lower <- read_csv('scripts/LC-MS-MS/3.Cleanning_Annotation/input/manually_curated_tables/RP_level_2_lower_mzCloud_ready.csv')


L3 <- read_csv('scripts/LC-MS-MS/3.Cleanning_Annotation/input/manually_curated_tables/RP_level_3_CD_Sirius_CMMRT_ready.csv')


L3_lower <- read_csv('scripts/LC-MS-MS/3.Cleanning_Annotation/input/manually_curated_tables/RP_level_3_lower_CD_Sirius_CMMRT_ready.csv')


## Data wrangling

L1_ready <- L1 %>% 
  filter(Final_decision != 'delete') %>% 
  select(FeatureID, Name, Formula, Final_decision) %>% 
  rename(Annotation_level = Final_decision) %>% 
  mutate(use_Sirius_class = NA)

L2_ready <- L2 %>%
  filter(!Final_decision == 'delete') %>% 
  select(FeatureID, Final_name, Final_formula, Final_decision) %>% 
  rename(Annotation_level = Final_decision,
         Name = Final_name,
         Formula = Final_formula) %>% 
  mutate(use_Sirius_class = NA)

L2_lower_ready <- L2_lower %>% 
  filter(!Final_decision == 'delete') %>% 
  select(FeatureID, Final_names, Final_formula, Final_decision) %>% 
  rename(Annotation_level = Final_decision,
         Name = Final_names,
         Formula = Final_formula) %>% 
  mutate(use_Sirius_class = NA)

## Data wrangling L3 levels

L3_ready <- rbind(L3, L3_lower) %>% 
  filter(!Final_annotation_comment == 'delete') %>% 
  select(FeatureID, Final_name, Keep_formula, Final_annotation_comment, use_Sirius_class) %>% 
  rename(Annotation_level = Final_annotation_comment,
         Name = Final_name,
         Formula = Keep_formula) %>% 
  distinct()

## Joining in a single table


final_table <- rbind(L1_ready, L2_ready, L2_lower_ready, L3_ready)

write_csv(final_table, "scripts/LC-MS-MS/3.Cleanning_Annotation/output/RP_final_annotation_updated.csv")


## Parsing with checked file

rp_check <- read_csv("scripts/LC-MS-MS/3.Cleanning_Annotation/output/RP_final_annotation_checked.csv")



rp_checkft <- rp_check %>% 
  filter(Annotation_level == 'L4')


final_ready <- final_table %>% 
  mutate(Formula = if_else(FeatureID %in% rp_checkft$FeatureID, 
                           'no formula',
                           Formula),
         Annotation_level = if_else(FeatureID %in% rp_checkft$FeatureID, 
                                    'L5',
                                    Annotation_level),
         use_Sirius_class = if_else(is.na(use_Sirius_class), 'na',
                                    use_Sirius_class)) %>% 
  distinct()

write_csv(final_ready, "scripts/LC-MS-MS/3.Cleanning_Annotation/output/RP_final_annotation_11_20_25.csv")




## Check number of features/ annot level


per_annot <- final_ready %>% 
  count(Annotation_level) %>% 
  mutate(perc = (n/sum(n)*100))

plot_hilic <- per_annot %>% 
  ggplot(aes(x = Annotation_level,
             y = n))+
  geom_col()+
  labs(title = 'RP')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

plot_hilic  
  