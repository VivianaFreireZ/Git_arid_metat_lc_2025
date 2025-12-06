
library(tidyverse)

# Loading cleaned annotation files


L1 <- read_csv('scripts/LC-MS-MS/3.Cleanning_Annotation/input/manually_curated_tables/HILIC_level_1_mzvault.csv')


L2 <- read_csv('scripts/LC-MS-MS/3.Cleanning_Annotation/input/manually_curated_tables/HILIC_level_2_flat_mzCloud_ready.csv')


L2_lower <- read_csv('scripts/LC-MS-MS/3.Cleanning_Annotation/input/manually_curated_tables/HILIC_level_2_lower_mzCloud_ready.csv')


L3 <- read_csv('scripts/LC-MS-MS/3.Cleanning_Annotation/input/manually_curated_tables/HILIC_level_3_CD_Sirius.csv')


## Data wranglig

L1_ready <- L1 %>% 
  select(FeatureID, Name, Formula) %>% 
  mutate(Annotation_level = "L1",
         use_Sirius_class = NA)

L2_ready <- L2 %>% 
  filter(!Final_decision == 'delete') %>% 
  filter(!Final_decision == 'DELETE') %>% 
  select(FeatureID, Name, Formula, Final_decision) %>% 
  rename(Annotation_level = Final_decision) %>% 
  mutate(use_Sirius_class = NA)

L2_lower_ready <- L2_lower %>% 
  filter(!Final_decision == 'DELETE') %>% 
  select(FeatureID, Name, Formula, Final_decision) %>% 
  rename(Annotation_level = Final_decision) %>% 
  mutate(use_Sirius_class = NA)

## Data wrangling L3 levels

L3_ready <- L3 %>% 
  filter(!Final_annotation_comment == 'delete') %>% 
  select(FeatureID, Final_name, Keep_formula, Final_annotation_comment, use_Sirius_class) %>% 
  rename(Annotation_level = Final_annotation_comment,
         Name = Final_name,
         Formula = Keep_formula) %>% 
  distinct()

  

### Sanity check - check if formulas to keep are identical to the source 
  
all(L3_ready$Formula_sirius == L3_ready$Keep_formula)

## TRUE


## Joining in a single table


final_table <- rbind(L1_ready, L2_ready, L2_lower_ready, L3_ready)

write_csv(final_table, "scripts/LC-MS-MS/3.Cleanning_Annotation/output/HILIC_final_annotation_updated.csv")


## Parsing with checked file

hilic_check <- read_csv("scripts/LC-MS-MS/3.Cleanning_Annotation/output/HILIC_final_annotation_checked.csv")



hilic_checkft <- hilic_check %>% 
  filter(Annotation_level == 'L4',
         Name == 'no name')


final_ready <- final_table %>% 
  mutate(Formula = if_else(FeatureID %in% hilic_checkft$FeatureID, 
                            'no formula',
                            Formula),
         Annotation_level = if_else(FeatureID %in% hilic_checkft$FeatureID, 
                           'L5',
                           Annotation_level),
         use_Sirius_class = if_else(is.na(use_Sirius_class), 'na',
                                    use_Sirius_class),
         use_Sirius_class = str_replace(use_Sirius_class, 'YES', 'yes'))

write_csv(final_ready, "scripts/LC-MS-MS/3.Cleanning_Annotation/output/HILIC_final_annotation_11_20_25.csv")


## Check number of features/ annot level


per_annot <- final_ready %>% 
  count(Annotation_level) %>% 
  mutate(perc = (n/sum(n)*100))

plot_hilic <- per_annot %>% 
  ggplot(aes(x = Annotation_level,
             y = n))+
  geom_col()+
  labs(title = 'HILIC')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

plot_hilic  
  