Date: 1/9/2023

Checking HILIC spectra, I found a duplicated feature (see below)

compounds_table <- cd_results_table %>% 
  select(FeatureID, contains('Area:')) %>% 
  mutate(duplicate = duplicated(FeatureID)) %>% 
  filter(duplicate == "TRUE")
  
### Note: I found a duplicated feature * FT_1.664_349.09316 *

### Similar intensities, different adduct and different predicted formula
### These are duplicate features, delete the one with lower annotation error
### Manually editing the LC-MS report and deleting the feature
### Readme created

Deleting the feature with reference ion [M-H]-1