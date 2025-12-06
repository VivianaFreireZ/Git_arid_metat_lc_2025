Date: 1-8-2023

## Important note ##

# Running the analysis in R with RP results, I found there is a duplicated Feature ID that's the reason I got an error in the
# previous line of code (column_rownames(var = FeatureID).

## Finding duplicate feature
  
temp <- compounds_table %>% 
  mutate(duplicate = duplicated(FeatureID)) %>% 
  filter(duplicate == "TRUE")

## The duplicated feature is FT_8.445_223.05767
## I manually check the intensities and CD annotation of these features and
## they are identical. 

## Solution: manually delete one

The file within this folder, CD_report_modified, includes a copy of the RP raw report where the duplicated feature was manually deleted

