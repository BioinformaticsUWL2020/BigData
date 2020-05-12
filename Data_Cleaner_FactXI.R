### Script to trim and chart relevant gene data from the API returns
### 12.05.2020

# Clearing Working Environment
rm(list = ls())

# Package Import
library(httr)
library(tidyverse)
library(jsonlite)

# Target disease phenotypes
fact_xi_geneids <- 'F11'

# Setting wd variables
git_wd <- 'C:/Users/zacha/Documents/BigData/'
fact_xi_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/fact_xi_pheno_data/'
XCEL_dir <- 'H:/BigData/'

# Function to clean the data frames
df_cleaner <- function(dataset) {
  # read RDS files from wd
  readin_df <- jsonlite::fromJSON(
    content(readRDS(dataset), as = 'text')
  )
  return(readin_df)
}

# Final df colnames for both phenotypes
cleaned_col_names <- c('Database ID', 'Associated Gene', 'Location',
                       'Description', 'Variation', 'Risk Allele', 'P - Value',
                       'Odds Ratio', 'Beta Coefficient', 'Clinical Significance')

##############################################
### Factor XI Deficiency Gene Data Cleaner ###
##############################################

### Initial RDS read in
setwd(fact_xi_dir)

# Building empty list to append raw dfs to
fact_xi_gene_dfs <- list()
for (gene in fact_xi_geneids) {
  fact_xi_gene_dfs[[gene]] <- list()
}

# Appends raw df where phenotypes match desired
rds_files <- list.files(fact_xi_dir)
for(i in 1:length(rds_files)) {
  fact_xi_gene_dfs[[i]] <- filter(df_cleaner(rds_files[i]),
                                  grepl('factor xi', tolower(description)))
}

### Cleaning and compiling the data

# Step 1 - Build first temp df from the index
init_cols <- c('location', 'description', 'Gene',
               'Variation')
temp_df <- fact_xi_gene_dfs[[1]][, names(fact_xi_gene_dfs[[1]]) %in% init_cols]

# Step 2 - Build second temp df from the attributes nested in the indexed df
attr_subset <- fact_xi_gene_dfs[[1]]$attributes
attr_cols_removed <- c('external_reference', 'external_id', 'MIM')
attr_clean <- attr_subset[, !names(attr_subset) %in% attr_cols_removed]

# Step 3 - Merge the two together
final_temp_df <- data.frame(c(temp_df, attr_clean))

# Step 4 - Organise all the columns into uniform postions and rename them
final_df_organised <- final_temp_df[,c("Gene", "associated_gene", "location", "description",
                                       "Variation", "risk_allele", "p_value", "odds_ratio",
                                       "beta_coefficient", "clinical_significance")]
names(final_df_organised) <- cleaned_col_names

# Step 5 - Rename and write to local
fact_xi_cleaned_data <- final_df_organised # final dim 559, 10
setwd(XCEL_dir) 
write.csv(fact_xi_cleaned_data, 'factor_xi_deficiency_data.csv') # Data saved to local

