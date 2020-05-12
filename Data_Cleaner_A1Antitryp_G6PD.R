### Script to trim and chart relevant gene data from the AI returns
### 12.05.2020

# Clearing Working Environment
rm(list = ls())

# Package Import
library(httr)
library(tidyverse)
library(jsonlite)

# Target disease phenotypes
gluc6phos_dd_geneids <- c('G6PD', 'IKBKG')
alpha1_antitryp_geneids <- c('SERPINA1')

# Setting wd variables
git_wd <- 'C:/Users/zacha/Documents/BigData/'
gluc6phos_dd_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/gluc6phos_dehydr_defic_pheno_data/'
alpha1_antitryp_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/alpha1_antitryp_pheno_data/'
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

#########################################################
### Glucose-6-phosphate Dehydrogenase Gene Data Clean ###
#########################################################

### Initial RDS read in
setwd(gluc6phos_dd_dir)

# Building empty list to append raw dfs to
gluc6phos_dd_gene_dfs <- list()
for (gene in gluc6phos_dd_geneids) {
  gluc6phos_dd_gene_dfs[[gene]] <- list()
}

# Appends raw df where phenotypes match desired
rds_files <- list.files(gluc6phos_dd_dir)
for(i in 1:length(rds_files)) {
  gluc6phos_dd_gene_dfs[[i]] <- filter(df_cleaner(rds_files[i]),
                                       grepl('(g6pd|glucose 6 phosphate)', tolower(description)))
}

### Cleaning and compiling the data

# ALL GENE DATA ON IKBKG IS ASSOCIATED WITH G6PD AND IS THEREFOR DUPLICATE
# Only cleaning dataframe for gluc6phos_dd_gene_dfs[[1]]

# Step 1 - Build first temp df from the index
init_cols <- c('location', 'description', 'Gene',
               'Variation')
temp_df <- gluc6phos_dd_gene_dfs[[1]][, names(gluc6phos_dd_gene_dfs[[1]]) %in% init_cols]

# Step 2 - Build second temp df from the attributes nested in the indexed df
attr_subset <- gluc6phos_dd_gene_dfs[[1]]$attributes
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
gluc6phos_dd_cleaned_data <- final_df_organised # final dim 559, 10
setwd(XCEL_dir) 
write.csv(gluc6phos_dd_cleaned_data, 'gp6d_data.csv') # Data saved to local

######################################################
### Alpha-1 Antitrypsin Deficiency Gene Data Clean ###
######################################################

### Initial RDS read in
setwd(alpha1_antitryp_dir)

# Building empty list to append raw dfs to
alpha1_antitryp_genes_dfs <- list()
for (gene in alpha1_antitryp_geneids) {
  alpha1_antitryp_genes_dfs[[gene]] <- list()
}

# Appends raw df where phenotypes match desired
rds_files <- list.files(alpha1_antitryp_dir)
for (i in 1:length(rds_files)) {
  alpha1_antitryp_genes_dfs[[i]] <- filter(df_cleaner(rds_files[i]),
                                           grepl('(alpha-1-antitrypsin|gallstone)', tolower(description)))
}

### Cleaning and compling data

# Step 1 - Build first temp df from the index
init_cols <- c('location', 'description', 'Gene',
               'Variation')
temp_df <- alpha1_antitryp_genes_dfs[[1]][, names(alpha1_antitryp_genes_dfs[[1]]) %in% init_cols]

# Step 2 - Build second temp df from the attributes nested in the indexed df
attr_subset <- alpha1_antitryp_genes_dfs[[1]]$attributes
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
alpha1_antitryp_cleaned_data <- final_df_organised # final dim 996, 10
setwd(XCEL_dir) 
write.csv(alpha1_antitryp_cleaned_data, 'alpha1_antitrypsin_deficiency_data.csv') # Data saved to local