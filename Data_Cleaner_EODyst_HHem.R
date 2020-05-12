### Script to Clean and Chart relevant genids from the API returns
### 12.05.2020

# Clearing the Working Directory
rm(list = ls())

# Package Import
library(httr)
library(tidyverse)
library(jsonlite)

# Target Disease Phenotypes
h_hem_geneids <- c('HFE', 'SLC40A1', 'TFR2', 'HAMP')
eo_dyst_geneids <- c('TOR1A')

# Setting wd variables
git_wd <- 'C:/Users/zacha/Documents/BigData'
eo_dyst_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/eo_prim_dyst_pheno_data/'
h_hem_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/h_hem_pheno_data/'
XCEL_dir <- 'H:/BigData'

# Function to Clean Dataframes
df_cleaner <- function(dataset) {
  # Read RDS files from wd
  readin_df <- jsonlite::fromJSON(
    content(readRDS(dataset), as = 'text')
  )
  return(readin_df)
}

# Final df colnames for both disease phenotype
cleaned_col_names <- c('Database ID', 'Associated Gene', 'Location',
                       'Description', 'Variation', 'Risk Allele', 'P - Value',
                       'Odds Ration', 'Beta Coefficient', 'Clinical Significance')

###############################################################
### Early Onset Primary Dystonia Associated Gene Data Clean ###
###############################################################

### Initial list of dataframes
setwd(eo_dyst_dir)

# Starting a list to append raw dfs of corresponding geneid to
eo_dyst_genes_dfs <- list()
for (gene in eo_dyst_geneids) {
  eo_dyst_genes_dfs[[gene]] <- list()
}

# Creates opens and appends rds files with API resp to the empty list made above
rds_files <- list.files(eo_dyst_dir)
for (i in 1:length(rds_files)) {
  eo_dyst_genes_dfs[[i]] <- filter(df_cleaner(rds_files[i]),
                                   grepl('dystonia', tolower(description)) & !grepl('(late-onset|torsion)', tolower(description)))
}

### Cleaning and Compiling the Data 
# Opening a df with a dead row to bind the first iteration of the loop below to
cleaned_eo_dyst_data <- data.frame(t(rep('REMOVE', 10)))
names(cleaned_eo_dyst_data) <- cleaned_col_names

# Building one cleaned df from all the dfs in the list created above (only one for eo_dyst)
for (i in 1:length(eo_dyst_genes_dfs)) {
  
  # Step 1 - Build first temp df from the index
  init_cols <- c('location', 'description', 'Gene',
                 'Variation')
  temp_df <- eo_dyst_genes_dfs[[i]][, names(eo_dyst_genes_dfs[[i]]) %in% init_cols]
  
  # Step 2 - Build second temp df from nested attributes df in the index
  attr_subset <- eo_dyst_genes_dfs[[i]]$attributes
  attr_cols_removed <- c('external_reference', 'MIM', 'external_id')
  if (!('beta_coefficient' %in% names(attr_subset))) {
    beta_coefficient <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, beta_coefficient)
  }
  
  attr_clean <- attr_subset[, !names(attr_subset) %in% attr_cols_removed]
  
  # Step 3 - Merge the two temp dfs into one
  final_temp_df <- data.frame(temp_df, attr_clean)
  
  # Step 4 - Organise all the columns into uniform positions and rename them
  final_temp_df_organised <- final_temp_df[, c("Gene", "associated_gene", "location", "description",
                                               "Variation", "risk_allele", "p_value", "odds_ratio", 
                                               "beta_coefficient", "clinical_significance")]
  names(final_temp_df_organised) <- cleaned_col_names
  
  # Step 5 -Bind all newly cleaned df into one
  cleaned_eo_dyst_data <- rbind(cleaned_eo_dyst_data, final_temp_df_organised) # final dim 100, 10
}

# Step 6 - Remove the dead row and write to a CSV
cleaned_eo_dyst_data <- cleaned_eo_dyst_data[2:101,]
setwd(XCEL_dir)
write.csv(cleaned_eo_dyst_data, 'eo_dystonia_data.csv') # Data saved to local

#############################################################
### Hereditary Hemochromatosis Associated Gene Data Clean ###
#############################################################

### Initial list of dataframes
setwd(h_hem_dir)

# Starting a list to append raw dfs of corresponding geneid to
h_hem_gene_dfs <- list()
for (gene in h_hem_geneids) {
  h_hem_gene_dfs[[gene]] <- list()
}

# Creates opens and appends rds files with API resp to the empty list made above
rds_files <- list.files(h_hem_dir) 
for (i in 1:length(rds_files)) {
  h_hem_gene_dfs[[i]] <- filter(df_cleaner(rds_files[i]),
                                grepl('(hereditary hemochromatosis|arthritis)', tolower(description)))
}


### Cleaning and Compiling the Data 
# Opening a df with a dead row to bind the first iteration of the loop below to
cleaned_h_hem_data <- data.frame(t(rep('REMOVE', 10)))
names(cleaned_h_hem_data) <- cleaned_col_names

# Building one cleaned df from all the dfs in the list created above
for (i in 1:length(h_hem_gene_dfs)) {
  
  # Step 1 - Build first temp df from the index
  init_cols <- c('location', 'description', 'Gene',
                 'Variation')
  temp_df <- h_hem_gene_dfs[[i]][, names(h_hem_gene_dfs[[i]]) %in% init_cols]
  
  # Step 2 - build second temp df from the nested attr df
  attr_subset <- h_hem_gene_dfs[[i]]$attributes
  attr_cols_removed <- c('external_reference', 'MIM', 'external_id')
  if (!('odds_ratio' %in% names(attr_subset))) {
    odds_ratio <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, odds_ratio)
  }
  attr_clean <- attr_subset[, !names(attr_subset) %in% attr_cols_removed]
  
  # Step 3 - Merge two temporary dfs into one
  final_temp_df <- data.frame(temp_df, attr_clean)
  
  # Step 4 - Organise all the columns into uniform positions and rename them
  final_temp_df_organised <- final_temp_df[, c("Gene", "associated_gene", "location", "description",
                                               "Variation", "risk_allele", "p_value", "odds_ratio",
                                               "beta_coefficient", "clinical_significance")]
  names(final_temp_df_organised) <- cleaned_col_names
  
  # Step 5 - Bind all the newly cleaned df into one
  cleaned_h_hem_data <- rbind(cleaned_h_hem_data, final_temp_df_organised) # final dim 379, 10
}

# Step 6 - Remove dead row and write to local
cleaned_h_hem_data <-cleaned_h_hem_data[2:379, ]
setwd(XCEL_dir)
write.csv(cleaned_h_hem_data, 'hereditary_hemochromatosis_data') # Data saved to local