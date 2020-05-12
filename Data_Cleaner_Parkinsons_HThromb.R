### Script to Clean and Chart relevant genids from the API returns
### 12.05.2020

# Clearing the Working Directory
rm(list = ls())

# Package Import
library(httr)
library(tidyverse)
library(jsonlite)

# Target Disease Phenotypes
pd_geneids <- c('GLUD2', 'TBP', 'LRRK2', 'MAPT', 'SNCA', 'ATP13A2',
                'VPS35', 'SYNJ1', 'UCHL1', 'SNCAIP', 'NR4A2', 'CHCHD2',
                'PLA2G6', 'GBA', 'DNAJC6', 'ATXN2', 'ADH1C')
h_thromb_geneids <- c('F9', 'HRG', 'PROC', 'PROS1', 'SERPINC1')

# Setting wd variables
git_wd <- 'C:/Users/zacha/Documents/BigData'
pd_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/pd_pheno_data/'
h_thromb_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/h_thromb_pheno_data/'
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

#####################################################
### Parkinsons Disease Associated Gene Data Clean ###
#####################################################

### Initial RDS read in 
setwd(pd_dir)

# Starting a list to append raw dfs of corresponding geneid to
pd_genes_dfs <- list()
for (gene in pd_geneids) {
  pd_genes_dfs[[gene]] <- list()
}

# Creates opens and appends rds files with API resp to the empty list made above
rds_files <- list.files(pd_dir)
for (i in 1:length(rds_files)) {
  pd_genes_dfs[[i]] <- filter(df_cleaner(rds_files[i]),
                              grepl('parkinson', tolower(description)))
}

### Cleaning and Compiling the Data 
# Opening a df with a dead row to bind the first iteration of the loop below to
cleaned_pd_data <- data.frame(t(rep('REMOVE', 10)))
names(cleaned_pd_data) <- cleaned_col_names

# Building one cleaned df from all the dfs in the list created above (only one for eo_dyst)
for (i in 1:length(pd_genes_dfs)) {
  
  # Step 1 - Build first temp df from the index
  init_cols <- c('location', 'description', 'Gene',
                 'Variation')
  temp_df <- pd_genes_dfs[[i]][, names(pd_genes_dfs[[i]]) %in% init_cols]
  
  # Step 2 - Build second temp df from nested attributes df in the index
  attr_subset <- pd_genes_dfs[[i]]$attributes
  attr_cols_removed <- c('external_reference', 'MIM', 'external_id')
  
  if (!('beta_coefficient' %in% names(attr_subset))) {
    beta_coefficient <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, beta_coefficient)
  }
  
  if (!('clinical_significance' %in% names(attr_subset))) {
    clinical_significance <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, clinical_significance)
  }
  
  if (!('odds_ratio' %in% names(attr_subset))) {
    odds_ratio <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, odds_ratio)
  }
  
  if (!('p_value' %in% names(attr_subset))) {
    p_value <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, p_value)
  }
  # next if statement if necessary
  
  attr_clean <- attr_subset[, !names(attr_subset) %in% attr_cols_removed]
  
  # Step 3 - Merge the two temp dfs into one
  final_temp_df <- data.frame(temp_df, attr_clean)
  
  # Step 4 - Organise all the columns into uniform positions and rename them
  final_temp_df_organised <- final_temp_df[, c("Gene", "associated_gene", "location", "description",
                                               "Variation", "risk_allele", "p_value", "odds_ratio", 
                                               "beta_coefficient", "clinical_significance")]
  names(final_temp_df_organised) <- cleaned_col_names
  
  # Step 5 -Bind all newly cleaned df into one
  cleaned_pd_data <- rbind(cleaned_pd_data, final_temp_df_organised) # final dim 1838 , 10
}

# Step 6 - Remove dead row and write to a CSV
cleaned_pd_data <- cleaned_pd_data[2:1838, ]
setwd(XCEL_dir)
write.csv(cleaned_pd_data, 'parkinsons_data.csv') # Data saved to local

###########################################################
### Hereditary Thrombophilia Associated Gene Data Clean ###
###########################################################

### Initial list of dataframes
setwd(h_thromb_dir)

# Starting a list to append raw dfs of corresponding geneid to
h_thromb_gene_dfs <- list()
for (gene in h_thromb_geneids) {
  h_thromb_gene_dfs[[gene]] <- list()
}

# Creates opens and appends rds files with API resp to the empty list made above
rds_files <- list.files(h_thromb_dir) 
for (i in 1:length(rds_files)) {
  h_thromb_gene_dfs[[i]] <- filter(df_cleaner(rds_files[i]),
                                   grepl('hereditary thrombophilia|thrombophilia|hemophilia)', tolower(description)))
}

### Cleaning and Compiling the Data 
# Opening a df with a dead row to bind the first iteration of the loop below to
cleaned_h_thromb_data <- data.frame(t(rep('REMOVE', 10)))
names(cleaned_h_thromb_data) <- cleaned_col_names

# Building one cleaned df from all the dfs in the list created above (only one for eo_dyst)
for (i in 1:length(h_thromb_gene_dfs)) {
  
  # Step 1 - Build first temp df from the index
  init_cols <- c('location', 'description', 'Gene',
                 'Variation')
  temp_df <- h_thromb_gene_dfs[[i]][, names(h_thromb_gene_dfs[[i]]) %in% init_cols]
  
  print(dim(temp_df))
  
  # Step 2 - Build second temp df from nested attributes df in the index
  attr_subset <- h_thromb_gene_dfs[[i]]$attributes
  attr_cols_removed <- c('external_reference', 'MIM', 'external_id')
  
  if (!('beta_coefficient' %in% names(attr_subset))) {
    beta_coefficient <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, beta_coefficient)
  }
  
  if (!('clinical_significance' %in% names(attr_subset))) {
    clinical_significance <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, clinical_significance)
  }
  
  if (!('odds_ratio' %in% names(attr_subset))) {
    odds_ratio <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, odds_ratio)
  }
  
  if (!('p_value' %in% names(attr_subset))) {
    p_value <- rep(NA, length(temp_df$description))
    attr_subset <- cbind(attr_subset, p_value)
  }
  # next if statement if necessary
  
  attr_clean <- attr_subset[, !names(attr_subset) %in% attr_cols_removed]
  
  # Step 3 - Merge the two temp dfs into one
  final_temp_df <- data.frame(temp_df, attr_clean)
  
  # Step 4 - Organise all the columns into uniform positions and rename them
  final_temp_df_organised <- final_temp_df[, c("Gene", "associated_gene", "location", "description",
                                               "Variation", "risk_allele", "p_value", "odds_ratio", 
                                               "beta_coefficient", "clinical_significance")]
  names(final_temp_df_organised) <- cleaned_col_names
  
  # Step 5 -Bind all newly cleaned df into one
  cleaned_h_thromb_data <- rbind(cleaned_h_thromb_data, final_temp_df_organised) # final dim 560 , 10
}

# Step 6 - Remove dead row and write to local
cleaned_h_thromb_data <-cleaned_h_thromb_data[2:560, ]
setwd(XCEL_dir)
write.csv(cleaned_h_thromb_data, 'hereditary_thrombophilia_data') # Data saved to local