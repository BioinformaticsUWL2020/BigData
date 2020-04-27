### Script to trim and chart relevant gene data from the API returns
### 12.04.2020 | Happy Easter

# Clearing the working env prior to work
rm(list = ls())

# Package import
library(httr)
library(tidyverse)
library(jsonlite)

# Target disease phenotypes
celiac_geneids <- c('HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'MYO9B', 'CTLA4')
lo_alz_geneids_1 <- c('SCAPER', 'INPPD5', 'EPHA1', 'AK128216', 'STAU2',
                      'ZCWPW1', 'SORL1', 'ATA', 'HLA-DRB5', 'HLA-DRB1',
                      'HDAC9', 'MS4A4A', 'CR1', 'FST', 'MPZL1', 'COL25A1',
                      'UGT1A8', 'UGT1A10', 'AP2A2', 'TMCO4', 'NDUFAF6',
                      'FBXL7', 'COBL', 'HS3ST1', 'PICALM', 'NEGR1', 'NEGR1-IT1',
                      'CD2AP', 'ABCA7', 'BIN1', 'GAB2', 'AHNAK', 'IGH',
                      'CASS4', 'APOE', 'SCIMP', 'GRIN3B', 'TIMP2', 'ABI3',
                      'MS4A6A', 'SELP', 'PCDH7', 'MIR5582', 'F2', 'CKAP5',
                      'TAS2R5', 'CEACAM16', 'BCL3', 'PVRL2', 'TOMM40',
                      'PPP1R37', 'AOX1', 'ADAMTS20', 'HMCN1', 'MTHFD1L',
                      'FERMT2', 'CLU', 'PTK2B', 'OZD4', 'HS3ST1', 'CD33')
lo_alz_geneids_2 <- c('BSG', 'ACE', 'SPPL2A', 'TECTA', 'FAM19A5', 'ECHDC3',
                      'TRPM1', 'MEF2C', 'SELO', 'NME8', 'SQSTM1', 'JPH3', 
                      'SESTD1', 'SEC24B', 'PLCG2', 'SLC24A4', 'RIN3',
                      'LUZP2', 'DLC1', 'BAALC', 'TREML2', 'NEK10', 'TRIP4',
                      'ECHDC3', 'CDH1', 'FAT3', 'ATP5C1', 'GLTSCR1', 'CCZ1B',
                      'SYNGAP1', 'KSR2', 'EGLN3', 'C1orf88', 'IRAK1BP1', 'DNAH6',
                      'SLMAP', 'TBXAS1', 'SUCLG2P4', 'CNTNAP2', 'THSD4', 'CEP63',
                      'RHBDF1', 'TRPM1', 'WWOX', 'MAF', 'ANKRD22', 'GPR180', 
                      'BZRAP1-AS1', 'MECOM', 'TFEB', 'GTF2H3', 'IQGAP2', 'SUZ12P1',
                      'IGSF23', 'TREM2', 'PCNX', 'CUGBP2', 'OARD1', 'ABCA8', 
                      'RASSF5', 'AP2A2', 'NHLRC3', 'PROSER1', 'ATXN7L1',
                      'AKR7A3', 'CELF1', 'ADAM10', 'NDUFAF6', 'CACNA2D3', 'HMGA2', 
                      'ISYNA1', 'ELL', 'SSBP4', 'MIR142', 'TSPOAP1-AS1', 'SLC10A2',
                      'CDON')

# setting wd variables
git_wd <- 'C:/Users/zacha/Documents/BigData/'
lo_alz_1_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/lo_alz_1_pheno_data/'
lo_alz_2_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/lo_alz_2_pheno_data/'
celiac_dir <- 'C:/Users/zacha/Documents/BigData/ProjectAPISaves/celiac_pheno_data/'
XCEL_dir <- 'H:/BigData/'

# Function to  clean data frames 
df_cleaner <- function(dataset) {
  # Read rds files from wd
  print(dataset)
  readin_df <- jsonlite::fromJSON(
    content(readRDS(dataset), as = 'text')
  ) 
  return(readin_df)
}

# Final df colnames for both phenotypes
cleaned_col_names <- c('Database ID', 'Associated Gene', 'Location',
                       'Description', 'Variation', 'Risk Allele', 'P - Value',
                       'Odds Ratio', 'Beta Coefficient', 'Clinical Significance')

#########################################
### Celiac Associated Gene Data Clean ###  # Data Successfully Cleaned 
#########################################

# Setting wd to clean celiac gene data
setwd(celiac_dir)


# Starting a list to append cleaned dataframes of corresponding gene to 
celiac_genes_dfs <- list()
for (gene in celiac_geneids) {
  celiac_genes_dfs[[gene]] <- list()
}

# Creates an iterable list of .rds files to read in 
rds_files <- list.files(celiac_dir)
for (i in 1:length(rds_files)) {
  celiac_genes_dfs[[i]] <- filter(df_cleaner(rds_files[i]), 
                                  grepl('celiac', tolower(description)))
}

# Opening a df with a dead row to bind the first iteration return from the loop below
cleaned_celiac_data <- data.frame(t(rep('REMOVE', 10)))
names(cleaned_celiac_data) <- cleaned_col_names

# Proper implementation of the code above to build one cleaned df
# from all the dfs created from the RDS files
for (i in 1:length(celiac_genes_dfs)) {

  # Step 1 - Build first temp df from the index
  init_cols <- c('location', 'description', 'Gene',
                 'Variation')
  temp_df <- celiac_genes_dfs[[i]][,names(celiac_genes_dfs[[i]]) %in% init_cols]
  
  # Step 2 - Build second temp df from the attributes df nested in the indexed df
  attr_subset <- celiac_genes_dfs[[i]]$attributes
  attr_cols_removed <- c('external_reference', 'external_id', 'MIM')
  if ('clinical_significance' %in% names(attr_subset)){
    attr_clean <- attr_subset[,!names(attr_subset) %in% attr_cols_removed]
  } else {
    clinical_significance <- rep(NA, length(temp_df$description))
    attr_clean <- attr_subset[,!names(attr_subset) %in% attr_cols_removed]
    attr_clean <- cbind(attr_clean, clinical_significance)
  }
  
  # Step 3 - Merge the two df into one 
  final_temp_df <- data.frame(c(temp_df, attr_clean))
  
  # Step 4 - Organise all the columns into uniform postions and rename them
  final_df_organised <- final_temp_df[,c("Gene", "associated_gene", "location", "description",
                                         "Variation", "risk_allele", "p_value", "odds_ratio",
                                         "beta_coefficient", "clinical_significance")]
  names(final_df_organised) <- cleaned_col_names
  
  # Step 5 - Bind all the newly cleaned df into one
  cleaned_celiac_data <- rbind(cleaned_celiac_data, final_df_organised) # final dim 237, 10
}

# Step 6 - Removing dead row (row 1)
cleaned_celiac_data <- cleaned_celiac_data[2:237,]

setwd(XCEL_dir)
write.csv(cleaned_celiac_data, 'celiac_data.xlsx') # Data saved to local

#########################################################
### Late-Onset Alzheimer's Associated Gene Data Clean ###
#########################################################

lo_alz_1_genes_dfs <- list()
for (gene in lo_alz_geneids_1) {
  lo_alz_1_genes_dfs[[gene]] <- list()
}

# Setting wd to clean loalz_1 gene data
setwd(lo_alz_1_dir)
rds_files <- list.files(lo_alz_1_dir)
for (i in 1:length(lo_alz_1_genes_dfs)) {
  lo_alz_df <- df_cleaner(rds_files[i])
  if (class(lo_alz_df) == 'data.frame') {
    lo_alz_df <- filter(lo_alz_df,
                        grepl('(late-onset|alzheimer\'s)', tolower(description)))
    if (class(lo_alz_df) != 'list') {
      lo_alz_1_genes_dfs[[i]] <- lo_alz_df
    }
  } else {
    next
  }
}

for (index in lo_alz_1_genes_dfs) {
  if (length(lo_alz_1_genes_dfs[[index]]) == 0){
    remove(lo_alz_1_genes_dfs[[index]])
  }
}













lo_alz_2_gene_dfs <- list()
for (gene in lo_alz_geneids_2) {
  lo_alz_2_gene_dfs[[gene]] <- list()
}

setwd(lo_alz_2_dir)
rds_files <- list.files(lo_alz_2_dir)
for (i in 1:length(lo_alz_2_gene_dfs)) {
  lo_alz_df <- df_cleaner(rds_files[i])
  if (class(lo_alz_df) == 'data.frame') {
    lo_alz_2_gene_dfs[[i]] <- filter(lo_alz_df,
                                     grepl('(late-onset|alzheimer\'s)', tolower(description)))
  } else {
    next
  }
}

























### TEST CODE ###
step_1 <- c('location', 'description', 'Gene', 
            'Variation', 'SupportingStructuralVariation',
            'StructuralVariation') 
test_df <- celiac_genes_dfs[[1]][,names(celiac_genes_dfs[[1]]) %in% step_1]

attr_subset <- celiac_genes_dfs[[1]]$attributes
step_2 <- c('external_reference', 'external_id', 'MIM')
attr_clean <- attr_subset[,!names(attr_subset) %in% step_2]

final_test_df_1 <- data.frame(c(test_df, attr_clean))

step_1 <- c('location', 'description', 'Gene', 
            'Variation', 'SupportingStructuralVariation',
            'StructuralVariation') 
test_df <- celiac_genes_dfs[[2]][,names(celiac_genes_dfs[[2]]) %in% step_1]
########## second df for bind testing ################
attr_subset <- celiac_genes_dfs[[2]]$attributes
step_2 <- c('external_reference', 'external_id', 'MIM')
if (isTRUE(attr_subset$clinical_significance)) {
  attr_clean <- attr_subset[,!names(attr_subset) %in% step_2]
  final_test_df_2 <- data.frame(c(test_df, attr_clean))
} else {
  clinical_significance <- rep(NA, length(test_df$description))
  print(length(clinical_significance))
  attr_clean <- attr_subset[,!names(attr_subset) %in% step_2]
  attr_clean_clinSig <- cbind(attr_clean, clinical_significance)
  final_test_df_2 <- data.frame(test_df, attr_clean_clinSig)
}
test <- rbind(final_test_df_1, final_test_df_2) # this works

# Testing moving columns
test_organised_cols <- test[,c("Gene", "associated_gene", "location",
                               "description", "Variation", "risk_allele",
                               "p_value", "odds_ratio", "beta_coefficient",
                               "clinical_significance",
                               "StructuralVariation", "SupportingStructuralVariation")]

### TEST CODE ###