### openSNP API Construct to Pull Genotype Data of Users with Target Phenotype Associated SNPs
### 26.05.2020

# Library import
library(httr)
library(jsonlite)
library(xml2)
library(fs)
library(tidyverse)
library(curl)
library(rlist)

#####################
##### IMPORTANT #####
#####################

# This script uses the output of openSNP_TargetPhenoUsers.R
# Be sure to run that script first before running this

# Mildly obnoxious disclaimer... That script will take several hours to run

#####################
##### IMPORTANT #####
#####################

# Building empty df to append to avoiding overwriting data
rsID <- as.character('Remove')
Chromosome <- as.character('Remove')
Position <- as.character('Remove')
Genotype <- as.character('Remove')
UserID <- as.integer(000)

binding_df <- tibble(rsID, Chromosome, Position, Genotype, UserID)

### Step 5 - Run curl() on the genotype data returned for each user data response from Step 4

### a1a users genotype table build

# a1a sample geno extraction
a1a_sample_genotypes <- openSNP_users %>%
  filter(id %in% a1a_user_sample)

# List to store genotype data frames
a1a_pos_user_genotypes <- list()

for (i in 100:length(a1a_sample_genotypes$id)) {
  
  # checks for multiple download urls
  for (x in 1:length(a1a_sample_genotypes[[3]][[i]]$download_url)) {
    a1a_pos_user_genotypes <- readLines(curl(a1a_sample_genotypes[[3]][[i]]$download_url[x]))
    a1a_pos_user_genotypes <- data.frame(a1a_pos_user_genotypes)
    a1a_pos_user_genotypes <- a1a_pos_user_genotypes %>%
      filter(!(grepl('#', a1a_pos_user_genotypes))) %>%
      separate(a1a_pos_user_genotypes, 
               into = c('rsID', 'Chromosome', 'Position', 'Genotype'),
               sep = '\t')
    
    # adding the col with the current user id
    user_id_col_name <- 'UserID'
    id_vector <- as.data.frame(rep(a1a_sample_genotypes[[2]][[i]], length(a1a_pos_user_genotypes$rsID)))
    names(id_vector) <- user_id_col_name
    a1a_pos_user_genotypes <- data.frame(a1a_pos_user_genotypes, id_vector)
    print(length(a1a_pos_user_genotypes$rsID))
    
    # binding all genotype tables together 
    binding_df <- rbind(binding_df, a1a_pos_user_genotypes)
    
  }
  
  # storing genotype df into list 
  a1a_pos_user_genotypes[[a1a_sample_genotypes[[2]][[i]]]] <- list(binding_df)
  
}

binding_df <- tibble(rsID, Chromosome, Position, Genotype, UserID)


# g6pd sample geno extraction
g6pd_sample_genotypes <- openSNP_users %>% 
  filter(id %in% g6pd_user_sample)

# f11 sample geno extraction
f11_sample_genotypes <- openSNP_users %>% 
  filter(id %in% f11_user_sample)

# eo dyst sample geno extraction
eo_dyst_sample_genotypes <- openSNP_users %>%
  filter(id %in% eo_dyst_user_sample)

# celiac sample geno extraction
celiac_sample_genotypes <- openSNP_users %>%
  filter(id %in% celiac_user_sample)

# hhem sample geno extraction
hhem_sample_genotypes <- openSNP_users %>%
  filter(id %in% hhem_user_sample)

# hthromb sample geno extraction
hthromb_sample_genotypes <- openSNP_users %>%
  filter(id %in% hthromb_user_sample)

# pd sample geno extraction
pd_sample_genotypes <- openSNP_users %>%
  filter(id %in% pd_user_sample)

# lo alz sample geno extraction
lo_alz_sample_genotypes <- openSNP_users %>% 
  filter(id %in% lo_alz_user_sample)


### Step 6 - Create column with just user id to span the genotype data pulled

### Step 7 - Bind all genotype data for each 100 users with a specific disease together



