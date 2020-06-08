### openSNP API Construct to Pull Genotype Data of Users with Target Phenotype Associated SNPs
### 26.05.2020

# Clear wd
rm(list = ls())

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

# 08.06.2020 - Note | The above disclaimer and need for the openSNP_TargetPhenoUsers.R is no longer necessary
# Make sure the ID_samples_TextFiles folder is pulled from the GitHub and setw it to your wd
# This script now imports the ids from those txt files

#####################
##### IMPORTANT #####
#####################

# PATHs 
top_level_path <- 'C:/Users/zacha/Documents/BigData/' # Change to match your file structure
ext_path <- 'ID_samples_TextFiles/'

# Seting ID vectors to subset openSNP_users
setwd(paste0(top_level_path, ext_path))

pd_user_sample <- sample(
  as.vector(scan('pd_pos_userids_sample.txt',character(), sep = ' ')),
  50
)

lo_alz_user_sample <- sample(
  as.vector(scan('lo_alz_pos_userids_sample.txt', character(), sep = '')),
  50
)

# API to bring in OpenSNP total user table
openSNP_users <- GET('http://opensnp.org/users.json',
                     content_type('application/json'),
                     user_agent('zachary.h.sykes@gmail.com -- Genetic risk prediction project'))
openSNP_users <- fromJSON(
  content(openSNP_users, as = 'text')
)


# Building empty df to append to avoiding overwriting data
rsID <- as.character('Remove')
Chromosome <- as.character('Remove')
Position <- as.character('Remove')
Genotype <- as.character('Remove')
UserID <- as.integer(000)

binding_df <- tibble(rsID, Chromosome, Position, Genotype, UserID)

### Step 5 - Run curl() on the genotype data returned for each user data response from Step 4

### PD users genotype table build

# pd sample geno extraction
pd_sample_genotypes <- openSNP_users %>%
  filter(id %in% pd_user_sample)

# List to store genotype data frames
pd_pos_user_genotypes_list <- list()

for (i in 1:length(pd_sample_genotypes$id)) {
  
  # checks for multiple download urls
  for (x in 1:length(pd_sample_genotypes[[3]][[i]]$download_url)) {
    pd_pos_user_genotypes <- readLines(curl(pd_sample_genotypes[[3]][[i]]$download_url[x]))
    pd_pos_user_genotypes <- data.frame(pd_pos_user_genotypes)
    pd_pos_user_genotypes <- pd_pos_user_genotypes %>%
      filter(!(grepl('#', pd_pos_user_genotypes))) %>%
      separate(pd_pos_user_genotypes, 
               into = c('rsID', 'Chromosome', 'Position', 'Genotype'),
               sep = '\t')
    
    ### Step 6 - Create column with just user id to span the genotype data pulled
    user_id_col_name <- 'UserID'
    id_vector <- as.data.frame(rep(pd_sample_genotypes[[2]][[i]], length(pd_pos_user_genotypes$rsID)))
    names(id_vector) <- user_id_col_name
    pd_pos_user_genotypes <- data.frame(pd_pos_user_genotypes, id_vector)
    print(length(pd_pos_user_genotypes$rsID))
    
    # binding all genotype tables together 
    binding_df <- rbind(binding_df, pd_pos_user_genotypes)
    
  }
  
  # storing genotype df into list 
  pd_pos_user_genotypes_list[[pd_sample_genotypes[[2]][[i]]]] <- list(binding_df)
  
}

binding_df <- tibble(rsID, Chromosome, Position, Genotype, UserID)

### LOAD users genotype table build

# lo alz sample geno extraction
lo_alz_sample_genotypes <- openSNP_users %>% 
  filter(id %in% lo_alz_user_sample)

# List to store genotype data frames
lo_alz_pos_user_genotypes_list <- list()

for (i in 1:length(lo_alz_sample_genotypes$id)) {
  
  # checks for multiple download urls
  for (x in 1:length(lo_alz_sample_genotypes[[3]][[i]]$download_url)) {
    lo_alz_pos_user_genotypes <- readLines(curl(lo_alz_sample_genotypes[[3]][[i]]$download_url[x]))
    lo_alz_pos_user_genotypes <- data.frame(lo_alz_pos_user_genotypes)
    lo_alz_pos_user_genotypes <- lo_alz_pos_user_genotypes %>%
      filter(!(grepl('#', lo_alz_pos_user_genotypes))) %>%
      separate(lo_alz_pos_user_genotypes, 
               into = c('rsID', 'Chromosome', 'Position', 'Genotype'),
               sep = '\t')
    
    ### Step 6 - Create column with just user id to span the genotype data pulled
    user_id_col_name <- 'UserID'
    id_vector <- as.data.frame(rep(lo_alz_sample_genotypes[[2]][[i]], length(lo_alz_pos_user_genotypes$rsID)))
    names(id_vector) <- user_id_col_name
    lo_alz_pos_user_genotypes <- data.frame(lo_alz_pos_user_genotypes, id_vector)
    print(length(lo_alz_pos_user_genotypes$rsID))
    
    # binding all genotype tables together 
    binding_df <- rbind(binding_df, lo_alz_pos_user_genotypes)
    
  }
  
  # storing genotype df into list 
  lo_alz_pos_user_genotypes_list[[lo_alz_sample_genotypes[[2]][[i]]]] <- list(binding_df)
  
}

binding_df <- tibble(rsID, Chromosome, Position, Genotype, UserID)