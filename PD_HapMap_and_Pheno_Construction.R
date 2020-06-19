### Preparing HapMap and Phenotype Data for PD and Control
### 09.06.2020

# Lib import
library(tidyverse)
library(curl)
library(httr)
library(jsonlite)
library(sparklyr)
library(hablar)
library(plyr)

# PATHs 
top_level_path <- 'C:/Users/zacha/Documents/BigData/' # Change to match your file structure
ext1_path      <- 'ID_samples_TextFiles/'
ext2_path      <- 'HapMap_and_Pheno_Files/' 
XCEL_dump_top  <- 'L:/UWL/BigDataProjectXcelFiles/'

# Seting ID vectors to subset openSNP_users
setwd(paste0(top_level_path, ext1_path))

# PD
set.seed(1000)
pd_user_sample <- sample(
  as.vector(scan('pd_pos_userids_sample.txt',character(), sep = ' ')),
  5
)

# Control
set.seed(333)
control_1 <- sample(as.vector(scan('a1a_pos_userids_sample.txt', character(), sep = ' ')), 10)
set.seed(333)
control_2 <- sample(as.vector(scan('g6pd_pos_userids_sample.txt', character(), sep = ' ')), 10)
set.seed(333)
control_3 <- sample(as.vector(scan('hhem_pos_userids_sample.txt', character(), sep = ' ')), 10)
set.seed(333)
control_4 <- sample(as.vector(scan('f11_pos_userids_sample.txt', character(), sep = ' ')), 10)
set.seed(333)
control_5 <- sample(as.vector(scan('celiac_pos_userids_sample.txt', character(), sep = ' ')), 10)

set.seed(333)
control_user_sample <- c(sample(control_1, 1),
                         sample(control_2, 1),
                         sample(control_3, 1),
                         sample(control_4, 1),
                         sample(control_5, 1))

# API to bring in OpenSNP total user table
openSNP_users <- GET('http://opensnp.org/users.json',
                     content_type('application/json'),
                     user_agent('zachary.h.sykes@gmail.com -- Genetic risk prediction project'))
openSNP_users <- jsonlite::fromJSON(
  httr::content(openSNP_users, as = 'text')
)


##### PD SAMPLE GENO EXTRACTION #####

# binding df to keep the genotypes in mem
rs       <- as.character(NA)
chrom    <- as.character(NA)
pos      <- as.character(NA)
Genotype <- as.character(NA)
UserID   <- as.factor(NA)

binding_df <- tibble(rs, chrom, pos, Genotype, UserID)

# PD user sample filter from API pull
pd_sample_genotypes <- openSNP_users %>%
  filter(id %in% pd_user_sample)

for (i in 1:length(pd_sample_genotypes$id)) {
  
  # bringing in singular genotype table
  pd_case_genotype <- readLines(curl(pd_sample_genotypes[[3]][[i]]$download_url[1]))
  pd_case_genotype <- data.frame(pd_case_genotype)
  pd_case_genotype <- pd_case_genotype %>%
    filter(!(grepl('#', pd_case_genotype))) %>%
    separate(pd_case_genotype, 
             into = c('rs', 'chrom', 'pos', 'Genotype'),
             sep = '\t') 
  
  # add the userid col
  UserID <- rep(as.character(pd_sample_genotypes[[2]][[i]]), length(pd_case_genotype$rs))
  pd_case_genotype <- cbind(pd_case_genotype, UserID)
  
  # output test prior to binding all 5
  print(colnames(pd_case_genotype))
  print(paste0('Number of cols: ', length(pd_case_genotype)))
  print(paste0('Number of obs: ', length(pd_case_genotype$rs)))
  
  # bind the df into one
  binding_df <- rbind(binding_df, pd_case_genotype)
  
}
pd_case_genotypes <- binding_df[2:length(binding_df$rs), ]

# Removing '--'
genotype_col      <- pd_case_genotypes$Genotype
index_w_blanks    <- grep('--', genotype_col)
genotype_col      <- replace(genotype_col, index_w_blanks, 'NN')
pd_case_genotypes <- pd_case_genotypes %>%
  mutate(Genotype = genotype_col)


##### CONTROL SAMPLE GENO EXTRACTION #####

# binding df to keep the genotypes in mem
rs       <- as.character(NA)
chrom    <- as.character(NA)
pos      <- as.character(NA)
Genotype <- as.character(NA)
UserID   <- as.factor(NA)

binding_df <- tibble(rs, chrom, pos, Genotype, UserID)

# Control user sample filter from API pull
control_sample_genotypes <- openSNP_users %>%
  filter(id %in% control_user_sample)

for (i in 1:length(control_sample_genotypes$id)) {

  # bringing in control genotype table
  control_genotype <- readLines(curl(control_sample_genotypes[[3]][[i]]$download_url[1]))
  control_genotype <- data.frame(control_genotype)
  control_genotype <- control_genotype %>%
    filter(!(grepl('#', control_genotype))) %>%
    separate(control_genotype, 
             into = c('rs', 'chrom', 'pos', 'Genotype'),
             sep = '\t')
  
  # add the userid col
  UserID <- rep(as.character(control_sample_genotypes[[2]][[i]]), length(control_genotype$rs))
  control_genotype <- cbind(control_genotype, UserID)
  
  # output test prior to binding all 5
  print(colnames(control_genotype))
  print(paste0('Number of cols: ', length(control_genotype)))
  print(paste0('Number of obs: ', length(control_genotype$rs)))
  
  # bind the df into one
  binding_df <- rbind(binding_df, control_genotype)
  
}
control_genotype <- binding_df[2:length(binding_df$rs), ]

# Removing '--'
genotype_col     <- control_genotype$Genotype
index_w_blanks   <- grep('--', genotype_col)
genotype_col     <- replace(genotype_col, index_w_blanks, 'NN')
control_genotype <- control_genotype %>%
  mutate(Genotype = genotype_col)


##### HAPMAP CONSTRUCTION #####
setwd(top_level_path)

pd_hapmap <- rbind(pd_case_genotypes, control_genotype)
pd_hapmap <- pd_hapmap %>%
  spread(UserID, Genotype, fill = 'NN') %>%
  convert(num(pos)) %>%
  arrange(chrom, pos)

pd_hapmap$chrom <- factor(pd_hapmap$chrom, levels = c('1', '2', '3', '4', '5', '6', '7',
                                                      '8', '9', '10', '11', '12', '13', '14', '15', '16',
                                                      '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT'))
pd_hapmap$chrom <- revalue(pd_hapmap$chrom, c('X' = '23', 'Y' = '24', 'MT' = '25'))
pd_hapmap$chrom <- as.numeric(pd_hapmap$chrom)
pd_hapmap$pos   <- as.numeric(pd_hapmap$pos)
pd_hapmap$rs    <- as.factor(pd_hapmap$rs)


pd_hapmap$alleles   <- as.character(rep(NA, length(pd_hapmap$rs)))
pd_hapmap$strand    <- as.character(rep(NA, length(pd_hapmap$rs)))
pd_hapmap$assembly  <- as.character(rep(NA, length(pd_hapmap$rs)))
pd_hapmap$center    <- as.character(rep(NA, length(pd_hapmap$rs)))
pd_hapmap$protLSID  <- as.character(rep(NA, length(pd_hapmap$rs)))
pd_hapmap$assayLSID <- as.character(rep(NA, length(pd_hapmap$rs)))
pd_hapmap$panel     <- as.character(rep(NA, length(pd_hapmap$rs)))
pd_hapmap$QCcode    <- as.character(rep(NA, length(pd_hapmap$rs)))

pd_hapmap <- pd_hapmap[, c(1, 14, 2, 3, 15, 16, 17, 18, 19,
                           20, 21, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)]

### PHENOTYPE CONSTRUCTION ###
pheno_col_names         <- c('UserID', 'Pheno')
pd_pheno_df             <- tibble(pd_user_sample,
                                  rep(1, length(pd_user_sample)))
control_pheno_df        <- tibble(control_user_sample,
                                  rep(0, length(control_user_sample)))
names(pd_pheno_df)      <- pheno_col_names
names(control_pheno_df) <- pheno_col_names

# Final pd phenotype df
pheno_df <- rbind(pd_pheno_df, control_pheno_df)
  
setwd(paste0(top_level_path, ext2_path))

write.table(pd_hapmap, 'pd_hapmap.txt', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(pheno_df, 'pd_pheno.txt', quote = FALSE, sep = '\t', row.names = FALSE)

### NOTE 18.06.2020 -- Needs spark cluster connection to show competence