### PD and LOAD GWAS and plot construction

# Fresh global env
rm(list = ls())

# lib import
library(rrBLUP)
library(GWASTools)
library(tidyverse)
library(sparklyr)

# PATHs
top_level_path <- 'C:/Users/zacha/Documents/BigData/' # Change to match your file structure
ext1_path      <- 'ID_samples_TextFiles/'
ext2_path      <- 'HapMap_and_Pheno_Files/' 
ext3_path      <- 'PD_Geno_by_Chrom/'
ext4_path      <- 'LOAD_Geno_by_Chrom/'

##### PD GWAS RUN #####

# set dir to PD pheno an bring in traits
setwd(paste0(top_level_path, ext2_path))
pd_pheno <- read.table('pd_pheno.txt', header = TRUE)

# set to dir with PD geno for GWAS
setwd(paste0(top_level_path, ext3_path))

# Building binding df
SNP        <- as.factor(000)
Chromosome <- as.integer(000)
Position   <- as.integer(000)
Pheno      <- as.numeric(000)
binding_df <- tibble(SNP, Chromosome, Position, Pheno)

for(file in list.files(paste0(top_level_path, ext3_path))) {
  
  # Bring in PD genotype
  geno <- read.table(file, header = TRUE, check.names = FALSE)
 
  # Run GWAS per chrom - rrBLUP
  gwas <- GWAS(
    pheno = pd_pheno,
    geno = geno,
    min.MAF = 0.05,
    plot = FALSE
  )
  
  # Output marker 
  print(paste0('Output for: ', file, ', complete'))
  
  # Bind together
  binding_df <- rbind(binding_df, gwas)
  
}
PD_GWAS_RESULTS <- binding_df[2:length(binding_df$SNP),]
PD_GWAS_RESULTS <- PD_GWAS_RESULTS %>%
  arrange(Chromosome) %>%
  filter(Pheno != 0.0000000)

# PD Plots
manhattanPlot(
  p = PD_GWAS_RESULTS$Pheno,
  chromosome = PD_GWAS_RESULTS$Chromosome,
  signif = 0.05
)

qqPlot(
  pval = PD_GWAS_RESULTS$Pheno,
  truncate = TRUE
)

##### LOAD GWAS RUN #####

# set dir to PD pheno an bring in traits
setwd(paste0(top_level_path, ext2_path))
load_pheno <- read.table('load_pheno.txt', header = TRUE)

# set to dir with LOAD geno for GWAS
setwd(paste0(top_level_path, ext4_path))

# Building binding df
SNP        <- as.factor(000)
Chromosome <- as.integer(000)
Position   <- as.integer(000)
Pheno      <- as.numeric(000)
binding_df <- tibble(SNP, Chromosome, Position, Pheno)

for(file in list.files(paste0(top_level_path, ext4_path))) {
  
  # Bring in LOAD genotype
  geno <- read.table(file, header = TRUE, check.names = FALSE)
  
  # Run GWAS per chrom - rrBLUP
  gwas <- GWAS(
    pheno = load_pheno,
    geno = geno,
    min.MAF = 0.05,
    plot = FALSE
  )
  
  # Output marker 
  print(paste0('Output for: ', file, ', complete'))
  
  # Bind together
  binding_df <- rbind(binding_df, gwas)
  
}
LOAD_GWAS_RESULTS <- binding_df[2:length(binding_df$SNP), ]
LOAD_GWAS_RESULTS <- LOAD_GWAS_RESULTS %>%
  arrange(Chromosome) %>%
  filter(Pheno != 0.0000000)

# LOAD Plots
manhattanPlot(
  p = LOAD_GWAS_RESULTS$Pheno,
  chromosome = LOAD_GWAS_RESULTS$Chromosome,
  signif = 0.05
)

qqPlot(
  pval = LOAD_GWAS_RESULTS$Pheno,
  truncate = TRUE
)