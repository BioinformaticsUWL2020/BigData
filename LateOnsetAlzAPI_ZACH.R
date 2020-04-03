### API GENE DATA PULL | Late-Onset Alzheimer's Disease | Zachary Sykes

# Library load
library(httr)
library(jsonlite)
library(xml2)
library(fs)

lo_alz_geneids_1 <- c('SCAPER', 'INPPD5', 'EPHA1', 'AK128216', 'STAU2',
                      'ZCWPW1', 'SORL1', 'ATM', 'HLA-DRB5', 'HLA-DRB1',
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
                      'RASSF5', 'AP2A2', 'NHLRC3', 'PROSER1', 'IQCK', 'ATXN7L1',
                      'AKR7A3', 'CELF1', 'ADAM10', 'NDUFAF6', 'CACNA2D3', 'HMGA2', 
                      'ISYNA1', 'ELL', 'SSBP4', 'MIR142', 'TSPOAP1-AS1', 'SLC10A2',
                      'CDON')

target_pheno_list <- list(
  lateonset_alzheimers_disease_1 = lo_alz_geneids_1,
  lateonset_alzheimers_disease_2 = lo_alz_geneids_2
)

######################
### VERY IMPORTANT ###
######################

# Set directory below to a folder OTHER than the git env as I have done below
# If you upload the .rds files saved below to the repository I will kill you
start_dir <- getwd()
# CHANGE THIS PATH TO SOMETHING OTHER THAN YOUR GIT FOLDER
# rds_save_dir <- setwd('H:/R_Scripts/ProjectAPISaves/') # Zach's Windows PATH
rds_save_dir <- setwd('/media/sykes/BLUE/R_Scripts/ProjectAPISaves/') # Zach's Linux PATH

######################
### VERY IMPORTANT ###
######################

# Function to shift files from working directory
file_shift <- function(namedir, filename) {
  
  dir.create(paste0(rds_save_dir, namedir)) # Create new directory to store the .rds files in
  
  files_to_move <- list.files(rds_save_dir,
                              pattern = '\\.rds')
  save_path <- paste0(rds_save_dir, namedir)
  
  for (file in files_to_move) {
    file_move(file, save_path)
  }
}

##### GET request function to pull in disease data #####

# SOURCE: Ensembl
disease_pheno_json_pull <- function(geneid, filename) {
  
  # URL set
  core_url <- 'https://rest.ensembl.org'
  extension_url <- paste0('/phenotype/gene/homo_sapiens/', geneid, '?include_overlap=1;include_associated=1')
  
  pheno_json <- GET(paste0(core_url, extension_url),
                    content_type('application/json'),
                    user_agent('zachary.h.sykes@gmail.com -- Genetic risk prediction project')) # DO NOT CHANGE MY EMAIL
  
  # Save and move files into disease specific directory
  if(http_type(pheno_json) != 'application/json') {
    warning('No information returned for this gene id')
  } else {
    saveRDS(pheno_json, filename)
  }
} # API COMPLETE

### NOTE: Add path to Ophranet xml file here and parse data out in each disease section below

# Handling response JSON files for each target phenotype requested
### BE PATIENT IF YOU HAVE TO RUN THE API IT TAKES ABOUT 10 - 40 MINUTES


################################################
##### Late-Onset Alzheimer's JSON response #####
################################################

# first gene set
lo_alz_geneids_1_resp <- list()
for (gene in target_pheno_list[[1]]) {
  disease_pheno_json_pull(
    geneid = gene,
    filename = paste0(tolower(gene), '_data_loalz1_json.rds')
  )
  lo_alz_geneids_1_resp[[gene]] <- list(
    content(readRDS(paste0(tolower(gene), '_data_loalz1_json.rds')))
  )
  Sys.sleep(15)
}

# Run this after the GET request above to clear the working directory
file_shift(
  namedir = '/lo_alz_1_pheno_data/',
  filename = paste0(tolower(target_pheno_list[[1]]), '_data_loalz1_json.rds')
)

##################
### Second Set ###
##################

# second gene set
lo_alz_geneids_2_resp <- list()
for (gene in target_pheno_list[[2]]) {
  disease_pheno_json_pull(
    geneid = gene,
    filename = paste0(tolower(gene), '_data_loalz2_json.rds')
  )
  lo_alz_geneids_2_resp[[gene]] <- list(
    content(readRDS(paste0(tolower(gene), '_data_loalz2_json.rds')))
  )
  Sys.sleep(15)
}

# Run this after the GET request above to clear the working directory
file_shift(
  namedir = '/lo_alz_2_pheno_data/',
  filename = paste0(tolower(target_pheno_list[[2]]), '_data_loalz2_json.rds')
)