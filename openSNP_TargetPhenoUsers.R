### openSNP API Construct to Pull User IDs with Target Phenotype Associated SNPs
### 11.05.2020

# Clear environment prior to work
rm(list = ls())

# Package import for API
library(httr)
library(jsonlite)
library(xml2)
library(fs)
library(tidyverse)
library(curl)
library(rlist)

# PATHs - change for your machine
top_level_dir <- 'C:/Users/zacha/Documents/BigData'
snp_rds_save_dir <- 'C:/Users/zacha/Documents/BigData/SNP_RDS'
id_samples_dir <- 'C:/Users/zacha/Documents/BigData/ID_samples_TextFiles/'

# Creating Vectors of rsIDs with target disease association
# I am pulling all I can find from studies so there may be duplicates
# to fix this I am applying unique() over each by vector
g6pd_rsIDs <- c('rs1059828', 'rs1050829', 'rs2230037', 'rs2071429', 'rs5030868')

a1a_rsIDs <- c('rs1243168', 'rs361525', 'rs8034191')

celiac_rsIDs <- c('rs2241880', 'rs11209026', 'rs10883365', 'rs888208', 'rs1248696', 'rs601338',
                  'rs2298428', 'rs4374642', 'rs6974491', 'rs3087456', 'rs4774', 'rs1805087',
                  'rs10892279', 'rs1800795', 'rs458046', 'rs7259292', 'rs2305767', 'rs2305765',
                  'rs1004819', 'rs7517847', 'rs10754558', 'rs6822844', 'rs1799969', 'rs3184504',
                  'rs1801274', 'rs7040561', 'rs11552708', 'rs11203203', 'rs1545620', 'rs1801274',
                  'rs6840978', 'rs1801394', 'rs1800795', 'rs11209026', 'rs13119723', 'rs3184504',
                  'rs1799969', 'rs2305767', 'rs7259292', 'rs10754558', 'rs17810546', 'rs10892279',
                  'rs11209026', 'rs1800795', 'rs802734') 
celiac_rsIDs <- unique(celiac_rsIDs)

eo_dystonia_rsIDs <- c('rs10483639', 'rs1559510', 'rs1182', 'rs1801968')

F11_rsIDs <- c('rs925453', 'rs710446', 'rs925451', 'rs925451', 'rs35927125', 
               'rs6427340', 'rs945508', 'rs10505346', 'rs1884613', 'rs1333040',
               'rs9298506', 'rs7927894', 'rs699947') 
F11_rsIDs <- unique(F11_rsIDs)

hhem_rsIDs <- c('rs807212', 'rs738409', 'rs1800562', 'rs2304704', 'rs855791', 'rs1695', 'rs738409',
                'rs1799945', 'rs855791', 'rs1800730', 'rs1799963', 'rs1801133', 'rs6025',
                'rs1801131') 
hhem_rsIDs <- unique(hhem_rsIDs)

hthromb_rsIDs <- c('rs1799963', 'rs2731672', 'rs13146272', 'rs2036914', 'rs1800872', 'rs6589488', 'rs2066865',
                   'rs3862019', 'rs4680', 'rs6025', 'rs7542281', 'rs1063856', 'rs9898', 'rs710446', 'rs1799963', 
                   'rs1045642', 'rs6050', 'rs10974944', 'rs2036914', 'rs4680', 'rs5370', 'rs8176719', 'rs2066865',
                   'rs1800790', 'rs1042031', 'rs3917643', 'rs4680', 'rs8176719', 'rs1045642', 'rs1143634', 'rs4220',
                   'rs6048', 'rs1799963', 'rs1041296', 'rs1613662', 'rs1063856', 'rs253061', 'rs2036914', 
                   'rs2731672', 'rs10974944','rs2269648', 'rs333948', 'rs8176719') 
hthromb_rsIDs <- unique(hthromb_rsIDs)

lo_alz_rsIDs <- c('rs4420638', 'rs16934131', 'rs668387', 'rs641120', 'rs6330', 'rs1768208', 'rs2072374', 'rs1042522',
                  'rs727153', 'rs5984894', 'rs727153', 'rs5984894', 'rs63750847', 'rs10884402', 'rs7078098', 'rs2070045',
                  'rs11610206', 'rs744373', 'rs3851179', 'rs1524107', 'rs4343', 'rs3731211', 'rs1927907', 'rs2075650',
                  'rs242557', 'rs6584777', 'rs3818361', 'rs7294919', 'rs3851179', 'rs9349407', 'rs8063', 'rs11136000', 'rs2230806',
                  'rs760678', 'rs6656401', 'rs12344615', 'rs17070145', 'rs10793294', 'rs1937', 'rs5984894', 'rs1800587',
                  'rs2075650', 'rs561655', 'rs10972300', 'rs11767557', 'rs2075650', 'rs165932', 'rs1385600', 'rs1160985',
                  'rs1800587', 'rs17070145', 'rs12344615', 'rs2075650', 'rs11136000', 'rs3818361', 'rs5984894',
                  'rs1937', 'rs2986017', 'rs7561528', 'rs6656401', 'rs17070145', 'rs3851179', 'rs760678', 'rs1801270',
                  'rs610932', 'rs2306604', 'rs2027432', 'rs17070145', 'rs7561528', 'rs4680', 'rs1990622', 'rs11136000', 
                  'rs2373115', 'rs597668', 'rs2772677', 'rs11515', 'rs5984894', 'rs5984894', 'rs892086', 'rs2075650', 'rs2177369',
                  'rs2571598', 'rs689021', 'rs2282649') 
lo_alz_rsIDs <- unique(lo_alz_rsIDs)

park_rsIDs <- c('rs11931532', 'rs823156', 'rs823128', 'rs4538475', 'rs334558', 'rs823144', 'rs708730', 'rs421016',
                'rs3775444', 'rs1572931', 'rs2736990', 'rs356220', 'rs3822086', 'rs356219', 'rs11931074', 
                'rs4947342', 'rs4248166', 'rs2858324', 'rs3129299', 'rs2858880', 'rs3129888', 'rs9277489',
                'rs7769979', 'rs2395173', 'rs4434496', 'rs3129859', 'rs4418214', 'rs3763349', 'rs169494', 'rs7773756',
                'rs2040410', 'rs2395148', 'rs3129882', 'rs2395185', 'rs3129763', 'rs7454108', 'rs1150758', 'rs646984',
                'rs3134942', 'rs1612904', 'rs3763313', 'rs2395175', 'rs6936204', 'rs2187818', 'rs399604', 'rs2075800',
                'rs1060619', 'rs6489721', 'rs3754775', 'rs6740826', 'rs10847864', 'rs6723108', 'rs11724635', 'rs4538475',
                'rs11248060', 'rs823156', 'rs6532194', 'rs7077361', 'rs1491942', 'rs947211', 'rs2390669', 'rs34778348',
                'rs356219', 'rs156429', 'rs6812193', 'rs2301134', 'rs708723', 'rs2245801', 'rs4889603', 
                'rs12185268', 'rs7412', 'rs10410544', 'rs11716740', 'rs4837628', 'rs12185268', 'rs1994090',
                'rs34637584', 'rs12431733', 'rs1223271', 'rs13312', 'rs797906', 'rs823128', 'rs10464059', 
                'rs1799836', 'rs2205108', 'rs872606', 'rs6280', 'rs356220', 'rs1043424', 'rs7617877', 'rs1876828',
                'rs4954218', 'rs12063142', 'rs2010795', 'rs6812193', 'rs1079597', 'rs1801133', 'rs1805874', 'rs10737170',
                'rs11868035', 'rs3129882', 'rs1801582', 'rs4880', 'rs11248060', 'rs4698412', 'rs823156', 'rs17115100',
                'rs11175655', 'rs11030104', 'rs7077361', 'rs10200894', 'rs4130047', 'rs10513789', 'rs9917256', 
                'rs2282048', 'rs6599389', 'rs17329669', 'rs1491942', 'rs5174', 'rs7412', 'rs11060112', 'rs2823357',
                'rs11649804', 'rs11724635', 'rs1799836', 'rs878396', 'rs4149589', 'rs10784293', 'rs6734184', 'rs13014473',
                'rs2059198', 'rs13006838', 'rs7315790', 'rs6714092', 'rs10450989', 'rs2548278', 'rs3775444', 'rs3822086',
                'rs11701', 'rs6518956', 'rs7308720', 'rs10878371', 'rs1427263', 'rs33949390', 'rs11176013',
                'rs3761863', 'rs7133914', 'rs7966550', 'rs11564148', 'rs10878405', 'rs34778348', 'rs11175964')
park_rsIDs <- unique(park_rsIDs)

# Create the save location for the snp data rds files
if (!(dir.exists(snp_rds_save_dir))) {
  dir.create(snp_rds_save_dir)
}

setwd(snp_rds_save_dir)

# Function to employ API that pulls genotype data of users with SNPs associated with specific diseases  
rsID_based_data_extractor <- function(rsID) {
  
  # URL setup
  core_url <- 'http://opensnp.org/snps/'
  ext_url <- paste0(rsID, '.json')
  
  # API execute - Don't forget to delay the API runs
  rsID_users_resp <- GET(
    paste0(core_url, ext_url), content_type('application/json'),
    user_agent('zachary.h.sykes@gmail.com -- Genetic risk prediction project')
  )
  
  # Save the SNP response as a RDS file
  if(http_type(rsID_users_resp) != 'application/json') {
    warning('No information returned for this gene id')
  } else {
    setwd(snp_rds_save_dir)
    saveRDS(rsID_users_resp, paste0(snp, '_json.rds'))
    print(paste('RDS file saved in', getwd(), sep = ' '))
  }
  
  # DF return 
  return(
    jsonlite::fromJSON(content(rsID_users_resp, as = 'text'))
  )
  
}


### Step 1 - Build df of all users with variants associated with specific disease (only those with genotype data)

# Building a df to bind each iteration to to avoid overwriting
name <- as.character('REMOVE')
id <- as.integer(000)
genotypes <- list()

binding_df <- tibble(name, id, genotypes) # Proper structure should avoid unwanted coercion

##########################################################################################################
# SNP disease association is based on publications surrounding each variant sourced on OpenSNP's website #
##########################################################################################################


##### Users with G6PD assoc SNP #####
# Timer start
loop_start <- Sys.time()

count <- 0

for (snp in g6pd_rsIDs) {
  
  # Bring in data from API
  g6pd_pos_users <- rsID_based_data_extractor(snp) # g6pd pos means users we assume have g6pd
  g6pd_pos_users <- g6pd_pos_users$user # remove snp df
  indicies_w_genotypes <- c()
  
  # Remove users without genotype data
  for (i in 1:length(g6pd_pos_users$genotypes)) {
    
    if (!(is.null(unlist(g6pd_pos_users$genotypes[i])))) {
      indicies_w_genotypes <- c(indicies_w_genotypes, i)
    } # End if statement
    
  } # End for loop checking indicies
  
  # bind everything together
  g6pd_pos_users <- g6pd_pos_users[indicies_w_genotypes, ]
  binding_df <- rbind(binding_df, g6pd_pos_users)
  
  count <- count + 1
  
  Sys.sleep(10)
  print(paste(count, 'loop complete', sep = ' '))
  print(Sys.time())

} # End g6pd rsID for loop

g6pd_pos_users <- binding_df # Stores data pull allowing binding df to be overwritten
binding_df <- tibble(name, id, genotypes)

# Stop timer
end_time <- Sys.time()

print(end_time - loop_start)


##### Users with Alpha-1-Antitrypsin Deficiency assoc SNP #####
# Timer start
loop_start <- Sys.time()

count <- 0

for (snp in a1a_rsIDs) {
  
  # Bring in data from API
  a1a_pos_users <- rsID_based_data_extractor(snp) # a1a pos means users we assume have a1a
  a1a_pos_users <- a1a_pos_users$user # remove snp df
  indicies_w_genotypes <- c()
  
  # Remove users without genotype data
  for (i in 1:length(a1a_pos_users$genotypes)) {
    
    if (!(is.null(unlist(a1a_pos_users$genotypes[i])))) {
      indicies_w_genotypes <- c(indicies_w_genotypes, i)
    } # End if statement
    
  } # End for loop checking indicies
  
  # Bind everything together
  a1a_pos_users <- a1a_pos_users[indicies_w_genotypes, ]
  binding_df <- rbind(binding_df, a1a_pos_users)
  
  count <- count + 1
  
  Sys.sleep(10)
  print(paste(count, 'loop complete', sep = ' '))
  print(Sys.time())
  
} # End a1a rsID for loop

a1a_pos_users <- binding_df
binding_df <- tibble(name, id, genotypes)

# Stop timer
end_time <- Sys.time()

print(end_time - loop_start)


##### Users with Factor XI assoc SNP #####
# Timer start
loop_start <- Sys.time()

count <- 0

for (snp in F11_rsIDs) {
  
  # Bring in data from API 
  f11_pos_users <- rsID_based_data_extractor(snp) # see reasoning in two data pulls above
  f11_pos_users <- f11_pos_users$user # see reasoning in two data pulls above
  indicies_w_genotypes <- c()
  
  # Remove users without genotype data
  for (i in 1:length(f11_pos_users$genotypes)) {
    
    if(!(is.null(unlist(f11_pos_users$genotypes[i])))) {
      indicies_w_genotypes <- c(indicies_w_genotypes, i)
    } # End if statement
    
  } # End for loop checking indicies
  
  # Bind everything together
  f11_pos_users <- f11_pos_users[indicies_w_genotypes, ]
  binding_df <- rbind(binding_df, f11_pos_users)
  
  count <- count + 1
  
  Sys.sleep(10)
  print(paste(count, 'loop complete', sep = ' '))
  print(Sys.time())
  
} # End f11 rsID for loop

f11_pos_users <- binding_df
binding_df <- tibble(name, id, genotypes)

# Stop timer
end_time <- Sys.time()

print(end_time - loop_start)

##### Users with Celiac assoc SNP #####
# Timer start
loop_start <- Sys.time()

count <- 0

for (snp in celiac_rsIDs) {
  
  # Bring in data from API 
  celiac_pos_users <- rsID_based_data_extractor(snp)
  celiac_pos_users <- celiac_pos_users$user
  indicies_w_genotypes <- c()
  
  # Remove users without genotype data
  for (i in 1:length(celiac_pos_users$genotypes)) {
    
    if(!(is.null(unlist(celiac_pos_users$genotypes[i])))) {
      indicies_w_genotypes <- c(indicies_w_genotypes, i)
    } # End if statement
    
  } # End for loop checking indicies
  
  # Bind everything together
  celiac_pos_users <- celiac_pos_users[indicies_w_genotypes, ]
  binding_df <- rbind(binding_df, celiac_pos_users)
  
  count <- count + 1
  
  Sys.sleep(10)
  print(paste(count, 'loop complete', sep = ' '))
  print(Sys.time())
  
} # End celiac rsID for loop

celiac_pos_users <- binding_df
binding_df <- tibble(name, id, genotypes)

# Stop timer
end_time <- Sys.time()

print(end_time - loop_start)


##### Users with Hereditary Hemochromatosis assoc SNP #####
# Timer start
loop_start <- Sys.time()

count <- 0

for (snp in hhem_rsIDs) {
  
  # Bring in data from API 
  hhem_pos_users <- rsID_based_data_extractor(snp)
  hhem_pos_users <- hhem_pos_users$user
  indicies_w_genotypes <- c()
  
  # Remove users without genotype data
  for (i in 1:length(hhem_pos_users$genotypes)) {
    
    if(!(is.null(unlist(hhem_pos_users$genotypes[i])))) {
      indicies_w_genotypes <- c(indicies_w_genotypes, i)
    } # End if statement
    
  } # End for loop checking indicies
  
  # Bind everything together
  hhem_pos_users <- hhem_pos_users[indicies_w_genotypes, ]
  binding_df <- rbind(binding_df, hhem_pos_users)
  
  count <- count + 1
  
  Sys.sleep(10)
  print(paste(count, 'loop complete', sep = ' '))
  print(Sys.time())
  
} # End hhem rsID for loop

hhem_pos_users <- binding_df
binding_df <- tibble(name, id, genotypes)

# Stop timer
end_time <- Sys.time()

print(end_time - loop_start)


##### Users with Hereditary Thrombophilia assoc SNP #####
# Timer start
loop_start <- Sys.time()

count <- 0

for (snp in hthromb_rsIDs) {
  
  # Bring in data from API 
  hthromb_pos_users <- rsID_based_data_extractor(snp)
  hthromb_pos_users <- hthromb_pos_users$user
  indicies_w_genotypes <- c()
  
  # Remove users without genotype data
  for (i in 1:length(hthromb_pos_users$genotypes)) {
    
    if(!(is.null(unlist(hthromb_pos_users$genotypes[i])))) {
      indicies_w_genotypes <- c(indicies_w_genotypes, i)
    } # End if statement
    
  } # End for loop checking indicies
  
  # Bind everything together
  hthromb_pos_users <- hthromb_pos_users[indicies_w_genotypes, ]
  binding_df <- rbind(binding_df, hthromb_pos_users)
  
  count <- count + 1
  
  Sys.sleep(10)
  print(paste(count, 'loop complete', sep = ' '))
  print(Sys.time())
  
} # End hthromb rsID for loop

hthromb_pos_users <- binding_df
binding_df <- tibble(name, id, genotypes)

# Stop timer
end_time <- Sys.time()

print(end_time - loop_start)


##### Users with Early Onset Primary Dystonia assoc SNP #####
# Timer start
loop_start <- Sys.time()

count <- 0

for (snp in eo_dystonia_rsIDs) {
  
  # Timer start
  loop_start <- Sys.time()
  
  # Bring in data from API 
  eo_dyst_pos_users <- rsID_based_data_extractor(snp)
  eo_dyst_pos_users <- eo_dyst_pos_users$user
  indicies_w_genotypes <- c()
  
  # Remove users without genotype data
  for (i in 1:length(eo_dyst_pos_users$genotypes)) {
    
    if(!(is.null(unlist(eo_dyst_pos_users$genotypes[i])))) {
      indicies_w_genotypes <- c(indicies_w_genotypes, i)
    } # End if statement
    
  } # End for loop checking indicies
  
  # Bind everything together
  eo_dyst_pos_users <- eo_dyst_pos_users[indicies_w_genotypes, ]
  binding_df <- rbind(binding_df, eo_dyst_pos_users)
  
  count <- count + 1
  
  Sys.sleep(10)
  print(paste(count, 'loop complete', sep = ' '))
  print(Sys.time())
  
} # End eo dyst rsID for loop

eo_dyst_pos_users <- binding_df
binding_df <- tibble(name, id, genotypes)

# Stop timer
end_time <- Sys.time()

print(end_time - loop_start)

##### Users with Parkinson's assoc SNP #####
# Timer start
loop_start <- Sys.time()

count <- 0

for (snp in park_rsIDs) {
  
  # Bring in data from API 
  pd_pos_users <- rsID_based_data_extractor(snp)
  pd_pos_users <- pd_pos_users$user
  indicies_w_genotypes <- c()
  
  # Remove users without genotype data
  for (i in 1:length(pd_pos_users$genotypes)) {
    
    if(!(is.null(unlist(pd_pos_users$genotypes[i])))) {
      indicies_w_genotypes <- c(indicies_w_genotypes, i)
    } # End if statement
    
  } # End for loop checking indicies
  
  # Bind everything together
  pd_pos_users <- pd_pos_users[indicies_w_genotypes, ]
  binding_df <- rbind(binding_df, pd_pos_users)
  
  count <- count + 1
  
  Sys.sleep(10)
  print(paste(count, 'loop complete', sep = ' '))
  print(Sys.time())
  
} # End pd rsID for loop

pd_pos_users <- binding_df
binding_df <- tibble(name, id, genotypes)

# Stop timer
end_time <- Sys.time()

print(end_time - loop_start)


##### Users with Late Onset Alzheimer's assoc SNP #####
# Timer start
loop_start <- Sys.time()

count <- 0

for (snp in lo_alz_rsIDs) {
  
  # Bring in data from API 
  lo_alz_pos_users <- rsID_based_data_extractor(snp)
  lo_alz_pos_users <- lo_alz_pos_users$user
  indicies_w_genotypes <- c()
  
  # Remove users without genotype data
  for (i in 1:length(lo_alz_pos_users$genotypes)) {
    
    if(!(is.null(unlist(lo_alz_pos_users$genotypes[i])))) {
      indicies_w_genotypes <- c(indicies_w_genotypes, i)
    } # End if statement
    
  } # End for loop checking indicies
  
  # Bind everything together
  lo_alz_pos_users <- lo_alz_pos_users[indicies_w_genotypes, ]
  binding_df <- rbind(binding_df, lo_alz_pos_users)
  
  count <- count + 1
  
  Sys.sleep(10)
  print(paste(count, 'loop complete', sep = ' '))
  print(Sys.time())
  
} # End lo alz rsID for loop

lo_alz_pos_users <- binding_df
binding_df <- tibble(name, id, genotypes)

# Stop timer
end_time <- Sys.time()

print(end_time - loop_start)


### Step 2 - Filter df of users by unique ids

# Filtering all data by unique users
g6pd_pos_users <- unique(g6pd_pos_users)
a1a_pos_users <- unique(a1a_pos_users)
f11_pos_users <- unique(f11_pos_users)
celiac_pos_users <- unique(celiac_pos_users)
hhem_pos_users <- unique(hhem_pos_users)
hthromb_pos_users <- unique(hthromb_pos_users)
eo_dyst_pos_users <- unique(eo_dyst_pos_users)
pd_pos_users <- unique(pd_pos_users)
lo_alz_pos_users <- unique(lo_alz_pos_users)

### Step 3 - Randomly sample 100 of the remaining users
setwd(id_samples_dir)

g6pd_user_sample <- g6pd_pos_users[sample(nrow(g6pd_pos_users), 100), ]
g6pd_user_sample <- as.vector(g6pd_user_sample$id) # creates a vector to be fed through the next API
write(g6pd_user_sample, file = 'g6pd_pos_userids_sample.txt')

a1a_user_sample <- a1a_pos_users[sample(nrow(a1a_pos_users), 100), ]
a1a_user_sample <- as.vector(a1a_user_sample$id)
write(a1a_user_sample, file = 'a1a_pos_userids_sample.txt')


f11_user_sample <- f11_pos_users[sample(nrow(f11_pos_users), 100), ]
f11_user_sample <- as.vector(f11_user_sample$id)
write(f11_user_sample, file = 'f11_pos_userids_sample.txt')

celiac_user_sample <- celiac_pos_users[sample(nrow(celiac_pos_users), 100), ]
celiac_user_sample <- as.vector(celiac_user_sample$id)
write(celiac_user_sample, file = 'celiac_pos_userids_sample.txt')

hhem_user_sample <- hhem_pos_users[sample(nrow(hhem_pos_users), 100), ]
hhem_user_sample <- as.vector(hhem_user_sample$id)
write(hhem_user_sample, file = 'hhem_pos_userids_sample.txt')

hthromb_user_sample <- hthromb_pos_users[sample(nrow(hthromb_pos_users), 100), ]
hthromb_user_sample <- as.vector(hthromb_user_sample$id)
write(hthromb_user_sample, file = 'hthromb_pos_userids_sample.txt')

eo_dyst_user_sample <- eo_dyst_pos_users[sample(nrow(eo_dyst_pos_users), 100), ]
eo_dyst_user_sample <- as.vector(eo_dyst_user_sample$id)
write(eo_dyst_user_sample, file = 'eo_dyst_pos_userids_sample.txt')

pd_user_sample <- pd_pos_users[sample(nrow(pd_pos_users), 100), ]
pd_user_sample <- as.vector(pd_user_sample$id)
write(pd_user_sample, file = 'pd_pos_userids_sample.txt')

lo_alz_user_sample <- lo_alz_pos_users[sample(nrow(lo_alz_pos_users), 100), ]
lo_alz_user_sample <- as.vector(lo_alz_user_sample$id)
write(lo_alz_user_sample, file = 'lo_alz_pos_userids_sample.txt')

### REMEMBER ###
# Because there is no formal diagnosis, we are taking patients with published rsIDs associated with a target pheno as a substitute