### Numericalisation and Separation of LOAD Genotypes
### 17.06.2020

#### VERY IMPORTANT ####
# Run LOAD_HapMap_and_Pheno_Construction.R first

#lib import
## GAPIT REQUIREMENTS
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
## GAPIT REQUIREMENTS

library(plyr)
library(hablar)

top_level_path <- 'C:/Users/zacha/Documents/BigData/' # Change to match your file structure
ext1_path      <- 'LOAD_Geno_by_Chrom/'
ext2_path      <- 'HapMap_and_Pheno_Files/' 

setwd(paste0(top_level_path, ext2_path))

# Genotype Numericalisation
myY  <- read.table('load_pheno.txt', head = TRUE)
myG  <- read.table('load_hapmap.txt', head = FALSE)


x <- GAPIT.HapMap(myG)

numeric_colanmes        <- as.vector(t(x$GT))
numeric_genotype        <- as.data.frame(t(x$GD))
names(numeric_genotype) <- numeric_colanmes

rs_chrom_pos <- x$GI

# Convert 0,1,2 to -1,0,1
numeric_genotype$`2890` <- factor(numeric_genotype$`2890`, levels = c(0, 1, 2))
numeric_genotype$`4038` <- factor(numeric_genotype$`4038`, levels = c(0, 1, 2))
numeric_genotype$`5003` <- factor(numeric_genotype$`5003`, levels = c(0, 1, 2))
numeric_genotype$`7056` <- factor(numeric_genotype$`7056`, levels = c(0, 1, 2))
numeric_genotype$`1639` <- factor(numeric_genotype$`1639`, levels = c(0, 1, 2))
numeric_genotype$`441`  <- factor(numeric_genotype$`441`, levels = c(0, 1, 2))
numeric_genotype$`6573` <- factor(numeric_genotype$`6573`, levels = c(0, 1, 2))
numeric_genotype$`4280` <- factor(numeric_genotype$`4280`, levels = c(0, 1, 2))
numeric_genotype$`7678` <- factor(numeric_genotype$`7678`, levels = c(0, 1, 2))
numeric_genotype$`9170` <- factor(numeric_genotype$`9170`, levels = c(0, 1, 2))

numeric_genotype$`2890` <- mapvalues(numeric_genotype$`2890`,
                                     from = c(0, 1, 2), to = c(-1, -0, 1))
numeric_genotype$`4038` <- mapvalues(numeric_genotype$`4038`, 
                                     from = c(0, 1, 2), to = c(-1, 0, 1))
numeric_genotype$`5003` <- mapvalues(numeric_genotype$`5003`, 
                                     from = c(0, 1, 2), to = c(-1, 0, 1))
numeric_genotype$`7056` <- mapvalues(numeric_genotype$`7056`, 
                                     from = c(0, 1, 2), to = c(-1, 0, 1))
numeric_genotype$`1639` <- mapvalues(numeric_genotype$`1639`, 
                                     from = c(0, 1, 2), to = c(-1, 0, 1))
numeric_genotype$`441`  <- mapvalues(numeric_genotype$`441`, 
                                     from = c(0, 1, 2), to = c(-1, 0, 1))
numeric_genotype$`6573` <- mapvalues(numeric_genotype$`6573`, 
                                     from = c(0, 1, 2), to = c(-1, 0, 1))
numeric_genotype$`4280` <- mapvalues(numeric_genotype$`4280`, 
                                     from = c(0, 1, 2), to = c(-1, 0, 1))
numeric_genotype$`7678` <- mapvalues(numeric_genotype$`7678`,
                                     from = c(0, 1, 2), to = c(-1, 0, 1))
numeric_genotype$`9170` <- mapvalues(numeric_genotype$`9170`, 
                                     from = c(0, 1, 2), to = c(-1, 0, 1))

numeric_swap <- numeric_genotype[, c('2890', '4038', '5003', '7056', '1639',  
                                     '441', '6573', '4280', '7678', '9170')]

numeric_swap <- numeric_genotype %>%
  convert(num(c(`2890`, `4038`, `5003`, `7056`, `1639`,
                `441`, `6573`, `4280`, `7678`, `9170`))) 
numeric_load_genotype <- cbind(rs_chrom_pos, numeric_swap)
numeric_load_genotype <- numeric_load_genotype %>%
  mutate(as.character(Chromosome))
numeric_load_genotype$Chromosome <- factor(numeric_load_genotype$Chromosome,
                                           levels = c('1', '2', '3', '4', '5', '6', '7',
                                                      '8', '9', '10', '11', '12', '13', '14', '15', '16',
                                                      '17', '18', '19', '20', '21', '22', '23', '24', '25'))
numeric_load_genotype <- numeric_load_genotype %>%
  convert(num(c(Chromosome, Position))) %>%
  mutate(as.character(SNP))
numeric_load_genotype <- as.data.frame(numeric_load_genotype)

# Separate data by chrom and write to disk
numeric_load_genotype_chrom1 <- numeric_load_genotype %>%
  filter(grepl('^1$', Chromosome))
numeric_load_genotype_chrom2 <- numeric_load_genotype %>%
  filter(grepl('^2$', Chromosome))
numeric_load_genotype_chrom3 <- numeric_load_genotype %>%
  filter(grepl('^3$', Chromosome))
numeric_load_genotype_chrom4 <- numeric_load_genotype %>%
  filter(grepl('^4$', Chromosome))
numeric_load_genotype_chrom5 <- numeric_load_genotype %>%
  filter(grepl('^5$', Chromosome))
numeric_load_genotype_chrom6 <- numeric_load_genotype %>%
  filter(grepl('^6$', Chromosome))
numeric_load_genotype_chrom7 <- numeric_load_genotype %>%
  filter(grepl('^7$', Chromosome))
numeric_load_genotype_chrom8 <- numeric_load_genotype %>%
  filter(grepl('^8$', Chromosome))
numeric_load_genotype_chrom9 <- numeric_load_genotype %>%
  filter(grepl('^9$', Chromosome))
numeric_load_genotype_chrom10 <- numeric_load_genotype %>%
  filter(grepl('^10$', Chromosome))
numeric_load_genotype_chrom11 <- numeric_load_genotype %>%
  filter(grepl('^11$', Chromosome))
numeric_load_genotype_chrom12 <- numeric_load_genotype %>%
  filter(grepl('^12$', Chromosome))
numeric_load_genotype_chrom13 <- numeric_load_genotype %>%
  filter(grepl('^13$', Chromosome))
numeric_load_genotype_chrom14 <- numeric_load_genotype %>%
  filter(grepl('^14$', Chromosome))
numeric_load_genotype_chrom15 <- numeric_load_genotype %>%
  filter(grepl('^15$', Chromosome))
numeric_load_genotype_chrom16 <- numeric_load_genotype %>%
  filter(grepl('^16$', Chromosome))
numeric_load_genotype_chrom17 <- numeric_load_genotype %>%
  filter(grepl('^17$', Chromosome))
numeric_load_genotype_chrom18 <- numeric_load_genotype %>%
  filter(grepl('^18$', Chromosome))
numeric_load_genotype_chrom19 <- numeric_load_genotype %>%
  filter(grepl('^19$', Chromosome))
numeric_load_genotype_chrom20 <- numeric_load_genotype %>%
  filter(grepl('^20$', Chromosome))
numeric_load_genotype_chrom21 <- numeric_load_genotype %>%
  filter(grepl('^21$', Chromosome))
numeric_load_genotype_chrom22 <- numeric_load_genotype %>%
  filter(grepl('^22$', Chromosome))
numeric_load_genotype_chromX <- numeric_load_genotype %>%
  filter(grepl('^23$', Chromosome))
numeric_load_genotype_chromY <- numeric_load_genotype %>%
  filter(grepl('^24$', Chromosome))
numeric_load_genotype_chromMT <- numeric_load_genotype %>%
  filter(grepl('^25$', Chromosome))

chrom_sep_genotypes <- list(numeric_load_genotype_chrom1, numeric_load_genotype_chrom2, 
                            numeric_load_genotype_chrom3, numeric_load_genotype_chrom4,
                            numeric_load_genotype_chrom5, numeric_load_genotype_chrom6,
                            numeric_load_genotype_chrom7, numeric_load_genotype_chrom8,
                            numeric_load_genotype_chrom9, numeric_load_genotype_chrom10,
                            numeric_load_genotype_chrom11, numeric_load_genotype_chrom12,
                            numeric_load_genotype_chrom13, numeric_load_genotype_chrom14,
                            numeric_load_genotype_chrom15, numeric_load_genotype_chrom16,
                            numeric_load_genotype_chrom17, numeric_load_genotype_chrom18,
                            numeric_load_genotype_chrom19, numeric_load_genotype_chrom20,
                            numeric_load_genotype_chrom21, numeric_load_genotype_chrom22,
                            numeric_load_genotype_chromX, numeric_load_genotype_chromY,
                            numeric_load_genotype_chromMT)

list_names <- as.character(c(1:22, 'X', 'Y', 'MT'))
names(chrom_sep_genotypes) <- list_names

if(!(dir.exists('LOAD_Geno_by_Chrom'))) {
  dir.create('LOAD_Geno_by_Chrom')
}

# Where the .txt files will be written
setwd(paste0(top_level_path, ext_path))

for(i in 1:length(chrom_sep_genotypes)) {
  write.table(
    chrom_sep_genotypes[[i]],
    paste0('LOAD_num_gen_chrom_', i, '.txt'),
    quote = FALSE,
    sep = '\t',
    row.names = FALSE
  )
}