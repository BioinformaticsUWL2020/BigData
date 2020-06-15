### GAPIT Tutorials
### 11.06.2020

# GAPIT - Genomic Association and Prediction Integrated Tool
# Designed by Zhiwu Zhang
# Written by Zhiwu Zhang, Alex Lipka, Feng Tian and You Tang
# Last update: September 15, 2015

#Install packages (Do this section only for new installation of R)
#-------------------------------------------------------------------------------
#source("http://www.bioconductor.org/biocLite.R") 
#biocLite("multtest")
#install.packages("gplots")
#install.packages("scatterplot3d")#The downloaded link at: http://cran.r-project.org/package=scatterplot3d

#Step 0: Import library and GAPIT functions run this section each time to start R)
#######################################################################################
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

#source("/Users/Zhiwu/Dropbox/Current/revolutionr/gapit/gapit_functions.txt")
#############################################################################################

#download tutorial data and save them in myGAPIT directory under C drive and  run tutorials
setwd("C:/Users/zacha/Documents/myGAPIT")


#Tutorial 1: Basic Scenario of Compressed MLM by Zhang and et. al. (Nature Genetics, 2010) 
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.delim("mdp_genotype_test.hmp.txt", head = FALSE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
  PCA.total=3,
)

#Tutorial 2: Using ECMLM by Li and et. al. (BMC Biology, 2014)
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  G=myG,
  PCA.total=3,
  kinship.cluster=c("average", "complete", "ward"),
  kinship.group=c("Mean", "Max"),
  group.from=200,
  group.to=1000000,
  group.by=10
)

#Tutorial 3: User defined Kinship and PCs
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)
myKI <- read.table("KSN.txt", head = FALSE)
myCV <- read.table("Copy of Q_First_Three_Principal_Components.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  G=myG,
  KI=myKI,
  CV=myCV,
)

#Tutorial 4: Genome Prediction
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myKI <- read.table("KSN.txt", head = FALSE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  G=myG,
  KI=myKI,
  PCA.total=3,
  SNP.test=FALSE
)

#Tutorial 5: Work with big data by spliting genotype Files
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  PCA.total=3,
  file.G="mdp_genotype_chr",
  file.Ext.G="hmp.txt",
  file.from=1,
  file.to=10
)

#Tutorial 6: Numeric Genotype Format
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myGD <- read.table("mdp_numeric.txt", head = TRUE)
myGM <- read.table("mdp_SNP_information.txt" , head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  GD=myGD,
  GM=myGM,
  PCA.total=3,
)

#Tutorial 7: Numerical Multiple Files
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  PCA.total=3,
  file.GD="mdp_numeric",
  file.GM="mdp_SNP_information",
  file.Ext.GD="txt",
  file.Ext.GM="txt",
  file.from=1,
  file.to=3,
  
)


#Tutorial 8: Improve speed of computing PC and kinship by using fractional snps 
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  PCA.total=3,
  file.GD="mdp_numeric",
  file.GM="mdp_SNP_information",
  file.Ext.GD="txt",
  file.Ext.GM="txt",
  file.from=1,
  file.to=3,
  SNP.fraction=0.6
)

#Tutorial 9: Reduce memory usage by loading fragment of file, one at a time, by defining fragment size (number of SNPs)
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  PCA.total=3,
  file.GD="mdp_numeric",
  file.GM="mdp_SNP_information",
  file.Ext.GD="txt",
  file.Ext.GM="txt",
  file.from=1,
  file.to=3,
  SNP.fraction=0.6,
  file.fragment = 128
)

#Tutorial 10: Optimization for number of PCs based on BIC
#The result is saved in GAPIT.TraitName.BIC.Model.Selection.Results.csv
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
  Y=myY[,1:2],
  G=myG,
  PCA.total=3,
  Model.selection = TRUE
)

#Tutorial 11: SUPER GWAS method by Wang and et. al. (PLoS One, 2014)
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myCV <- read.table("Copy of Q_First_Three_Principal_Components.txt", head = TRUE)
myY  <- read.table("mdp_traits.txt", head = TRUE)
myG <- read.table("mdp_genotype_test.hmp.txt" , head = FALSE)

#Step 2: Run GAPIT
myGAPIT_SUPER <- GAPIT(
  Y=myY[,1:2],			
  G=myG,				
  CV=myCV,
  #PCA.total=3,				
  sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1,
)


#Tutorial 12: Compare to Power against FDR for GLM,MLM,CMLM,ECMLM,SUPER(PLINK)
#Hint:Program runing time is more than 24 hours for repetition 100 times.
#Run description:Please refer to page 34,35 of the User manual on GAPIT website.
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myGD <-read.table("mdp_numeric.txt", head = TRUE)
myGM <-read.table("mdp_SNP_information.txt", head = TRUE)
myKI <- read.table("KSN.txt", head = FALSE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)

#Step 2: Run GAPIT
#GAPIT.Power.compare.plink
GAPIT.Power.compare(
  myG=myG,
  myGD=myGD,
  myGM=myGM,
  myKI=myKI,
  rel=100,
  h2=0.9,
  NQTN=5
)


#Tutorial 13: Genetic Prediction one time by cross validation
#-----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY<-read.table("mdp_traits.txt", head = TRUE)
myK<-read.table("KSN.txt", head = FALSE)

#Step 2: Run GAPIT
OnePre<-GAPIT.Prediction(
  myK=myK,
  y<-myY[,c(1,3)],
  ##y=y[,1:2],
  num=5
)


#Tutorial 14: Compare accuracy to different folds by replicate times
#-----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY<- read.table("mdp_traits.txt", head = TRUE)
myGD <-read.table("mdp_numeric.txt", head = TRUE)

#Step 2: Run GAPIT
myCross<-GAPIT.cross_validation.compare(
  myGD=myGD,
  y=myY,
  #y<-y[,c(1,3)],
  rel=100,
  tc<-c(2,5,10,20)  ##input compare to folds num
)


#Tutorial 15: Marker density and decade of linkage disequilibrium over distance
#-----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files

myGM <-read.table("mdp_SNP_information.txt", head = TRUE)
myGD <- read.table("mdp_numeric.txt", head = TRUE)

#Step 2: Run GAPIT

myGenotype<-GAPIT.Genotype.View(
  myGI=myGM,
  myGD=myGD[,-1],
  #chr=1,
  #w1_start=10,
  #w1_end=110,
  #mav1=10,
)

#Tutorial 16: Statistical distributions of phenotype
#-----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

myPhenotype<-GAPIT.Phenotype.View(
  myY=myY
)