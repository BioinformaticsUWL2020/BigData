# BigData

2020 Group project to analyse gene data and determine likelihood of a mutant phenotype developing based predictive modelling of known collected data.
Will use APIs to pull known genetic data from databanks such as Ensembl to train a model and test data from openSNP.org to feed into the model and predict genetic predisposition.

Members: Ellie Sparling, Zach Sykes, Andrew Gyesi, Josh Russell, Judith Mbuyi

## Project Outline:

Predicticting likelihood of having/developing nine target diseases:
Parkinson's Diseasse, Late-Onset Alzheimer's Disease, Hereditary Thrombophilia, Celiac's Disease, Alpha-1-Antitrypsin Deficiency, Early-Onset Primary Dystonia, Glucose-6-Phosphate Dehydrogenase Deficiency (G6PD), Hereditary Hemochromatosis, and Factor XI Deficiency

This study plans to do this by completing the 5 phases below.

## Phase 1:

- Create an API using R to pull data from Ensembl database (Training data)


- Create an API to pull genotype data from OpenSNP users

  - First set of data pulled will be genotype data from users with OpenSNP designated phenotype ids for one of the nine target diseases
  - This initial round of data will be used to further test the model can accurately predict disease when passing data with unknown      phenotypes

  - Second set of data will be genotype data from randomly selected OpenSNP users

## Phase 2:

- Clean data and ensure the observations for each set matches (except for data labels) to prepare for Phase 

## Phase 3:

- Begin construction of a supervised learning model (either KNN or Classification Trees)
- Train the model with preped labelled data
- Feed preped OpenSNP genotype data with known phenotype corresponding to each target disease to test predictive accuray of model
  - If predicitive accuracy from the Ensembl data not successful (will just use known phenotype OpenSNP data to train the model)
- Feed preped OpenSNP genotype data with unknown phenotypes to test the model's ability to predict disease likelihood
- Pass the results onto Phase 5 

## Phase 4:

- Re-analyse the data but instead of passing data straight to a predictive model, data will be processed using the R package sparklyr (making use of an Apache Spark cluster) before being passed through the model
- Data output will be sent to Phase 5

## Phase 5:

- Predictions from both Phase 3 and 4 will be plotted and graphed (using R packages GGplot and Shiny) to see where certain users were placed based on genotype data fed into model
- All outputs from both models (just ML and ML with Apache Spark processing) will be plotted against each other to highlight the benefits of processing the data with Apache Spark