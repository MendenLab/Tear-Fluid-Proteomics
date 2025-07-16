# Analysis of tear fluid proteomics

Scripts to analyse proteomics data of ALS and control patients with unsupervised and supervised machine learning. Includes plotting of the results.

Analysis includes
* PCA
* UMAP
* tSNE
* heatmap
* linear regression model
* random forest model
* support vector machine model
* gene enrichment analysis

## Scripts
* unsupervised.R **includes PCA, UMAP, tSNE and heatmap analysis and visualisations**
* fgsea.R **script for Gene Set Enrichment Analysis**
* functions.R **script containing the functions that are used in fgsea.R and ml_github.R**
* ml_github.R **script containing the machine learning analysis with plotting of the results**
* aucCombinations.R ** script containing AUC calculations based on previous run linear regression models (see analysis.md for more informations)**
* data_completeness.R **script to analyse data completeness and calculate coefficient of variation of the raw proteomics data**

## Data
The data folder contains mockdata and a file with annotated KEGG pathways.
