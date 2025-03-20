# Analysis

## Data preparation
Data is in the form of rows corresponding to samples and columns to proteins. The values are protein expression values. You will also need a column named "status", which encodes (with 0 and 1) the groups that the samples belong to. 0 corresponds to control, 1 to diseased patient.
Make sure the protein names are either already entrezIDs or can all be converted to entrezIDs using the **org.Hs.eg.db** R package with the following code.
```
ids = AnnotationDbi::select(org.Hs.eg.db, 
                 keys = VECTOR_OF_YOUR_PROTEIN_NAMES,
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = 'SYMBOL')
```
then check if any are unmapped with `ids$SYMBOL[which(is.na(ids$ENTREZID))]`

## Unsupervised Visualisations
For a quick overview 
* [PCA](/unsupervised.R?plain=1#L38)
<img src="/plots/pca.png" width="400">

* [UMAP](/unsupervised.R?plain=1#L83)
 <img src="/plots/umap.png" width="400">
 
* [tSNE](/unsupervised.R?plain=1#L65)
<img src="/plots/tsne.png" width="400">

## Heatmap
The code for the heatmap can also be found in the script for [unsupervised](/unsupervised.R) analysis. The vector *annotation* and *annotation_colours* can be changed to display the desired attributes. The option *cutree_cols* in pheatmap() can be changed to better display groupings in your own data.
<img src="/plots/heatmap.png" width="500">

## Supervised Machine Learning
All 3 machine learning alogrithms can be run with the runML() function from the [functions](/functions.R) script. Usage of the function and analysis including plotting of results and extracting feature weight/importance where possible, can be found in the [ml_github.R](/ml_github.R) script.
The functions takes the following arguments
* data_frame: data as data frame with a status column,
* algorithm: 'lm', 'rf', 'svm lin' or 'svm rad',
* cv: cross validation folds (default = 10),
* BS_number: number of bootstrap runs (default = 1)
 and returns a list of models, importance, predictions (probability), predictions (raw), indices (which samples had been used for training)
 
### Linear Regression (lasso)
Running rumML() with *algorithm = 'lm'* runs a linear regression using the lasso algorithm (to change this change the [responsible variable](/functions.R?plain=1#L79) in functions.R).
the results can be parsed into [calculateROC()](functions.R?plain=1#L133) to plot the averaged ROC curve. 

<img src="/plots/rocc_lm.png" width="400">

For the linear model the weights can be extracted and plotted using [plotWeights()](functions.R?plain=1#L192).

### Random Forest
For the random forest analysis the runML() function should be run as indicated in [ml_girhub.R](/ml_github.R?plain=1#L65). Using [calculateROC()](functions.R?plain=1#L133) and [plotWeights()](functions.R?plain=1#L192) the results can be analysed by plotting feature importance and averaged ROC curve.

<img src="/plots/importance_rf.png" width="400">

### Support Vector Machine
The [ml_github.R](/ml_github.R) script also includes analysis with the SVM algorithm both with the radial and linear kernel, however this can be extended by changing the [responsible variable](/functions.R?plain=1#L104) in functions.R.
Besided plotting the ROC curve, the feature weights from the linear kernel can also be plotted using [plotWeights()](functions.R?plain=1#L192).

## Gene Set Enrichment Analysis
The here documented gene enrichment analysis uses the **fgsea** package in R. It can be used to run a GSEA with or without adjusted background. The relevant genes have to be given as a numeric, named vector. The numbers can either be weights from previous analysis or kept at 1 for every gene to not have a weighted GSEA. The function uses the KEGG pathways but this can be changed in the [functions](/functions.R?plain=1#L331) script.
Pathways can be downloaded from https://www.pathwaycommons.org/pc2/datasources

The analysis is as easy as running [gsea()](fgsea.R) and also plots the results, if a plot path is given.

<img src="/plots/fgsea.png" width="500">

## AUC combinations
Since we know the equation that is the basis for the linear regression model we can calculate which group a sample belongs to with any chosen combination of features, if we already know their weights from previous analysis. For this there is a function in [functions](/functions.R?plain=1#L422) script. [auccurve()](aucCombinations.R?plain=1#L67) returns a data frame which lets you plot

* The best combinations out of all

<img src="/plots/AUCcomb.png" width="500">

* The change in AUC when adding the proteins of interes one by one

<img src="/plots/AUCadd.png" width="500">

* The "impact"/auc of single proteins

<img src="/plots/AUCbar.png" width="500">
