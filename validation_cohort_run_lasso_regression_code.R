library('ggplot2')
library('tidyverse')
library("viridis")
library("dplyr")
library('caret')
library('glmnet')
library('pROC')
library('matrixStats') # row standard deviation
library('naniar')
library('impute')
library('readxl')
library('Hmisc')#install.packages('acepack')
library('shapr')
library('ComplexHeatmap')
library(ggnewscale)
library(GGally)
library(imputeLCMD)
library(ggpubr)

##Setup
# if you are using Rstudio run the following command, otherwise, set the working directory to the folder where this script is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Create folders if they do not exist
folders <- c("validation_plots", "validation_results")
for (folder in folders) {
  if (!dir.exists(folder)) {
    dir.create(folder)
    message(paste("Created folder:", folder))
  } else {
    message(paste("Folder already exists:", folder))
  }
}

##Load dataset---------------------------
datasets = list()
datasets[[1]] = read_excel("data/240801-WB results VC.xlsx")
names(datasets) = c("RD_raw")
#Change column names
names(datasets[[1]])[names(datasets[[1]]) == "Kriterium"] = "status"
names(datasets[[1]])[names(datasets[[1]]) == "Probe"] = "patid"

datasets[[1]]$status = as.factor(datasets[[1]]$status)
#Convert protein values into numeric 
datasets[[1]][, 3:ncol(datasets[[1]])] = apply(datasets[[1]][, 3:ncol(datasets[[1]])], 2, as.numeric)

#--transform----
protein_start_col = 3
RD_log10 = datasets[["RD_raw"]]
RD_log10[, protein_start_col:ncol(RD_log10)] = log10(RD_log10[, protein_start_col:ncol(RD_log10)])
RD_log10[RD_log10 == -Inf] <- 0
datasets = c(datasets, list(RD_log10))
names(datasets)[length(datasets)] = c("RD_log10")

#--imputation----
d = datasets[[2]][, protein_start_col:ncol(datasets[[2]])]
set.seed(1234)  # Adjust the formula to create distinct seeds
imp = impute.MinDet(d, q = 0.01)
k = length(datasets) + 1
datasets[[k]] = cbind(datasets[[1]][, 1:2], imp)

# Rename the new dataset with "_mindet" suffix and iteration number
names(datasets)[k] = paste0(names(datasets)[1], "_mindet")


##Lasso regression---------------------------
#use log10_mindet data
source("Lasso_regression_code.R")
outputFolder = "validation_plots"
d = datasets[[3]]
d_name = names(datasets)[3]
d$status = (as.numeric(d$status))-1 
d = d[,!names(d) == "patid"]
original_row_names <- rownames(d)
drug = as.vector(d[,"status"]) 
feat = as.matrix(d[,!names(d) == "status"]) 
names(drug) <- original_row_names
rownames(feat) <- original_row_names
checkRandom(feat, drug, "lasso", paste0(outputFolder, "/009_lasso_check_random_",d_name,".pdf"))
cvLasso <- runCV(feat, drug, "lasso", paste0(outputFolder,"/010_lasso_cross_validation_",d_name, ".pdf"))
wLasso <- reTrainModel(feat, drug, cvLasso, paste0(outputFolder,"/011_lasso_plot_weights_",d_name, ".pdf"))
plotWeights(feat, drug, wLasso, paste0(outputFolder,"/012_lasso_plot_weights_zoomIn_",d_name,".pdf"))

