# if you are using Rstudio run the following command, otherwise, set the working directory to the folder where this script is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# create directory for results
dir.create(file.path(getwd(),'results'), showWarnings = FALSE)
# create directory for plots
dir.create(file.path(getwd(),'plots'), showWarnings = FALSE)

# load the functions needed for the analysis
source('functions.R')

# LOAD DATA
# data frame 1: row names = patientID, protein expression values (columns)
# data frame 2: row names = patientID, patientinfo (columns), status column values: 1 and 0
data_frame = read.csv('data/mockdata.csv', row.names = 1)
df_patient_info = read.csv('data/patientinfo.csv', row.names = 1)

# normalise expression values (optional)
data_frame = apply(data_frame, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

# MACHINE LEARNING MODELS

# data frame for ML needs to contain the Protein expression values and patient status
df_ml = transform(merge(data_frame, df_patient_info['status'], by = 'row.names', all.x = T), row.names = Row.names, Row.names = NULL)

# LINEAR REGRESSION MODEL (LASSO)
# arguments:  data_frame: data as data frame with a status column,
#             algorithm: 'lm','rf','svm lin' or 'svm rad',
#             cv: cross validation folds (default = 10),
#             BS_number: number of bootstrap runs (default = 1)  
# returns: list of models, importance, predictions (probability), predictions (raw), indices (which samples had been used for training)
lm = runML(df_ml,'lm', BS_number = 500)

# save results for later use
saveRDS(lm, file = 'results/linearModel.rds')

# plot averaged ROC curve
# takes result object from runML(), the original data frame and path to save plot (optional)
# returns data frame with values for plot
ROC_curve = calculateROC(lm,df_ml,'plots/rocc_lm.pdf')

# save results for later use
write.csv(ROC_curve, file = 'results/rocc_lm.csv')

# extract weights + plot averaged
# arguments:  list_from_ML: results from runML()
#             plot_path:    path to save plot(optional)
#             number:       number of proteins on plot
# returns data frame with values for plot
lm_weights = plotWeights(lm,"plots/weights_lm.pdf")
# save results for later use
write.csv(lm_weights, file = 'results/weights_lm.csv')





# RANDOM FOREST
# arguments:  data_frame: data as data frame with a status column,
#             algorithm: 'lm','rf','svm lin' or 'svm rad',
#             cv: cross validation folds (default = 10),
#             BS_number: number of bootstrap runs (default = 1)
# returns: list of models, importance, predictions (probability), predictions (raw), indices (which samples had been used for training)
rf = runML(df_ml,'rf',BS_number = 500)

# save results for later use
saveRDS(rf, file = 'results/randomForest.rds')

# plot averaged ROC curve
# takes result object from runML(), the original data frame and path to save plot (optional)
# returns data frame with values for plot
ROC_curve_rf = calculateROC(rf,df_ml,'plots/rocc_rf.pdf')

# save results for later use
write.csv(ROC_curve_rf, file = 'results/rocc_rf.csv')

# extract weights + plot averaged
# arguments:  list_from_ML: results from runML()
#             plot_path:    path to save plot(optional)
#             number:       number of proteins on plot
# returns data frame with values for plot
rf_importance = plotWeights(rf,'plots/importance_rf.pdf')
# save results for later use
write.csv(rf_importance, file = 'results/importance_rf.csv')





# SUPORT VECTOR MACHINE
# arguments:  data_frame: data as data frame with a status column,
#             algorithm: 'lm','rf','svm lin' or 'svm rad',
#             cv: cross validation folds (default = 10),
#             BS_number: number of bootstrap runs (default = 1)  
# returns: list of models, importance, predictions (probability), predictions (raw), indices (which samples had been used for training)
svm_l = runML(df_ml,'svm lin', BS_number = 500)
svm_r = runML(df_ml,'svm rad', BS_number = 500)
# save results for alter use
saveRDS(svm_l, file = 'results/svm_lin.rds')
saveRDS(svm_r, file = 'results/svm_rad.rds')

# plot averaged ROC curve
# takes result object from runML(), the original data frame and path to save plot (optional)
# returns data frame with values for plot
ROC_curve_l = calculateROC(svm_l,df_ml,'plots/rocc_svml.pdf')
ROC_curve_r = calculateROC(svm_r,df_ml,'plots/rocc_svmr.pdf')

# save results for alter use
write.csv(ROC_curve_l, file = 'results/rocc_svml.csv')
write.csv(ROC_curve_r, file = 'results/rocc_svmr.csv')

# extract weights + plot averaged
# arguments:  list_from_ML: results from runML()
#             plot_path:    path to save plot(optional)
#             number:       number of proteins on plot
# returns data frame with values for plot
svm_weights = plotWeights(svm_l,"plots/weights_svm.pdf")
# save results for later use
write.csv(svm_weights, file = 'results/weights_svm.csv')

# 