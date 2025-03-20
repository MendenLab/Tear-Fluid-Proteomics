# if you are using Rstudio run the following command, otherwise, set the working directory to the folder where this script is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# create directory for results
dir.create(file.path(getwd(),'results'), showWarnings = FALSE)
# create directory for plots
dir.create(file.path(getwd(),'plots'), showWarnings = FALSE)

# load the functions needed for the analysis
source('functions.R')


# LOAD DATA
# data frame 1: row names = patientID, protein expression values (columns)             here the data frame is of the data you want to run the model on (could be validation data)
# data frame 2: row names = patientID, patientinfo (columns), status column values: 1 and 0
data_frame = read.csv('data/mockdata.csv', row.names = 1)
df_patient_info = read.csv('data/patientinfo.csv', row.names = 1)

# normalise expression values (optional)
data_frame = apply(data_frame, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))


# bring data in same sample order
all_data = merge(data_frame, df_patient_info[ ,'status', drop = FALSE], by='row.names')
row.names(all_data) = all_data$Row.names 
all_data$Row.names = NULL

###
# create matrix with average proteins as weight
###

# load averaged weights
avgweights = read.csv('results/weights_lm.csv', row.names = 1)

# filter out proteins that were only picked in 1 run
# remove the ones where error > avg & picked once
avgweights = avgweights[which(avgweights$picks > 1), ]
avgweights = avgweights[which(abs(avgweights$avg)>abs(avgweights$error)), ]
avgweights = as.data.frame(avgweights)

# create data.frame with rows corresponding to all available proteins. Where possible average weights of proteins are stored, all others are 0
averageweights_m = matrix(nrow=length(names(all_data[ ,!names(all_data) == "status"])), ncol = 1)
averageweights_m = as.data.frame(averageweights_m)
averageweights_m[is.na(averageweights_m)] = 0
rownames(averageweights_m) = names(all_data[ ,!names(all_data) == "status"])
for(i in 1:nrow(averageweights_m)){
  for(j in 1:nrow(avgweights)){
    if(row.names(averageweights_m)[i] == row.names(avgweights)[j]){
      averageweights_m[i,"V1"] = avgweights[j, "avg"]
    }
  }
}



#### genes of interest
genes = c("GPS1", "ADRM1", "LAMA2", "NAP1L1", "CCT8", "PPP1R7", "RAP2A", "SBF1", "PGP")

# Get data frame for plotting with auc for different combinations
# arguments:    vectornames: character vector of genes of interest
#               weight_data: data frame as explained above (averageweights_m)
#               protein_data: data to calculate the auc on (needs protein columns and status column)
#               maxn: maximal length of combinations (default = length of vectornames)
#               add: boolean value if combinations should be built, or just cummulative addition of the genes in vectornames
# returns: data frame with combinations and auc column
curvy = auccurve(genes, averageweights_m,all_data)

# select colours for bar plot
farbe = c("#D17480","#C95D6B", "#C14655", "#BA2F40", "#B2182B", "#A01627",	"#8E1322", "#7D111E","#6B0E1A")

# plot overall best combinations
# number of combinations to plot
n = 9
ggplot(curvy[order(curvy$auc, decreasing = T)[1:n], ], aes(x = reorder(paste0(combinations),auc) , y = auc,group = 1)) + geom_line() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, face = 'bold', size = 11))+
  scale_x_discrete(labels = wrap_format(10))+
  ylab("AUC") +
  xlab("Genes") +
  ggtitle("AUC of combinations")
# save AUC bar  plot
ggsave("plots/AUCcomb.pdf", width = 11, height = 8, units = "in")


# bar plot to show single protein "IMPACT"
# number of combinations to plot
n = 9
ggplot(curvy[c(1:n), ], aes(x = reorder(paste0(combinations),auc), y = auc, fill = reorder(paste0(combinations),auc))) + geom_bar(stat = "identity") +
  ylab("AUC")  +theme(legend.position = "none", axis.text.x = element_text(face = 'bold', size = 11)) +
  scale_fill_manual(values = farbe) +
  xlab("Genes") + scale_y_continuous(breaks = c(0,0.2,0.4,seq(0.5,0.7,0.05)),limits = c(0,0.7))+
  ggtitle("AUC of single protein models")
# save AUC bar  plot
ggsave("plots/AUCbar.pdf", width = 11, height = 8, units = "in")


#plot of cumulative
# run function with 'add = TRUE'
curvy_add = auccurve(genes, averageweights_m,all_data, add = TRUE)

ggplot(curvy_add, aes(x = paste0(combinations), y = auc,group = 1)) + geom_line() +
  theme(axis.text.x = element_text(face = 'bold', size = 11))+
  scale_x_discrete(labels = wrap_format(10))+
  #scale_y_continuous(breaks = seq(0.55,0.61,0.005), limits = c(0.55,0.61) )+ # eeeh
  ylab("AUC") +
  xlab("Genes") + 
  ggtitle("AUC of models")
# save adding proteins plot
ggsave("plots/AUCadd.pdf", width = 11, height = 8, units = "in")

