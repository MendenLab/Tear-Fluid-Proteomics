# UNSUPERVISED PLOTS
library('ggplot2')
library('ggthemr')
ggthemr('greyscale', layout = "scientific")
library('tidyverse')
library('Rtsne')
library('umap')
library('scales')
library('pheatmap')
library('RColorBrewer')
library('ComplexHeatmap')

# if you are using Rstudio run the following command, otherwise, set the working directory to the folder where this script is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# create directory for results
dir.create(file.path(getwd(),'results'), showWarnings = FALSE)
# create directory for plots
dir.create(file.path(getwd(),'plots'), showWarnings = FALSE)

# LOAD DATA
# data frame 1: row names = patientID, protein expression values (columns)
# data frame 2: row names = patientID, patientinfo (columns), status column values: 1 and 0
data_frame = read.csv('data/mockdata.csv', row.names = 1)
df_patient_info = read.csv('data/patientinfo.csv', row.names = 1)

# normalise expression values (optional)
data_frame = apply(data_frame, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

# data frame needs to contain the protein expression values and patient status
df_ml = transform(merge(data_frame, df_patient_info['status'], by = 'row.names', all.x = T), row.names = Row.names, Row.names = NULL)

####
####
# define colours for plots
farben = c("#000000", "#B2182B")


## PCA --------
# set seed for reproducible results
set.seed(9)

# run pca function from base R
pca_patients = prcomp(df_ml[ , !names(df_ml) == "status"], scale = T)

# extract values to calculate explained Variance of each PC
pca_summary = summary(pca_patients)$importance
varExp = round(pca_summary[2,]*100,2)
# extract values and create data frame for plot
pca = data.frame(pca_patients$x)

# add status to data frame for plotting
if(all(rownames(pca) == rownames(df_ml))){
  pca$group = df_ml$status
}

# plot PCA
ggplot(pca, aes(x=PC1,y=PC2, colour = as.factor(group))) + 
  geom_point(size = 2) +
  scale_color_manual(name = "Group", labels = c("ctrl", "ALS"), values = farben) +
  xlab(paste("PC1 (",varExp[1],"%)")) + ylab(paste("PC2 (",varExp[2],"%)")) + ggtitle("PCA")
# save plot
ggsave("plots/pca.pdf", width = 11, height = 8, units = "in")


## t-SNE ------
# set seed for reproducible results
set.seed(9)

# run tsne function from Rtsne package
tsne_out = Rtsne(df_ml[ , !names(df_ml) == "status"])

# extract values to data frame and add status for plotting
tsne_plot = as.data.frame(tsne_out$Y)
tsne_plot$group = df_ml$status

# plot tSNE
ggplot(tsne_plot) + geom_point(aes(x=V1, y=V2, color = as.factor(group))) + ggtitle("tSNE Patients") +
  scale_color_manual(name = "Group", labels = c("ctrl", "ALS"), values = farben)
# save plot
ggsave("plots/tsne.pdf", width = 11, height = 8, units = "in")


## UMAP ---------
# set seed for reproducible results
set.seed(9)

# run umap function from umap package
umap_out = umap(df_ml[,!names(df_ml) == "status"])
# extract values for plot
umap_plot = as.data.frame(umap_out$layout)

# add status to data frame for plotting
if(all(rownames(umap_plot) == rownames(df_ml))){
  umap_plot$group = df_ml$status
}

# plot umap
ggplot(umap_plot) + geom_point(aes(x=V1, y=V2, color = as.factor(group))) +
  ggtitle("umap Patients") +
  scale_color_manual(name = "Group", labels = c("ctrl", "ALS"), values = farben)
# save umap
ggsave("plots/umap.pdf", width = 11, height = 8, units = "in")


## heatmap ---------
# change values of patient info to display control/ALS instead of binary values
df_patient_info[which(df_patient_info$status == 0), "status"] = "control"
df_patient_info[which(df_patient_info$status == 1), "status"] = "ALS"

# check same order of data and patient info
df_patient_info = df_patient_info[match(row.names(df_ml), row.names(df_patient_info)), ]

# get annotations ready
annotation = data.frame(group = as.factor(df_patient_info$status), sex = as.factor(df_patient_info$sex), age = as.numeric(df_patient_info$age))
rownames(annotation) = rownames(df_patient_info)
annotation_colours <- list(group = c(control = "#000000", ALS = "#B2182B"), sex = c(w="lightpink1",m="skyblue1"), age =c("white", "darkgreen"))


# save and plot heatmap
pdf("plots/heatmap.pdf", width = 11, height = 8)

ComplexHeatmap::pheatmap(as.matrix(t(df_ml[, !names(df_ml) == "status"])), name = "expression",cutree_cols = 2,
         show_colnames = T,
         show_rownames = FALSE,
         fontsize = 6,
         annotation_col = annotation,
         annotation_colors = annotation_colours,
         color = colorRampPalette(brewer_pal(palette = 'YlGnBu')(9))(100),
         main = "Heatmap")
dev.off()














