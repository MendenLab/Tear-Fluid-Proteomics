### Gene Set Enrichment Analysis, with background filter

# if you are using Rstudio run the following command, otherwise, set the working directory to the folder where this script is in
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# create directory for results
dir.create(file.path(getwd(),'results'), showWarnings = FALSE)
# create directory for plots
dir.create(file.path(getwd(),'plots'), showWarnings = FALSE)

# load the functions needed for the analysis
source('functions.R')

# LOAD DATA
# load data for background
background_data = read.csv('data/mockdata.csv', row.names = 1)
# get background
background_data = names(background_data)

## get gene set to 
lm_genes = read.csv('results/weights_lm.csv', row.names = 1)
lm_genes = setNames(lm_genes$avg,row.names(lm_genes))

# run fgsea
# arguments:  genes: named (gene name) vector with values (either geneXexpression, or weights from ML models)
#             min_size: min number of genes to be found in a pathway for it to show up in results (default = 3)
#             plot_path: path to save plot (optional)
#             number: number of pathways included in the plot (default = 20)
#             entrezID: logical, if input data (gene names and background) is entrezid or gene name
#             background: vector with background genes (optional)
# returns # returns data frame with values for plot
gsea_result = gsea(lm_genes,background = background_data, min_size = 1, plot_path = 'plots/fgsea.pdf')



