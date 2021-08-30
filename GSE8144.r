################################################################################
## koot_GSE8144.r
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

# Dataset info avaliable in:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58144

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

# LOAD PACKAGES -----------------------------------------------------------
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

functions_path = '/Users/francisco/Desktop/TFM/functions'
dataset_path = '/Users/francisco/Desktop/TFM/datasets/GSE8144_koot'
data_path = paste0(dataset_path, '/data')
results_path = '/Users/francisco/Desktop/TFM/datasets/results'

source(paste0(functions_path, '/function_plot_PCAscores.r'))
source(paste0(functions_path, '/function_optimal_nclusters.r'))
source(paste0(functions_path, '/function_unsupervised_clustering.r'))

getwd()
setwd(dataset_path)

################################################################################
# load(paste0(data_path, '/koot_data.rda'), verbose = T)
# Loading objects:
# ed
# koot_data_A
# koot_data_M
# data_batchcorrect_A
# data_batchcorrect_M
# data_time_batch_corrected_A
# data_time_batch_corrected
################################################################################

# READ PLATFORM FILE ------------------------------------------------------
# Platform = A-UMCU-HS44K-2.0
gpl_name = '/GPL15789.txt'

gpl = read.delim(file = paste0(data_path, gpl_name),
                  header = T,
                  sep = '\t',
                  comment.char = '#',
                  skip = 1,
                  quote = '',
                 stringsAsFactors = F)

dim(gpl) # [1] 45220    15

# READ EXPERIMENTAL DESSIGN -----------------------------------------------
# series_name = '/GSE58144_series_matrix.txt'
# series_file = paste0(data_path, series_name)
#
# # Read characteristics
# con = file(series_file, 'r') # open file
# characteristics = c() # prepare empty vector to save data
# while(TRUE) {
#   line = readLines(con, n=1)
#   if(length(line) == 0) {
#     break
#   } else if(startsWith(line, '!Sample_title')) {
#     titles = unlist(strsplit(line, '\t'))[-1]
#     titles = gsub('\\\'', '', titles)
#   } else if(startsWith(line, '!Sample_characteristics')) {
#     characteristics = c(characteristics, line)
#   } else if(startsWith(line, '!Sample_geo_accession')) {
#     accession = unlist(strsplit(line, '\t'))[-1]
#     accession = gsub('\\\'', '', accession)
#   }
# }
# close(con) # closes file
# # Now we parse the info:
# ed = data.frame(lapply(characteristics, function(x) {
#   values = unlist(strsplit(x, '\t'))[-1]
#   values = gsub('\\\'', '', values)
#   parts = strsplit(values, ': ')
#
#   name = parts[[1]][[1]]
#   values = sapply(parts, function(x) x[2])
#
#   out = list()
#   out[[name]] = values
#   return(out)
# }))
#
# ed = data.table(sample = accession, title = titles, ed)
# rownames(ed) = ed$sample
#
# type_column = gsub(pattern = 'control_.*',
#                    replacement = 'control',
#                    x = ed$title)
#
# type_column = gsub(pattern = 'recurrent implantation failure after IVF_.*',
#                    replacement = 'RIF',
#                    x = type_column)
#
# ed = data.frame(sample = ed$sample, type = type_column, tissue = ed$tissue,
#                 batch = ed$batch, biopsy_time = ed$time.of.biopsy,
#                 age = ed$age)

################################################################################
# OBJECT ED IS ALREADY CREATED & LOADED

# To Homogenize between dataframes:
ed$condition[ed$Condition == 'CONTROL'] = 'control'
ed$condition[ed$Condition == 'RIF'] = 'RIF'
ed$sample = rownames(ed)

ed = data.frame(sample = ed$sample, condition = ed$condition,
                tissue = ed$Tissue, time = ed$Time, batch = ed$Batch)
rownames(ed) = ed$sample
################################################################################

# DOWNLOAD RAW DATA FOR KOOT (GSE8144) ----------------------------------
# getGEOSuppFiles('GSE58144', baseDir = './data/')
# if 'error: Timeout of 60 seconds was reached' then 'options(timeout = 1000)'.
# Downloaded in '/raw_datasets/altmae_GSE71331/data/GSE71331'
# (Remember to unzip it and prepare paths)

# READ EXPRESSION DATA ----------------------------------------------------
# my_files = list.files(path = 'data/GSE58144_RAW', full.names = T)
# target_info = data.frame(FileName = my_files,
#                           RIF_CONTROL = ed[, 2],
#                           stringsAsFactors = F)

# ...

# SAVE RAW DATA, ED AND GPL ------------------------------------------------
save(ed, file = paste0(data_path, '/ed_koot.rda'), version = 2)
save(gpl, file = paste0(data_path, '/gpl4133_koot.rda'), version = 2)
# save(gse8144raw, file = paste0(data_path, '/gse8144raw.rda'), version = 2)

# This will be our dat (M & A after different channels in chip) and ed:-
dat = data_time_batch_corrected_M
save(dat, ed, file = paste0(results_path, '/koot.rda'), version = 2)

load(paste0(results_path, '/koot.rda'), verbose = T)

# NORMALIZED PLOTS --------------------------------------------------------
plot_PCAscores(dat = dat, ed = ed, components = c(1, 2),
               condition1 = 'condition', colors = c('#1B9E77', '#D95F02'),
               title = 'PCA: Koot (GSE8144)')

plot_PCAscores(dat = dat, ed = ed, components = c(3, 4),
               condition1 = 'condition', colors = c('#1B9E77', '#D95F02'),
               title = 'PCA: Koot (GSE8144)')

# Previously processed dataset, but just in case we perform aqm:
# OUTLIER DETECTION (ARRAY QUALITY METRICS) -------------------------------
outdir = paste0(dataset_path, '/arrayQuality_report')

# We need to create an koot_eset objet:
koot_eset = ExpressionSet(assayData = assayDataNew(exprs = dat))

# Now we can use aqm:
arrayQualityMetrics(expressionset = koot_eset,
                    outdir = outdir,
                    force = TRUE)

# Check index.html file in outdir for results.

################################################################################
# # Plot normalized data
# toplot = melt(data_time_batch_corrected_M)
# ed$sample = rownames(ed)
# toplot2 = merge(x = toplot,y = ed[,c('sample','condition')],
#                  by.x = 'Var2', by.y = 'sample')
# ggplot(toplot2, aes(x = Var2, y = value, fill = condition)) +
#        geom_boxplot() +
#        ggtitle('Quantile normalized expression data') +
#        scale_fill_manual(values = c(brewer.pal(8, 'Dark2'))) +
#        theme_light() +
#        xlab('') +
#        ylab('') +
#        theme(plot.title = element_text(size=35), legend.position = 'bottom',
#              axis.text.y = element_text(size = 25, color = 'darkgrey'),
#              axis.text.x = element_text(angle = 90),
#              axis.ticks.x = element_blank(),
#              legend.text = element_text(size = 12),
#              legend.title = element_blank())
#
# plot_PCAscores(dat = data_batchcorrect_M, ed = ed,
#                components = c(1,2),condition1 = 'time',
#                colors = brewer.pal(8,'Dark2'))
#
# plot_PCAscores(dat = data_time_batch_corrected_M, ed = ed,
#                components = c(1,2),condition1 = 'type',
#                colors = brewer.pal(8,'Dark2'))
