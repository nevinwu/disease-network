################################################################################
## bastu_GSE111974.r
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

# Dataset info avaliable in:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111974

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

# LOAD PACKAGES -----------------------------------------------------------
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(sva)

functions_path = '/Users/francisco/Desktop/TFM/functions'
dataset_path = '/Users/francisco/Desktop/TFM/datasets/GSE111974_bastu'
data_path = paste0(dataset_path, '/data')
results_path = '/Users/francisco/Desktop/TFM/datasets/results_sva_2'

source(paste0(functions_path, '/function_plot_PCAscores.r'))
source(paste0(functions_path, '/function_optimal_nclusters.r'))
source(paste0(functions_path, '/function_unsupervised_clustering.r'))

getwd()
setwd(dataset_path)

# READ PLATFORM FILE ------------------------------------------------------
# Platform = Agilent-014850
gpl_name = '/GPL17077.txt'

gpl = read.delim(file = paste0(data_path, gpl_name),
                  header = T,
                  sep = '\t',
                  comment.char = '#',
                  skip = 1,
                  quote = '',
                  stringsAsFactors = F)

dim(gpl) # [1] 50739     16

# READ EXPERIMENTAL DESSIGN -----------------------------------------------
series_name = '/GSE111974_series_matrix.txt'
series_file = paste0(data_path, series_name)

# Read characteristics
con = file(series_file, 'r') # open file
characteristics = c() # prepare empty vector to save data
while(TRUE) {
  line = readLines(con, n=1)
  if(length(line) == 0) {
    break
  } else if(startsWith(line, '!Sample_title')) {
    titles = unlist(strsplit(line, '\t'))[-1]
    titles = gsub('\\\'', '', titles)
  } else if(startsWith(line, '!Sample_characteristics')) {
    characteristics = c(characteristics, line)
  } else if(startsWith(line, '!Sample_geo_accession')) {
    accession = unlist(strsplit(line, '\t'))[-1]
    accession = gsub('\\\'', '', accession)
  }
}
close(con) # closes file

# Now we parse the info:
ed = data.frame(lapply(characteristics, function(x) {
  values = unlist(strsplit(x, '\t'))[-1]
  values = gsub('\\\'', '', values)
  parts = strsplit(values, ': ')

  name = parts[[1]][[1]]
  values = sapply(parts, function(x) x[2])

  out = list()
  out[[name]] = values
  return(out)
}))

ed = data.table(sample = accession, title = titles, ed)

# To Homogenize between dataframes:
ed = data.frame(sample = ed$sample, condition = ed$title, tissue = ed$tissue)

ed$tissue[ed$tissue == 'Endometrial tissue'] = 'endometrial tissue'
ed$condition[1:24] = 'control'
ed$condition[25:48] = 'RIF'
rownames(ed) = ed$sample

# DOWNLOAD RAW DATA FOR BASTU (GSE111974) ----------------------------------
getGEOSuppFiles('GSE111974', baseDir = './data/')
# if 'error: Timeout of 60 seconds was reached' then 'options(timeout = 300)'.
# Downloaded in '/raw_datasets/altmae_GSE71331/data/GSE71331'
# (Remember to unzip it and prepare paths)

# READ EXPRESSION DATA ----------------------------------------------------
my_files = list.files(path = 'data/GSE111974_RAW', full.names = T)
target_info = data.frame(FileName = my_files,
                          RIF_CONTROL = ed[,2],
                          stringsAsFactors = F)

gse111974raw = read.maimages(files = target_info,
                      source = 'agilent.median',
                      green.only = T)
# Use View(gse111974raw$E) and check with files that parsing is correct.

# SAVE RAW DATA, ED AND GPL ------------------------------------------------
save(ed, file = paste0(data_path, '/ed_bastu.rda'), version = 2)
save(gpl, file = paste0(data_path, '/gpl4133_bastu.rda'), version = 2)
save(gse111974raw, file = paste0(data_path, '/gse111974raw.rda'), version = 2)

# PRE-NORMALIZATION ANALYSIS ----------------------------------------------
load(paste0(data_path, '/ed_bastu.rda'), verbose = T)
load(paste0(data_path, '/GPL4133_bastu.rda'), verbose = T)
load(paste0(data_path, '/gse111974raw.rda'), verbose = T)

# RAW PLOTS ---------------------------------------------------------------
#### MDplot ####
plotMD(gse111974raw, column = 1,
       main = 'MD plot Bastu (GSM3045867): raw control-1')
plotMD(gse111974raw, column = 43,
       main = 'MD plot Bastu (GSM3045909): RIF')
# Warning message:
#In plotMD.EListRaw(gse111974raw, column = 1, main = 'MD plot Altmae (GSE414976):
# raw control-1') : NaNs produced

#### Boxplot ####
# boxplot(data.frame(log2(gse111974raw$Eb)), main = 'Green background')
boxplot(data.frame(log2(gse111974raw$E)), main = 'Raw data')

# RAW PCA -----------------------------------------------------------------
pca_raw = prcomp(t(log2(gse111974raw$E) + 1))
var_raw = round(summary(pca_raw)$importance[2, c(1,2)] * 100, 1)
toplot_raw = data.frame(pca_raw$x[, c(1,2)], stringsAsFactors = F)
lim_raw = max(abs(c(min(toplot_raw), max(toplot_raw))))

axis_limits_raw = c(-lim_raw, lim_raw)
# toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
toplot_raw$color = ed$condition

ggplot(data = toplot_raw,
       aes_string(x = colnames(toplot_raw)[1], y = colnames(toplot_raw)[2])) +
       geom_point(aes(color = color), size = 3) +
#       geom_text_repel(label = ed$sample, size = 3) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
       xlab(paste0('PC1', ': ', var_raw[1], '%')) +
       ylab(paste0('PC2', ': ', var_raw[2], '%')) +
       ggtitle('PCA: Bastu (GSE111974)') +
       xlim(axis_limits_raw) + ylim(axis_limits_raw) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))

# OUTLIER DETECTION (ARRAY QUALITY METRICS) -------------------------------
outdir = paste0(dataset_path, '/arrayQuality_report')

# We need to create an bastu_eset objet:
bastu_eset = ExpressionSet(assayData = assayDataNew(exprs = gse111974raw$E))

# Now we can use aqm:
arrayQualityMetrics(expressionset = bastu_eset,
                    outdir = outdir,
                    force = TRUE,
                    do.logtransform = TRUE) # Since data is not processed yet.

# Check index.html file in outdir for results.

# REMOVE OUTLIERS: --------------------------------------------------------
# We decide to remove 4 samples as they are probably outliers: GSM3045908,
# GSM3045909, GSM30459010, GSM3045911

# Remove GSM3045908:
gse111974raw$E = gse111974raw$E[, colnames(gse111974raw$E) != 'data/GSE111974_RAW/GSM3045908_SG12324209_253949427345_S001_GE1_1100_Jul11_1_2']
gse111974raw$Eb = gse111974raw$Eb[, colnames(gse111974raw$Eb) != 'data/GSE111974_RAW/GSM3045908_SG12324209_253949427345_S001_GE1_1100_Jul11_1_2']
gse111974raw$targets = gse111974raw$targets[, colnames(gse111974raw$targets) != 'data/GSE111974_RAW/GSM3045908_SG12324209_253949427345_S001_GE1_1100_Jul11_1_2']
ed = ed[rownames(ed) != 'GSM3045908',]

# Remove GSM3045909:
gse111974raw$E = gse111974raw$E[, colnames(gse111974raw$E) != 'data/GSE111974_RAW/GSM3045909_SG12324209_253949427345_S001_GE1_1100_Jul11_1_3']
gse111974raw$Eb = gse111974raw$Eb[, colnames(gse111974raw$Eb) != 'data/GSE111974_RAW/GSM3045909_SG12324209_253949427345_S001_GE1_1100_Jul11_1_3']
gse111974raw$targets = gse111974raw$targets[, colnames(gse111974raw$targets) != 'data/GSE111974_RAW/GSM3045909_SG12324209_253949427345_S001_GE1_1100_Jul11_1_3']
ed = ed[rownames(ed) != 'GSM3045909',]

# Remove GSM3045910:
gse111974raw$E = gse111974raw$E[, colnames(gse111974raw$E) != 'data/GSE111974_RAW/GSM3045910_SG12324209_253949427345_S001_GE1_1100_Jul11_1_4']
gse111974raw$Eb = gse111974raw$Eb[, colnames(gse111974raw$Eb) != 'data/GSE111974_RAW/GSM3045910_SG12324209_253949427345_S001_GE1_1100_Jul11_1_4']
gse111974raw$targets = gse111974raw$targets[, colnames(gse111974raw$targets) != 'data/GSE111974_RAW/GSM3045910_SG12324209_253949427345_S001_GE1_1100_Jul11_1_4']
ed = ed[rownames(ed) != 'GSM3045910',]

# Remove GSM3045911:
gse111974raw$E = gse111974raw$E[, colnames(gse111974raw$E) != 'data/GSE111974_RAW/GSM3045911_SG12324209_253949427345_S001_GE1_1100_Jul11_2_1']
gse111974raw$Eb = gse111974raw$Eb[, colnames(gse111974raw$Eb) != 'data/GSE111974_RAW/GSM3045911_SG12324209_253949427345_S001_GE1_1100_Jul11_2_1']
gse111974raw$targets = gse111974raw$targets[, colnames(gse111974raw$targets) != 'data/GSE111974_RAW/GSM3045911_SG12324209_253949427345_S001_GE1_1100_Jul11_2_1']
ed = ed[rownames(ed) != 'GSM3045911',]

# Check raw PCA:
pca_raw = prcomp(t(log2(gse111974raw$E) + 1))
var_raw = round(summary(pca_raw)$importance[2, c(1,2)] * 100, 1)
toplot_raw = data.frame(pca_raw$x[, c(1,2)], stringsAsFactors = F)
lim_raw = max(abs(c(min(toplot_raw), max(toplot_raw))))

axis_limits_raw = c(-lim_raw, lim_raw)
# toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
toplot_raw$color = ed$condition

ggplot(data = toplot_raw,
       aes_string(x = colnames(toplot_raw)[1], y = colnames(toplot_raw)[2])) +
       geom_point(aes(color = color), size = 3) +
       geom_text_repel(label = ed$sample, size = 3) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
       xlab(paste0('PC1', ': ', var_raw[1], '%')) +
       ylab(paste0('PC2', ': ', var_raw[2], '%')) +
       ggtitle('PCA: Bastu (GSE111974)') +
       xlim(axis_limits_raw) + ylim(axis_limits_raw) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))


# CORRECTING BACKGROUND ---------------------------------------------------
gse111974 = backgroundCorrect(gse111974raw, method = 'normexp', offset = 50)
# check if offset needed, we use 50 as default

# ANNOTATION --------------------------------------------------------------
eset = as.matrix(gse111974$E) # expression info
rownames(eset) = gse111974$genes$ProbeName

#### Filtering probes: controls and NAs ####
probesInfo = data.frame(gse111974$genes$ProbeName,
                        gse111974$genes$SystematicName,
                        gse111974$genes$ControlType,
                        stringsAsFactors = F)

probesInfo_aux = left_join(x = probesInfo,
                           y = gpl,
                           by = c('gse111974.genes.ProbeName' ='ID'))

dim(probesInfo) # 62976     3
dim(probesInfo_aux) # 62976    18
all(probesInfo$gse111974.genes.ProbeName == probesInfo_aux$gse111974.genes.ProbeName)

# Remove controls
indexestorm2 = which(probesInfo_aux$CONTROL_TYPE != FALSE)
# indexestorm2 = probesInfo_aux[probesInfo_aux$CONTROL_TYPE == FALSE,
#                               'gse111974.genes.ProbeName']
exprFilt = eset[-indexestorm2, ]
# exprFilt2 = eset[indexestorm2, ]


# Check number of genes
probesInfo2 = probesInfo_aux[probesInfo_aux$CONTROL_TYPE == FALSE,]
all(probesInfo2$gse111974.genes.ProbeName == rownames(exprFilt))

# GROUP PROBESETS INFORMATION ---------------------------------------------
# Condense replicate probes by their average
exprbyprobe = avereps(exprFilt,ID=probesInfo2$gse111974.genes.ProbeName)

which(rownames(exprbyprobe) == '') # 0
which(is.na(rownames(exprbyprobe))) # 0

# GROUP GENES BY PROBESET ID ----------------------------------------------
exprbygene = avereps(exprFilt, ID=probesInfo2$GENE_SYMBOL)
dim(exprbygene) # 32.081 genes
colnames(exprbygene) = ed$sample

all(rownames(ed) == colnames(exprbygene))

# PLOT EXPRESSION BY GENE -------------------------------------------------
#### Raw expression ####
toplot1 = melt(exprbygene)
p1 = ggplot(toplot1, aes(x = Var2, y = value)) +
            geom_boxplot() + ggtitle('Gene raw expression data') + xlab('') +
            ylab('') + theme_bw() +
            theme(plot.title = element_text(size = 35),
                  legend.position = 'none',
                  axis.text.y = element_text(size = 25, color = 'darkgrey'),
                  axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                  legend.text = element_text(size = 25),
                  legend.title = element_blank())

# APPLY LOG2 TRANSFORMATION -----------------------------------------------
expr_log2 = log2(exprbygene)
colnames(expr_log2) = ed$sample

#### Plot log2 data ####
toplot2 = melt(expr_log2)
p2 = ggplot(toplot2, aes(x = Var2, y = value)) + geom_boxplot() +
     ggtitle('Log2 expression data') + xlab('') + ylab('') + theme_bw() +
     theme(plot.title = element_text(size = 35), legend.position = 'none',
           axis.text.y = element_text(size = 25, color = 'darkgrey'),
           axis.text.x = element_blank(), axis.ticks.x = element_blank(),
           legend.text = element_text(size = 25),
           legend.title = element_blank())

# QUANTILE NORMALIZATION --------------------------------------------------
# Final step of normalization
dat_time = normalizeBetweenArrays(expr_log2, method = 'quantile')
sum(is.na(rownames(dat_time))) # 0
sum(rownames(dat_time) == '') # 1

which(rownames(dat_time) == '') # 18th row
dat_time = dat_time[-18,] # We remove that row

dim(dat_time) # 32.080 genes

# NORMALIZED PLOTS --------------------------------------------------------
toplot3 = melt(dat_time)
p3 = ggplot(toplot3, aes(x = Var2, y = value)) + geom_boxplot() +
     ggtitle('Quantile normalized expression data') + xlab('') + ylab('') +
     theme_bw() +
     theme(plot.title = element_text(size = 35),
           legend.position = 'none',
           axis.text.y = element_text(size = 25, color = 'darkgrey'),
           axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           legend.text = element_text(size = 25),
           legend.title = element_blank())

#### Normalized MD-plots ####
plotMD(dat_time, column = 1, main = 'MD plot Bastu (GSE111974): RIF-1')

# plotMD(gse111974, column = 43,
#        main = 'MD plot Bastu (GSM3045909): RIF')

#### Boxplot ####
boxplot(dat_time, main = 'Normalized data')

# NORMALIZED PCA ----------------------------------------------------------
pca_norm = prcomp(t(dat_time))

#### Components 1 & 2 ####
var_norm = round(summary(pca_norm)$importance[2, c(1,2)] * 100, 1)
toplot_norm = data.frame(pca_norm$x[, c(1,2)], stringsAsFactors = F)
lim_norm = max(abs(c(min(toplot_norm), max(toplot_norm))))

axis_limits_norm = c(-lim_norm, lim_norm)
# toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
toplot_norm$color = ed$condition

ggplot(data = toplot_norm,
       aes_string(x = colnames(toplot_norm)[1], y = colnames(toplot_norm)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values = brewer.pal(n = 5, 'Dark2')) +
       xlab(paste0('PC1', ': ', var_norm[1], '%')) +
       ylab(paste0('PC2', ': ', var_norm[2], '%')) +
       ggtitle('PCA: Bastu (GSE111974)') +
       xlim(axis_limits_norm) + ylim(axis_limits_norm) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))

# If want labels add:
# geom_text(label = ed$sample, size = 3)

# SAVE .RDA ---------------------------------------------------------------
save(dat_time, ed, file = paste0(data_path, '/bastu_time.rda'), version = 2)

# CORRECT BATCH EFFECT AND TIME EFFECT ------------------------------------
load(paste0(data_path, '/bastu_time.rda'), verbose = T)

#### Show normalized PCA: ####
plot_PCAscores(dat = dat_time, ed = ed, condition1 = 'condition',
               colors = c('#D95F02', '#1B9E77'),
               title = 'PCA: Bastu (GSE111974)')

# We'want to remove the unknow batch effect that seems to affect our data
# (probably time-effect)t:

# At the beginning, we used non-supervised clustering as follows:
# (but finally we will use sva):

# Seems a good separation in PCAplot

# pca = prcomp(t(dat))
#
# var = round(summary(pca)$importance[2, c(1,2)] * 100, 1)
# toplot = data.frame(pca$x[, c(1,2)], stringsAsFactors = F)
# lim = max(abs(c(min(toplot), max(toplot))))
#
# axis_limits = c(-lim, lim)
# # toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
# toplot$color = c(rep('control', 24), rep('RIF', 24))
#
# ggplot(data = toplot,
#        aes_string(x = colnames(toplot)[1], y = colnames(toplot)[2])) +
#        geom_point(aes(color = color), size = 3) +
#        scale_color_manual(name = 'RIF', values = brewer.pal(n = 5, 'Dark2')) +
#        xlab(paste0('PC1', ': ', var[1], '%')) +
#        ylab(paste0('PC2', ': ', var[2], '%')) +
#        ggtitle('PCA: Bastu (GSE111974)') +
#        xlim(axis_limits) + ylim(axis_limits) +
#        theme_light() +
#        theme(legend.position = 'bottom',
#              axis.title = element_text(size = 18),
#              axis.text = element_text(size = 15),
#              plot.title = element_text(size = 22, hjust = 0.5),
#              legend.title = element_blank(),
#              legend.text = element_text(size = 13))

# #### Calculate how many clusters to use ####
# optimal_nclusters(dat = as.matrix(dat_time),
#                   method = 'silhouette',
#                   max_nclusters = 5)
#
# # optimal number of clusters = 3
#
# #### Unsupervised clustering ####
# clustering = unsupervised_clustering(dat = dat_time,
#                                      method = 'kmeans',
#                                      nclusters = 3)
# # clustering$cluster_results
# # clustering$cluster_group
# clustering$cluster_plot
#
# ed$cluster = as.character(clustering$cluster_group)
#
# #### Plot normalized PCA with calculated clusters ###
# plot_PCAscores(dat = dat_time, ed = ed,
#                condition1 = 'condition', condition2 = 'cluster',
#                components = c(1,2), colors = brewer.pal(n = 6, 'Dark2'),
#                title = 'PCA: Bastu (GSE111974)')
#
# # pca = prcomp(t(dat))
# #
# # var = round(summary(pca)$importance[2, c(1,2)] * 100, 1)
# # toplot = data.frame(pca$x[, c(1,2)], stringsAsFactors = F)
# # lim = max(abs(c(min(toplot), max(toplot))))
# #
# # axis_limits = c(-lim, lim)
# # # toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
# # toplot$color = c(rep('control', 24), rep('RIF', 24))
# # toplot$cluster = as.character(clustering$cluster_group)
# #
# # ggplot(data = toplot,
# #        aes_string(x = colnames(toplot)[1], y = colnames(toplot)[2])) +
# #        geom_point(aes(color = color, shape = cluster),  size = 3) +
# #        scale_color_manual(name = 'RIF', values = brewer.pal(n = 5, 'Dark2')) +
# #        xlab(paste0('PC1', ': ', var[1], '%')) +
# #        ylab(paste0('PC2', ': ', var[2], '%')) +
# #        ggtitle('PCA: Bastu (GSE111974)') +
# #        xlim(axis_limits) + ylim(axis_limits) +
# #        theme_light() +
# #        theme(legend.position = 'bottom',
# #              axis.title = element_text(size = 18),
# #              axis.text = element_text(size = 15),
# #              plot.title = element_text(size = 22, hjust = 0.5),
# #              legend.title = element_blank(),
# #              legend.text = element_text(size = 13))
#
# #### Remove time effect ####
# mdesign = model.matrix(~ed$condition)
#
# dat = removeBatchEffect(x = dat_time,
#                         design = mdesign,
#                         batch = ed$cluster)
#
# plot_PCAscores(dat = dat, ed = ed,
#                condition1 = 'condition',
#                components = c(1,2), colors = brewer.pal(n = 6, 'Dark2'),
#                title = 'PCA: Bastu (GSE111974)')
#
# # pca = prcomp(t(dat))
# #
# # var = round(summary(pca)$importance[2, c(1,2)] * 100, 1)
# # toplot = data.frame(pca$x[, c(1,2)], stringsAsFactors = F)
# # lim = max(abs(c(min(toplot), max(toplot))))
# #
# # axis_limits = c(-lim, lim)
# # # toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
# # toplot$color = c(rep('control', 24), rep('RIF', 24))
# #
# # ggplot(data = toplot,
# #        aes_string(x = colnames(toplot)[1], y = colnames(toplot)[2])) +
# #        geom_point(aes(color = color), size = 3) +
# #        scale_color_manual(name = 'RIF', values = brewer.pal(n = 5, 'Dark2')) +
# #        xlab(paste0('PC1', ': ', var[1], '%')) +
# #        ylab(paste0('PC2', ': ', var[2], '%')) +
# #        ggtitle('PCA: Bastu (GSE111974)') +
# #        xlim(axis_limits) + ylim(axis_limits) +
# #        theme_light() +
# #        theme(legend.position = 'bottom',
# #              axis.title = element_text(size = 18),
# #              axis.text = element_text(size = 15),
# #              plot.title = element_text(size = 22, hjust = 0.5),
# #              legend.title = element_blank(),
# #              legend.text = element_text(size = 13))

# SVA FOR UNKNOWN BATCH-EFFECT CORRECTION ---------------------------------
mod = model.matrix( ~ 0 + ed$condition)
n = num.sv(dat = dat_time, mod = mod, method = 'leek') # 2

mod0 = model.matrix( ~ 1, data = ed)
svobj = sva(dat = dat_time, mod = mod, mod0 = mod0, n.sv = n)

dat = removeBatchEffect(x = dat_time, covariates = svobj$sv)

plot_PCAscores(dat = dat, ed = ed, condition1 = 'condition',
               components = c(1, 2),
               colors = c('#D95F02', '#1B9E77'),
               title = 'PCA: Bastu (GSE111974)')

# SAVE .RDA ---------------------------------------------------------------
save(dat, ed, file = paste0(results_path, '/bastu.rda'), version = 2)

load(paste0(results_path, '/bastu.rda'), verbose = T)

head(dat)
head(ed)

plot_PCAscores(dat = dat, ed = ed, condition1 = 'condition',
               components = c(1, 2),
               colors = c('#D95F02', '#1B9E77'),
               title = 'PCA: Bastu (GSE111974)')

################################################################################
pca_raw = prcomp(t(dat))
var_raw = round(summary(pca_raw)$importance[2, c(1,2)] * 100, 1)
toplot_raw = data.frame(pca_raw$x[, c(1,2)], stringsAsFactors = F)
lim_raw = max(abs(c(min(toplot_raw), max(toplot_raw))))

axis_limits_raw = c(-lim_raw, lim_raw)
# toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
toplot_raw$color = ed$condition

ggplot(data = toplot_raw,
       aes_string(x = colnames(toplot_raw)[1], y = colnames(toplot_raw)[2])) +
       geom_point(aes(color = color), size = 3) +
       geom_text_repel(label = ed$sample, size = 3) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
       xlab(paste0('PC1', ': ', var_raw[1], '%')) +
       ylab(paste0('PC2', ': ', var_raw[2], '%')) +
       ggtitle('PCA: Bastu (GSE111974)') +
       xlim(axis_limits_raw) + ylim(axis_limits_raw) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))
