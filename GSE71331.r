################################################################################
## shi_GSE71331.R
## 2021-02-15
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

# Dataset info avaliable in:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71331

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

# LOAD PACKAGES -----------------------------------------------------------
library(Biobase)
library(GEOquery)
library(ggplot2)
library(RColorBrewer)
library(AnnotationDbi)
library(limma)
library(dplyr)
library(data.table)
library(foreach)
library(hgug4112a.db)
library(cluster)
library(ggrepel)
library(plyr)

functions_path = '/Users/francisco/Desktop/TFM/functions'
dataset_path = '/Users/francisco/Desktop/TFM/datasets/GSE71331_shi'
data_path = paste0(dataset_path, '/data')
raw_plots_path = paste0(dataset_path, '/raw_plots')
normal_plots_path = paste0(dataset_path, '/normal_plots')
results_path = '/Users/francisco/Desktop/TFM/datasets/results/results_sva_2'

source(paste0(functions_path, '/function_arbol.r'))
source(paste0(functions_path, '/function_plot_PCAscores.r'))
source(paste0(functions_path, '/function_optimal_nclusters.r'))
source(paste0(functions_path, '/function_unsupervised_clustering.r'))

getwd()
setwd(dataset_path)

# READ PLATFORM FILE ------------------------------------------------------
# Platform = Agilent-052909
gpl_name = '/GPL19072.txt'

gpl = read.delim(file = paste0(data_path, gpl_name),
                 header = T,
                 sep = '\t',
                 comment.char = '#',
                 skip = 42,
                 stringsAsFactors = F)

dim(gpl) # [1] 180881     20
gpl = gpl[!gpl$SPOT_ID == '', ]

head(gpl)

# But this file only contains probe names and not gene names, symbols... We use
# hgug4112a.db to annotate these

#### Adding gene symbol to gpl ####
keytypes(hgug4112a.db)

annotation_info = AnnotationDbi::select(hgug4112a.db,
                                  keys = gpl$SPOT_ID,
                                  columns = c('ENTREZID', 'ENSEMBL', 'SYMBOL'),
                                  keytype = 'PROBEID')
# 'select()' returned many:many mapping between keys and columns

# We keep first match
posiciones = match(unique(annotation_info[, 1]), annotation_info[, 1]) # PROBEID
annotation_info = annotation_info[posiciones, ]
posiciones = match(unique(annotation_info[, 2]), annotation_info[, 2]) # ENTREZ
annotation_info = annotation_info[posiciones, ]

dim(annotation_info) # [1] 14725     4
head(annotation_info)
colnames(annotation_info) = c('SPOT_ID','ENTREZID', 'ENSEMBL', 'SYMBOL')

#### Remove NAs ####
annotation_info = annotation_info[!is.na(annotation_info[, 'ENTREZID']), ]
dim(annotation_info) # [1] 14724     4

#### Now, we join gpl with genes info in annotation_info ####
new_gpl = join(x = gpl, y = annotation_info, by = 'SPOT_ID')
dim(new_gpl) # [1] 180880     23
head(new_gpl)

#### Create a correct gpl ####
gpl = data.frame('ID' = new_gpl$ID, 'COL' = new_gpl$COL, 'ROW' = new_gpl$ROW,
                 'PROBE_NAME' = new_gpl$NAME, 'SPOT_ID' = new_gpl$SPOT_ID,
                 'CONTROL_TYPE' = new_gpl$CONTROL_TYPE,
                 'GENE_SYMBOL' = new_gpl$SYMBOL,
                 'ENSEMBL_ID' = new_gpl$ENSEMBL,
                 'ENTREZ_ID' = new_gpl$ENTREZID,
                 'CHROMOSOMAL_LOCATION' = new_gpl$CHROMOSOMAL_LOCATION,
                 'SEQUENCE' = new_gpl$SEQUENCE,
                 stringsAsFactors = F)

dim(gpl) # 180.880 rows 11 columns
head(gpl)

# READ EXPERIMENTAL DESSIGN -----------------------------------------------
series_name = '/GSE71331_series_matrix.txt'
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
rownames(ed) = ed$sample

# To Homogenize between dataframes:
ed = data.frame(sample = ed$sample, condition = ed$patients, title = ed$title,
                tissue = ed$tissue, collect.time = ed$collect.time)

ed$condition[ed$condition == 'RIF (recurrent implantation failure)'] = 'RIF'
ed$condition[ed$condition == 'Control'] = 'control'
ed$tissue[ed$tissue == 'Endometrium'] = 'endometrium'
ed$collect.time[ed$collect.time == 'Window of implantation'] = 'woi'

rownames(ed) = ed$sample

# DOWNLOAD RAW DATA FOR SHI (GSE26787) ----------------------------------
getGEOSuppFiles('GSE71331', baseDir = './data/')
# if 'error: Timeout of 60 seconds was reached' then 'options(timeout = 300)'.
# Downloaded in '/raw_datasets/shi_GSE71331/data/GSE71331'
# (Remember to unzip it and prepare paths)

# READ EXPRESSION DATA ----------------------------------------------------
my_files = list.files(path = 'data/GSE71331_RAW', full.names = T)
target_info = data.frame(FileName = my_files,
                          RIF_CONTROL = ed[,2],
                          stringsAsFactors = F)
# View(target_info)

## Read and normalize data
# We get warnings. Problems with enconding.
# > Sys.getlocale()
# [1] 'en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8'
# Then do:
# Sys.setlocale(locale='C')
# > Sys.getlocale()
# [1] 'C/C/C/C/C/en_US.UTF-8'
# To read.maimages correctly.
# You can change it back with:
# > Sys.setlocale(locale='en_US.UTF-8')
# [1] 'en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8'
# More info:
# https://support.bioconductor.org/p/32640/
# https://www.rdocumentation.org/packages/base/versions/3.6.1/topics/locales
gse71331raw = read.maimages(files = target_info, source = 'agilent.median',
                    green.only = T)
# Use View(gse71331raw$E) and check with files that parsing is correct.

# SAVE RAW DATA, ED AND GPL ------------------------------------------------
save(ed, file = paste0(data_path, '/ed_shi.rda'), version = 2)
save(gpl, file = paste0(data_path, '/gpl19072_shi.rda'), version = 2)
save(gse71331raw, file = paste0(data_path, '/gse71331raw.rda'), version = 2)

# PRE-NORMALIZATION ANALYSIS ----------------------------------------------
load(paste0(data_path, '/ed_shi.rda'), verbose = T)
load(paste0(data_path, '/gpl19072_shi.rda'), verbose = T)
load(paste0(data_path, '/gse71331raw.rda'), verbose = T)

# RAW PLOTS ---------------------------------------------------------------
#### MDplot ####
plotMD(gse71331raw, column = 1, main = 'MD plot Shi (GSE71331): raw RIF-1')
plotMD(gse71331raw, column = 2, main = 'MD plot Shi (GSE71331): raw RIF-2')
plotMD(gse71331raw, column = 3, main = 'MD plot Shi (GSE71331): raw RIF-3')
plotMD(gse71331raw, column = 4, main = 'MD plot Shi (GSE71331): raw RIF-4')
plotMD(gse71331raw, column = 5, main = 'MD plot Shi (GSE71331): raw RIF-5')
plotMD(gse71331raw, column = 6, main = 'MD plot Shi (GSE71331): raw RIF-6')
plotMD(gse71331raw, column = 7, main = 'MD plot Shi (GSE71331): raw RIF-7')
plotMD(gse71331raw, column = 8, main = 'MD plot Shi (GSE71331): raw Control-1')
plotMD(gse71331raw, column = 9, main = 'MD plot Shi (GSE71331): raw Control-2')
plotMD(gse71331raw, column = 10, main = 'MD plot Shi (GSE71331): raw Control-3')
plotMD(gse71331raw, column = 11, main = 'MD plot Shi (GSE71331): raw Control-4')
plotMD(gse71331raw, column = 12, main = 'MD plot Shi (GSE71331): raw Control-5')
# Warning message: In plotMD.EListRaw(gse71331raw) : NaNs produced

#### Boxplot ####
# boxplot(data.frame(log2(gse71331raw$Eb)), main = 'Green background')
boxplot(data.frame(log2(gse71331raw$E)), main = 'Raw data')
# imageplot(log2(gse71331raw$Eb[, 4]), gse71331raw$printer) # optional

# RAW PCA -----------------------------------------------------------------
pca_raw = prcomp(t(log2(gse71331raw$E) + 1)) # t cause we need samples in rows

var_raw = round(summary(pca_raw)$importance[2, c(1,2)] * 100, 1)
toplot_raw = data.frame(pca_raw$x[, c(1,2)], stringsAsFactors = F)
lim_raw= max(abs(c(min(toplot_raw), max(toplot_raw))))

axis_limits_raw = c(-lim_raw, lim_raw)
# toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
toplot_raw$color = c(rep('RIF', 7), rep('control', 5))

ggplot(data = toplot_raw,
       aes_string(x = colnames(toplot_raw)[1], y = colnames(toplot_raw)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
#       geom_text_repel(label = rownames(pData(gse26787raw)), size = 3) +
       xlab(paste0('PC1', ': ', var_raw[1], '%')) +
       ylab(paste0('PC2', ': ', var_raw[2], '%')) +
       ggtitle('PCA: Shi (GSE71331)') +
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

# We need to create an shi_eset objet:
shi_eset = ExpressionSet(assayData = assayDataNew(exprs = gse71331raw$E))

# Now we can use aqm:
arrayQualityMetrics(expressionset = shi_eset,
                    outdir = outdir,
                    force = TRUE,
                    do.logtransform = TRUE) # Since data is not processed yet.

# Check index.html file in outdir for results.

# CORRECTING BACKGROUND ---------------------------------------------------
gse71331 = backgroundCorrect(gse71331raw, method = 'normexp', offset = 50)
# check if offset needed, we use 50 as default

# ANNOTATION --------------------------------------------------------------
eset = as.matrix(gse71331$E) # expression info
dim(eset) # [1] 180880     12

#### Filtering probes: controls and NAs ####
probesInfo = data.frame('ProbeName' = gse71331$genes$ProbeName, # Probes
                        'GeneSymbol' = gse71331$genes$GeneName,
                        'Control' = gse71331$genes$ControlType, # Control?
                        stringsAsFactors = F)

# Remove possible NAs
probesInfo_noNA = probesInfo[!is.na(probesInfo$ProbeName), ]
gpl_noNA = gpl[!is.na(gpl$SPOT_ID), ]

all(gpl_noNA$SPOT_ID == probesInfo_noNA$ProbeName) # TRUE

eset_noNA = eset[!is.na(probesInfo$ProbeName),]
dim(eset_noNA) # 180.880 rows
dim(probesInfo_noNA) # 180.880 rows

# Remove controls
probesInfo_noctrl = probesInfo_noNA[probesInfo_noNA$Control == '0',]
gpl_noctrl = gpl_noNA[probesInfo_noNA$Control == '0',]
eset_noctrl = eset_noNA[probesInfo_noNA$Control == '0',]

#### Check if they have the same number of rows ####
dim(probesInfo_noctrl) # Probe info
dim(gpl_noctrl) # Annotationinfo
dim(eset_noctrl) # Expression info
all(gpl_noctrl$SPOT_ID == probesInfo_noctrl$ProbeName) # TRUE

# GROUP PROBESETS INFORMATION ---------------------------------------------
# Condense replicate probes by their average
exprbyprobe = avereps(eset_noctrl, ID = probesInfo_noctrl$ProbeName)

which(rownames(exprbyprobe) == '') # 0
which(is.na(rownames(exprbyprobe))) # 0

# GROUP GENES BY PROBESET ID ----------------------------------------------
indexNA_symbol = which(is.na(gpl_noctrl$GENE_SYMBOL))
gpl_noctrl_notNA = gpl_noctrl[-indexNA_symbol,]
eset_noctrl_notNA = eset_noctrl[-indexNA_symbol,]
dim(eset_noctrl_notNA)
dim(gpl_noctrl_notNA)

exprbygene = avereps(eset_noctrl_notNA, ID = gpl_noctrl_notNA$GENE_SYMBOL)

dim(exprbygene)
colnames(exprbygene) = ed$sample

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
sum(rownames(dat_time) == '') # 0

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
plotMD(dat_time, column = 1, main = 'MD plot Shi (GSE71331): RIF-1')
plotMD(dat_time, column = 2, main = 'MD plot Shi (GSE71331): RIF-2')
plotMD(dat_time, column = 3, main = 'MD plot Shi (GSE71331): RIF-3')
plotMD(dat_time, column = 4, main = 'MD plot Shi (GSE71331): RIF-4')
plotMD(dat_time, column = 5, main = 'MD plot Shi (GSE71331): RIF-5')
plotMD(dat_time, column = 6, main = 'MD plot Shi (GSE71331): RIF-6')
plotMD(dat_time, column = 7, main = 'MD plot Shi (GSE71331): RIF-7')
plotMD(dat_time, column = 8, main = 'MD plot Shi (GSE71331): Control-1')
plotMD(dat_time, column = 9, main = 'MD plot Shi (GSE71331): Control-2')
plotMD(dat_time, column = 10, main = 'MD plot Shi (GSE71331): Control-3')
plotMD(dat_time, column = 11, main = 'MD plot Shi (GSE71331): Control-4')
plotMD(dat_time, column = 12, main = 'MD plot Shi (GSE71331): Control-5')

#### Boxplot ####
boxplot(dat_time, main = 'Normalized data')

# NORMALIZED PCA ----------------------------------------------------------
pca_norm = prcomp(t(dat_time))

#### Components 1 & 2 ####
var_norm = round(summary(pca_norm)$importance[2, c(1,2)] * 100, 1)
toplot_norm = data.frame(pca_norm$x[, c(1,2)], stringsAsFactors = F)
lim_norm = max(abs(c(min(toplot_norm), max(toplot_norm))))

axis_limits_norm = c(-lim_norm, lim_norm)
toplot_norm$color = c(rep('RIF', 7), rep('control', 5))

ggplot(data = toplot_norm,
       aes_string(x = colnames(toplot_norm)[1], y = colnames(toplot_norm)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values =  c('#D95F02', '#1B9E77')) +
#       geom_text_repel(label = ed$sample, size = 3) +
       xlab(paste0('PC1', ': ', var_norm[1], '%')) +
       ylab(paste0('PC2', ': ', var_norm[2], '%')) +
       ggtitle('PCA: Shi (GSE71331)') +
       xlim(axis_limits_norm) + ylim(axis_limits_norm) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))


#### Components 3 & 4 ####
var_norm2 = round(summary(pca_norm)$importance[2, c(3 ,4)] * 100, 1)
toplot_norm2 = data.frame(pca_norm$x[, c(3, 4)], stringsAsFactors = F)
lim_norm2 = max(abs(c(min(toplot_norm2), max(toplot_norm2))))

axis_limits_norm2 = c(-lim_norm2, lim_norm2)
# toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
toplot_norm2$color = c(rep('RIF', 7), rep('control', 5))

ggplot(data = toplot_norm2,
       aes_string(x = colnames(toplot_norm2)[1],
                  y = colnames(toplot_norm2)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values =  c('#D95F02', '#1B9E77')) +
       xlab(paste0('PC3', ': ', var_norm2[1], '%')) +
       ylab(paste0('PC4', ': ', var_norm2[2], '%')) +
       ggtitle('PCA: Shi (GSE71331)') +
       xlim(axis_limits_norm2) + ylim(axis_limits_norm2) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))

# boxplot(data.frame(dat), main = 'Normalized data')

# SAVE .RDA ---------------------------------------------------------------
save(dat_time, ed, file = paste0(data_path, '/shi_time.rda'), version = 2)

# CORRECT BATCH EFFECT AND TIME EFFECT ------------------------------------
load(paste0(data_path, '/shi_time.rda'), verbose = T)

#### Show normalized PCA: ####
plot_PCAscores(dat = dat_time, ed = ed, condition1 = 'condition',
               title = 'PCA: Shi (GSE71331)',
               colors = c('#D95F02', '#1B9E77'))

# We'want to remove the unknow batch effect that seems to affect our data
# (probably time-effect)t:

# At the beginning, we used non-supervised clustering as follows:
# (but finally we will use sva):

# #### Calculate how many clusters to use ####
# optimal_nclusters(dat = as.matrix(dat_time),
#                   method = 'silhouette',
#                   max_nclusters = 5)
#
# # optimal number of clusters = 2
#
# #### Unsupervised clustering ####
# # clustering = unsupervised_clustering(dat = dat_time,
# #                                      method = 'kmeans',
# #                                      nclusters = 2, groups = ed$condition,
# #                                      group_colors = c(RIF='red',
# #                                                       control='green'))
#
# clustering = unsupervised_clustering(dat = dat_time,
#                                      method = 'kmeans',
#                                      nclusters = 2)
# clustering$cluster_plot
# # clustering$cluster_results
# # clustering$cluster_group
#
# ed$cluster = as.character(clustering$cluster_group)
#
# #### Plot normalized PCA with calculated clusters ###
# plot_PCAscores(dat = dat_time, ed = ed,
#                 condition1 = 'condition', condition2 = 'cluster',
#                 components = c(1,2), colors = brewer.pal(n = 6, 'Dark2'),
#                 title = 'PCA: Shi (GSE26787)')
#
# #### Remove time effect ####
# mdesign = model.matrix(~ed$condition)
#
# dat_corrected = removeBatchEffect(x = dat_time,
#                                   design = mdesign,
#                                   batch = ed$cluster)
#
# plot_PCAscores(dat = dat_corrected, ed = ed,
#                condition1 = 'condition',
#                components = c(1,2),
#                colors = brewer.pal(n = 6, 'Dark2'),
#                title = 'PCA: Shi(GSE26787)')
#
# # pca_correct = prcomp(t(dat))
# #
# # var_correct = round(summary(pca_correct)$importance[2, c(1,2)] * 100, 1)
# # toplot_correct = data.frame(pca_correct$x[, c(1,2)], stringsAsFactors = F)
# # lim_correct = max(abs(c(min(toplot_correct), max(toplot_correct))))
# #
# # axis_limits_correct = c(-lim_correct, lim_correct)
# # # toplot$color = c(paste0(rep('control-'), 1:2), paste0(rep('RIF-'), 1:3))
# # toplot_correct$color = c(rep('RIF', 7), rep('control', 5))
# #
# # ggplot(data = toplot_correct,
# #        aes_string(x = colnames(toplot_correct)[1],
# #                   y = colnames(toplot_correct)[2])) +
# #        geom_point(aes(color = color), size = 3) +
# #        geom_text_repel(label = ed$sample, size = 3) +
# #        scale_color_manual(name = 'RIF', values = brewer.pal(n = 5, 'Dark2')) +
# #        xlab(paste0('PC1', ': ', var_correct[1], '%')) +
# #        ylab(paste0('PC2', ': ', var_correct[2], '%')) +
# #        ggtitle('PCA: Shi (GSE71331)') +
# #        xlim(axis_limits_correct) + ylim(axis_limits_correct) +
# #        stat_ellipse(aes(fill=color), condition = 't', level = 0.90)+
# #        theme_light() +
# #        theme(legend.position = 'bottom',
# #              axis.title = element_text(size = 18),
# #              axis.text = element_text(size = 15),
# #              plot.title = element_text(size = 22, hjust = 0.5),
# #              legend.title = element_blank(),
# #              legend.text = element_text(size = 13))

# Decide to remove GSM1832698 according to ellipses.
# dat = dat_corrected[, colnames(dat_corrected)!='GSM1832698']
# ed = ed[rownames(ed)!='GSM1832698',]

# SVA FOR UNKNOWN BATCH-EFFECT CORRECTION ---------------------------------
# mod = model.matrix( ~ 0 + ed$condition)
# n = num.sv(dat = dat_time, mod = mod, method = 'leek') # 6
#
# mod0 = model.matrix( ~ 1, data = ed)
# svobj = sva(dat = dat_time, mod = mod, mod0 = mod0, n.sv = n)
#
# dat = removeBatchEffect(x = dat_time, covariates = svobj$sv)

plot_PCAscores(dat = dat, ed = ed, condition1 = 'condition',
               components = c(1, 2),
               title = 'PCA: Shi (GSE71331)',
               colors = c('#D95F02', '#1B9E77'))

# SAVE .RDA ---------------------------------------------------------------
save(dat, ed, file = paste0(results_path, '/shi.rda'), version = 2)

load(paste0(results_path, '/shi.rda'), verbose = T)
head(dat)
head(ed)

plot_PCAscores(dat = dat, ed = ed,
               condition1 = 'condition',
               components = c(1,2),
               colors = c('#D95F02', '#1B9E77'),
               title = 'PCA: Shi (GSE26787)')
