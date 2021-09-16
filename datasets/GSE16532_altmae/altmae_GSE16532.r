################################################################################
## altmae_GSE16532.r
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

# Dataset info avaliable in:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16532

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

# LOAD PACKAGES -----------------------------------------------------------
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(hgug4112a.db)

functions_path = '/Users/francisco/Desktop/TFM/functions'
dataset_path = '/Users/francisco/Desktop/TFM/datasets/GSE16532_altmae'
data_path = paste0(dataset_path, '/data')
results_path = '/Users/francisco/Desktop/TFM/datasets/results/results_sva_2'

source(paste0(functions_path, '/function_plot_PCAscores.r'))

getwd()
setwd(dataset_path)

# READ PLATFORM FILE ------------------------------------------------------
# Platform = Agilent-014850
gpl_name = '/GPL4133.txt'

gpl = read.delim(file = paste0(data_path, gpl_name),
                  header = T,
                  sep = '\t',
                  comment.char = '#',
                  skip = 1,
                  quote = '',
                  stringsAsFactors = F)

dim(gpl) # [1] 45220     22

# READ EXPERIMENTAL DESSIGN -----------------------------------------------
series_name = '/GSE16532_series_matrix.txt'
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
ed = data.frame(sample = ed$sample, condition = ed$sample.type, title = ed$title,
                tissue = ed$tissue, phase = ed$phase, group = ed$group)

ed$condition[ed$condition == 'experimental'] = 'RIF'
rownames(ed) = ed$sample

# DOWNLOAD RAW DATA FOR ALTMAE (GSE16532) ----------------------------------
getGEOSuppFiles('GSE16532', baseDir = './data/')
# if 'error: Timeout of 60 seconds was reached' then 'options(timeout = 300)'.
# Downloaded in '/raw_datasets/altmae_GSE71331/data/GSE71331'
# (Remember to unzip it and prepare paths)

# READ EXPRESSION DATA ----------------------------------------------------
my_files = list.files(path = 'data/GSE16532_RAW', full.names = T)
target_info = data.frame(FileName = my_files,
                          RIF_CONTROL = ed[,2],
                          stringsAsFactors = F)

# View(target_info)

# Check colnames in file:
scan('data/GSE16532_RAW/GSM414976.gpr', nlines = 1, what = 'c', sep = '\t')

columns = list(E = 'F532 Median', Eb = 'B532 Median')

gse16532raw = read.maimages(files = target_info, columns = columns,
                      source = 'agilent',
                      annotation = c('Ref','GeneName','ControlType'))
# Use View(gse16532raw$E) and check with files that parsing is correct.

# SAVE RAW DATA, ED AND GPL ------------------------------------------------
save(ed, file = paste0(data_path, '/ed_altmae.rda'), version = 2)
save(gpl, file = paste0(data_path, '/GPL4133_altmae.rda'), version = 2)
save(gse16532raw, file = paste0(data_path, '/gse16532raw.rda'), version = 2)

# PRE-NORMALIZATION ANALYSIS ----------------------------------------------
load(paste0(data_path, '/ed_altmae.rda'), verbose = T)
load(paste0(data_path, '/gpl4133_altmae.rda'), verbose = T)
load(paste0(data_path, '/gse16532raw.rda'), verbose = T)

# RAW PLOTS ---------------------------------------------------------------
#### MDplot ####
plotMD(gse16532raw, column = 1,
       main = 'MD plot Altmae (GSE16532): raw control-1')
plotMD(gse16532raw, column = 7,
       main = 'MD plot Altmae (GSE16532): raw rif-2')
# Warning message:
# In plotMD.EListRaw(gse16532raw,
# column = 1, main = 'MD plot Altmae (GSE414976):
# raw control-1') : NaNs produced

#### Boxplot ####
# boxplot(data.frame(log2(gse16532raw$Eb)), main = 'Green background')
boxplot(data.frame(log2(gse16532raw$E)), main = 'Raw data')

# RAW PCA -----------------------------------------------------------------
pca_raw = prcomp(t(log2(gse16532raw$E) + 1))
var_raw = round(summary(pca_raw)$importance[2, c(1,2)] * 100, 1)
toplot_raw = data.frame(pca_raw$x[, c(1,2)], stringsAsFactors = F)
lim_raw = max(abs(c(min(toplot_raw), max(toplot_raw))))

axis_limits_raw = c(-lim_raw, lim_raw)
toplot_raw$color = c(rep('control', 5), rep('RIF', 4))

ggplot(data = toplot_raw,
       aes_string(x = colnames(toplot_raw)[1], y = colnames(toplot_raw)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
       xlab(paste0('PC1', ': ', var_raw[1], '%')) +
       ylab(paste0('PC2', ': ', var_raw[2], '%')) +
       ggtitle('PCA: Altmae (GSE16532)') +
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

# We need to create an altmae_eset objet:
altmae_eset = ExpressionSet(assayData = assayDataNew(exprs = gse16532raw$E))

# Now we can use aqm:
arrayQualityMetrics(expressionset = altmae_eset,
                    outdir = outdir,
                    force = TRUE,
                    do.logtransform = TRUE) # Since data is not processed yet.

# Check index.html file in outdir for results.

# CORRECTING BACKGROUND ---------------------------------------------------
gse16532 = backgroundCorrect(gse16532raw, method = 'normexp', offset = 50)
# check if offset needed, we use 50 as default

# ANNOTATION --------------------------------------------------------------
eset = as.matrix(gse16532$E) # expression info
dim(eset) # [1] 45.220     9

#### Filtering probes: controls and NAs ####
probesInfo = data.frame('ProbeName' = gse16532$genes$Ref, # Probes
                        'GeneSymbol' = gse16532$genes$GeneName,
                        'Control' = gse16532$genes$ControlType, # Control?
                        stringsAsFactors = F)

# Remove possible NAs
probesInfo_noNA = probesInfo[!is.na(probesInfo$ProbeName), ]
gpl_noNA = gpl[!is.na(gpl$SPOT_ID), ]

all(gpl_noNA$SPOT_ID == probesInfo_noNA$ProbeName) # TRUE

eset_noNA = eset[!is.na(probesInfo$ProbeName),]
dim(eset_noNA) # 45.015 rows
dim(probesInfo_noNA) # 45.015 rows

# Remove controls
probesInfo_noctrl = probesInfo_noNA[probesInfo_noNA$Control == 'false',]
gpl_noctrl = gpl_noNA[probesInfo_noNA$Control == 'false',]
eset_noctrl = eset_noNA[probesInfo_noNA$Control == 'false',]

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
# gpl_noctrl does not have NAs but '' symbols. We need to remove those empty
# values.
indexNA_symbol = which(gpl_noctrl$GENE_SYMBOL == '')
gpl_noctrl_notNA = gpl_noctrl[-indexNA_symbol,]
eset_noctrl_notNA = eset_noctrl[-indexNA_symbol,]
dim(eset_noctrl_notNA)
dim(gpl_noctrl_notNA)
# Still we have 32.696 genes

exprbygene = avereps(eset_noctrl_notNA, ID = gpl_noctrl_notNA$GENE_SYMBOL)

dim(exprbygene) # 19.749 genes
colnames(exprbygene) = ed$sample

# PLOT EXPRESSION BY GENE -------------------------------------------------
#### Raw expression ####
toplot = melt(exprbygene)
p1 = ggplot(toplot, aes(x = Var2, y = value)) +
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
toplot = melt(expr_log2)
p2 = ggplot(toplot, aes(x = Var2, y = value)) + geom_boxplot() +
     ggtitle('Log2 expression data') + xlab('') + ylab('') + theme_bw() +
     theme(plot.title = element_text(size = 35), legend.position = 'none',
           axis.text.y = element_text(size = 25, color = 'darkgrey'),
           axis.text.x = element_blank(), axis.ticks.x = element_blank(),
           legend.text = element_text(size = 25),
           legend.title = element_blank())

# QUANTILE NORMALIZATION --------------------------------------------------
# Final step of normalization
dat = normalizeBetweenArrays(expr_log2, method = 'quantile')
sum(is.na(rownames(dat))) # 0
sum(rownames(dat) == '') # 0
dim(dat) # 19.749

# NORMALIZED PLOTS --------------------------------------------------------
toplot = melt(dat)
p3 = ggplot(toplot, aes(x = Var2, y = value)) + geom_boxplot() +
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
plotMD(dat, column = 1, main = 'MD plot Altmae (GSE16532): RIF-1')

#### Boxplot ####
boxplot(dat, main = 'Normalized data')

# NORMALIZED PCA ----------------------------------------------------------
pca_norm = prcomp(t(dat))

#### Components 1 & 2 ####
var_norm = round(summary(pca_norm)$importance[2, c(1,2)] * 100, 1)
toplot_norm = data.frame(pca_norm$x[, c(1,2)], stringsAsFactors = F)
lim_norm = max(abs(c(min(toplot_norm), max(toplot_norm))))

axis_limits_norm = c(-lim_norm, lim_norm)
# toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
toplot_norm$color = c(rep('RIF', 4), rep('control', 5))

ggplot(data = toplot_norm,
       aes_string(x = colnames(toplot_norm)[1], y = colnames(toplot_norm)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
#      geom_text_repel(label = rownames(pData(gse26787raw)), size = 3) +
       xlab(paste0('PC1', ': ', var_norm[1], '%')) +
       ylab(paste0('PC2', ': ', var_norm[2], '%')) +
       ggtitle('PCA: Altmae (GSE16532)') +
       xlim(axis_limits_norm) + ylim(axis_limits_norm) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))

#### Components 3 & 4 ####
var = round(summary(pca)$importance[2, c(3 ,4)] * 100, 1) # Proportion of varian
toplot = data.frame(pca$x[, c(3, 4)], stringsAsFactors = F)
lim = max(abs(c(min(toplot), max(toplot))))

axis_limits = c(-lim, lim)
# toplot$color = c(paste0(rep('Control-'), 1:2), paste0(rep('RIF-'), 1:3))
toplot$color = c(rep('RIF', 7), rep('Control', 5))

ggplot(data = toplot,
       aes_string(x = colnames(toplot)[1], y = colnames(toplot)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values = brewer.pal(n = 5, 'Dark2')) +
       xlab(paste0('PC3', ': ', var[1], '%')) +
       ylab(paste0('PC4', ': ', var[2], '%')) +
       ggtitle('PCA: Altmae (GSE16532)') +
       xlim(axis_limits) + ylim(axis_limits) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))

# boxplot(data.frame(dat), main = 'Normalized data')

# SAVE .RDA ---------------------------------------------------------------
save(dat, ed, file = paste0(results_path, '/altmae.rda'), version = 2)

load(paste0(results_path, '/altmae.rda'), verbose = T)
head(dat)
head(ed)

plot_PCAscores(dat = dat, ed = ed,
               condition1 = 'condition',
               components = c(1,2),
               colors = c('#D95F02', '#1B9E77'),
               title = 'PCA: Altmae(GSE16532)')
