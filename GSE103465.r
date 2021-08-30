################################################################################
## guo_GSE103465.r
## 2021-02-15
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

# Dataset info avaliable in:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103465

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

# LOAD PACKAGES -----------------------------------------------------------
library(Biobase)
library(GEOquery)
library(tseries)
library(affy)
library(affyPLM)
library(ggplot2)
library(RColorBrewer)
library(AnnotationDbi)
library(limma)
library(dplyr)
library(sva)
library(arrayQualityMetrics)

functions_path = '/Users/francisco/Desktop/TFM/functions'
dataset_path = '/Users/francisco/Desktop/TFM/datasets/GSE103465_guo'
data_path = paste0(dataset_path, '/data')
raw_plots_path = paste0(dataset_path, '/raw_plots')
normal_plots_path = paste0(data_path, '/normal_plots')
results_path = '/Users/francisco/Desktop/TFM/datasets/results/results_sva_2' # modify

source(paste0(functions_path, '/function_arbol.r'))
source(paste0(functions_path, '/function_plot_PCAscores.r'))

getwd()
setwd(dataset_path)

# DOWNLOAD RAW DATA FOR GUO (GSE103465) -------------------------------
getGEOSuppFiles('GSE103465') # Raw dataset provided as supplementary file
system('tar xvf GSE103465/GSE103465_RAW.tar')

gse103465raw = ReadAffy()

system('rm -fr GSE103465')
system('rm *.gz')

# SAVE RAW DATA -----------------------------------------------------------
save(gse103465raw, file = paste0(data_path, '/gse103465raw.rda'), version = 2)

load(paste0(data_path, '/gse103465raw.rda'), verbose = T) # raw data

# PRE-NORMALIZATION ANALYSIS ----------------------------------------------
# GSE = getGEO('GSE103465') # series_matrix

# dim(exprs(gse103465raw)) # 535.824 probes y 6 samples
# annotation(gse103465raw) # chip 'primeview'
head(exprs(gse103465raw))
# colnames(gse103465raw) = gsub(x = colnames(gse103465raw),
#                               pattern = '_.*.CEL.gz', replacement = '')

#### MDplot ####
# Change names to get better plots:
old_rownames_raw = rownames(pData(gse103465raw))
fancy_rownames_raw = gsub(x = colnames(gse103465raw),
                   pattern = '_.*.CEL.gz', replacement = '')

rownames(pData(gse103465raw)) = fancy_rownames_raw

affy::MAplot(gse103465raw, type = 'pm', plot.method = 'smoothScatter')

#### Density estimator ####
affy::hist(gse103465raw)
legend(x = 'topright',          # Position
       legend = rownames(pData(gse103465raw)),  # Legend
       lty = c(1:6),           # Line types
       col = c(1:6),           # Line colors
       lwd = 2)                # line width

#### Boxplot ####
affy::boxplot(gse103465raw)

#### Clustering ####
datos_raw = exprs(gse103465raw)
colnames(datos_raw) = rownames(pData(gse103465raw))

# Make a data.frame to colour our groups
group.raw = c(rep('control', 3), rep('RIF', 3))
sample.name = colnames(datos_raw)
sinfo_raw = as.data.frame(cbind(sample.name, group.raw), as.is = T)
char.group = as.character(sinfo_raw[, 'group.raw'])
tan = unique(char.group)
micolor = rainbow(length(tan))
names(micolor) = tan
sinfo_raw$color.group = micolor[char.group]
sinfo_raw

## Clustering by correlation distance
correlacion_raw = cor(datos_raw)
distancia_raw = as.dist((1 - correlacion_raw) / 2)
hc_raw = hclust(distancia_raw)

table(hc_raw$labels == colnames(datos_raw)) # same layout?

hc_raw$clase = sinfo_raw$group # Colour

## Plot using 'function_arbol.r' ##
arbol(cluster = hc_raw, main = 'Clustering by Correlation Distance')

# RAW PCA -----------------------------------------------------------------
pca_raw = prcomp(t(log2(exprs(gse103465raw)) + 1))
# t cause we need samples in rows

# We can check PCA results
# see:
# https://towardsdatascience.com/principal-component-analysis-pca-101-using-r
# consider 'center = TRUE, scale = TRUE'
# screeplot(pca, type = 'l', npcs = 6, main = 'Screeplot of the first 6 PCs')
# abline(h = 1, col = 'red', lty = 5)
# legend('topright', legend = c('Eigenvalue = 1'),
#        col = c('red'), lty = 5, cex = 0.6)

var_raw = round(summary(pca_raw)$importance[2, c(1,2)] * 100, 1) # rounds 1st d
toplot_raw = data.frame(pca_raw$x[, c(1,2)], stringsAsFactors = F)
lim_raw = max(abs(c(min(toplot_raw), max(toplot_raw)))) # axis values limit
axis_limits_raw = c(-lim_raw, lim_raw)
# toplot$color = c(paste0('Control-', 1:3), paste0('RIF-', 1:3))
toplot_raw$color = c(rep('control', 3), rep('RIF', 3))

# png(filename = paste0(raw_plots_path, '/PCAplot_group.png'), # File name
#     width = 600, height = 600, # Width and height in pixels
#     unit = 'px', # Width and height units, pixels
#     bg = 'white') # Background color

ggplot(data = toplot_raw,
       aes_string(x = colnames(toplot_raw)[1], y = colnames(toplot_raw)[2])) +
       geom_point(aes(color = color), size = 3) +
#      scale_color_manual(name = 'RIF', values = brewer.pal(n = 3, 'Dark2')) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
       geom_text_repel(label = rownames(pData(gse103465raw)), size = 3) +
       xlab(paste0('PC1', ': ', var_raw[1], '%')) +
       ylab(paste0('PC2', ': ', var_raw[2], '%')) +
       ggtitle('PCA: Guo (GSE103465)') +
       xlim(axis_limits_raw) + ylim(axis_limits_raw) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))
# dev.off()

# OUTLIER DETECTION (ARRAY QUALITY METRICS) -------------------------------
outdir = paste0(dataset_path, '/arrayQuality_report')

arrayQualityMetrics(expressionset = gse103465raw,
                    outdir = outdir,
                    force = TRUE,
                    do.logtransform = TRUE) # Since data is not processed yet.

# Check index.html file in outdir for results.

# EXCLUDING GSM2771797 (CONTROL-3) ----------------------------------------
rownames(pData(gse103465raw)) = old_rownames_raw # needed for rma work properly

# According to all previous analysis we decide to exclude control-3 (GSM2771797)
gse103465raw_no3 = gse103465raw[, c(1, 2, 4, 5, 6)]

# NORMALIZATION -----------------------------------------------------------
gse103465 = affy::rma(gse103465raw_no3)

# POST-NORMALIZATION ANALYSIS ---------------------------------------------
#### MDplot ####
old_rownames_norm = rownames(pData(gse103465))
fancy_rownames_norm = gsub(x = colnames(gse103465),
                   pattern = '_.*.CEL.gz', replacement = '')

rownames(pData(gse103465)) = fancy_rownames_norm

affy::MAplot(gse103465, plot.method = 'smoothScatter')

#### Density estimator ####
affy::hist(gse103465)
legend(x = 'topright',          # Position
       legend = rownames(pData(gse103465)),  # Legend texts
       lty = c(1:6),           # Line types
       col = c(1:6),           # Line colors
       lwd = 2)                # line width

#### Boxplot ####
affy::boxplot(gse103465)

#### Clustering ####
datos_norm = exprs(gse103465)
colnames(datos_norm) = rownames(pData(gse103465))

# Make a data.frame to colour our groups
group.norm = c(rep('control', 2), rep('RIF', 3))
sample.name = colnames(datos_norm)
sinfo_norm = as.data.frame(cbind(sample.name, group.norm), as.is = T)
char.group = as.character(sinfo_norm[, 'group.norm'])
tan = unique(char.group)
micolor = rainbow(length(tan))
names(micolor) = tan
sinfo_norm$color.group = micolor[char.group]
sinfo_norm

# Clustering by correlation distance
correlacion_norm = cor(datos_norm)
distancia_norm = as.dist((1 - correlacion_norm) / 2)
hc_norm = hclust(distancia_norm)

table(hc_norm$labels ==  colnames(datos_norm)) # same layout?

hc_norm$clase = sinfo_norm$group # Colour

# Plot using 'function_arbol.r'
arbol(cluster = hc_norm, main = 'Clustering by Correlation Distance')

# NORMALIZED PCA ----------------------------------------------------------
pca_norm = prcomp(t(exprs(gse103465)))

# screeplot(pca, type = 'l', npcs = 6, main = 'Screeplot of the first 6 PCs')
# abline(h = 1, col = 'red', lty = 5)
# legend('topright', legend = c('Eigenvalue = 1'),
#        col = c('red'), lty = 5, cex = 0.6)

var_norm = round(summary(pca_norm)$importance[2, c(1,2)] * 100, 1)
toplot_norm = data.frame(pca_norm$x[, c(1,2)], stringsAsFactors = F)
lim_norm = max(abs(c(min(toplot_norm), max(toplot_norm))))

axis_limits_norm = c(-lim_norm, lim_norm)
toplot_norm$color = c(rep('control', 2), rep('RIF', 3))

# png(filename = paste0(normal_plots_path, '/PCAplot_group.png'), # File name
#     width = 650, height = 650, # Width and height in pixels
#     unit = 'px', # Width and height units, pixels
#     bg = 'white') # Background color

ggplot(data = toplot_norm,
       aes_string(x = colnames(toplot_norm)[1], y = colnames(toplot_norm)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
#       geom_text_repel(label = rownames(pData(gse103465)), size = 3) +
       xlab(paste0('PC1', ': ', var_norm[1], '%')) +
       ylab(paste0('PC2', ': ', var_norm[2], '%')) +
       ggtitle('PCA: Guo (GSE103465)') +
       xlim(axis_limits_norm) + ylim(axis_limits_norm) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))

# dev.off()

# ANNOTATION --------------------------------------------------------------
#### pData ####
factores = c(0, 0, 1, 1, 1)

condition = factor(factores, levels = 0:1, labels = c('control','RIF'))

infor_fenotip = data.frame(pData(gse103465), condition)
pData(gse103465) = infor_fenotip
rownames(pData(gse103465)) = fancy_rownames_norm
pData(gse103465)$sample = fancy_rownames_norm
pData(gse103465)

#### fData ####
annotation(gse103465) # primeview
# We need to create primeview.db since this annotation info does not exist.
# More info in:
# /Estudios_InSilico/Trabajo/RIF_gse103465.Rmd
# https://support.bioconductor.org/p/130727/
# require(primeview.db)

# Or alternatively, we read GPL (chip/platform) info:
gpl = read.delim(file = paste0(data_path, '/GPL16043-14507.txt'), header = T,
                  comment.char = '#', stringsAsFactors = F, sep = '\t')

# Expression info:
exprinfo = exprs(gse103465)

# head(gpl)
dim(gpl) # 49.495 rows
dim(exprinfo)# 49.495 rows

# Are rows ordered?
gpl$ID[(!gpl$ID == rownames(exprs(gse103465)))] # No order, show not-ordered

# When multiple gene symbols with '///', we keep the first:
gpl$Gene.Symbol = unlist(lapply(gpl$Gene.Symbol, function(my_symbol){
        unlist(strsplit(x = my_symbol, split = ' /// '))[1]
}))

# Which positions in gpl are not in order?:
which(gpl$ID%in%gpl$ID[(!gpl$ID == rownames(exprs(gse103465)))])
# [1] 49301 49302 49303 49304 49305 49306 49307 49308 49309 49310 49311 49312
# 49313 49314 49315 49316 49317 49318 49319 49320 49321 49322
# [23] 49323 49324 49325 49326 49327 49328 49329 49330 49331 49332 49333 49334
# 49335 49336 49337 49338 49339 49340 49341 49342 49343 49344
# [45] 49345 49346 49347 49348 49349 49350 49351 49352 49353 49354 49355 49356
# 49357 49358 49359 49360 49361 49362 49363 49364 49365 49366
# [67] 49367 49368 49369 49370 49371 49372 49373 49374 49375 49376 49377 49378
# 49379 49380 49381 49382 49383 49384 49385 49386 49387 49388
# [89] 49389 49390 49391 49392 49393 49394 49395

gpl_filtered = gpl[-c(49301:49395), ] # Ordered rows in gpl.
not_ordered_rows_gpl = gpl[c(49301:49395), ]
dim(gpl_filtered) # [1] 49400    24

exprinfo_filtered = exprinfo[-c(49301:49395), ]
not_ordered_rows_exprinfo = as.data.frame(exprinfo[c(49301:49395),])
dim(exprinfo_filtered) # [1] 49400     5

# Ordered positions are the same:
all(gpl_filtered$ID == rownames(exprinfo_filtered))

# Add a column with probe id to order the disordered rows:
not_ordered_rows_exprinfo$ID = rownames(not_ordered_rows_exprinfo)

rows_to_unit = full_join(x = not_ordered_rows_gpl,
                     y = not_ordered_rows_exprinfo,
                     by = 'ID')
rownames(rows_to_unit) = rows_to_unit$ID
dim(rows_to_unit) # 95 rows to add in order

gpl_filtered = rbind(gpl_filtered, rows_to_unit[, 1:24])
dim(gpl_filtered) # [1] 49495    24

exprinfo_filtered = rbind(exprinfo_filtered, rows_to_unit[, 25:29])
dim(exprinfo_filtered) # [1] 49495     5

table(gpl_filtered$ID == rownames(exprinfo_filtered))

#### experimentData ####
infoData = new('MIAME',
        name = 'Guo et al.',
        lab = 'Ruijin Hospital, Shanghai Jiaotong University',
        contact ='Mingjie Wang <huzai920621@126.com>',
        title = ' 	Characterization of gene expression profile during the
        implantation phase in IVF cycles in endometrium of women with recurrent
        implantation failure',
        abstract = 'Our study was designed with an aim to compare the gene
        expression profile of RIF women with that of healthy proven fertile
        controls to explore the influence factors of endometrial receptivity and
        the mechanisms of RIF.',
        url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103465')

experimentData(gse103465) = infoData
experimentData((gse103465))

# SAVE NORMALIZED AND ANNOTATED DATA --------------------------------------
save(gse103465, file = paste0(data_path, '/gse103465.rda'), version = 2)

load(paste0(data_path, '/gse103465.rda'), verbose = T)

# CREATE ED AND DAT -------------------------------------------------------
#### Removing NAs from annotation ####
gpl_filtered_notNA = gpl_filtered[!is.na(gpl_filtered$Gene.Symbol), ]
exprinfo_filtered_notNA = exprinfo_filtered[!is.na(gpl_filtered$Gene.Symbol), ]
dim(gpl_filtered_notNA) # [1] 49372    24
dim(exprinfo_filtered_notNA) # [1] 49372     5

#### Group expression info by gene symbol ####
dat = avereps(exprinfo_filtered_notNA, ID = gpl_filtered_notNA$Gene.Symbol)

# Correct sample name from dat:
colnames(dat) = gsub(pattern = '_.*', replacement = '', x = colnames(dat))
head(dat)

#### We create the ed ####
pData(gse103465)

ed = pData(gse103465) # 'ed' after experimental dessign
# rownames(ed) = colnames(dat)

# SAVE .RDA ---------------------------------------------------------------
save(dat, ed, file = paste0(results_path, '/guo.rda'), version = 2)

# CHECK DAT AND ED ----------------------------------------------------------
load(paste0(results_path, '/guo.rda'), verbose = T)
head(dat)
head(ed)

plot_PCAscores(dat = dat, ed = ed,
               condition1 = 'condition',
               components = c(1,2),
               colors = c('#D95F02', '#1B9E77'),
               title = 'PCA: Guo (GSE103465)')
