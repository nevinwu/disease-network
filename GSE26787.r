################################################################################
## ledee_GSE26787.r
## 2021-02-15
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

# Dataset info avaliable in:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26787

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

# LOAD PACKAGES -----------------------------------------------------------
library(Biobase)
library(GEOquery)
library(affy)
library(affyPLM)
library(ggplot2)
library(RColorBrewer)
library(AnnotationDbi)
library(cluster)
library(ggrepel)
library(plyr)
library(hgu133plus2.db)
library(sva)
library(limma)

functions_path = '/Users/francisco/Desktop/TFM/functions'
dataset_path = '/Users/francisco/Desktop/TFM/datasets/GSE26787_ledee'
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

# DOWNLOAD RAW DATA FOR LEDEE (GSE26787) ----------------------------------
getGEOSuppFiles('GSE26787')
system('tar xvf GSE26787/GSE26787_RAW.tar')

gse26787raw = ReadAffy()

system('rm -fr GSE26787')
system('rm *.gz')

# SAVE RAW DATA -----------------------------------------------------------
save(gse26787raw, file = './gse26787raw.rda', version = 2)

load(paste0(data_path, '/gse26787raw.rda'), verbose = T) # raw data

# PRE-NORMALIZATION ANALYSIS ----------------------------------------------
# GSE = getGEO('GSE26787') # series_matrix

# Remove RSA since its not our targetted condition
gse26787raw = gse26787raw[, 1:10]

# dim(exprs(gse26787raw)) # 1.354.896 probes y 10 samples
# annotation(gse26787raw) # chip 'hgu133plus2'
head(exprs(gse26787raw))
# colnames(gse26787raw) = gsub(x = colnames(gse26787raw),
#                               pattern = '.CEL.gz', replacement = '')

#### MDplot ####
old_rownames = rownames(pData(gse26787raw))
fancy_rownames = gsub(x = colnames(gse26787raw),
                       pattern = '.CEL.gz', replacement = '')

rownames(pData(gse26787raw)) = fancy_rownames

affy::MAplot(gse26787raw, type = 'pm', plot.method = 'smoothScatter')

#### Density estimator ####
affy::hist(gse26787raw)
legend(x = 'topright',          # Position
       legend = rownames(pData(gse26787raw)),
       lty = c(1:15),           # Line types
       col = c(1:15),           # Line colors
       lwd = 2)                # line width

#### Boxplot ####
affy::boxplot(gse26787raw)

#### Clustering ####
datos_raw = exprs(gse26787raw)
colnames(datos_raw) = rownames(pData(gse26787raw))

# Make a data.frame to colour our groups
group.raw = c(rep('control', 5), # Fertile
              rep('RIF', 5)) # Implantatio failure
sample.name = colnames(datos_raw)
sinfo_raw = as.data.frame(cbind(sample.name, group.raw), as.is = T)
char.group = as.character(sinfo_raw[, 'group.raw'])
tan = unique(char.group)
micolor = rainbow(length(tan))
names(micolor) = tan
sinfo_raw$color.group = micolor[char.group]
sinfo_raw

# Clustering by correlation distance
correlacion_raw = cor(datos_raw)
distancia_raw = as.dist((1 - correlacion_raw) / 2)
hc_raw = hclust(distancia_raw)

table(hc_raw$labels ==  colnames(datos_raw)) # same layout?

hc_raw$clase = sinfo_raw$group # Colour

# Plot using 'function_arbol.r'
arbol(cluster = hc_raw, main = 'Clustering by Correlation Distance')

# RAW PCA -----------------------------------------------------------------
pca_raw = prcomp(t(log2(exprs(gse26787raw)) + 1))
# t cause we need samples in rows

# screeplot(pca, type = 'l', npcs = 6, main = 'Screeplot of the first 6 PCs')
# abline(h = 1, col = 'red', lty = 5)
# legend('topright', legend = c('Eigenvalue = 1'),
#        col = c('red'), lty = 5, cex = 0.6)

var_raw = round(summary(pca_raw)$importance[2, c(1,2)] * 100, 1) # rounds 1st d
toplot_raw = data.frame(pca_raw$x[, c(1,2)], stringsAsFactors = F) ###
lim_raw = max(abs(c(min(toplot_raw), max(toplot_raw)))) # axis values limit
axis_limits_raw = c(-lim_raw, lim_raw)
toplot_raw$color = c(rep('control', 5), rep('RIF', 5))

# png(filename = paste0(raw_plots_path, '/PCAplot_group.png'), # File name
#     width = 600, height = 600, # Width and height in pixels
#     unit = 'px', # Width and height units, pixels
#     bg = 'white') # Background color

ggplot(data = toplot_raw,
       aes_string(x = colnames(toplot_raw)[1], y = colnames(toplot_raw)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
#       geom_text_repel(label = rownames(pData(gse26787raw)), size = 3) +
       xlab(paste0('PC1', ': ', var_raw[1], '%')) +
       ylab(paste0('PC2', ': ', var_raw[2], '%')) +
       ggtitle('PCA: Ledee (GSE26787)') +
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

arrayQualityMetrics(expressionset = gse26787raw,
                    outdir = outdir,
                    force = TRUE,
                    do.logtransform = TRUE) # Since data is not processed yet.

# Check index.html file in outdir for results.

# NORMALIZATION -----------------------------------------------------------
rownames(pData(gse26787raw)) = old_rownames # needed for rma work properly

gse26787 = affy::rma(gse26787raw)

# load(paste0(data_path, '/gse26787.rda'), verbose = T)

# POST-NORMALIZATION ANALYSIS ---------------------------------------------
# colnames(exprs(gse26787)) = c(rep(paste0('control-', 1:5)),
#                               rep(paste0('RIF-', 1:5)))
# length(colnames((exprs(gse26787))))

#### MDplot ####
old_rownames = rownames(pData(gse26787))
fancy_rownames = gsub(x = colnames(gse26787),
                       pattern = '.CEL.gz', replacement = '')

rownames(pData(gse26787)) = fancy_rownames

affy::MAplot(gse26787, plot.method = 'smoothScatter')

#### Density estimator ####
affy::hist(gse26787)
legend(x = 'topright',          # Position
       legend = rownames(pData(gse26787)),  # Legend texts
       lty = c(1:6),           # Line types
       col = c(1:6),           # Line colors
       lwd = 2)                # line width

#### Boxplot ####
affy::boxplot(gse26787)

#### Clustering ####
datos_norm = exprs(gse26787)
colnames(datos_norm) = rownames(pData(gse26787))

# Make a data.frame to colour our groups
group.norm = c(rep('control', 5), rep('RIF', 5))
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
pca_norm = prcomp(t(exprs(gse26787)))

# https://towardsdatascience.com/principal-component-analysis-pca-101-using-r
# consider 'center = TRUE, scale = TRUE'
# screeplot(pca, type = 'l', npcs = 6, main = 'Screeplot of the first 6 PCs')
# abline(h = 1, col = 'red', lty = 5)
# legend('topright', legend = c('Eigenvalue = 1'),
#        col = c('red'), lty = 5, cex = 0.6)

var_norm = round(summary(pca_norm)$importance[2, c(1,2)] * 100, 1)
toplot_norm = data.frame(pca_norm$x[, c(1,2)], stringsAsFactors = F)
lim_norm = max(abs(c(min(toplot_norm), max(toplot_norm))))
axis_limits_norm = c(-lim_norm, lim_norm)
toplot_norm$color = c(rep('control', 5), rep('RIF', 5))

# png(filename = paste0(normal_plots_path, '/PCAplot_group.png'), # File name
#     width = 600, height = 600, # Width and height in pixels
#     unit = 'px', # Width and height units, pixels
#     bg = 'white') # Background color

ggplot(data = toplot_norm,
       aes_string(x = colnames(toplot_norm)[1], y = colnames(toplot_norm)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values = c('#D95F02', '#1B9E77')) +
#       geom_text_repel(label = rownames(pData(gse26787)), size = 3) +
       xlab(paste0('PC1', ': ', var_norm[1], '%')) +
       ylab(paste0('PC2', ': ', var_norm[2], '%')) +
       ggtitle('PCA: Ledee (GSE26787)') +
       xlim(axis_limits_norm ) + ylim(axis_limits_norm) +
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
factores = c(rep(0, 5), rep(1, 5))

condition = factor(factores, levels = 0:1, labels = c('control','RIF'))

infor_fenotip = data.frame(pData(gse26787), condition)
pData(gse26787) = infor_fenotip
pData(gse26787)$sample = rownames(pData(gse26787))
pData(gse26787)

#### fData ####
annotation(gse26787) # hgu133plus2

probeInfo = AnnotationDbi::select(hgu133plus2.db, keys = featureNames(gse26787),
                                  columns = c('ENTREZID', 'ENSEMBL', 'SYMBOL'),
                                  keytype = 'PROBEID')
# Warning: 'select()' returned 1:many mapping between keys and columns

dim(probeInfo) # [1] 64522     4
dim(exprs(gse26787)) # [1] 54675    10

# When multiple probe-gen matches, we keep only the first match:
posiciones = match(unique(probeInfo[,1]), probeInfo[,1]) # PROBEID
probeInfo = probeInfo[posiciones,]
posiciones = match(unique(probeInfo[,2]), probeInfo[,2]) # ENTREZ
probeInfo = probeInfo[posiciones,]

probeInfo = probeInfo[!is.na(probeInfo[,'ENTREZID']),]
dim(probeInfo) # [1] 20947     4
head(probeInfo, n = 100)

# To keep the order between eset and annotation
probes_in_eset = data.frame('PROBEID' = rownames(exprs(gse26787)),
                             stringsAsFactors = F)
sum(duplicated(probes_in_eset))

annotation = join(x = probes_in_eset, y = probeInfo, by = 'PROBEID')
dim(annotation) # [1] 54675     4
dim(exprs(gse26787)) # [1] 54675    10

fData(gse26787) = annotation
head(fData(gse26787))

# Care since yet we have NAs in annotation... Removed later before ed & dat.

#### experimentData ####
infoData = new('MIAME',
        name = 'Ledee et al.',
        lab = 'implantation et dialogue mother-conceptus',
        contact ='Nathalie Ledee <nathalie-ledee@orange.fr>',
        title = 'Specific and extensive endometrial deregulation is present
        before conception in IVF/ICSI repeated implantation failures (IF) or
        recurrent miscarriages',
        abstract = 'Summary: In order to identify pre-conceptional endometrial
        dysregulations, we compared the endometrial expression between fertile
        and IF and RM patients.',
        url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26787')

experimentData(gse26787) = infoData
experimentData((gse26787))

# SAVE NORMALIZED AND ANNOTATED DATA --------------------------------------
save(gse26787, file = paste0(data_path, '/gse26787.rda'), version = 2)

load(paste0(data_path, '/gse26787.rda'), verbose = T)

# CREATE ED AND DAT -------------------------------------------------------
exprinfo = exprs(gse26787)

#### Removing NAs ####
indexNA = which(is.na(annotation$SYMBOL))
annotation[is.na(annotation$SYMBOL), ]
exprinfo = exprinfo[-indexNA, ]
annotation = annotation[-indexNA, ]
dim(annotation)
dim(exprinfo)

#### Group expression info by gene symbol ####
dat_time = avereps(exprinfo, ID = annotation$SYMBOL)
head(dat_time)

#### We create ed ####
pData(gse26787)
colnames(dat_time) = gsub(pattern = '.CEL.gz',replacement = '',
                             x = colnames(dat_time))
ed = pData(gse26787)
rownames(ed) = colnames(dat_time)

ed$sample = rownames(ed)

# SAVE .RDA ---------------------------------------------------------------
save(dat_time, ed, file = paste0(data_path, '/ledee_time.rda'), version = 2)

load(paste0(data_path, '/ledee_time.rda'), verbose = T)

# CORRECT BATCH EFFECT AND TIME EFFECT ------------------------------------
# At the beginningn, we removed RIF-1 (GSM659108) sample according to
# normalized PCA.
# But, finally, we won't do it since it doesn't seem to be an outlier.

#### Remove outliers and clustering ####
ed
head(dat_time)

# GSM659108 = RIF-1
# dat_time = dat_time[, colnames(dat_time) != 'GSM659108']
# ed = ed[rownames(ed) !='GSM659108',]

#### Show normalized PCA: ####
# Function 'plot_PCAscores' obtained from FIVI:
plot_PCAscores(dat = dat_time, ed = ed, condition1 = 'condition',
               colors = c('#D95F02', '#1B9E77'),
               title = 'PCA: Ledee (GSE26787)')

# We'want to remove the unknow batch effect that seems to affect our data
# (probably time-effect)t:

# At the beginning, we used non-supervised clustering as follows:
# (but finally we will use sva):

#### Calculate how many clusters to use ####
# We want to know how many groups should we use in clustering. We can use this
# function made in FIVI:
# optimal_nclusters(dat = as.matrix(dat_time),
#                   method = 'silhouette',
#                   max_nclusters = 5)
#  Error: number of cluster centres must lie between 1 and nrow(x)
# Doesn't work, error.

# sil = rep(0, 8)
# for (i in 2:8) {
#         km.res = kmeans(t(dat_time), centers = i, nstart = 25)
#         ss = silhouette(km.res$cluster, dist(t(dat_time)))
#         sil[i] = mean(ss[, 3])
#         print(i)
# }
# silPlot = data.frame(sil = sil, clusters = 1:8)
# ggplot(silPlot, aes(x = factor(clusters), y = sil, fill = '1')) +
#         geom_point(shape = 1, size = 4) +
#         geom_line(data = silPlot, aes(x = clusters, y = sil), size = 1,
#                   colour = 'azure4') +
#         geom_vline(xintercept = silPlot$clusters[which.max(silPlot$sil)],
#                    linetype = 2, color = 'steelblue') + theme_bw() +
#         xlab('Number of clusters (k)') +
#         ylab('Average silhouette coefficient') +
#         ggtitle('Silhoutte method') +
#         theme(legend.position = 'none',
#               plot.title = element_text(size = 22, hjust = 0.5),
#               axis.text = element_text(size = 15, color = 'dimgray'),
#               axis.title = element_text(size = 18, hjust = 0.5))
#
# # optimal number of clusters = 2
#
# #### Unsupervised clustering ####
# clustering = unsupervised_clustering(dat = dat_time,
#                                      method = 'kmeans',
#                                      nclusters = 2)
# clustering$cluster_results
# clustering$cluster_group
# clustering$cluster_plot
# ## Problem with cluster 1:
# # Too few points to calculate an ellipse
# # Warning message:
# # Use of `toplot$color` is discouraged. Use `color` instead.
#
# #### Add the clustering classfication to ed ####
# # all(names(clustering$cluster_group) == ed$sample)
# # [1] TRUE
# ed$cluster = as.character(clustering$cluster_group)
#
# #### Plot PCA with new condition ####
# # Colours obtained from: brewer.pal(n = 6, 'Dark2')
# plot_PCAscores(dat = dat_time, ed = ed, condition1 = 'condition',
#                condition2 = 'cluster', components = c(1,2),
#                colors = c('#D95F02', '#1B9E77'),
#                title = 'PCA: Ledee (GSE26787)')
# # check clusters visually
#
# #### Remove batch effect ####
# mdesign = model.matrix(~ed$condition)
#
# # CARE: we modify dat to remove time effect.
# dat = removeBatchEffect(x = dat_time,
#                         design = mdesign,
#                         batch = ed$cluster)
#
# plot_PCAscores(dat = dat, ed = ed, condition1 = 'condition',
#                components = c(1,2),
#                colors = c('#D95F02', '#1B9E77'),
#                title = 'PCA: Ledee (GSE26787)')


# SVA FOR UNKNOWN BATCH-EFFECT CORRECTION ---------------------------------
# mod = model.matrix( ~ 0 + ed$condition)
# n = num.sv(dat = dat_time, mod = mod, method = 'leek') # 4
#
# mod0 = model.matrix( ~ 1, data = ed)
# svobj = sva(dat = dat_time, mod = mod, mod0 = mod0, n.sv = n)
#
# dat = removeBatchEffect(x = dat_time, covariates = svobj$sv)
#
# plot_PCAscores(dat = dat, ed = ed, condition1 = 'condition',
#                components = c(1, 2),
#                title = 'PCA: Ledee (GSE26787)',
#                colors =  c('#D95F02', '#1B9E77'))

# SAVE DAT & ED -----------------------------------------------------------
save(dat, ed, file = paste0(results_path, '/ledee.rda'), version = 2)

# CHECK DAT & ED ----------------------------------------------------------
load(paste0(results_path, '/ledee.rda'), verbose = T)
head(dat)
head(ed)
