################################################################################
## join_datasets.r
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

# KOOT INCLUDED:

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

# LOAD PACKAGES -----------------------------------------------------------
library(limma)
library(venn)
# library(gsrmUtils)

functions_path = '/Users/francisco/Desktop/TFM/functions'
join_path = '/Users/francisco/Desktop/TFM/datasets/join_datasets'
results_path = '/Users/francisco/Desktop/TFM/datasets/results/results_sva_2'

getwd()
setwd(join_path)

source(paste0(functions_path, '/function_plot_PCAscores.r'))
source(paste0(functions_path, '/function_diffExprAnalysis.r'))

# LOAD DATASETS (DAT & ED) ------------------------------------------------
files_koot = dir(results_path, full.names = T)

## DAT OBJECTS (expression data):
dats2_koot = lapply(files_koot, function(f){
  cat(f, '\n')
  load(f,verbose = T)
  return(dat)
})

names(dats2_koot) = gsub(pattern = '.rda', replacement = '',
                         x = basename(files_koot))

str(dats2_koot)

# How many genes in common between datasets?
common_genes_koot = Reduce('intersect', lapply(dats2_koot, rownames))
length(common_genes_koot) # [1] 12156

dats2_koot = do.call('cbind', lapply(dats2_koot, function(x){
  x[common_genes_koot, ]
}))

dim(dats2_koot) # [1] 12156   195
# 12.156 genes (in common) in all datasets and 195 samples.

## ED OBJECTS (experimental dessigns):
eds2_koot = do.call('rbind', lapply(files_koot, function(f){
  load(f)
  ed = data.frame(sample = rownames(ed),
                  condition = ed$condition,
                  experiment = gsub(pattern = '.rda',replacement = '',
                                   x = basename(f)),
                  stringsAsFactors = F)
  return(ed)
}))

rownames(eds2_koot) = eds2_koot$sample
str(eds2_koot)

# SAMPLES AMONG DATASETS --------------------------------------------------
# We expect samples to be clustered in datasets. We should remove samples that
# are not clustered in datasets.
plot_PCAscores(dat = dats2_koot, ed = eds2_koot, condition1 = 'experiment',
               title = 'PCA: muestras agrupadas por datasets')

# REMOVE EXPERIMENT EFFECT ------------------------------------------------
design_koot = model.matrix(~condition, data = eds2_koot)
dats3_koot = removeBatchEffect(dats2_koot, batch = eds2_koot$experiment)

plot_PCAscores(dat = dats3_koot, ed = eds2_koot, condition1 = 'experiment',
               title = 'PCA: eliminado el efecto del experimento')

# CHECK DEGs AMONG CONTROLS -----------------------------------------------
control_dat = dats3_koot[, eds2_koot$condition == 'control']
control_ed = eds2_koot[eds2_koot$condition == 'control', ]

DEGs_control = diffExprAnalysis(dat = control_dat, ed = control_ed,
                                condition = 'experiment')

sigs = lapply(DEGs_control, function(deg){
  rownames(deg)[deg$adj.P.Val <= 0.05]
})

str(sigs)
to_remove_control = unique(unlist(sigs))

# Significantly differentially expressed genes among controls:
length(to_remove_control) # [1] 2628

# CHECK DEGs AMONG RIF ----------------------------------------------------
rif_dat = dats3_koot[, eds2_koot$condition == 'RIF']
rif_ed = eds2_koot[eds2_koot$condition == 'RIF', ]

DEGs_rif = diffExprAnalysis(dat = rif_dat, ed = rif_ed,
                            condition = 'experiment')

sigs = lapply(DEGs_rif, function(deg){
  rownames(deg)[deg$adj.P.Val <= 0.05]
})

str(sigs)
to_remove_rif = unique(unlist(sigs))

# Significantly differentially expressed genes among RIF:
length(to_remove_rif) # [1] 2154

# VENN DIAGRAM (DEGs CONTROL & DEGs RIF) ----------------------------------
venn(list(Control = to_remove_control, RIF = to_remove_rif))
# 1632 in common

# REMOVE DEGs IN CONTROLS & RIF -------------------------------------------
to_remove = unique(c(to_remove_control, to_remove_rif))
length(to_remove) # [1] 3150

dats3_koot = dats3_koot[!rownames(dats3_koot)%in%to_remove, ]
dim(dats3_koot) # 9006 genes & 195 samples

plot_PCAscores(dat = dats3_koot, ed = eds2_koot, condition1 = 'condition',
               components = c(1,2),
               title = 'PCA: control vs RIF')

# DEGs BETWEEN CONDITIONS -------------------------------------------------
DEGs_network_koot = diffExprAnalysis(dat = dats3_koot, ed = eds2_koot,
                        condition = 'condition')

head(DEGs_network_koot$`control-RIF`)

sigs_DEGs_network_koot = lapply(DEGs_network_koot, function(deg){
  rownames(deg)[deg$adj.P.Val <= 0.05]
})

length(sigs_DEGs_network_koot$`control-RIF`) # 436 DEGs

# CHECK GENES AMONG EXPERIMENTS -------------------------------------------
DEGs = diffExprAnalysis(dat = dats3_koot, ed = eds2_koot,
                        condition = 'experiment')

str(DEGs)

sigs = lapply(DEGs, function(deg){
  rownames(deg)[deg$adj.P.Val <= 0.05]
})
str(sigs)

to_remove = unique(unlist(sigs))
length(to_remove) # 0 --> No genes to be removed


################################################################################
################################################################################

# KOOT NOT INCLUDED:

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

# LOAD PACKAGES -----------------------------------------------------------
library(limma)
library(venn)
library(ggplot2)
# library(gsrmUtils)

functions_path = '/Users/francisco/Desktop/TFM/functions'
join_path = '/Users/francisco/Desktop/TFM/datasets/join_datasets'
results_path = '/Users/francisco/Desktop/TFM/datasets/results/results_sva_2'

getwd()
setwd(join_path)

source(paste0(functions_path, '/function_plot_PCAscores.r'))
source(paste0(functions_path, '/function_diffExprAnalysis.r'))

# LOAD DATASETS (DAT & ED) ------------------------------------------------
files = dir(results_path, full.names = T)
files = files[-4] # Remove koot since it provides much noise

## DAT OBJECTS (expression data):
dats2 = lapply(files, function(f){
  cat(f, '\n')
  load(f,verbose = T)
  return(dat)
})

names(dats2) = gsub(pattern = '.rda', replacement = '', x = basename(files))

str(dats2)

# How many genes in common between datasets?
common_genes = Reduce('intersect', lapply(dats2, rownames))
length(common_genes) # [1] 12358

dats2 = do.call('cbind', lapply(dats2, function(x){
  x[common_genes, ]
}))

dim(dats2) # [1] 12358   80
# 12.358 genes (in common) in all datasets and 80 samples.

## ED OBJECTS (experimental dessigns):
eds2 = do.call('rbind', lapply(files, function(f){
  load(f)
  ed = data.frame(sample = rownames(ed),
                  condition = ed$condition,
                  experiment = gsub(pattern = '.rda',replacement = '',
                                   x = basename(f)),
                  stringsAsFactors = F)
  return(ed)
}))

rownames(eds2) = eds2$sample
str(eds2)
table(eds2$condition)
table(eds2$experiment)

# SAMPLES AMONG DATASETS --------------------------------------------------
# We expect samples to be clustered in datasets. We should remove samples that
# are not clustered in datasets.
plot_PCAscores(dat = dats2, ed = eds2, condition1 = 'experiment',
               title = 'PCA: muestras agrupadas por datasets (sin koot)')

#############
pca = prcomp(t(dats2))

var_norm = round(summary(pca)$importance[2, c(1,2)] * 100, 1)
toplot_norm = data.frame(pca$x[, c(1,2)], stringsAsFactors = F)
lim_norm = max(abs(c(min(toplot_norm), max(toplot_norm))))

axis_limits_norm = c(-lim_norm, lim_norm)
toplot_norm$color = eds2$experiment

# To fit the title:
wrapper = function(x, ...)
{
paste(strwrap(x, ...), collapse = '\n')
}

my_title = 'PCA: muestras agrupadas por datasets (sin koot)'
# + ggtitle(wrapper(my_title, width = 20))

ggplot(data = toplot_norm,
       aes_string(x = colnames(toplot_norm)[1], y = colnames(toplot_norm)[2])) +
       geom_point(aes(color = color), size = 3) +
       scale_color_manual(name = 'RIF', values = c('#E69F00', '#56B4E9',
                                                   '#009E73', '#F0E442',
                                                   '#0072B2')) +
#       geom_text_repel(label = ed$sample, size = 3) +
       xlab(paste0('PC1', ': ', var_norm[1], '%')) +
       ylab(paste0('PC2', ': ', var_norm[2], '%')) +
       ggtitle(wrapper(my_title, width = 35)) +
       xlim(axis_limits_norm) + ylim(axis_limits_norm) +
       theme_light() +
       theme(legend.position = 'bottom',
             axis.title = element_text(size = 18),
             axis.text = element_text(size = 15),
             plot.title = element_text(size = 22, hjust = 0.5),
             legend.title = element_blank(),
             legend.text = element_text(size = 13))

# REMOVE EXPERIMENT EFFECT ------------------------------------------------
design = model.matrix(~condition, data = eds2)
dats3 = removeBatchEffect(dats2, batch = eds2$experiment)

plot_PCAscores(dat = dats3, ed = eds2, condition1 = 'experiment',
               title = 'PCA: eliminado el efecto del experimento')

# CHECK DEGs AMONG CONTROLS -----------------------------------------------
control_dat = dats3[, eds2$condition == 'control']
control_ed = eds2[eds2$condition == 'control', ]

DEGs_control = diffExprAnalysis(dat = control_dat, ed = control_ed,
                                condition = 'experiment')

sigs = lapply(DEGs_control, function(deg){
  rownames(deg)[deg$adj.P.Val <= 0.05]
})

str(sigs)
to_remove_control = unique(unlist(sigs))

# Significantly differentially expressed genes among controls:
length(to_remove_control) # [1] 88

# CHECK DEGs AMONG RIF ----------------------------------------------------
rif_dat = dats3[, eds2$condition == 'RIF']
rif_ed = eds2[eds2$condition == 'RIF', ]

DEGs_rif = diffExprAnalysis(dat = rif_dat, ed = rif_ed,
                            condition = 'experiment')

sigs = lapply(DEGs_rif, function(deg){
  rownames(deg)[deg$adj.P.Val <= 0.05]
})

str(sigs)
to_remove_rif = unique(unlist(sigs))

# Significantly differentially expressed genes among RIF:
length(to_remove_rif) # [1] 129

# VENN DIAGRAM (DEGs CONTROL & DEGs RIF) ----------------------------------
venn(list(Control = to_remove_control, RIF = to_remove_rif))
# 36 in common

# REMOVE DEGs IN CONTROLS & RIF -------------------------------------------
to_remove = unique(c(to_remove_control, to_remove_rif))
length(to_remove) # [1] 181

dats3 = dats3[!rownames(dats3)%in%to_remove, ]
dim(dats3) # 12177 genes & 80 samples

plot_PCAscores(dat = dats3, ed = eds2, condition1 = 'condition',
               components = c(1,2),
               title = 'PCA: control vs RIF')

# DEGs BETWEEN CONDITIONS -------------------------------------------------
DEGs_network = diffExprAnalysis(dat = dats3, ed = eds2, condition = 'condition')

head(DEGs_network$`control-RIF`)

sigs_DEGs_network = lapply(DEGs_network, function(deg){
  rownames(deg)[deg$adj.P.Val <= 0.05]
})

length(sigs_DEGs_network$`control-RIF`) # 2152 DEGs

# Differential analysis
DEGs = diffExprAnalysis(dat = dats3, ed = eds2, condition = 'experiment')

str(DEGs)

sigs = lapply(DEGs, function(deg){
  rownames(deg)[deg$adj.P.Val <= 0.05]
})
str(sigs)

to_remove = unique(unlist(sigs))
length(to_remove) # 0 --> No genes to be removed


################################################################################

# sigs_DEGs_network_koot includes sig DEGs when Koot included
network_koot_v2 = sigs_DEGs_network_koot$`control-RIF`
length(network_koot_v2) # 436 DEGs when Koot included

# sigs_DEGs_network includes sig DEGS when koot NOT included
network_v2 = sigs_DEGs_network$`control-RIF`
length(network_v2) # 2152 DEGs when koot NOT included

network_koot_v2%in%network_v2 # returns TRUE/FALSE if genes in network koot are in network

core_DEGs_v2 = network_koot_v2[network_koot_v2%in%network_v2]
core_DEGs_v2
length(core_DEGs_v2) # 368

venn(list(Koot_included = network_koot_v2, Koot_NOT_included = network_v2))

save(network_koot_v2, network_v2, core_DEGs_v2,
     file = '/Users/francisco/Desktop/network_genes_v2.rda', version = 2)

load('/Users/francisco/Desktop/network_genes_v2.rda', verbose = T)

venn(list(core_genes_v1 = core_DEGs_v1, core_genes_v2 = core_DEGs_v2))

save(DEGs_network, )


# SAVING FOR FUNCTIONAL ANALYSIS ------------------------------------------

save(DEGs_network, sigs_DEGs_network,
     file = '/Users/francisco/Desktop/DEGs_network_v2.rda', version = 2)

# DEGs_network es el resultado de la expresión diferencial (tanto aquellos que
# la tienen como los que no), y sigs_DEGs_network contiene sólo aquellos genes
#con expresión diferencial significativa.
