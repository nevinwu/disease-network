################################################################################
## node_analysis.r
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

analysis_path = '/Users/francisco/Desktop/TFM/network/analysis'

# LOAD PACKAGES -----------------------------------------------------------

# ANALYSIS ----------------------------------------------------------------
network_analysis = read.csv(paste0(analysis_path, '/network_analysis.csv'))

# BETWEENNESS CENTRALITY
boxplot(network_analysis[, 'BetweennessCentrality'])

summary(network_analysis[, 'BetweennessCentrality'])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000000 0.0000000 0.0002836 0.0048863 0.0039543 0.2364035

IQR(network_analysis[, 'BetweennessCentrality'])
# 0.003954346

# So, our relative maximum would be:
bc_rel_max = 0.0039543 + 1.5 * 0.003954346
# 0.009885819

table(network_analysis[, 'BetweennessCentrality'] > bc_rel_max)
# 81 genes

# DEGREE
boxplot(network_analysis[, 'Degree'])

summary(network_analysis[, 'Degree'])
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 1.000   1.000   2.000   3.304   3.000  72.000

IQR(network_analysis[, 'Degree'])
# 2

# So, our relative maximum would be:
d_rel_max = 3 + 1.5 * 2
# 6

table(network_analysis[, 'Degree'] > d_rel_max)
# 80 genes

# IN COMMON
bc_above = network_analysis[network_analysis[, 'BetweennessCentrality'] > bc_rel_max,]
d_above = network_analysis[network_analysis[, 'Degree'] > d_rel_max,]

bc_names = bc_above[, 'name']
d_names = d_above[, 'name']

table(bc_names%in%d_names)
# FALSE  TRUE
#    25    56

common_names = bc_names[bc_names%in%d_names]
# We have 56 common genes above relative maximum in betweenness centrality and
# degree.

positions = network_analysis[, 'name']%in%common_names
target_genes = network_analysis[positions,]


targets = target_genes[order(target_genes$BetweennessCentrality, target_genes$Degree, decreasing = T),]


# OUTPUT ------------------------------------------------------------------
file = paste0(analysis_path, '/target_genes.csv')

write.csv(targets, file = file)
