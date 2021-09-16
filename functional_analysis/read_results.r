################################################################################
## read_results.r
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

# Biological Process - Overrepresented:
bp_over = read.delim(file = '/Users/francisco/Desktop/TFM/network/functional_analysis/network/bp_overrep_disease_network.bgo',
           skip = 21, header = T)

# Biological Process - Underrepresented:
bp_under = read.delim(file = '/Users/francisco/Desktop/TFM/network/functional_analysis/network/bp_underrep_disease_network.bgo',
           skip = 21, header = T)
# x_<cluster name> : the number of genes in your cluster annotated to a certain GO class
# X_<cluster name> : the total number of genes in your cluster. This number may be different from the number of genes you selected in the graph or put into the text field, since genes without any annotation are discarded (see FAQ)
# n_<cluster name> : the number of genes in the reference set (graph or annotation) annotated to a certain GO class
# N_<cluster name> : the total number of genes in your reference set. This number may be different from the number of genes in your reference graph, since genes without any annotation are discarded (see FAQ)

file = '/Users/francisco/Desktop/TFM/network/functional_analysis/network/bp_over_table.csv'
write.csv(bp_over, file = file)

# Read REVIGO:

revigo = read.delim(file = '/Users/francisco/Desktop/TFM/network/functional_analysis/network/revigo.csv',
                    header = T, sep = ',')
dim(revigo)
# [1] 227  11

condensed_revigo = revigo[revigo[,'Eliminated'] == ' False',]
dim(condensed_revigo)
# [1] 102  11
