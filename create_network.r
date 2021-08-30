################################################################################
## create_network.r
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

functions_path = '/Users/francisco/Desktop/TFM/functions'
results_path = '/Users/francisco/Desktop/TFM/datasets/results_sva_2'
network_path = '/Users/francisco/Desktop/TFM/network'

# LOAD PACKAGES -----------------------------------------------------------
# library(gsrmUtils)
library(data.table)
library(venn)

# LOAD DATA ---------------------------------------------------------------
load(paste0(network_path, '/network_genes_v2.rda'), verbose = T)
# Loading objects:
#   network_koot_v2
#   network_v2
#   core_DEGs_v2
# We will use object network_v2 since it is the network WITHOUT including Koot.

degs_koot = network_koot_v2
degs_nokoot = network_v2

# LOAD INTERACTOME
load(paste0(network_path, '/interactome_final.rda'),verbose = T)
# Loading objects:
# interactome
# interactome_genename

# CHECK GENES/DEGS:
length(degs_nokoot) # 2152
length(degs_koot) # 436

venn(list('With Koot' = degs_koot, 'Without_Koot' = degs_nokoot))
# There are 368 genes in common

# INTERACTOME -- -----------------------------------------------------------

# WITH KOOT
venn(list('Interactome' = unique(c(interactome_genename$node1, interactome_genename$node2)),
          'With Koot' = degs_koot))
# We lose 138 genes, maybe we should seek synonyms.

network_koot = interactome_genename[interactome_genename$node1%in%degs_koot &
                                       interactome_genename$node2%in%degs_koot,]

dim(network_koot)
# [1] 38  2
length(unique(c(network_koot$node1, network_koot$node2))) # 61

# WITHOUT KOOT
venn(list('Interactome' = unique(c(interactome_genename$node1, interactome_genename$node2)),
          'Without Koot' = degs_nokoot))
# We lose 693 genes, maybe we should seek synonyms.

network_nokoot = interactome_genename[interactome_genename$node1%in%degs_nokoot &
                                         interactome_genename$node2%in%degs_nokoot,]

dim(network_nokoot)
# [1] 1296    2
length(unique(c(network_nokoot$node1, network_nokoot$node2))) # 824

write.table(x = network_nokoot,
            file = paste0(network_path,
                          '/network_noKoot_v2.tsv'),
            quote = F,sep = '\t',
            row.names = F)

# CREATE CORE GENES ATTRIBUTE ---------------------------------------------
# Create column genes with all genes in network
length(unique(c(network_nokoot$node1, network_nokoot$node2))) # 824 genes
genes = unique(c(network_nokoot$node1, network_nokoot$node2))
genes

# check how many core genes are in network
table(core_DEGs_v2%in%genes) # 143 core genes are in network
length(core_DEGs_v2)# 368

# Prepare column core (TRUE if gene is core, FALSE otherwise)
core = rep(FALSE, length(genes))
length(core)

# Create data.frame
core_attribute = data.frame(genes, core, row.names = genes)
dim(core_attribute)

core_attribute[,'core'] = rownames(core_attribute)%in%core_DEGs_v2
table(core_attribute[,'core']) # check
all(genes%in%core_DEGs_v2 == core_attribute[,'core'])

core_attribute

# Write core attribute table:
write.table(x = core_attribute,
            file = paste0(network_path,
                          '/core_attribute.tsv'),
            quote = F,sep = '\t',
            row.names = F)
