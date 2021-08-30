################################################################################
## gene_analysis.r
## Francisco Martínez Picó - francisco9896@gmail.com
################################################################################

Sys.info()[c('nodename', 'user')]
rm(list = ls())
R.version.string # 'R version 4.0.3 (2020-10-10)'

promiscuity_path = '/Users/francisco/Desktop/TFM/drugs/promiscuity'
network_path = '/Users/francisco/Desktop/TFM/network'

# LOAD PACKAGES -----------------------------------------------------------
pacman::p_load(stringr)

# LOAD DATA ---------------------------------------------------------------
# Drug Bank info:
drug_bank = read.delim(file = paste0(promiscuity_path,
                                     '/drugbank_info_2020_07.txt'),
           header = T)

db = drug_bank[, c('X.ID', 'name', 'targets', 'groups')]

# Genes from central network:
network_genes = read.csv(file = paste0(network_path, '/central_network.csv'))

# Selected drugs:
# Care since copper appears twice us a drug in our data, we only write it here
# once:
my_drugs_id = c('DB01136', 'DB00898', 'DB11338', 'DB00945', 'DB01593',
                'DB14487', 'DB14533', 'DB14548', 'DB12010', 'DB08889',
                'DB00188', 'DB00759', 'DB09130', 'DB00615', 'DB00716',
                'DB03904', 'DB00123', 'DB00129', 'DB00895')

my_drugs_name = c('Carvedilol', 'Ethanol', 'Clove oil', 'Acetylsalicylic acid',
                  'Zinc', 'Zinc acetate', 'Zinc chloride',
                  'Zinc sulfate, unspecified form', 'Fostamatinib',
                  'Carfilzomib', 'Bortezomib', 'Tetracycline', 'Copper',
                  'Rifabutin', 'Nedocromil', 'Urea', 'L-Lysine',
                  'Ornithine', 'Benzylpenicilloyl polylysine')

my_drugs = data.frame( id = my_drugs_id, name = my_drugs_name)

# So we have DrugBank info in 'drug_bank', our central network genes in
# 'network_genes' and our drugs info in 'my_drugs'. We aim to know how many
# genes of our network are targeted by our drugs.

# PREPARE DATA.FRAME ------------------------------------------------------
df_final = data.frame()

for (e in my_drugs[, 'id']) {
     df <- data.frame(db[(db[, 'X.ID'] == e),])
     df_final <- rbind(df_final,df)
}

# In 'df_final' we have our drug info collected from 'db' = 'drug_bank' =
# DrugBank database.

# For some reason, 'Copper' has two elements in its 'targets' column (same
# element duplicated). We fix it manually:

# Then we create data.frame:
targets = data.frame()

for (e in df_final[, 'name']){
    # Get gene/target string
    genes_string = df_final[df_final[,'name'] == e, 'targets']

    # Separate in ';'
    genes = str_split(genes_string, pattern = ';')

    genes_final = data.frame(genes)
    colnames(genes_final) = c('targets')
    genes_final$drug = e

    targets = rbind(targets, genes_final)
}

# In 'targets' we have a row for each target. Each target is associated with the
# drug that targets it. Next step is to check if these targets appear in our
# network (beside our selected ones).
my_genes = network_genes[,'name']

positions = targets[,'targets']%in%my_genes

promiscous_drugs = targets[positions,]

# Finally, we get 'promiscous_drugs', which shows genes that are targeted in our
# network and by which drug.
View(promiscous_drugs)

for (e in my_drugs[,'name']) {
    print(e)
    print(table(promiscous_drugs[, 'drug'] == e))
}

# Note that each drug must appear at least 1 time. We consider promiscous those
# drugs that appear more than once.
table(promiscous_drugs[, 'drug'] == 'Carvedilol')

# [1] "Carvedilol"
# FALSE  TRUE
#    63     1

# [1] "Ethanol"
# FALSE  TRUE
#    63     1

# [1] "Clove oil"
# FALSE  TRUE
#    63     1

# [1] "Acetylsalicylic acid"
# FALSE  TRUE
#    62     2

# [1] "Zinc"
# FALSE  TRUE
#    58     6

# [1] "Zinc acetate"
# FALSE  TRUE
#    58     6

# [1] "Zinc chloride"
# FALSE  TRUE
#    59     5

# [1] "Zinc sulfate, unspecified form"
# FALSE  TRUE
#    59     5

# [1] "Fostamatinib"
# FALSE  TRUE
#    47    17

# [1] "Carfilzomib"
# FALSE  TRUE
#    62     2

# [1] "Bortezomib"
# FALSE  TRUE
#    63     1

# [1] "Tetracycline"
# FALSE  TRUE
#    63     1

# [1] "Copper"
# FALSE  TRUE
#    54    10

# [1] "Rifabutin"
# FALSE  TRUE
#    63     1

# [1] "Nedocromil"
# FALSE  TRUE
#    63     1

# [1] "Urea"
# FALSE  TRUE
#    63     1

# [1] "L-Lysine"
# FALSE  TRUE
#    63     1

# [1] "Ornithine"
# FALSE  TRUE
#    63     1

# [1] "Benzylpenicilloyl polylysine"
# FALSE  TRUE
#    63     1
