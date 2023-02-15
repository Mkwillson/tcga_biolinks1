# Install biomart
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

# Connect to human genes biomart
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Build a biomart query

listEnsemblGenomes()

filters = listFilters(ensembl)
filters[1:5,]

attributes = listAttributes(ensembl)
attributes[1:5,]

searchFilters(mart = ensembl, pattern = "ensembl.*id")

mart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
head(keytypes(mart), n = 3)

brca.data$mutational$Entrez_Gene_Id

ensembl = getBM(attributes = c('entrezgene_id', 'ensembl_gene_id'), 
              filters = 'entrezgene_id', 
              values = brca.data$mutational$Entrez_Gene_Id, 
              mart = ensembl)

###################################
library(dplyr)

# Select for specific column headings
head(brca.data$mutational)
BRCAmutations <- dplyr::select(brca.data$mutational, Hugo_Symbol, Entrez_Gene_Id, Gene, Variant_Classification, SIFT, PolyPhen)
BRCAmutations
brca.data$mutational$Consequence
