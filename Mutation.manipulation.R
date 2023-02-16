# Install biomart
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

# Connect to human genes biomart
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Build a biomart query
brca.data$mutational$Entrez_Gene_Id

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes <- getBM(
  filters="entrezgene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=brca.data$mutational$Entrez_Gene_Id,
  mart=mart)
genes

###################################

library(dplyr)

# Select for specific column headings
head(brca.data$mutational)
BRCAmutations <- dplyr::select(brca.data$mutational, Hugo_Symbol, Entrez_Gene_Id, Gene, Variant_Classification, SIFT, PolyPhen)
BRCAmutations
table(BRCAmutations$Gene)

# Calculate Frequency of Mutation



###################################

# Compare 2 lists of ensemble genes
genes$ensembl_gene_id
BRCAmutations$Gene
genes$ensembl_gene_id %in% BRCAmutations$Gene

NumberTRUE <- table(genes$ensembl_gene_id %in% BRCAmutations$Gene)["TRUE"]
NumberFALSE <- table(genes$ensembl_gene_id %in% BRCAmutations$Gene)["FALSE"]
(NumberFALSE/NumberTRUE)*100
