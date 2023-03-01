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
###################################

# List frequencies of all mutations
brca.data$mutational %>% dplyr::count(Hugo_Symbol, "MTOR")
dplyr::count(brca.data$mutational$Variant_Classification) 
brca.data$mutational %>% dplyr::count(Hugo_Symbol, Variant_Classification)
brca.data$mutational %>% dplyr::count(case_id, Hugo_Symbol)
brca.data$mutational %>% dplyr::count(Variant_Classification != "Silent")
brca.data$mutational %>% dplyr::filter(Variant_Classification != "Silent")
dplyr::select(brca.data$mutational, Variant_Classification)

test <- as.data.frame(brca.data$mutational$Variant_Classification, brca.data$mutational$Hugo_Symbol)
test
##############
prad.data$mutational

##
one <- brca.data$mutational %>% dplyr::filter(Variant_Classification != "Silent")
one
help("group_by")
two <- dplyr::group_by(one, Hugo_Symbol)
two
dplyr::select(two, Hugo_Symbol, Entrez_Gene_Id, Gene, Variant_Classification, SIFT, PolyPhen)
dplyr::count(two, Hugo_Symbol, Entrez_Gene_Id, Gene, Variant_Classification, SIFT, PolyPhen)

# Filtering:
no.silent <- brca.data$mutational %>% dplyr::filter(Variant_Classification != "Silent")
no.silent %>% dplyr::count(Hugo_Symbol)



##
BRCAmutations
BRCAmutations <- dplyr::count(brca.data$mutational, Hugo_Symbol, Entrez_Gene_Id, Gene, Variant_Classification)
BRCAmutations

#################################

# Compare 2 lists of ensemble genes
genes$ensembl_gene_id
BRCAmutations$Gene
genes$ensembl_gene_id %in% BRCAmutations$Gene
head(data.frame(genes$ensembl_gene_id, BRCAmutations$Gene))

NumberTRUE <- table(genes$ensembl_gene_id %in% BRCAmutations$Gene)["TRUE"]
NumberFALSE <- table(genes$ensembl_gene_id %in% BRCAmutations$Gene)["FALSE"]
(NumberFALSE/NumberTRUE)*100
