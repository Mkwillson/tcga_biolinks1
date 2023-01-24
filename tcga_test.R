TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

# Creating a singular function for data extraction of:
# Transcriptomic/expression data
# Mutational data
# Clinical data

# Downloading expression data
query.exp <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  legacy = FALSE,
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query.exp)

expression <- GDCprepare(query = query.exp,
                         save = TRUE,
                         save.filename = "TCGA-BRCA Expression")

# Downloading clinical data
query.clin <- GDCquery(
  project = "TCGA-BRCA",
  file.type = "xml",
  data.category = "Clinical",
  legacy = FALSE,
  data.type = "Clinical Supplement"
)

GDCdownload(query.clin)
clinical <- GDCprepare_clinic(query = query.clin, clinical.info = "patient")

# Downloading mutational data following helper code
library(maftools)

# GDCquery_Maf function was removed in last update so method does not work!!!
mutational <- query.maf(tumor = tumor, pipelines = "mutect2")

#Downloading mutational data alternative

query.maf <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )

GDCdownload(query.maf)

# Prepare into R object
maf <- GDCprepare(query.maf)

# View data in datatable
datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# All 3 data types - ERROR: NO FUNCTION TO RETURN FROM!!!
return(list(expression = expression,
            clinical = clinical,
            mutational = mutational))

# Return all 3 data types
multi_return <- function() {
  my_list <- list("expression" = expression, "clinical" = clinical, "mutational" = maf)
  return(my_list) 
}

# real.pull.data isn't a function so what is this supposed to be doing?
brca.data <- real.pull.data(project = "TCGA-BRCA",
                           tumor = "BRCA",
                           save.filename = "BRCA.Exp.rda"
)

