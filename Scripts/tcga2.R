# Downloading expression data
query.exp.prad <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  legacy = FALSE,
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query.exp.prad)

expression.prad <- GDCprepare(query = query.exp.prad,
                         save = TRUE,
                         save.filename = "TCGA-PRAD Expression")

# Downloading clinical data
query.clin.prad <- GDCquery(
  project = "TCGA-PRAD",
  file.type = "xml",
  data.category = "Clinical",
  legacy = FALSE,
  data.type = "Clinical Supplement"
)

GDCdownload(query.clin.prad)
clinical.prad <- GDCprepare_clinic(query = query.clin.prad, clinical.info = "patient")

# Downloading mutational data following helper code
library(maftools)

# GDCquery_Maf function was removed in last update so method does not work!!!
mutational <- query.maf(tumor = tumor, pipelines = "mutect2")

#Downloading mutational data alternative

query.maf.prad <- GDCquery(
  project = "TCGA-PRAD", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query.maf.prad)

# Prepare into R object
maf.prad <- GDCprepare(query.maf.prad)

# View data in datatable
datatable(maf.prad[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# All 3 data types - ERROR: NO FUNCTION TO RETURN FROM!!!
return(list(expression = expression,
            clinical = clinical,
            mutational = mutational))

# Return all 3 data types
multi_return <- function() {
  my_list <- list("expression" = expression.prad, "clinical" = clinical.prad, "mutational" = maf.prad)
  return(my_list) 
}