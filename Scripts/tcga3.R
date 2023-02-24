# Downloading expression data
query.exp.aml <- GDCquery(
  project = "TARGET-AML",
  data.category = "Transcriptome Profiling",
  legacy = FALSE,
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query.exp.aml)

expression <- GDCprepare(query = query.exp.aml,
                         save = TRUE,
                         save.filename = "TARGET-AML Expression")

# Downloading clinical data
query.clin.aml <- GDCquery(
  project = "TARGET-AML",
  file.type = "xml",
  data.category = "Clinical",
  legacy = FALSE,
  data.type = "Clinical Supplement"
)

GDCdownload(query.clin.aml)
clinical.aml <- GDCprepare_clinic(query = query.clin.aml, clinical.info = "patient")

# Downloading mutational data following helper code
library(maftools)

# GDCquery_Maf function was removed in last update so method does not work!!!
mutational.aml <- query.maf(tumor = tumor, pipelines = "mutect2")

#Downloading mutational data alternative

query.maf.aml <- GDCquery(
  project = "TARGET-AML", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query.maf.aml)

# Prepare into R object
maf.aml <- GDCprepare(query.maf.aml)

# View data in datatable
datatable(maf.aml[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

# All 3 data types - ERROR: NO FUNCTION TO RETURN FROM!!!
return(list(expression = expression,
            clinical = clinical,
            mutational = mutational))

# Return all 3 data types
multi_return <- function() {
  my_list <- list("expression" = expression.aml, "clinical" = clinical.aml, "mutational" = maf.aml)
  return(my_list) 
}
