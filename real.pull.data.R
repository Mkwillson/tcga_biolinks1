# Creating a singular function for data extraction of:
# Transcriptomic/expression data
# Clinical data
# Mutational data

real.pull.data <- function(
    project = "",
    workflow.type = "STAR - Counts",
    save.filename = ""
){
  
  # Downloading expression data
  query.exp <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    legacy = FALSE,
    data.type = "Gene Expression Quantification",
    workflow.type = workflow.type
  )
  
  GDCdownload(query.exp)
  expression <- GDCprepare(query = query.exp,
                           save = TRUE,
                           save.filename = save.filename)
  
  
  # Downloading clinical data
  query.clin <- GDCquery(
    project = project,
    file.type = "xml",
    data.category = "Clinical",
    legacy = FALSE,
    data.type = "Clinical Supplement"
  )
  
  GDCdownload(query.clin)
  clinical <- GDCprepare_clinic(query = query.clin, clinical.info = "patient")
  
  
  # Downloading mutational data
  query.maf <- GDCquery(
    project = project, 
    data.category = "Simple Nucleotide Variation", 
    access = "open", 
    legacy = FALSE, 
    data.type = "Masked Somatic Mutation",
    workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
  )
  
  GDCdownload(query.maf)
  maf <- GDCprepare(query.maf)
  
  
  return(list(expression = expression,
              clinical = clinical,
              maf = maf))
}

# Function
brca.data <- real.pull.data(project = "TCGA-BRCA",
                           save.filename = "BRCA.Exp.rda"
)
