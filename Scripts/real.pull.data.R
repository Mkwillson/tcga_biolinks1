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
  mutational <- GDCprepare(query.maf)
  
  
  return(list(expression = expression,
              clinical = clinical,
              mutational = mutational))
}

# Extracting TCGA-BRCA data
brca.data <- real.pull.data(project = "TCGA-BRCA",
                           save.filename = "BRCA.Exp.rda"
)

# Print
brca.data$expression
brca.data$clinical
brca.data$mutational

# Extracting TCGA-PRAD data
prad.data <- real.pull.data(project = "TCGA-PRAD",
                            save.filename = "PRAD.Exp.rda"
)

# Print
prad.data$expression
prad.data$clinical
prad.data$mutational

#Extracting TCGA-LUAD data
luad.data <- real.pull.data(project = "TCGA-LUAD",
                              save.filename = "LUAD.Exp.rda"
)

# Print
luad.data$expression
luad.data$clinical
luad.data$mutational

# Extracting TCGA-KIRC data
kirc.data <- real.pull.data(project = "TCGA-KIRC",
                              save.filename = "KIRC.Exp.rda"
)

# Print
kirc.data$expression
kirc.data$clinical
kirc.data$mutational

# Extracting TCGA-STAD data
stad.data <- real.pull.data(project = "TCGA-STAD",
                            save.filename = "STAD.Exp.rda"
)

# Print
stad.data$expression
stad.data$clinical
stad.data$mutational

# Extracting TCGA-KIRP data
kirp.data <- real.pull.data(project = "TCGA-KIRP",
                            save.filename = "KIRP.Exp.rda"
)

# Print
kirp.data$expression
kirp.data$clinical
kirp.data$mutational