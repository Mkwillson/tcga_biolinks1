#### Downloading and preparing TARGET-AML data ####
projects <- getGDCprojects()$project_id
as.data.frame(grep("TARGET", projects,value = TRUE))

query_Target <- GDCquery(project = "TARGET-AML", 
                         data.category = "Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification",
                         workflow.type = "STAR - Counts")

GDCdownload(query_Target)

#### Summarized experiment = TRUE ####
query.tp.cases <- getResults(query_Target, cols="cases")
length(query.tp.cases)
query.tp.cases.dups <- query.tp.cases[duplicated(query.tp.cases)]
length(query.tp.cases.dups)
query.tp.cases.unique <- unique(query.tp.cases)
length(query.tp.cases.unique)
query.tp.cases.nodups <- setdiff(query.tp.cases.unique,query.tp.cases.dups)
length(query.tp.cases.nodups)


query.tp.nodups <- GDCquery(project = "TARGET-AML", 
                         data.category = "Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification",
                         workflow.type = "STAR - Counts",
                         barcode = query.tp.cases.nodups)
aml.data <- GDCprepare(query.tp.nodups)
aml.data
query.tp.nodups$results

#### Downloading and preparing TARGET-ALL-P2 data ####
# Expression data
projects <- getGDCprojects()$project_id
as.data.frame(grep("TARGET", projects,value = TRUE))

query_Target2 <- GDCquery(project = "TARGET-ALL-P2", 
                         data.category = "Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification",
                         workflow.type = "STAR - Counts")

GDCdownload(query_Target2)

### SummarizedExperiment = TRUE
allp2.exp.data <- GDCprepare(query_Target2)

# Mutational data
query_mutationalall2 <- GDCquery(project = "TARGET-ALL-P2", 
                          data.category = "Simple Nucleotide Variation", 
                          data.type = "Masked Somatic Mutation",
                          workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")


GDCdownload(query_mutationalall2)