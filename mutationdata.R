library(dplyr)

#### PRAD ######################################################################
head(prad.data$mutational)

colnames(prad.data$mutational)

# Count number of mutations per gene
prad.data$mutational %>%
  count(Gene) -> Total.Number.Mutations
Total.Number.Mutations

# Sort Variant_Classification into main categories
# OPTION 1 - BASED ON ?
table(prad.data$mutational$Variant_Classification)

#1
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation" | Variant_Classification=="Nonsense_Mutation" | Variant_Classification=="Frame_Shift_Del" |  Variant_Classification=="Frame_Shift_Ins" | Variant_Classification=="Nonstop_Mutation" ) %>%
  count(Gene) -> High.Importance
All.Mutations

#2
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="3'UTR" | Variant_Classification=="5'UTR" | Variant_Classification=="3'Flank" |  Variant_Classification=="5'Flank" | Variant_Classification=="In_Frame_Del" | Variant_Classification=="In_Frame_Ins" | Variant_Classification=="Intron" | Variant_Classification=="Splice_Region" | Variant_Classification=="Splice_Site" | Variant_Classification=="Translation_Start_Site") %>%
  count(Gene) -> Medium.Importance
Medium.Importance

#3
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Silent" | Variant_Classification=="IGR" | Variant_Classification=="RNA") %>%
  count(Gene) -> Low.Importance
Low.Importance

##################################################################################################################

# OPTION 2 - BASED ON IMPACT RATING
prad.data$mutational %>% dplyr::count(prad.data$mutational$Variant_Classification, prad.data$mutational$IMPACT)

# High impact
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="HIGH") %>% 
  count(Gene) -> High.Impact.Mutations
High.Impact.Mutations

# Moderate impact
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODERATE") %>% 
  count(Gene) -> Moderate.Impact.Mutations
Moderate.Impact.Mutations

# Low impact
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="LOW") %>% 
  count(Gene) -> Low.Impact.Mutations
Low.Impact.Mutations

# Modifier impact
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODIFIER") %>% 
  count(Gene) -> Modifier.Mutations
Modifier.Mutations

# Filter for missense mutations, find number of missense mutations per gene
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
count(Gene) -> no.missense
no.missense

# Filter for deleterious missense mutations, find number per gene
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  filter(grepl("deleterious", SIFT))%>%
  count(Gene) -> no.sift.del
no.sift.del

# Filter for possibly damaging missense mutations, find number per gene
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  filter(grepl("possibly", PolyPhen))%>%
  count(Gene) -> no.possible.polyphen
no.possible.polyphen

# Filter for probably damaging missense mutations, find number per gene
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  filter(grepl("probably", PolyPhen))%>%
  count(Gene) -> no.probable.polyphen
no.probable.polyphen

# Separate SIFT data into category and numerical
# SIFT category
library(qdapRegex)
rm_round(text.var =  prad.data$mutational$SIFT, extract = F) -> prad.data$mutational$SIFT.cat

prad.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(SIFT.cat) %>%
  select(Gene,SIFT.cat) -> sift.categories

# Count by category type
table(sift.categories) -> sift.categories.tab
as.data.frame.matrix(sift.categories.tab) -> sift.categories.df

colnames(sift.categories.df)[1]  <- "SIFT.Deleterious"
colnames(sift.categories.df)[2]  <- "SIFT.Deleterious_low_confidence"
colnames(sift.categories.df)[3]  <- "SIFT.Tolerated"
colnames(sift.categories.df)[4]  <- "SIFT.Tolerated_low_confidence"
sift.categories.df

# SIFT numerical values
unlist(rm_round(text.var =  prad.data$mutational$SIFT, extract = T)) -> prad.data$mutational$SIFT.num

prad.data$mutational$SIFT.num <- as.numeric(prad.data$mutational$SIFT.num)

prad.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(SIFT.num) %>%
  select(Gene,SIFT.num) -> sift.number

sift.number

# Summarise sift.number - ignore NAs
# Mean
sift.mean <- sift.number %>%
  group_by(Gene) %>%
  summarise(mean = mean(SIFT.num, na.rm = TRUE))

sift.mean

# Median
sift.median <- sift.number %>%
  group_by(Gene) %>%
  summarise(median = median(SIFT.num, na.rm = TRUE))

sift.median

# Maximum
sift.maximum <- sift.number %>%
  group_by(Gene) %>%
  summarise(max = max(SIFT.num, na.rm = TRUE))

sift.maximum

# Standard deviation
sift.sd <- sift.number %>%
  group_by(Gene) %>%
  summarise(sd = sd(SIFT.num, na.rm = TRUE))
sift.sd

# Separate PolyPhen data into category and numerical
# polyphen category
library(qdapRegex)
rm_round(text.var =  prad.data$mutational$PolyPhen, extract = F) -> prad.data$mutational$PolyPhen.cat

prad.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(PolyPhen.cat) %>%
  select(Gene,PolyPhen.cat) -> PolyPhen.categories

PolyPhen.categories

# Count categories
table(PolyPhen.categories) -> PolyPhen.categories.tab
as.data.frame.matrix(PolyPhen.categories.tab) -> PolyPhen.categories.df
colnames(PolyPhen.categories.df)[1]  <- "PolyPhen.Benign"
colnames(PolyPhen.categories.df)[2]  <- "PolyPhen.Possibly.Damaging"
colnames(PolyPhen.categories.df)[3]  <- "PolyPhen.Probably.Damaging"
colnames(PolyPhen.categories.df)[4]  <- "PolyPhen.Unknown"
PolyPhen.categories.df

# PolyPhen numerical values
unlist(rm_round(text.var =  prad.data$mutational$PolyPhen, extract = T)) -> prad.data$mutational$PolyPhen.num

prad.data$mutational$PolyPhen.num <- as.numeric(prad.data$mutational$PolyPhen.num)

prad.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(PolyPhen.num) %>%
  select(Gene,PolyPhen.num) -> PolyPhen.number

PolyPhen.number

# Summarise PolyPhen.number - ignore NAs
# Mean
PolyPhen.mean <- PolyPhen.number %>%
  group_by(Gene) %>%
  summarise(mean = mean(PolyPhen.num, na.rm = TRUE))

PolyPhen.mean

# Median
PolyPhen.median <- PolyPhen.number %>%
  group_by(Gene) %>%
  summarise(median = median(PolyPhen.num, na.rm = TRUE))

PolyPhen.median

# Maximum
PolyPhen.maximum <- PolyPhen.number %>%
  group_by(Gene) %>%
  summarise(max = max(PolyPhen.num, na.rm = TRUE))

PolyPhen.maximum

# Standard deviation
PolyPhen.sd <- PolyPhen.number %>%
  group_by(Gene) %>%
  summarise(sd = sd(PolyPhen.num, na.rm = TRUE))
PolyPhen.sd

# Case_id
prad.data$mutational %>% count(Gene,case_id) -> test
table(test$Gene) -> test2
as.data.frame(test2) -> Patients.With.Mutation
colnames(Patients.With.Mutation)[1] <- "Gene"
colnames(Patients.With.Mutation)[2] <- "Patients.With.Mutation"
Patients.With.Mutation

# Annotate with GO terms using Biomart
library(biomaRt)

annotate.HTseq.IDs <- function(HTseq.IDs){
  ENSEMBL_DB_HOST = "uswest.ensembl.org" # Set back to default, once they are up and running again
  ENSEMBL_VERSION = "Ensembl Genes 105"  # Try to fix https://support.bioconductor.org/p/104454/
  
  # mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
  # host = ENSEMBL_DB_HOST,
  #                version = ENSEMBL_VERSION)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  ensemblID <- gsub("\\..*", '', HTseq.IDs)
  #fetch the annotation information from the Ensembl database
  #ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
  symbols <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'go_id'), filters='ensembl_gene_id', ensemblID, mart=mart)
  #combine the annotation with the RNA-Seq data
  annotatedix <- match(ensemblID,symbols$ensembl_gene_id)
  symbols[annotatedix,] -> annotatedGenes
  return(cbind(HTseq.IDs,annotatedGenes))
}

retrieved.data <- annotate.HTseq.IDs(PRAD.final2$Gene)
retrieved.data
dim(retrieved.data)
colnames(retrieved.data)[2] <- "Gene"
retrieved.data %>% dplyr::select(Gene, go_id) -> GO
Gene.Ontology <- GO                                  # Duplicate data frame
Gene.Ontology[Gene.Ontology == ""] <- NA                     # Replace blank by NA
Gene.Ontology

class(retrieved.data)
class(Gene.Ontology)
class(Number.Patients.With.Mutation)
  
# Construct function for matching different columns of data to same genes and length
match.select <- function(x, y, label){
  left_join(x, y, by = "Gene") -> temp.df
  temp.df[,c(1,3)] -> temp.df
  colnames(temp.df)[2] <- label
  return(temp.df)
}

# Match all columns
Gene.Ontology.Match <- match.select(x = Total.Number.Mutations, y = Gene.Ontology, label = "Gene.Ontology")
Number.Patients.With.Mutation <- match.select(x = Total.Number.Mutations, y = Patients.With.Mutation, label = "Number.Patients.With.Mutation")
Low.Impact.Mutations <- match.select(x = Total.Number.Mutations, y = Low.Impact.Mutations, label = "Low.Impact.Mutations")
Moderate.Impact.Mutations <- match.select(x = Total.Number.Mutations, y = Moderate.Impact.Mutations, label = "Moderate.Impact.Mutations")
High.Impact.Mutations <- match.select(x = Total.Number.Mutations, y = High.Impact.Mutations, label = "High.Impact.Mutations")
Modifier.Mutations <- match.select(x = Total.Number.Mutations, y = Modifier.Mutations, label = "Modifier.Mutations")
Deleterious.Missense <- match.select(x = Total.Number.Mutations, y = no.sift.del, label = "no.sift.del")
Missense.Mutation <- match.select(x = Total.Number.Mutations, y = no.missense, label = "no.missense")
Possibly.Damaging.Missense <- match.select(x = Total.Number.Mutations, y = no.possible.polyphen, label = "no.possible.polyphen")
Probably.Damaging.Missense <- match.select(x = Total.Number.Mutations, y = no.probable.polyphen, label = "no.probable.polyphen")
SIFT.Mean <- match.select(x = Total.Number.Mutations, y = sift.mean, label = "sift.mean")
SIFT.Median <- match.select(x = Total.Number.Mutations, y = sift.median, label = "sift.median")
SIFT.Maximum <- match.select(x = Total.Number.Mutations, y = sift.maximum, label = "sift.maximum")
SIFT.StandardDeviation <- match.select(x = Total.Number.Mutations, y = sift.sd, label = "sift.sd")
PolyPhen.Mean <- match.select(x = Total.Number.Mutations, y = PolyPhen.mean, label = "PolyPhen.mean")
PolyPhen.Median <- match.select(x = Total.Number.Mutations, y = PolyPhen.median, label = "PolyPhen.median")
PolyPhen.Maximum <- match.select(x = Total.Number.Mutations, y = PolyPhen.maximum, label = "PolyPhen.maximum")
PolyPhen.StandardDeviation <- match.select(x = Total.Number.Mutations, y = PolyPhen.sd, label = "PolyPhen.sd")

# Create dataframe
data.frame(Total.Number.Mutations,
           Gene.Ontology[,2],
           Number.Patients.With.Mutation[,2],
           Low.Impact.Mutations[,2],
           Moderate.Impact.Mutations[,2],
           High.Impact.Mutations[,2],
           Modifier.Mutations[,2],
           Missense.Mutation[,2],
           Deleterious.Missense[,2],
           Possibly.Damaging.Missense[,2],
           Probably.Damaging.Missense[,2],
           SIFT.Mean[,2],
           SIFT.Median[,2],
           SIFT.Maximum [,2],
           SIFT.StandardDeviation[,2],
           PolyPhen.Mean[,2],
           PolyPhen.Median[,2],
           PolyPhen.Maximum[,2],
           PolyPhen.StandardDeviation[,2]     
) -> df
View(df)
# Join together df, with all other data frames 
# Add another "Gene" column to allow matching
data.frame(Gene = rownames(sift.categories.df),sift.categories.df) -> sift.categories.df.first.col
data.frame(Gene = rownames(PolyPhen.categories.df),PolyPhen.categories.df) -> PolyPhen.categories.df.first.col


# Full outer join of multiple data frames
full_join(df, sift.categories.df.first.col, by = "Gene") %>%
  full_join(., PolyPhen.categories.df.first.col, by = "Gene") -> PRAD.final
View(PRAD.final)

colnames(PRAD.final)[2] <- "Frequency.Gene.Is.Mutated"

# Percentages
install.packages("formattable")
library(formattable)
PRAD.final$Percent.Low.Impact.Mutations = percent(PRAD.final$Low.Impact.Mutations/PRAD.final$Frequency.Gene.Is.Mutated) -> Percent.Low.Impact.Mutations
PRAD.final$Percent.Moderate.Impact.Mutations = percent(PRAD.final$Moderate.Impact.Mutations/PRAD.final$Frequency.Gene.Is.Mutated) -> Percent.Moderate.Impact.Mutations
PRAD.final$Percent.High.Impact.Mutations = percent(PRAD.final$High.Impact.Mutations/PRAD.final$Frequency.Gene.Is.Mutated) -> Percent.High.Impact.Mutations
PRAD.final$Percent.Modifier.Mutations = percent(PRAD.final$Modifier.Mutations/PRAD.final$Frequency.Gene.Is.Mutated) -> Percent.Modifier.Mutations
PRAD.final$Percent.Missense.Mutation = percent(PRAD.final$Missense.Mutation/PRAD.final$Frequency.Gene.Is.Mutated) -> Percent.Missense.Mutation

# Percent of patients with that mutation
# Find total number of patients!
prad.data$mutational %>% count(Gene,case_id) -> case.id.gene
data.frame(case.id.gene)
df_uniq <- unique(case.id.gene$case_id)
length(df_uniq)

# Calculate %
PRAD.final$Percent.Patients.With.Mutation = percent(PRAD.final$Number.Patients.With.Mutation/(length(df_uniq)))


# Reorder columns
PRAD.final2 = PRAD.final %>% select(Gene, 
                                   Gene.Ontology...2.,
                                   Frequency.Gene.Is.Mutated, 
                                   Number.Patients.With.Mutation...2., 
                                   Percent.Patients.With.Mutation, 
                                   Modifier.Mutations...2.,
                                   Percent.Modifier.Mutations,
                                   Low.Impact.Mutations...2.,
                                   Percent.Low.Impact.Mutations,
                                   Moderate.Impact.Mutations...2., 
                                   Percent.Moderate.Impact.Mutations,
                                   High.Impact.Mutations...2.,
                                   Percent.High.Impact.Mutations,
                                   Missense.Mutation...2., 
                                   Percent.Missense.Mutation,
                                   SIFT.Deleterious, 
                                   SIFT.Deleterious_low_confidence, 
                                   SIFT.Tolerated, 
                                   SIFT.Tolerated_low_confidence, 
                                   SIFT.Mean...2.,
                                   SIFT.Median...2., 
                                   SIFT.Maximum...2., 
                                   SIFT.StandardDeviation...2., 
                                   PolyPhen.Benign, 
                                   PolyPhen.Possibly.Damaging,
                                   PolyPhen.Probably.Damaging, 
                                   PolyPhen.Unknown,
                                   PolyPhen.Mean...2.,
                                   PolyPhen.Median...2.,
                                   PolyPhen.Maximum...2.,
                                   PolyPhen.StandardDeviation...2.)


View(PRAD.final2)


PRAD.final3 = PRAD.final %>% select(Gene, 
                                    Frequency.Gene.Is.Mutated, 
                                    Number.Patients.With.Mutation...2., 
                                    Percent.Patients.With.Mutation, 
                                    Modifier.Mutations...2.,
                                    Percent.Modifier.Mutations,
                                    Low.Impact.Mutations...2.,
                                    Percent.Low.Impact.Mutations,
                                    Moderate.Impact.Mutations...2., 
                                    Percent.Moderate.Impact.Mutations,
                                    High.Impact.Mutations...2.,
                                    Percent.High.Impact.Mutations,
                                    Missense.Mutation...2., 
                                    Percent.Missense.Mutation,
                                    SIFT.Deleterious, 
                                    SIFT.Deleterious_low_confidence, 
                                    SIFT.Tolerated, 
                                    SIFT.Tolerated_low_confidence, 
                                    SIFT.Mean...2.,
                                    SIFT.Median...2., 
                                    SIFT.Maximum...2., 
                                    SIFT.StandardDeviation...2., 
                                    PolyPhen.Benign, 
                                    PolyPhen.Possibly.Damaging,
                                    PolyPhen.Probably.Damaging, 
                                    PolyPhen.Unknown,
                                    PolyPhen.Mean...2.,
                                    PolyPhen.Median...2.,
                                    PolyPhen.Maximum...2.,
                                    PolyPhen.StandardDeviation...2.)

write.csv(PRAD.final2, "~/tcga_biolinks1/PRAD.final2.csv", row.names=TRUE)
########################################################################################################################
# Annotate with GO terms using Biomart
library(biomaRt)

annotate.HTseq.IDs <- function(HTseq.IDs){
  ENSEMBL_DB_HOST = "uswest.ensembl.org" # Set back to default, once they are up and running again
  ENSEMBL_VERSION = "Ensembl Genes 105"  # Try to fix https://support.bioconductor.org/p/104454/
  
  # mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
  # host = ENSEMBL_DB_HOST,
  #                version = ENSEMBL_VERSION)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  ensemblID <- gsub("\\..*", '', HTseq.IDs)
  #fetch the annotation information from the Ensembl database
  #ensembl <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
  symbols <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'go_id'), filters='ensembl_gene_id', ensemblID, mart=mart)
  #combine the annotation with the RNA-Seq data
  annotatedix <- match(ensemblID,symbols$ensembl_gene_id)
  symbols[annotatedix,] -> annotatedGenes
  return(cbind(HTseq.IDs,annotatedGenes))
}

retrieved.data <- annotate.HTseq.IDs(PRAD.final2$Gene)
retrieved.data
dim(retrieved.data)
colnames(retrieved.data)[2] <- "Gene"
retrieved.data %>% dplyr::select(Gene, go_id) -> GO
Gene.Ontology <- GO                                  # Duplicate data frame
Gene.Ontology[Gene.Ontology == ""] <- NA                     # Replace blank by NA
Gene.Ontology


listAttributes(ensembl)
listAttributes(mart)

