library(dplyr)

#### PRAD ######################################################################
head(prad.data$mutational)
colnames(prad.data$mutational)

# Count number of mutations per gene
prad.data$mutational %>%
  count(Gene) -> Total.Number.Mutations.prad
Total.Number.Mutations.prad

# Classify mutations into 3 groups based on IMPACT score
prad.data$mutational %>% dplyr::count(prad.data$mutational$Variant_Classification, prad.data$mutational$IMPACT)

# High impact
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="HIGH") %>% 
  count(Gene) -> High.Impact.Mutations.prad
High.Impact.Mutations.prad

# Moderate impact
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODERATE") %>% 
  count(Gene) -> Moderate.Impact.Mutations.prad
Moderate.Impact.Mutations.prad

# Low impact
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="LOW") %>% 
  count(Gene) -> Low.Impact.Mutations.prad
Low.Impact.Mutations.prad

# Modifier impact
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODIFIER") %>% 
  count(Gene) -> Modifier.Mutations.prad
Modifier.Mutations.prad

# Filter for missense mutations, find number of missense mutations per gene
prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  count(Gene) -> no.missense.prad
no.missense.prad

# Separate SIFT data into category and numerical
# SIFT category
library(qdapRegex)
rm_round(text.var =  prad.data$mutational$SIFT, extract = F) -> prad.data$mutational$SIFT.cat

prad.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(SIFT.cat) %>%
  select(Gene,SIFT.cat) -> sift.categories.prad

# Count by category type
table(sift.categories.prad) -> sift.categories.tab.prad
as.data.frame.matrix(sift.categories.tab.prad) -> sift.categories.df.prad

# Rename column headings
colnames(sift.categories.df.prad)[1]  <- "SIFT.Deleterious.prad"
colnames(sift.categories.df.prad)[2]  <- "SIFT.Deleterious_low_confidence.prad"
colnames(sift.categories.df.prad)[3]  <- "SIFT.Tolerated.prad"
colnames(sift.categories.df.prad)[4]  <- "SIFT.Tolerated_low_confidence.prad"
sift.categories.df.prad

# SIFT numerical values
unlist(rm_round(text.var =  prad.data$mutational$SIFT, extract = T)) -> prad.data$mutational$SIFT.num

prad.data$mutational$SIFT.num <- as.numeric(prad.data$mutational$SIFT.num)

prad.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(SIFT.num) %>%
  select(Gene,SIFT.num) -> sift.number.prad

sift.number.prad

# Summarise sift.number - ignore NAs
# Mean
sift.mean.prad <- sift.number.prad %>%
  group_by(Gene) %>%
  summarise(mean = mean(SIFT.num, na.rm = TRUE))

sift.mean.prad

# Median
sift.median.prad <- sift.number.prad %>%
  group_by(Gene) %>%
  summarise(median = median(SIFT.num, na.rm = TRUE))

sift.median.prad

# Maximum
sift.maximum.prad <- sift.number.prad %>%
  group_by(Gene) %>%
  summarise(max = max(SIFT.num, na.rm = TRUE))

sift.maximum.prad

# Standard deviation
sift.sd.prad <- sift.number.prad %>%
  group_by(Gene) %>%
  summarise(sd = sd(SIFT.num, na.rm = TRUE))
sift.sd.prad

# Separate PolyPhen data into category and numerical
# polyphen category
library(qdapRegex)
rm_round(text.var =  prad.data$mutational$PolyPhen, extract = F) -> prad.data$mutational$PolyPhen.cat

prad.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(PolyPhen.cat) %>%
  select(Gene,PolyPhen.cat) -> PolyPhen.categories.prad

PolyPhen.categories.prad

# Count categories
table(PolyPhen.categories.prad) -> PolyPhen.categories.tab.prad
as.data.frame.matrix(PolyPhen.categories.tab.prad) -> PolyPhen.categories.df.prad

# Rename column headings
colnames(PolyPhen.categories.df.prad)[1]  <- "PolyPhen.Benign.prad"
colnames(PolyPhen.categories.df.prad)[2]  <- "PolyPhen.Possibly.Damaging.prad"
colnames(PolyPhen.categories.df.prad)[3]  <- "PolyPhen.Probably.Damaging.prad"
colnames(PolyPhen.categories.df.prad)[4]  <- "PolyPhen.Unknown.prad"
PolyPhen.categories.df.prad

# PolyPhen numerical values
unlist(rm_round(text.var =  prad.data$mutational$PolyPhen, extract = T)) -> prad.data$mutational$PolyPhen.num

prad.data$mutational$PolyPhen.num <- as.numeric(prad.data$mutational$PolyPhen.num)

prad.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(PolyPhen.num) %>%
  select(Gene,PolyPhen.num) -> PolyPhen.number.prad

PolyPhen.number.prad

# Summarise PolyPhen.number - ignore NAs
# Mean
PolyPhen.mean.prad <- PolyPhen.number.prad %>%
  group_by(Gene) %>%
  summarise(mean = mean(PolyPhen.num, na.rm = TRUE))

PolyPhen.mean.prad

# Median
PolyPhen.median.prad <- PolyPhen.number.prad %>%
  group_by(Gene) %>%
  summarise(median = median(PolyPhen.num, na.rm = TRUE))

PolyPhen.median.prad

# Maximum
PolyPhen.maximum.prad <- PolyPhen.number.prad %>%
  group_by(Gene) %>%
  summarise(max = max(PolyPhen.num, na.rm = TRUE))

PolyPhen.maximum.prad

# Standard deviation
PolyPhen.sd.prad <- PolyPhen.number.prad %>%
  group_by(Gene) %>%
  summarise(sd = sd(PolyPhen.num, na.rm = TRUE))

PolyPhen.sd.prad

# Case_id - find number of patients with that gene mutated
prad.data$mutational %>% count(Gene,case_id) -> Gene.Patient.prad
table(Gene.Patient.prad$Gene) -> Gene.Patient.tab.prad
as.data.frame(Gene.Patient.tab.prad) -> Patients.With.Mutation.prad
colnames(Patients.With.Mutation.prad)[1] <- "Gene"
colnames(Patients.With.Mutation.prad)[2] <- "Patients.With.Mutation"
Patients.With.Mutation.prad

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

# Retrieve GO terms for mutatated genes
retrieved.data.prad <- annotate.HTseq.IDs(Total.Number.Mutations.prad$Gene)
retrieved.data.prad

# Rename column
colnames(retrieved.data.prad)[2] <- "Gene"

# Select columns needed
retrieved.data.prad %>% dplyr::select(Gene, go_id) -> GO.prad

# Duplicate data frame
Gene.Ontology.prad <- GO.prad  

# Replace blanks with NA
Gene.Ontology.prad[Gene.Ontology.prad == ""] <- NA   

# Print 
Gene.Ontology.prad

# Select columns needed
retrieved.data.prad %>% dplyr::select(Gene, hgnc_symbol) -> hgnc.prad

# Duplicate data frame
Hugo.Symbol.prad <- hgnc.prad

# Print
Hugo.Symbol.prad


# Construct function for matching different columns of data to same genes and length
match.select <- function(x, y, label){
  left_join(x, y, by = "Gene") -> temp.df
  temp.df[,c(1,3)] -> temp.df
  colnames(temp.df)[2] <- label
  return(temp.df)
}

# Match columns
Number.Patients.With.Mutation.prad <- match.select(x = Total.Number.Mutations.prad, y = Patients.With.Mutation.prad, label = "Number.Patients.With.Mutation")
Low.Impact.Mutations.prad <- match.select(x = Total.Number.Mutations.prad, y = Low.Impact.Mutations.prad, label = "Low.Impact.Mutations")
Moderate.Impact.Mutations.prad <- match.select(x = Total.Number.Mutations.prad, y = Moderate.Impact.Mutations.prad, label = "Moderate.Impact.Mutations")
High.Impact.Mutations.prad <- match.select(x = Total.Number.Mutations.prad, y = High.Impact.Mutations.prad, label = "High.Impact.Mutations")
Modifier.Mutations.prad <- match.select(x = Total.Number.Mutations.prad, y = Modifier.Mutations.prad, label = "Modifier.Mutations")
Missense.Mutation.prad <- match.select(x = Total.Number.Mutations.prad, y = no.missense.prad, label = "no.missense")
SIFT.Mean.prad <- match.select(x = Total.Number.Mutations.prad, y = sift.mean.prad, label = "sift.mean")
SIFT.Median.prad <- match.select(x = Total.Number.Mutations.prad, y = sift.median.prad, label = "sift.median")
SIFT.Maximum.prad <- match.select(x = Total.Number.Mutations.prad, y = sift.maximum.prad, label = "sift.maximum")
SIFT.StandardDeviation.prad <- match.select(x = Total.Number.Mutations.prad, y = sift.sd.prad, label = "sift.sd")
PolyPhen.Mean.prad <- match.select(x = Total.Number.Mutations.prad, y = PolyPhen.mean.prad, label = "PolyPhen.mean")
PolyPhen.Median.prad <- match.select(x = Total.Number.Mutations.prad, y = PolyPhen.median.prad, label = "PolyPhen.median")
PolyPhen.Maximum.prad <- match.select(x = Total.Number.Mutations.prad, y = PolyPhen.maximum.prad, label = "PolyPhen.maximum")
PolyPhen.StandardDeviation.prad <- match.select(x = Total.Number.Mutations.prad, y = PolyPhen.sd.prad, label = "PolyPhen.sd")

# Create dataframe
data.frame(Total.Number.Mutations.prad,
           Hugo.Symbol.prad[,2],
           Gene.Ontology.prad[,2],
           Number.Patients.With.Mutation.prad[,2],
           Low.Impact.Mutations.prad[,2],
           Moderate.Impact.Mutations.prad[,2],
           High.Impact.Mutations.prad[,2],
           Modifier.Mutations.prad[,2],
           Missense.Mutation.prad[,2],
           SIFT.Mean.prad[,2],
           SIFT.Median.prad[,2],
           SIFT.Maximum.prad [,2],
           SIFT.StandardDeviation.prad[,2],
           PolyPhen.Mean.prad[,2],
           PolyPhen.Median.prad[,2],
           PolyPhen.Maximum.prad[,2],
           PolyPhen.StandardDeviation.prad[,2]     
) -> df.prad

View(df.prad)
# Join together df, with all other data frames 
# Add another "Gene" column to allow matching
data.frame(Gene = rownames(sift.categories.df.prad),sift.categories.df.prad) -> sift.categories.df.first.col.prad
data.frame(Gene = rownames(PolyPhen.categories.df.prad),PolyPhen.categories.df.prad) -> PolyPhen.categories.df.first.col.prad


# Full outer join of multiple data frames
full_join(df.prad, sift.categories.df.first.col.prad, by = "Gene") %>%
  full_join(., PolyPhen.categories.df.first.col.prad, by = "Gene") -> PRAD.final
View(PRAD.final)

colnames(PRAD.final)[2] <- "Frequency.Gene.Is.Mutated.prad"

# Percentages
install.packages("formattable")
library(formattable)

# Percent Low Impact Mutations
PRAD.final$Percent.Low.Impact.Mutations.prad = percent(PRAD.final$Low.Impact.Mutations.prad/PRAD.final$Frequency.Gene.Is.Mutated.prad) -> Percent.Low.Impact.Mutations

# Remove % symbol
PRAD.final$Percent.Low.Impact.Mutations.prad = as.numeric(gsub("[\\%,]", "", PRAD.final$Percent.Low.Impact.Mutations.prad))

# Percent Moderate Impact Mutations
PRAD.final$Percent.Moderate.Impact.Mutations.prad = percent(PRAD.final$Moderate.Impact.Mutations.prad/PRAD.final$Frequency.Gene.Is.Mutated.prad) -> Percent.Moderate.Impact.Mutations

# Remove % symbol
PRAD.final$Percent.Moderate.Impact.Mutations.prad = as.numeric(gsub("[\\%,]", "", PRAD.final$Percent.Moderate.Impact.Mutations.prad))

# Percent High Impact Mutations
PRAD.final$Percent.High.Impact.Mutations.prad = percent(PRAD.final$High.Impact.Mutations.prad/PRAD.final$Frequency.Gene.Is.Mutated.prad) -> Percent.High.Impact.Mutations

# Remove % symbol
PRAD.final$Percent.High.Impact.Mutations.prad = as.numeric(gsub("[\\%,]", "", PRAD.final$Percent.High.Impact.Mutations.prad))

# Percent Modifier Mutations
PRAD.final$Percent.Modifier.Mutations.prad = percent(PRAD.final$Modifier.Mutations.prad/PRAD.final$Frequency.Gene.Is.Mutated.prad) -> Percent.Modifier.Mutations

# Remove % symbol
PRAD.final$Percent.Modifier.Mutation.prad = as.numeric(gsub("[\\%,]", "", PRAD.final$Percent.Modifier.Mutations.prad))

# Percent Missense Mutation
PRAD.final$Percent.Missense.Mutation.prad = percent(PRAD.final$Missense.Mutation.prad/PRAD.final$Frequency.Gene.Is.Mutated.prad) -> Percent.Missense.Mutation

# Remove % symbol
PRAD.final$Percent.Missense.Mutation.prad = as.numeric(gsub("[\\%,]", "", PRAD.final$Percent.Missense.Mutation.prad))

# Percent of patients with that mutation
# Find total number of patients!
prad.data$mutational %>% count(Gene,case_id) -> case.id.gene.prad
data.frame(case.id.gene.prad)
df_uniq.prad <- unique(case.id.gene.prad$case_id)
length(df_uniq.prad)

# Calculate %
PRAD.final$Percent.Patients.With.Mutation.prad = percent(PRAD.final$Number.Patients.With.Mutation.prad/(length(df_uniq.prad)))

# Remove % symbol
PRAD.final$Percent.Patients.With.Mutation.prad = as.numeric(gsub("[\\%,]", "",PRAD.final$Percent.Patients.With.Mutation.prad))


# Reorder columns
PRAD.final2 = PRAD.final %>% select(Gene,
                                    Hugo.Symbol.prad...2.,
                                    Gene.Ontology.prad...2.,
                                    Frequency.Gene.Is.Mutated.prad, 
                                    Number.Patients.With.Mutation.prad...2., 
                                    Percent.Patients.With.Mutation.prad, 
                                    Modifier.Mutations.prad...2.,
                                    Percent.Modifier.Mutations.prad,
                                    Low.Impact.Mutations.prad...2.,
                                    Percent.Low.Impact.Mutations.prad,
                                    Moderate.Impact.Mutations.prad...2., 
                                    Percent.Moderate.Impact.Mutations.prad,
                                    High.Impact.Mutations.prad...2.,
                                    Percent.High.Impact.Mutations.prad,
                                    Missense.Mutation.prad...2., 
                                    Percent.Missense.Mutation.prad,
                                    SIFT.Deleterious.prad, 
                                    SIFT.Deleterious_low_confidence.prad, 
                                    SIFT.Tolerated.prad, 
                                    SIFT.Tolerated_low_confidence.prad, 
                                    SIFT.Mean.prad...2.,
                                    SIFT.Median.prad...2., 
                                    SIFT.Maximum.prad...2., 
                                    SIFT.StandardDeviation.prad...2., 
                                    PolyPhen.Benign.prad, 
                                    PolyPhen.Possibly.Damaging.prad,
                                    PolyPhen.Probably.Damaging.prad, 
                                    PolyPhen.Unknown.prad,
                                    PolyPhen.Mean.prad...2.,
                                    PolyPhen.Median.prad...2.,
                                    PolyPhen.Maximum.prad...2.,
                                    PolyPhen.StandardDeviation.prad...2.)


View(PRAD.final2)

colnames(PRAD.final2)[1] <- "Gene"
colnames(PRAD.final2)[2] <- "Hugo Symbol"
colnames(PRAD.final2)[3] <- "Gene Ontology"
colnames(PRAD.final2)[4] <- "Frequency Gene is Mutated"
colnames(PRAD.final2)[5] <- "Number of Patients With Mutation"
colnames(PRAD.final2)[6] <- "Percent of Patients With Mutation"
colnames(PRAD.final2)[7] <- "Number of Modifier Mutations"
colnames(PRAD.final2)[8] <- "Percent of Mutations that are Modifiers"
colnames(PRAD.final2)[9] <- "Number of Low Impact Mutations"
colnames(PRAD.final2)[10] <- "Percent of Mutations that are Low Impact"
colnames(PRAD.final2)[11] <- "Number of Moderate Impact Mutations"
colnames(PRAD.final2)[12] <- "Percent of Mutations that are Moderate Impact"
colnames(PRAD.final2)[13] <- "Number of High Impact Mutations"
colnames(PRAD.final2)[14] <- "Percent of Mutations that are High Impact"
colnames(PRAD.final2)[15] <- "Number of Missense Mutations"
colnames(PRAD.final2)[16] <- "Percent of Mutations that are Missense"
colnames(PRAD.final2)[17] <- "SIFT Deleterious"
colnames(PRAD.final2)[18] <- "SIFT Deleterious Low Confidence"
colnames(PRAD.final2)[19] <- "SIFT Tolerated"
colnames(PRAD.final2)[20] <- "SIFT Tolerated Low Confidence"
colnames(PRAD.final2)[21] <- "SIFT Mean"
colnames(PRAD.final2)[22] <- "SIFT Median"
colnames(PRAD.final2)[23] <- "SIFT Maximum"
colnames(PRAD.final2)[24] <- "SIFT Standard Deviation"
colnames(PRAD.final2)[25] <- "PolyPhen Benign"
colnames(PRAD.final2)[26] <- "PolyPhen Possibly Damaging"
colnames(PRAD.final2)[27] <- "PolyPhen Probably Damaging"
colnames(PRAD.final2)[28] <- "PolyPhen Unknown"
colnames(PRAD.final2)[29] <- "PolyPhen Mean"
colnames(PRAD.final2)[30] <- "PolyPhen Median"
colnames(PRAD.final2)[31] <- "PolyPhen Maximum"
colnames(PRAD.final2)[32] <- "PolyPhen Standard Deviation"

write.csv(PRAD.final2, "~/tcga_biolinks1/PRAD.mutation.table.final.csv", row.names=TRUE, na = "")

################################################################################################

#### BRCA ######################################################################
head(brca.data$mutational)
colnames(brca.data$mutational)

# Count number of mutations per gene
brca.data$mutational %>%
  count(Gene) -> Total.Number.Mutations.brca
Total.Number.Mutations.brca

# Classify mutations into 3 groups based on IMPACT score
brca.data$mutational %>% dplyr::count(brca.data$mutational$Variant_Classification, brca.data$mutational$IMPACT)

# High impact
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="HIGH") %>% 
  count(Gene) -> High.Impact.Mutations.brca
High.Impact.Mutations.brca

# Moderate impact
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODERATE") %>% 
  count(Gene) -> Moderate.Impact.Mutations.brca
Moderate.Impact.Mutations.brca

# Low impact
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="LOW") %>% 
  count(Gene) -> Low.Impact.Mutations.brca
Low.Impact.Mutations.brca

# Modifier impact
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODIFIER") %>% 
  count(Gene) -> Modifier.Mutations.brca
Modifier.Mutations.brca

# Filter for missense mutations, find number of missense mutations per gene
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  count(Gene) -> no.missense.brca
no.missense.brca

# Separate SIFT data into category and numerical
# SIFT category
library(qdapRegex)
rm_round(text.var =  brca.data$mutational$SIFT, extract = F) -> brca.data$mutational$SIFT.cat

brca.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(SIFT.cat) %>%
  select(Gene,SIFT.cat) -> sift.categories.brca

# Count by category type
table(sift.categories.brca) -> sift.categories.tab.brca
as.data.frame.matrix(sift.categories.tab.brca) -> sift.categories.df.brca

# Rename column headings
colnames(sift.categories.df.brca)[1]  <- "SIFT.Deleterious.brca"
colnames(sift.categories.df.brca)[2]  <- "SIFT.Deleterious_low_confidence.brca"
colnames(sift.categories.df.brca)[3]  <- "SIFT.Tolerated.brca"
colnames(sift.categories.df.brca)[4]  <- "SIFT.Tolerated_low_confidence.brca"
sift.categories.df.brca

# SIFT numerical values
unlist(rm_round(text.var =  brca.data$mutational$SIFT, extract = T)) -> brca.data$mutational$SIFT.num

brca.data$mutational$SIFT.num <- as.numeric(brca.data$mutational$SIFT.num)

brca.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(SIFT.num) %>%
  select(Gene,SIFT.num) -> sift.number.brca

sift.number.brca

# Summarise sift.number - ignore NAs
# Mean
sift.mean.brca <- sift.number.brca %>%
  group_by(Gene) %>%
  summarise(mean = mean(SIFT.num, na.rm = TRUE))

sift.mean.brca

# Median
sift.median.brca <- sift.number.brca %>%
  group_by(Gene) %>%
  summarise(median = median(SIFT.num, na.rm = TRUE))

sift.median.brca

# Maximum
sift.maximum.brca <- sift.number.brca %>%
  group_by(Gene) %>%
  summarise(max = max(SIFT.num, na.rm = TRUE))

sift.maximum.brca

# Standard deviation
sift.sd.brca <- sift.number.brca %>%
  group_by(Gene) %>%
  summarise(sd = sd(SIFT.num, na.rm = TRUE))
sift.sd.brca

# Separate PolyPhen data into category and numerical
# polyphen category
library(qdapRegex)
rm_round(text.var =  brca.data$mutational$PolyPhen, extract = F) -> brca.data$mutational$PolyPhen.cat

brca.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(PolyPhen.cat) %>%
  select(Gene,PolyPhen.cat) -> PolyPhen.categories.brca

PolyPhen.categories.brca

# Count categories
table(PolyPhen.categories.brca) -> PolyPhen.categories.tab.brca
as.data.frame.matrix(PolyPhen.categories.tab.brca) -> PolyPhen.categories.df.brca

# Rename column headings
colnames(PolyPhen.categories.df.brca)[1]  <- "PolyPhen.Benign.brca"
colnames(PolyPhen.categories.df.brca)[2]  <- "PolyPhen.Possibly.Damaging.brca"
colnames(PolyPhen.categories.df.brca)[3]  <- "PolyPhen.Probably.Damaging.brca"
colnames(PolyPhen.categories.df.brca)[4]  <- "PolyPhen.Unknown.brca"
PolyPhen.categories.df.brca

# PolyPhen numerical values
unlist(rm_round(text.var =  brca.data$mutational$PolyPhen, extract = T)) -> brca.data$mutational$PolyPhen.num

brca.data$mutational$PolyPhen.num <- as.numeric(brca.data$mutational$PolyPhen.num)

brca.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(PolyPhen.num) %>%
  select(Gene,PolyPhen.num) -> PolyPhen.number.brca

PolyPhen.number.brca

# Summarise PolyPhen.number - ignore NAs
# Mean
PolyPhen.mean.brca <- PolyPhen.number.brca %>%
  group_by(Gene) %>%
  summarise(mean = mean(PolyPhen.num, na.rm = TRUE))

PolyPhen.mean.brca

# Median
PolyPhen.median.brca <- PolyPhen.number.brca %>%
  group_by(Gene) %>%
  summarise(median = median(PolyPhen.num, na.rm = TRUE))

PolyPhen.median.brca

# Maximum
PolyPhen.maximum.brca <- PolyPhen.number.brca %>%
  group_by(Gene) %>%
  summarise(max = max(PolyPhen.num, na.rm = TRUE))

PolyPhen.maximum.brca

# Standard deviation
PolyPhen.sd.brca <- PolyPhen.number.brca %>%
  group_by(Gene) %>%
  summarise(sd = sd(PolyPhen.num, na.rm = TRUE))

PolyPhen.sd.brca

# Case_id - find number of patients with that gene mutated
brca.data$mutational %>% count(Gene,case_id) -> Gene.Patient.brca
table(Gene.Patient.brca$Gene) -> Gene.Patient.tab.brca
as.data.frame(Gene.Patient.tab.brca) -> Patients.With.Mutation.brca
colnames(Patients.With.Mutation.brca)[1] <- "Gene"
colnames(Patients.With.Mutation.brca)[2] <- "Patients.With.Mutation"
Patients.With.Mutation.brca

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

# Retrieve GO terms for mutatated genes
retrieved.data.brca <- annotate.HTseq.IDs(Total.Number.Mutations.brca$Gene)
retrieved.data.brca

# Rename column
colnames(retrieved.data.brca)[2] <- "Gene"

# Select columns needed
retrieved.data.brca %>% dplyr::select(Gene, go_id) -> GO.brca

# Duplicate data frame
Gene.Ontology.brca <- GO.brca  

# Replace blanks with NA
Gene.Ontology.brca[Gene.Ontology.brca == ""] <- NA   

# Print 
Gene.Ontology.brca

# Select columns needed
retrieved.data.brca %>% dplyr::select(Gene, hgnc_symbol) -> hgnc.brca

# Duplicate data frame
Hugo.Symbol.brca <- hgnc.brca

# Print
Hugo.Symbol.brca


# Construct function for matching different columns of data to same genes and length
match.select <- function(x, y, label){
  left_join(x, y, by = "Gene") -> temp.df
  temp.df[,c(1,3)] -> temp.df
  colnames(temp.df)[2] <- label
  return(temp.df)
}

# Match columns
Number.Patients.With.Mutation.brca <- match.select(x = Total.Number.Mutations.brca, y = Patients.With.Mutation.brca, label = "Number.Patients.With.Mutation")
Low.Impact.Mutations.brca <- match.select(x = Total.Number.Mutations.brca, y = Low.Impact.Mutations.brca, label = "Low.Impact.Mutations")
Moderate.Impact.Mutations.brca <- match.select(x = Total.Number.Mutations.brca, y = Moderate.Impact.Mutations.brca, label = "Moderate.Impact.Mutations")
High.Impact.Mutations.brca <- match.select(x = Total.Number.Mutations.brca, y = High.Impact.Mutations.brca, label = "High.Impact.Mutations")
Modifier.Mutations.brca <- match.select(x = Total.Number.Mutations.brca, y = Modifier.Mutations.brca, label = "Modifier.Mutations")
Missense.Mutation.brca <- match.select(x = Total.Number.Mutations.brca, y = no.missense.brca, label = "no.missense")
SIFT.Mean.brca <- match.select(x = Total.Number.Mutations.brca, y = sift.mean.brca, label = "sift.mean")
SIFT.Median.brca <- match.select(x = Total.Number.Mutations.brca, y = sift.median.brca, label = "sift.median")
SIFT.Maximum.brca <- match.select(x = Total.Number.Mutations.brca, y = sift.maximum.brca, label = "sift.maximum")
SIFT.StandardDeviation.brca <- match.select(x = Total.Number.Mutations.brca, y = sift.sd.brca, label = "sift.sd")
PolyPhen.Mean.brca <- match.select(x = Total.Number.Mutations.brca, y = PolyPhen.mean.brca, label = "PolyPhen.mean")
PolyPhen.Median.brca <- match.select(x = Total.Number.Mutations.brca, y = PolyPhen.median.brca, label = "PolyPhen.median")
PolyPhen.Maximum.brca <- match.select(x = Total.Number.Mutations.brca, y = PolyPhen.maximum.brca, label = "PolyPhen.maximum")
PolyPhen.StandardDeviation.brca <- match.select(x = Total.Number.Mutations.brca, y = PolyPhen.sd.brca, label = "PolyPhen.sd")

# Create dataframe
data.frame(Total.Number.Mutations.brca,
           Hugo.Symbol.brca[,2],
           Gene.Ontology.brca[,2],
           Number.Patients.With.Mutation.brca[,2],
           Low.Impact.Mutations.brca[,2],
           Moderate.Impact.Mutations.brca[,2],
           High.Impact.Mutations.brca[,2],
           Modifier.Mutations.brca[,2],
           Missense.Mutation.brca[,2],
           SIFT.Mean.brca[,2],
           SIFT.Median.brca[,2],
           SIFT.Maximum.brca [,2],
           SIFT.StandardDeviation.brca[,2],
           PolyPhen.Mean.brca[,2],
           PolyPhen.Median.brca[,2],
           PolyPhen.Maximum.brca[,2],
           PolyPhen.StandardDeviation.brca[,2]     
) -> df.brca

View(df.brca)
# Join together df, with all other data frames 
# Add another "Gene" column to allow matching
data.frame(Gene = rownames(sift.categories.df.brca),sift.categories.df.brca) -> sift.categories.df.first.col.brca
data.frame(Gene = rownames(PolyPhen.categories.df.brca),PolyPhen.categories.df.brca) -> PolyPhen.categories.df.first.col.brca


# Full outer join of multiple data frames
full_join(df.brca, sift.categories.df.first.col.brca, by = "Gene") %>%
  full_join(., PolyPhen.categories.df.first.col.brca, by = "Gene") -> BRCA.final
View(BRCA.final)

colnames(BRCA.final)[2] <- "Frequency.Gene.Is.Mutated.brca"

# Percentages
install.packages("formattable")
library(formattable)

# Percent Low Impact Mutations
BRCA.final$Percent.Low.Impact.Mutations.brca = percent(BRCA.final$Low.Impact.Mutations.brca/BRCA.final$Frequency.Gene.Is.Mutated.brca) -> Percent.Low.Impact.Mutations

# Remove % symbol
BRCA.final$Percent.Low.Impact.Mutations.brca = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.Low.Impact.Mutations.brca))

# Percent Moderate Impact Mutations
BRCA.final$Percent.Moderate.Impact.Mutations.brca = percent(BRCA.final$Moderate.Impact.Mutations.brca/BRCA.final$Frequency.Gene.Is.Mutated.brca) -> Percent.Moderate.Impact.Mutations

# Remove % symbol
BRCA.final$Percent.Moderate.Impact.Mutations.brca= as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.Moderate.Impact.Mutations.brca))

# Percent High Impact Mutations
BRCA.final$Percent.High.Impact.Mutations.brca = percent(BRCA.final$High.Impact.Mutations.brca/BRCA.final$Frequency.Gene.Is.Mutated.brca) -> Percent.High.Impact.Mutations

# Remove % symbol
BRCA.final$Percent.High.Impact.Mutations.brca = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.High.Impact.Mutations.brca))

# Percent Modifier Mutations
BRCA.final$Percent.Modifier.Mutations.brca = percent(BRCA.final$Modifier.Mutations.brca/BRCA.final$Frequency.Gene.Is.Mutated.brca) -> Percent.Modifier.Mutations

# Remove % symbol
BRCA.final$Percent.Modifier.Mutation.brca = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.Modifier.Mutations.brca))

# Percent Missense Mutation
BRCA.final$Percent.Missense.Mutation.brca = percent(BRCA.final$Missense.Mutation.brca/BRCA.final$Frequency.Gene.Is.Mutated.brca) -> Percent.Missense.Mutation

# Remove % symbol
BRCA.final$Percent.Missense.Mutation.brca = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.Missense.Mutation.brca))

# Percent of patients with that mutation
# Find total number of patients!
brca.data$mutational %>% count(Gene,case_id) -> case.id.gene.brca
data.frame(case.id.gene.brca)
df_uniq.brca <- unique(case.id.gene.brca$case_id)
length(df_uniq.brca)

# Calculate %
BRCA.final$Percent.Patients.With.Mutation.brca = percent(BRCA.final$Number.Patients.With.Mutation.brca/(length(df_uniq.brca)))

# Remove % symbol
BRCA.final$Percent.Patients.With.Mutation.brca = as.numeric(gsub("[\\%,]", "",BRCA.final$Percent.Patients.With.Mutation.brca))


# Reorder columns
BRCA.final2 = BRCA.final %>% select(Gene,
                                    Hugo.Symbol.brca...2.,
                                    Gene.Ontology.brca...2.,
                                    Frequency.Gene.Is.Mutated.brca, 
                                    Number.Patients.With.Mutation.brca...2., 
                                    Percent.Patients.With.Mutation.brca, 
                                    Modifier.Mutations.brca...2.,
                                    Percent.Modifier.Mutations.brca,
                                    Low.Impact.Mutations.brca...2.,
                                    Percent.Low.Impact.Mutations.brca,
                                    Moderate.Impact.Mutations.brca...2., 
                                    Percent.Moderate.Impact.Mutations.brca,
                                    High.Impact.Mutations.brca...2.,
                                    Percent.High.Impact.Mutations.brca,
                                    Missense.Mutation.brca...2., 
                                    Percent.Missense.Mutation.brca,
                                    SIFT.Deleterious.brca, 
                                    SIFT.Deleterious_low_confidence.brca, 
                                    SIFT.Tolerated.brca, 
                                    SIFT.Tolerated_low_confidence.brca, 
                                    SIFT.Mean.brca...2.,
                                    SIFT.Median.brca...2., 
                                    SIFT.Maximum.brca...2., 
                                    SIFT.StandardDeviation.brca...2., 
                                    PolyPhen.Benign.brca, 
                                    PolyPhen.Possibly.Damaging.brca,
                                    PolyPhen.Probably.Damaging.brca, 
                                    PolyPhen.Unknown.brca,
                                    PolyPhen.Mean.brca...2.,
                                    PolyPhen.Median.brca...2.,
                                    PolyPhen.Maximum.brca...2.,
                                    PolyPhen.StandardDeviation.brca...2.)


View(BRCA.final2)


write.csv(BRCA.final2, "~/tcga_biolinks1/brca.mutation.table.final.csv", row.names=TRUE, na = "")

################################################################################################


library(igraph)

read.graph(file = "/home/ug_megan/tcga_biolinks1/brca5.network.rename.txt(1).graphml", format = "graphml") -> graph
graph














#### BRCA ######################################################################
head(brca.data$mutational)
colnames(brca.data$mutational)

# Count number of mutations per gene
brca.data$mutational %>%
  count(Gene) -> Total.Number.Mutations
Total.Number.Mutations

# Classify mutations into 3 groups based on IMPACT score
brca.data$mutational %>% dplyr::count(brca.data$mutational$Variant_Classification, brca.data$mutational$IMPACT)

# High impact
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="HIGH") %>% 
  count(Gene) -> High.Impact.Mutations
High.Impact.Mutations

# Moderate impact
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODERATE") %>% 
  count(Gene) -> Moderate.Impact.Mutations
Moderate.Impact.Mutations

# Low impact
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="LOW") %>% 
  count(Gene) -> Low.Impact.Mutations
Low.Impact.Mutations

# Modifier impact
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODIFIER") %>% 
  count(Gene) -> Modifier.Mutations
Modifier.Mutations

# Filter for missense mutations, find number of missense mutations per gene
brca.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  count(Gene) -> no.missense
no.missense

# Separate SIFT data into category and numerical
# SIFT category
library(qdapRegex)
rm_round(text.var =  brca.data$mutational$SIFT, extract = F) -> brca.data$mutational$SIFT.cat

brca.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(SIFT.cat) %>%
  select(Gene,SIFT.cat) -> sift.categories

# Count by category type
table(sift.categories) -> sift.categories.tab
as.data.frame.matrix(sift.categories.tab) -> sift.categories.df

# Rename column headings
colnames(sift.categories.df)[1]  <- "SIFT.Deleterious"
colnames(sift.categories.df)[2]  <- "SIFT.Deleterious_low_confidence"
colnames(sift.categories.df)[3]  <- "SIFT.Tolerated"
colnames(sift.categories.df)[4]  <- "SIFT.Tolerated_low_confidence"
sift.categories.df

# SIFT numerical values
unlist(rm_round(text.var =  brca.data$mutational$SIFT, extract = T)) -> brca.data$mutational$SIFT.num

brca.data$mutational$SIFT.num <- as.numeric(brca.data$mutational$SIFT.num)

brca.data$mutational %>% 
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
rm_round(text.var =  brca.data$mutational$PolyPhen, extract = F) -> brca.data$mutational$PolyPhen.cat

brca.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(PolyPhen.cat) %>%
  select(Gene,PolyPhen.cat) -> PolyPhen.categories

PolyPhen.categories

# Count categories
table(PolyPhen.categories) -> PolyPhen.categories.tab
as.data.frame.matrix(PolyPhen.categories.tab) -> PolyPhen.categories.df

# Rename column headings
colnames(PolyPhen.categories.df)[1]  <- "PolyPhen.Benign"
colnames(PolyPhen.categories.df)[2]  <- "PolyPhen.Possibly.Damaging"
colnames(PolyPhen.categories.df)[3]  <- "PolyPhen.Probably.Damaging"
colnames(PolyPhen.categories.df)[4]  <- "PolyPhen.Unknown"
PolyPhen.categories.df

# PolyPhen numerical values
unlist(rm_round(text.var =  brca.data$mutational$PolyPhen, extract = T)) -> brca.data$mutational$PolyPhen.num

brca.data$mutational$PolyPhen.num <- as.numeric(brca.data$mutational$PolyPhen.num)

brca.data$mutational %>% 
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

# Case_id - find number of patients with that gene mutated
brca.data$mutational %>% count(Gene,case_id) -> Gene.Patient
table(Gene.Patient$Gene) -> Gene.Patient.tab
as.data.frame(Gene.Patient.tab) -> Patients.With.Mutation
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

# Retrieve GO terms for mutatated genes
retrieved.data <- annotate.HTseq.IDs(Total.Number.Mutations$Gene)
retrieved.data

# Rename column
colnames(retrieved.data)[2] <- "Gene"

# Select columns needed
retrieved.data %>% dplyr::select(Gene, go_id) -> GO

# Duplicate data frame
Gene.Ontology <- GO  

# Replace blanks with NA
Gene.Ontology[Gene.Ontology == ""] <- NA   

# Print 
Gene.Ontology

# Select columns needed
retrieved.data %>% dplyr::select(Gene, hgnc_symbol) -> hgnc

# Duplicate data frame
Hugo.Symbol <- hgnc

# Print
Hugo.Symbol

# Construct function for matching different columns of data to same genes and length
match.select <- function(x, y, label){
  left_join(x, y, by = "Gene") -> temp.df
  temp.df[,c(1,3)] -> temp.df
  colnames(temp.df)[2] <- label
  return(temp.df)
}

# Match all columns
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
           Hugo.Symbol[,2],
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
  full_join(., PolyPhen.categories.df.first.col, by = "Gene") -> BRCA.final
View(BRCA.final)

colnames(BRCA.final)[2] <- "Frequency.Gene.Is.Mutated"

# Percentages
install.packages("formattable")
library(formattable)
# Percent Low Impact Mutations
BRCA.final$Percent.Low.Impact.Mutations = percent(BRCA.final$Low.Impact.Mutations/BRCA.final$Frequency.Gene.Is.Mutated) -> Percent.Low.Impact.Mutations

# Remove % symbol
BRCA.final$Percent.Low.Impact.Mutations = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.Low.Impact.Mutations))

# Percent Moderate Impact Mutations
BRCA.final$Percent.Moderate.Impact.Mutations = percent(BRCA.final$Moderate.Impact.Mutations/BRCA.final$Frequency.Gene.Is.Mutated) -> Percent.Moderate.Impact.Mutations

# Remove % symbol
BRCA.final$Percent.Moderate.Impact.Mutations = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.Moderate.Impact.Mutations))

# Percent High Impact Mutations
BRCA.final$Percent.High.Impact.Mutations = percent(BRCA.final$High.Impact.Mutations/BRCA.final$Frequency.Gene.Is.Mutated) -> Percent.High.Impact.Mutations

# Remove % symbol
BRCA.final$Percent.High.Impact.Mutations = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.High.Impact.Mutations))

# Percent Modifier Mutations
BRCA.final$Percent.Modifier.Mutations = percent(BRCA.final$Modifier.Mutations/BRCA.final$Frequency.Gene.Is.Mutated) -> Percent.Modifier.Mutations

# Remove % symbol
BRCA.final$Percent.Modifier.Mutations = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.Modifier.Mutations))

# Percent Missense Mutation
BRCA.final$Percent.Missense.Mutation = percent(BRCA.final$Missense.Mutation/BRCA.final$Frequency.Gene.Is.Mutated) -> Percent.Missense.Mutation

# Remove % symbol
BRCA.final$Percent.Missense.Mutation = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.Missense.Mutation))

# Percent of patients with that mutation
# Find total number of patients!
brca.data$mutational %>% count(Gene,case_id) -> case.id.gene
data.frame(case.id.gene)
df_uniq <- unique(case.id.gene$case_id)
length(df_uniq)

# Calculate %
BRCA.final$Percent.Patients.With.Mutation = percent(BRCA.final$Number.Patients.With.Mutation/(length(df_uniq)))

# Remove % symbol
BRCA.final$Percent.Patients.With.Mutation = as.numeric(gsub("[\\%,]", "", BRCA.final$Percent.Patients.With.Mutation))

# Annotate with GO terms using Biomart
library(biomaRt)

# Reorder columns
BRCA.final2 = BRCA.final %>% select(Gene, 
                                    Hugo.Symbol...2.,
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


View(BRCA.final2)

sapply(BRCA.final2, mode)
sapply(BRCA.final2, class)

write.csv(BRCA.final2, "~/tcga_biolinks1/BRCA.final2.csv", row.names=TRUE, na = "")

#### LUAD ######################################################################
head(luad.data$mutational)
colnames(luad.data$mutational)

# Count number of mutations per gene
luad.data$mutational %>%
  count(Gene) -> Total.Number.Mutations
Total.Number.Mutations

# Classify mutations into 3 groups based on IMPACT score
luad.data$mutational %>% dplyr::count(luad.data$mutational$Variant_Classification, luad.data$mutational$IMPACT)

# High impact
luad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="HIGH") %>% 
  count(Gene) -> High.Impact.Mutations
High.Impact.Mutations

# Moderate impact
luad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODERATE") %>% 
  count(Gene) -> Moderate.Impact.Mutations
Moderate.Impact.Mutations

# Low impact
luad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="LOW") %>% 
  count(Gene) -> Low.Impact.Mutations
Low.Impact.Mutations

# Modifier impact
luad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODIFIER") %>% 
  count(Gene) -> Modifier.Mutations
Modifier.Mutations

# Filter for missense mutations, find number of missense mutations per gene
luad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  count(Gene) -> no.missense
no.missense

# Separate SIFT data into category and numerical
# SIFT category
library(qdapRegex)
rm_round(text.var =  luad.data$mutational$SIFT, extract = F) -> luad.data$mutational$SIFT.cat

luad.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(SIFT.cat) %>%
  select(Gene,SIFT.cat) -> sift.categories

# Count by category type
table(sift.categories) -> sift.categories.tab
as.data.frame.matrix(sift.categories.tab) -> sift.categories.df

# Rename column headings
colnames(sift.categories.df)[1]  <- "SIFT.Deleterious"
colnames(sift.categories.df)[2]  <- "SIFT.Deleterious_low_confidence"
colnames(sift.categories.df)[3]  <- "SIFT.Tolerated"
colnames(sift.categories.df)[4]  <- "SIFT.Tolerated_low_confidence"
sift.categories.df

# SIFT numerical values
unlist(rm_round(text.var =  luad.data$mutational$SIFT, extract = T)) -> luad.data$mutational$SIFT.num

luad.data$mutational$SIFT.num <- as.numeric(luad.data$mutational$SIFT.num)

luad.data$mutational %>% 
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
rm_round(text.var =  luad.data$mutational$PolyPhen, extract = F) -> luad.data$mutational$PolyPhen.cat

luad.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(PolyPhen.cat) %>%
  select(Gene,PolyPhen.cat) -> PolyPhen.categories

PolyPhen.categories

# Count categories
table(PolyPhen.categories) -> PolyPhen.categories.tab
as.data.frame.matrix(PolyPhen.categories.tab) -> PolyPhen.categories.df

# Rename column headings
colnames(PolyPhen.categories.df)[1]  <- "PolyPhen.Benign"
colnames(PolyPhen.categories.df)[2]  <- "PolyPhen.Possibly.Damaging"
colnames(PolyPhen.categories.df)[3]  <- "PolyPhen.Probably.Damaging"
colnames(PolyPhen.categories.df)[4]  <- "PolyPhen.Unknown"
PolyPhen.categories.df

# PolyPhen numerical values
unlist(rm_round(text.var =  luad.data$mutational$PolyPhen, extract = T)) -> luad.data$mutational$PolyPhen.num

luad.data$mutational$PolyPhen.num <- as.numeric(luad.data$mutational$PolyPhen.num)

luad.data$mutational %>% 
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

# Case_id - find number of patients with that gene mutated
luad.data$mutational %>% count(Gene,case_id) -> Gene.Patient
table(Gene.Patient$Gene) -> Gene.Patient.tab
as.data.frame(Gene.Patient.tab) -> Patients.With.Mutation
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

# Retrieve GO terms for mutatated genes
retrieved.data <- annotate.HTseq.IDs(Total.Number.Mutations$Gene)
retrieved.data

# Rename column
colnames(retrieved.data)[2] <- "Gene"

# Select columns needed
retrieved.data %>% dplyr::select(Gene, go_id) -> GO

# Duplicate data frame
Gene.Ontology <- GO  

# Replace blanks with NA
Gene.Ontology[Gene.Ontology == ""] <- NA   

# Print 
Gene.Ontology

# Select columns needed
retrieved.data %>% dplyr::select(Gene, hgnc_symbol) -> hgnc

# Duplicate data frame
Hugo.Symbol <- hgnc

# Print
Hugo.Symbol

# Construct function for matching different columns of data to same genes and length
match.select <- function(x, y, label){
  left_join(x, y, by = "Gene") -> temp.df
  temp.df[,c(1,3)] -> temp.df
  colnames(temp.df)[2] <- label
  return(temp.df)
}

# Match all columns
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
           Hugo.Symbol[,2],
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
  full_join(., PolyPhen.categories.df.first.col, by = "Gene") -> LUAD.final
View(LUAD.final)

colnames(LUAD.final)[2] <- "Frequency.Gene.Is.Mutated"

# Percentages
install.packages("formattable")
library(formattable)
LUAD.final$Percent.Low.Impact.Mutations = percent(LUAD.final$Low.Impact.Mutations/LUAD.final$Frequency.Gene.Is.Mutated) -> Percent.Low.Impact.Mutations

# Remove % symbol
LUAD.final$Percent.Low.Impact.Mutations = as.numeric(gsub("[\\%,]", "", LUAD.final$Percent.Low.Impact.Mutations))

# Percent Moderate Impact Mutations
LUAD.final$Percent.Moderate.Impact.Mutations = percent(LUAD.final$Moderate.Impact.Mutations/LUAD.final$Frequency.Gene.Is.Mutated) -> Percent.Moderate.Impact.Mutations

# Remove % symbol
LUAD.final$Percent.Moderate.Impact.Mutations = as.numeric(gsub("[\\%,]", "", LUAD.final$Percent.Moderate.Impact.Mutations))

# Percent High Impact Mutations
LUAD.final$Percent.High.Impact.Mutations = percent(LUAD.final$High.Impact.Mutations/LUAD.final$Frequency.Gene.Is.Mutated) -> Percent.High.Impact.Mutations

# Remove % symbol
LUAD.final$Percent.High.Impact.Mutations = as.numeric(gsub("[\\%,]", "", LUAD.final$Percent.High.Impact.Mutations))

# Percent Modifier Mutations
LUAD.final$Percent.Modifier.Mutations = percent(LUAD.final$Modifier.Mutations/LUAD.final$Frequency.Gene.Is.Mutated) -> Percent.Modifier.Mutations

# Remove % symbol
LUAD.final$Percent.Modifier.Mutations = as.numeric(gsub("[\\%,]", "", LUAD.final$Percent.Modifier.Mutations))

# Percent Missense Mutation
LUAD.final$Percent.Missense.Mutation = percent(LUAD.final$Missense.Mutation/LUAD.final$Frequency.Gene.Is.Mutated) -> Percent.Missense.Mutation

# Remove % symbol
LUAD.final$Percent.Missense.Mutation = as.numeric(gsub("[\\%,]", "", LUAD.final$Percent.Missense.Mutation))


# Percent of patients with that mutation
# Find total number of patients!
luad.data$mutational %>% count(Gene,case_id) -> case.id.gene
data.frame(case.id.gene)
df_uniq <- unique(case.id.gene$case_id)
length(df_uniq)

# Calculate %
LUAD.final$Percent.Patients.With.Mutation = percent(LUAD.final$Number.Patients.With.Mutation/(length(df_uniq)))

# Remove % symbol
LUAD.final$Percent.Patients.With.Mutation = as.numeric(gsub("[\\%,]", "",LUAD.final$Percent.Patients.With.Mutation))

# Annotate with GO terms using Biomart
library(biomaRt)

# Reorder columns
LUAD.final2 = LUAD.final %>% select(Gene, 
                                    Hugo.Symbol...2.,
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


View(LUAD.final2)

write.csv(LUAD.final2, "~/tcga_biolinks1/LUAD.final2.csv", row.names=TRUE,na = "")

#### STAD ######################################################################
head(stad.data$mutational)
colnames(stad.data$mutational)

# Count number of mutations per gene
stad.data$mutational %>%
  count(Gene) -> Total.Number.Mutations
Total.Number.Mutations

# Classify mutations into 3 groups based on IMPACT score
stad.data$mutational %>% dplyr::count(stad.data$mutational$Variant_Classification, stad.data$mutational$IMPACT)

# High impact
stad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="HIGH") %>% 
  count(Gene) -> High.Impact.Mutations
High.Impact.Mutations

# Moderate impact
stad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODERATE") %>% 
  count(Gene) -> Moderate.Impact.Mutations
Moderate.Impact.Mutations

# Low impact
stad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="LOW") %>% 
  count(Gene) -> Low.Impact.Mutations
Low.Impact.Mutations

# Modifier impact
stad.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODIFIER") %>% 
  count(Gene) -> Modifier.Mutations
Modifier.Mutations

# Filter for missense mutations, find number of missense mutations per gene
stad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  count(Gene) -> no.missense
no.missense

# Separate SIFT data into category and numerical
# SIFT category
library(qdapRegex)
rm_round(text.var =  stad.data$mutational$SIFT, extract = F) -> stad.data$mutational$SIFT.cat

stad.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(SIFT.cat) %>%
  select(Gene,SIFT.cat) -> sift.categories

# Count by category type
table(sift.categories) -> sift.categories.tab
as.data.frame.matrix(sift.categories.tab) -> sift.categories.df

# Rename column headings
colnames(sift.categories.df)[1]  <- "SIFT.Deleterious"
colnames(sift.categories.df)[2]  <- "SIFT.Deleterious_low_confidence"
colnames(sift.categories.df)[3]  <- "SIFT.Tolerated"
colnames(sift.categories.df)[4]  <- "SIFT.Tolerated_low_confidence"
sift.categories.df

# SIFT numerical values
unlist(rm_round(text.var =  stad.data$mutational$SIFT, extract = T)) -> stad.data$mutational$SIFT.num

stad.data$mutational$SIFT.num <- as.numeric(stad.data$mutational$SIFT.num)

stad.data$mutational %>% 
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
rm_round(text.var =  stad.data$mutational$PolyPhen, extract = F) -> stad.data$mutational$PolyPhen.cat

stad.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(PolyPhen.cat) %>%
  select(Gene,PolyPhen.cat) -> PolyPhen.categories

PolyPhen.categories

# Count categories
table(PolyPhen.categories) -> PolyPhen.categories.tab
as.data.frame.matrix(PolyPhen.categories.tab) -> PolyPhen.categories.df

# Rename column headings
colnames(PolyPhen.categories.df)[1]  <- "PolyPhen.Benign"
colnames(PolyPhen.categories.df)[2]  <- "PolyPhen.Possibly.Damaging"
colnames(PolyPhen.categories.df)[3]  <- "PolyPhen.Probably.Damaging"
colnames(PolyPhen.categories.df)[4]  <- "PolyPhen.Unknown"
PolyPhen.categories.df

# PolyPhen numerical values
unlist(rm_round(text.var =  stad.data$mutational$PolyPhen, extract = T)) -> stad.data$mutational$PolyPhen.num

stad.data$mutational$PolyPhen.num <- as.numeric(stad.data$mutational$PolyPhen.num)

stad.data$mutational %>% 
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

# Case_id - find number of patients with that gene mutated
stad.data$mutational %>% count(Gene,case_id) -> Gene.Patient
table(Gene.Patient$Gene) -> Gene.Patient.tab
as.data.frame(Gene.Patient.tab) -> Patients.With.Mutation
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

# Retrieve GO terms for mutatated genes
retrieved.data <- annotate.HTseq.IDs(Total.Number.Mutations$Gene)
retrieved.data

# Rename column
colnames(retrieved.data)[2] <- "Gene"

# Select columns needed
retrieved.data %>% dplyr::select(Gene, go_id) -> GO

# Duplicate data frame
Gene.Ontology <- GO  

# Replace blanks with NA
Gene.Ontology[Gene.Ontology == ""] <- NA   

# Print 
Gene.Ontology

# Select columns needed
retrieved.data %>% dplyr::select(Gene, hgnc_symbol) -> hgnc

# Duplicate data frame
Hugo.Symbol <- hgnc

# Print
Hugo.Symbol

# Construct function for matching different columns of data to same genes and length
match.select <- function(x, y, label){
  left_join(x, y, by = "Gene") -> temp.df
  temp.df[,c(1,3)] -> temp.df
  colnames(temp.df)[2] <- label
  return(temp.df)
}

# Match all columns
Gene.Ontology.Match <- match.select(x = Total.Number.Mutations, y = Gene.Ontology, label = "Gene.Ontology")
Hugo.Symbol.Match <- match.select(x = Total.Number.Mutations, y = Hugo.Symbol, label = "Hugo.Symbol")
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
           Hugo.Symbol.Match[,2],
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
  full_join(., PolyPhen.categories.df.first.col, by = "Gene") -> STAD.final
View(STAD.final)

colnames(STAD.final)[2] <- "Frequency.Gene.Is.Mutated"

# Percentages
install.packages("formattable")
library(formattable)
STAD.final$Percent.Low.Impact.Mutations = percent(STAD.final$Low.Impact.Mutations/STAD.final$Frequency.Gene.Is.Mutated) -> Percent.Low.Impact.Mutations
STAD.final$Percent.Moderate.Impact.Mutations = percent(STAD.final$Moderate.Impact.Mutations/STAD.final$Frequency.Gene.Is.Mutated) -> Percent.Moderate.Impact.Mutations
STAD.final$Percent.High.Impact.Mutations = percent(STAD.final$High.Impact.Mutations/STAD.final$Frequency.Gene.Is.Mutated) -> Percent.High.Impact.Mutations
STAD.final$Percent.Modifier.Mutations = percent(STAD.final$Modifier.Mutations/STAD.final$Frequency.Gene.Is.Mutated) -> Percent.Modifier.Mutations
STAD.final$Percent.Missense.Mutation = percent(STAD.final$Missense.Mutation/STAD.final$Frequency.Gene.Is.Mutated) -> Percent.Missense.Mutation

# Percent of patients with that mutation
# Find total number of patients!
stad.data$mutational %>% count(Gene,case_id) -> case.id.gene
data.frame(case.id.gene)
df_uniq <- unique(case.id.gene$case_id)
length(df_uniq)

# Calculate %
STAD.final$Percent.Patients.With.Mutation = percent(STAD.final$Number.Patients.With.Mutation/(length(df_uniq)))

# Annotate with GO terms using Biomart
library(biomaRt)

# Reorder columns
STAD.final2 = STAD.final %>% select(Gene, 
                                    Hugo.Symbol.Match...2.,
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


View(STAD.final2)

write.csv(STAD.final2, "~/tcga_biolinks1/STAD.final2.csv", row.names=TRUE) 
#### KIRP ######################################################################
head(kirp.data$mutational)
colnames(kirp.data$mutational)

# Count number of mutations per gene
kirp.data$mutational %>%
  count(Gene) -> Total.Number.Mutations
Total.Number.Mutations

# Classify mutations into 3 groups based on IMPACT score
kirp.data$mutational %>% dplyr::count(kirp.data$mutational$Variant_Classification, kirp.data$mutational$IMPACT)

# High impact
kirp.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="HIGH") %>% 
  count(Gene) -> High.Impact.Mutations
High.Impact.Mutations

# Moderate impact
kirp.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODERATE") %>% 
  count(Gene) -> Moderate.Impact.Mutations
Moderate.Impact.Mutations

# Low impact
kirp.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="LOW") %>% 
  count(Gene) -> Low.Impact.Mutations
Low.Impact.Mutations

# Modifier impact
kirp.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODIFIER") %>% 
  count(Gene) -> Modifier.Mutations
Modifier.Mutations

# Filter for missense mutations, find number of missense mutations per gene
kirp.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  count(Gene) -> no.missense
no.missense

# Separate SIFT data into category and numerical
# SIFT category
library(qdapRegex)
rm_round(text.var =  kirp.data$mutational$SIFT, extract = F) -> kirp.data$mutational$SIFT.cat

kirp.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(SIFT.cat) %>%
  select(Gene,SIFT.cat) -> sift.categories

# Count by category type
table(sift.categories) -> sift.categories.tab
as.data.frame.matrix(sift.categories.tab) -> sift.categories.df

# Rename column headings
colnames(sift.categories.df)[1]  <- "SIFT.Deleterious"
colnames(sift.categories.df)[2]  <- "SIFT.Deleterious_low_confidence"
colnames(sift.categories.df)[3]  <- "SIFT.Tolerated"
colnames(sift.categories.df)[4]  <- "SIFT.Tolerated_low_confidence"
sift.categories.df

# SIFT numerical values
unlist(rm_round(text.var =  kirp.data$mutational$SIFT, extract = T)) -> kirp.data$mutational$SIFT.num

kirp.data$mutational$SIFT.num <- as.numeric(kirp.data$mutational$SIFT.num)

kirp.data$mutational %>% 
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
rm_round(text.var =  kirp.data$mutational$PolyPhen, extract = F) -> kirp.data$mutational$PolyPhen.cat

kirp.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(PolyPhen.cat) %>%
  select(Gene,PolyPhen.cat) -> PolyPhen.categories

PolyPhen.categories

# Count categories
table(PolyPhen.categories) -> PolyPhen.categories.tab
as.data.frame.matrix(PolyPhen.categories.tab) -> PolyPhen.categories.df

# Rename column headings
colnames(PolyPhen.categories.df)[1]  <- "PolyPhen.Benign"
colnames(PolyPhen.categories.df)[2]  <- "PolyPhen.Possibly.Damaging"
colnames(PolyPhen.categories.df)[3]  <- "PolyPhen.Probably.Damaging"
colnames(PolyPhen.categories.df)[4]  <- "PolyPhen.Unknown"
PolyPhen.categories.df

# PolyPhen numerical values
unlist(rm_round(text.var =  kirp.data$mutational$PolyPhen, extract = T)) -> kirp.data$mutational$PolyPhen.num

kirp.data$mutational$PolyPhen.num <- as.numeric(kirp.data$mutational$PolyPhen.num)

kirp.data$mutational %>% 
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

# Case_id - find number of patients with that gene mutated
kirp.data$mutational %>% count(Gene,case_id) -> Gene.Patient
table(Gene.Patient$Gene) -> Gene.Patient.tab
as.data.frame(Gene.Patient.tab) -> Patients.With.Mutation
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

# Retrieve GO terms for mutatated genes
retrieved.data <- annotate.HTseq.IDs(Total.Number.Mutations$Gene)
retrieved.data

# Rename column
colnames(retrieved.data)[2] <- "Gene"

# Select columns needed
retrieved.data %>% dplyr::select(Gene, go_id) -> GO

# Duplicate data frame
Gene.Ontology <- GO  

# Replace blanks with NA
Gene.Ontology[Gene.Ontology == ""] <- NA   

# Print 
Gene.Ontology

# Select columns needed
retrieved.data %>% dplyr::select(Gene, hgnc_symbol) -> hgnc

# Duplicate data frame
Hugo.Symbol <- hgnc

# Print
Hugo.Symbol

# Construct function for matching different columns of data to same genes and length
match.select <- function(x, y, label){
  left_join(x, y, by = "Gene") -> temp.df
  temp.df[,c(1,3)] -> temp.df
  colnames(temp.df)[2] <- label
  return(temp.df)
}

# Match all columns
Gene.Ontology.Match <- match.select(x = Total.Number.Mutations, y = Gene.Ontology, label = "Gene.Ontology")
Hugo.Symbol.Match <- match.select(x = Total.Number.Mutations, y = Hugo.Symbol, label = "Hugo.Symbol")
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
           Hugo.Symbol.Match[,2],
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
  full_join(., PolyPhen.categories.df.first.col, by = "Gene") -> KIRP.final
View(KIRP.final)

colnames(KIRP.final)[2] <- "Frequency.Gene.Is.Mutated"

# Percentages
install.packages("formattable")
library(formattable)
KIRP.final$Percent.Low.Impact.Mutations = percent(KIRP.final$Low.Impact.Mutations/KIRP.final$Frequency.Gene.Is.Mutated) -> Percent.Low.Impact.Mutations
KIRP.final$Percent.Moderate.Impact.Mutations = percent(KIRP.final$Moderate.Impact.Mutations/KIRP.final$Frequency.Gene.Is.Mutated) -> Percent.Moderate.Impact.Mutations
KIRP.final$Percent.High.Impact.Mutations = percent(KIRP.final$High.Impact.Mutations/KIRP.final$Frequency.Gene.Is.Mutated) -> Percent.High.Impact.Mutations
KIRP.final$Percent.Modifier.Mutations = percent(KIRP.final$Modifier.Mutations/KIRP.final$Frequency.Gene.Is.Mutated) -> Percent.Modifier.Mutations
KIRP.final$Percent.Missense.Mutation = percent(KIRP.final$Missense.Mutation/KIRP.final$Frequency.Gene.Is.Mutated) -> Percent.Missense.Mutation

# Percent of patients with that mutation
# Find total number of patients!
kirp.data$mutational %>% count(Gene,case_id) -> case.id.gene
data.frame(case.id.gene)
df_uniq <- unique(case.id.gene$case_id)
length(df_uniq)

# Calculate %
KIRP.final$Percent.Patients.With.Mutation = percent(KIRP.final$Number.Patients.With.Mutation/(length(df_uniq)))

# Annotate with GO terms using Biomart
library(biomaRt)

# Reorder columns
KIRP.final2 = KIRP.final %>% select(Gene,
                                    Hugo.Symbol.Match...2.,
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


View(KIRP.final2)

write.csv(KIRP.final2, "~/tcga_biolinks1/KIRP.final2.csv", row.names=TRUE) 
#### KIRC ######################################################################
head(kirc.data$mutational)
colnames(kirc.data$mutational)

# Count number of mutations per gene
kirc.data$mutational %>%
  count(Gene) -> Total.Number.Mutations
Total.Number.Mutations

# Classify mutations into 3 groups based on IMPACT score
kirc.data$mutational %>% dplyr::count(kirc.data$mutational$Variant_Classification, kirc.data$mutational$IMPACT)

# High impact
kirc.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="HIGH") %>% 
  count(Gene) -> High.Impact.Mutations
High.Impact.Mutations

# Moderate impact
kirc.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODERATE") %>% 
  count(Gene) -> Moderate.Impact.Mutations
Moderate.Impact.Mutations

# Low impact
kirc.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="LOW") %>% 
  count(Gene) -> Low.Impact.Mutations
Low.Impact.Mutations

# Modifier impact
kirc.data$mutational %>%
  #group_by(Gene) %>%
  filter(IMPACT=="MODIFIER") %>% 
  count(Gene) -> Modifier.Mutations
Modifier.Mutations

# Filter for missense mutations, find number of missense mutations per gene
kirc.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  count(Gene) -> no.missense
no.missense

# Separate SIFT data into category and numerical
# SIFT category
library(qdapRegex)
rm_round(text.var =  kirc.data$mutational$SIFT, extract = F) -> kirc.data$mutational$SIFT.cat

kirc.data$mutational %>% 
  group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  dplyr::count(SIFT.cat) %>%
  select(Gene,SIFT.cat) -> sift.categories

# Count by category type
table(sift.categories) -> sift.categories.tab
as.data.frame.matrix(sift.categories.tab) -> sift.categories.df

# Rename column headings
colnames(sift.categories.df)[1]  <- "SIFT.Deleterious"
colnames(sift.categories.df)[2]  <- "SIFT.Deleterious_low_confidence"
colnames(sift.categories.df)[3]  <- "SIFT.Tolerated"
colnames(sift.categories.df)[4]  <- "SIFT.Tolerated_low_confidence"
sift.categories.df

# SIFT numerical values
unlist(rm_round(text.var =  kirc.data$mutational$SIFT, extract = T)) -> kirc.data$mutational$SIFT.num

kirc.data$mutational$SIFT.num <- as.numeric(kirc.data$mutational$SIFT.num)

kirc.data$mutational %>% 
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
rm_round(text.var =  kirc.data$mutational$PolyPhen, extract = F) -> kirc.data$mutational$PolyPhen.cat

kirc.data$mutational %>% 
  filter(Variant_Classification=="Missense_Mutation") %>%
  group_by(Gene) %>%
  dplyr::count(PolyPhen.cat) %>%
  select(Gene,PolyPhen.cat) -> PolyPhen.categories

PolyPhen.categories

# Count categories
table(PolyPhen.categories) -> PolyPhen.categories.tab
as.data.frame.matrix(PolyPhen.categories.tab) -> PolyPhen.categories.df

# Rename column headings
colnames(PolyPhen.categories.df)[1]  <- "PolyPhen.Benign"
colnames(PolyPhen.categories.df)[2]  <- "PolyPhen.Possibly.Damaging"
colnames(PolyPhen.categories.df)[3]  <- "PolyPhen.Probably.Damaging"
colnames(PolyPhen.categories.df)[4]  <- "PolyPhen.Unknown"
PolyPhen.categories.df

# PolyPhen numerical values
unlist(rm_round(text.var =  kirc.data$mutational$PolyPhen, extract = T)) -> kirc.data$mutational$PolyPhen.num

kirc.data$mutational$PolyPhen.num <- as.numeric(kirc.data$mutational$PolyPhen.num)

kirc.data$mutational %>% 
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

# Case_id - find number of patients with that gene mutated
kirc.data$mutational %>% count(Gene,case_id) -> Gene.Patient
table(Gene.Patient$Gene) -> Gene.Patient.tab
as.data.frame(Gene.Patient.tab) -> Patients.With.Mutation
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

# Retrieve GO terms for mutatated genes
retrieved.data <- annotate.HTseq.IDs(Total.Number.Mutations$Gene)
retrieved.data

# Rename column
colnames(retrieved.data)[2] <- "Gene"

# Select columns needed
retrieved.data %>% dplyr::select(Gene, go_id) -> GO

# Duplicate data frame
Gene.Ontology <- GO  

# Replace blanks with NA
Gene.Ontology[Gene.Ontology == ""] <- NA   

# Print 
Gene.Ontology

# Select columns needed
retrieved.data %>% dplyr::select(Gene, hgnc_symbol) -> hgnc

# Duplicate data frame
Hugo.Symbol <- hgnc

# Print
Hugo.Symbol

# Construct function for matching different columns of data to same genes and length
match.select <- function(x, y, label){
  left_join(x, y, by = "Gene") -> temp.df
  temp.df[,c(1,3)] -> temp.df
  colnames(temp.df)[2] <- label
  return(temp.df)
}

# Match all columns
Gene.Ontology.Match <- match.select(x = Total.Number.Mutations, y = Gene.Ontology, label = "Gene.Ontology")
Hugo.Symbol.Match <- match.select(x = Total.Number.Mutations, y = Hugo.Symbol, label = "Hugo.Symbol")
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
           Hugo.Symbol.Match[,2],
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
  full_join(., PolyPhen.categories.df.first.col, by = "Gene") -> KIRC.final
View(KIRC.final)

colnames(KIRC.final)[2] <- "Frequency.Gene.Is.Mutated"

# Percentages
install.packages("formattable")
library(formattable)
KIRC.final$Percent.Low.Impact.Mutations = percent(KIRC.final$Low.Impact.Mutations/KIRC.final$Frequency.Gene.Is.Mutated) -> Percent.Low.Impact.Mutations
KIRC.final$Percent.Moderate.Impact.Mutations = percent(KIRC.final$Moderate.Impact.Mutations/KIRC.final$Frequency.Gene.Is.Mutated) -> Percent.Moderate.Impact.Mutations
KIRC.final$Percent.High.Impact.Mutations = percent(KIRC.final$High.Impact.Mutations/KIRC.final$Frequency.Gene.Is.Mutated) -> Percent.High.Impact.Mutations
KIRC.final$Percent.Modifier.Mutations = percent(KIRC.final$Modifier.Mutations/KIRC.final$Frequency.Gene.Is.Mutated) -> Percent.Modifier.Mutations
KIRC.final$Percent.Missense.Mutation = percent(KIRC.final$Missense.Mutation/KIRC.final$Frequency.Gene.Is.Mutated) -> Percent.Missense.Mutation

# Percent of patients with that mutation
# Find total number of patients!
kirc.data$mutational %>% count(Gene,case_id) -> case.id.gene
data.frame(case.id.gene)
df_uniq <- unique(case.id.gene$case_id)
length(df_uniq)

# Calculate %
KIRC.final$Percent.Patients.With.Mutation = percent(KIRC.final$Number.Patients.With.Mutation/(length(df_uniq)))

# Annotate with GO terms using Biomart
library(biomaRt)

# Reorder columns
KIRC.final2 = KIRC.final %>% select(Gene, 
                                    Hugo.Symbol.Match...2.,
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


View(KIRC.final2)

write.csv(KIRC.final2, "~/tcga_biolinks1/KIRC.final2.csv", row.names=TRUE) 








###############################################################################

# Change expression ensembl gene names into correct format (no .number) to be compatible with mutation data
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)

# Construct query ensembl_gene_id_version -> ensembl_gene_id
Ensembl.Gene <- getBM(attributes=c('ensembl_gene_id', 'ensemble_gene_id_version'), filters='ensembl_gene_id', ensemblID, mart=mart)

getBM(attributes = c('ensembl_gene_id', 'ensembl_gene_id_version'),
      filters = 'ensembl_gene_id_version',
      values = brca5.network$Regulator, 
      mart = mart) -> Gene.Regulator

# Rename columns for matching
colnames(Gene.Regulator)[1] <- "Gene.Regulator"
colnames(Gene.Regulator)[2] <- "Regulator"
Gene.Regulator

getBM(attributes = c('ensembl_gene_id', 'ensembl_gene_id_version'),
      filters = 'ensembl_gene_id_version',
      values = brca5.network$Target, 
      mart = mart) -> Gene.Target

colnames(Gene.Target)[1] <- "Gene.Target"
colnames(Gene.Target)[2] <- "Target"
Gene.Target

# Construct function for matching different columns of data to same genes and length
match.select <- function(x, y, label){
  left_join(x, y, by = "Gene") -> temp.df
  temp.df[,c(1,3)] -> temp.df
  colnames(temp.df)[2] <- label
  return(temp.df)
}

# Match all columns
Gene.Regulator.Match <- match.select(x = brca5.network, y = Gene.Regulator, label = "Gene.Regulator")


data.frame(brca5.network,
           Gene.Target[,1],
           Gene.Regulator[,1]
) -> brca.df.test

###############################################################################

write.csv(BRCA5.network.final, "~/tcga_biolinks1/BRCA5.network.final.csv", row.names=TRUE)
dim(brca5.network)
dim(BRCA5.network.final)

x <- "~/tcga_biolinks1/ARACNe.data/BRCA5/brca5.network.txt"


rename.network.file <- function(x){
temp.network <- read.delim(file = x)
gsub("\\..*", '', temp.network$Regulator) -> temp.network$Regulator
gsub("\\..*", '', temp.network$Target) -> temp.network$Target
write.table(temp.network, file = gsub("network.txt", "network.rename.txt", x), sep = "\t")
}

rename.network.file(x)
list.files(path = "~/tcga_biolinks1/ARACNe.data/", pattern = "*.network.txt", recursive = T, full.names = T) -> files.to.convert

sapply(files.to.convert, rename.network.file)
