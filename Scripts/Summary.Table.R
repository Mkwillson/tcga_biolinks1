### you are going to need some kind of summary statistics for your mutational tables design a simple summary table and then produce
### some dplyr style code to construct a data.frame with the relevant info to match your table

#### BRCA ####
readRDS(file = "~/tcga_biolinks1/RDS/BRCA.final2") -> BRCA.final2
readRDS(file = "~/tcga_biolinks1/RDS/brca.data") -> brca.data

"BRCA" -> cancer.type.brca

length(BRCA.final2$Gene) -> Gene.num.brca
sum(!is.na(BRCA.final2$`Frequency Gene is Mutated`)) -> Gene.freq.brca
sum(!is.na(BRCA.final2$`Number of Patients With Mutation`)) -> Mut.freq.brca
sum(!is.na(BRCA.final2$`Number of High Impact Mutations`)) -> High.freq.brca
sum(!is.na(BRCA.final2$`Number of Low Impact Mutations`)) -> Low.freq.brca
sum(!is.na(BRCA.final2$`Number of Moderate Impact Mutations`)) -> Moderate.freq.brca
sum(!is.na(BRCA.final2$`Number of Modifier Mutations`)) -> Modifier.freq.brca
sum(!is.na(BRCA.final2$`Number of Patients with a Missense Mutation`)) -> Miss.freq.brca
length(which(BRCA.final2$`PolyPhen Probably Damaging` != 0)) -> Prob.dam.brca
length(which(BRCA.final2$`PolyPhen Possibly Damaging` != 0)) -> Poss.dam.brca
length(which(BRCA.final2$`SIFT Deleterious` != 0)) -> sift.del.brca
length(which(BRCA.final2$`SIFT Deleterious Low Confidence` !=0)) -> sift.del.lc.brca

brca.data$mutational %>% count(Gene,case_id) -> case.id.gene.brca
data.frame(case.id.gene.brca)
df_uniq.brca <- unique(case.id.gene.brca$case_id)
length(df_uniq.brca) -> Number.Patients.brca



df.sum.brca <- data.frame(cancer.type.brca,
                          Number.Patients.brca,
                          Gene.num.brca,
                          Gene.freq.brca,
                          Mut.freq.brca,
                          Low.freq.brca,
                          Moderate.freq.brca,
                          High.freq.brca,
                          Modifier.freq.brca,
                          Miss.freq.brca,
                          Prob.dam.brca,
                          Poss.dam.brca,
                          sift.del.brca,
                          sift.del.lc.brca)
View(df.sum.brca)
# Calculate %
df.sum.brca$Percent.Gene.freq.brca = percent(df.sum.brca$Gene.freq.brca/df.sum.brca$Gene.num.brca) -> Percent.Gene.freq.brca
df.sum.brca$Percent.Mut.freq.brca = percent(df.sum.brca$Mut.freq.brca/df.sum.brca$Gene.num.brca) -> Percent.Mut.freq.brca
df.sum.brca$Percent.Low.freq.brca = percent(df.sum.brca$Low.freq.brca/df.sum.brca$Gene.num.brca) -> Percent.Low.freq.brca
df.sum.brca$Percent.Moderate.freq.brca = percent(df.sum.brca$Moderate.freq.brca/df.sum.brca$Gene.num.brca) -> Percent.Moderate.freq.brca
df.sum.brca$Percent.High.freq.brca = percent(df.sum.brca$High.freq.brca/df.sum.brca$Gene.num.brca) -> Percent.High.freq.brca
df.sum.brca$Percent.Modifier.freq.brca = percent(df.sum.brca$Modifier.freq.brca/df.sum.brca$Gene.num.brca) -> Percent.Modifier.freq.brca
df.sum.brca$Percent.Miss.freq.brca = percent(df.sum.brca$Miss.freq.brca/df.sum.brca$Gene.num.brca) -> Percent.Miss.freq.brca
df.sum.brca$Percent.Prob.dam.brca = percent(df.sum.brca$Prob.dam.brca/df.sum.brca$Gene.num.brca) -> Percent.Prob.dam.brca
df.sum.brca$Percent.Poss.dam.brca = percent(df.sum.brca$Poss.dam.brca/df.sum.brca$Gene.num.brca) -> Percent.Poss.dam.brca
df.sum.brca$Percent.sift.del.brca = percent(df.sum.brca$sift.del.brca/df.sum.brca$Gene.num.brca) -> Percent.sift.del.brca
df.sum.brca$Percent.sift.del.lc.brca = percent(df.sum.brca$sift.del.lc.brca/df.sum.brca$Gene.num.brca) -> Percent.sift.del.lc.brca

df.sum.brca2 <- df.sum.brca %>% dplyr::select(cancer.type.brca,
                                              Number.Patients.brca,
                                              Gene.num.brca,
                                              Percent.Mut.freq.brca,
                                              Percent.Miss.freq.brca,
                                              Percent.Prob.dam.brca,
                                              Percent.Poss.dam.brca,
                                              Percent.sift.del.brca,
                                              Percent.sift.del.lc.brca
                                              )
View(df.sum.brca2)

colnames(df.sum.brca2)[1] <- "Cancer Type"
colnames(df.sum.brca2)[2] <- "Number of Patients"
colnames(df.sum.brca2)[3] <- "Number of Genes"
colnames(df.sum.brca2)[4] <- "Percent of Genes that are Mutated in at least 1 Patient"
colnames(df.sum.brca2)[5] <- "Percent of Genes that are Missense Mutated in at least 1 Patient"
colnames(df.sum.brca2)[6] <- "Percent of Genes that are Probably Damaging"
colnames(df.sum.brca2)[7] <- "Percent of Genes that are Possibly Damaging"
colnames(df.sum.brca2)[8] <- "Percent of Genes that are Deleterious"
colnames(df.sum.brca2)[9] <- "Percent of Genes that are Deleterious (Low Confidence)"

View(df.sum.brca2)
saveRDS(df.sum.brca2, file = "~/tcga_biolinks1/RDS/df.sum.brca2")

#### PRAD ####
"PRAD" -> cancer.type.prad

readRDS(file = "~/tcga_biolinks1/RDS/PRAD.final2") -> PRAD.final2
readRDS(file = "~/tcga_biolinks1/RDS/prad.data") -> prad.data

length(PRAD.final2$Gene) -> Gene.num.prad
sum(!is.na(PRAD.final2$`Frequency Gene is Mutated`)) -> Gene.freq.prad
sum(!is.na(PRAD.final2$`Number of Patients With Mutation`)) -> Mut.freq.prad
sum(!is.na(PRAD.final2$`Number of High Impact Mutations`)) -> High.freq.prad
sum(!is.na(PRAD.final2$`Number of Low Impact Mutations`)) -> Low.freq.prad
sum(!is.na(PRAD.final2$`Number of Moderate Impact Mutations`)) -> Moderate.freq.prad
sum(!is.na(PRAD.final2$`Number of Modifier Mutations`)) -> Modifier.freq.prad
sum(!is.na(PRAD.final2$`Number of Patients with a Missense Mutation`)) -> Miss.freq.prad
length(which(PRAD.final2$`PolyPhen Probably Damaging` != 0)) -> Prob.dam.prad
length(which(PRAD.final2$`PolyPhen Possibly Damaging` != 0)) -> Poss.dam.prad
length(which(PRAD.final2$`SIFT Deleterious` != 0)) -> sift.del.prad
length(which(PRAD.final2$`SIFT Deleterious Low Confidence` !=0)) -> sift.del.lc.prad

prad.data$mutational %>% count(Gene,case_id) -> case.id.gene.prad
data.frame(case.id.gene.prad)
df_uniq.prad <- unique(case.id.gene.prad$case_id)
length(df_uniq.prad) -> Number.Patients.prad



df.sum.prad <- data.frame(cancer.type.prad,
                          Number.Patients.prad,
                          Gene.num.prad,
                          Gene.freq.prad,
                          Mut.freq.prad,
                          Low.freq.prad,
                          Moderate.freq.prad,
                          High.freq.prad,
                          Modifier.freq.prad,
                          Miss.freq.prad,
                          Prob.dam.prad,
                          Poss.dam.prad,
                          sift.del.prad,
                          sift.del.lc.prad)
View(df.sum.prad)
# Calculate %
df.sum.prad$Percent.Gene.freq.prad = percent(df.sum.prad$Gene.freq.prad/df.sum.prad$Gene.num.prad) -> Percent.Gene.freq.prad
df.sum.prad$Percent.Mut.freq.prad = percent(df.sum.prad$Mut.freq.prad/df.sum.prad$Gene.num.prad) -> Percent.Mut.freq.prad
df.sum.prad$Percent.Low.freq.prad = percent(df.sum.prad$Low.freq.prad/df.sum.prad$Gene.num.prad) -> Percent.Low.freq.prad
df.sum.prad$Percent.Moderate.freq.prad = percent(df.sum.prad$Moderate.freq.prad/df.sum.prad$Gene.num.prad) -> Percent.Moderate.freq.prad
df.sum.prad$Percent.High.freq.prad = percent(df.sum.prad$High.freq.prad/df.sum.prad$Gene.num.prad) -> Percent.High.freq.prad
df.sum.prad$Percent.Modifier.freq.prad = percent(df.sum.prad$Modifier.freq.prad/df.sum.prad$Gene.num.prad) -> Percent.Modifier.freq.prad
df.sum.prad$Percent.Miss.freq.prad = percent(df.sum.prad$Miss.freq.prad/df.sum.prad$Gene.num.prad) -> Percent.Miss.freq.prad
df.sum.prad$Percent.Prob.dam.prad = percent(df.sum.prad$Prob.dam.prad/df.sum.prad$Gene.num.prad) -> Percent.Prob.dam.prad
df.sum.prad$Percent.Poss.dam.prad = percent(df.sum.prad$Poss.dam.prad/df.sum.prad$Gene.num.prad) -> Percent.Poss.dam.prad
df.sum.prad$Percent.sift.del.prad = percent(df.sum.prad$sift.del.prad/df.sum.prad$Gene.num.prad) -> Percent.sift.del.prad
df.sum.prad$Percent.sift.del.lc.prad = percent(df.sum.prad$sift.del.lc.prad/df.sum.prad$Gene.num.prad) -> Percent.sift.del.lc.prad

df.sum.prad2 <- df.sum.prad %>% dplyr::select(cancer.type.prad,
                                              Number.Patients.prad,
                                              Gene.num.prad,
                                              Percent.Mut.freq.prad,
                                              Percent.Miss.freq.prad,
                                              Percent.Prob.dam.prad,
                                              Percent.Poss.dam.prad,
                                              Percent.sift.del.prad,
                                              Percent.sift.del.lc.prad)
View(df.sum.prad2)

colnames(df.sum.prad2)[1] <- "Cancer Type"
colnames(df.sum.prad2)[2] <- "Number of Patients"
colnames(df.sum.prad2)[3] <- "Number of Genes"
colnames(df.sum.prad2)[4] <- "Percent of Genes that are Mutated in at least 1 Patient"
colnames(df.sum.prad2)[5] <- "Percent of Genes that are Missense Mutated in at least 1 Patient"
colnames(df.sum.prad2)[6] <- "Percent of Genes that are Probably Damaging"
colnames(df.sum.prad2)[7] <- "Percent of Genes that are Possibly Damaging"
colnames(df.sum.prad2)[8] <- "Percent of Genes that are Deleterious"
colnames(df.sum.prad2)[9] <- "Percent of Genes that are Deleterious (Low Confidence)"

View(df.sum.prad2)
saveRDS(df.sum.prad2, file = "~/tcga_biolinks1/RDS/df.sum.prad2")

rbind(df.sum.brca2, df.sum.prad2) ->test


identical(names(df.sum.brca2[[1]]), names(df.sum.prad2[[2]]) )
