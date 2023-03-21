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

brca.data$mutational %>% dplyr::count(Gene,case_id) -> case.id.gene.brca
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

prad.data$mutational %>% dplyr::count(Gene,case_id) -> case.id.gene.prad
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

################################################################################

#### LUAD ####
readRDS(file = "~/tcga_biolinks1/RDS/LUAD.final2") -> LUAD.final2
readRDS(file = "~/tcga_biolinks1/RDS/luad.data") -> luad.data

"LUAD" -> cancer.type.luad

length(LUAD.final2$Gene) -> Gene.num.luad
sum(!is.na(LUAD.final2$`Frequency Gene is Mutated`)) -> Gene.freq.luad
sum(!is.na(LUAD.final2$`Number of Patients With Mutation`)) -> Mut.freq.luad
sum(!is.na(LUAD.final2$`Number of High Impact Mutations`)) -> High.freq.luad
sum(!is.na(LUAD.final2$`Number of Low Impact Mutations`)) -> Low.freq.luad
sum(!is.na(LUAD.final2$`Number of Moderate Impact Mutations`)) -> Moderate.freq.luad
sum(!is.na(LUAD.final2$`Number of Modifier Mutations`)) -> Modifier.freq.luad
sum(!is.na(LUAD.final2$`Number of Patients With a Missense Mutation`)) -> Miss.freq.luad
length(which(LUAD.final2$`PolyPhen Probably Damaging` != 0)) -> Prob.dam.luad
length(which(LUAD.final2$`PolyPhen Possibly Damaging` != 0)) -> Poss.dam.luad
length(which(LUAD.final2$`SIFT Deleterious` != 0)) -> sift.del.luad
length(which(LUAD.final2$`SIFT Deleterious Low Confidence` !=0)) -> sift.del.lc.luad
LUAD.final2$number
luad.data$mutational %>% dplyr::count(Gene,case_id) -> case.id.gene.luad
data.frame(case.id.gene.luad)
df_uniq.luad <- unique(case.id.gene.luad$case_id)
length(df_uniq.luad) -> Number.Patients.luad



df.sum.luad <- data.frame(cancer.type.luad,
                          Number.Patients.luad,
                          Gene.num.luad,
                          Gene.freq.luad,
                          Mut.freq.luad,
                          Low.freq.luad,
                          Moderate.freq.luad,
                          High.freq.luad,
                          Modifier.freq.luad,
                          Miss.freq.luad,
                          Prob.dam.luad,
                          Poss.dam.luad,
                          sift.del.luad,
                          sift.del.lc.luad)
View(df.sum.luad)
# Calculate %
df.sum.luad$Percent.Gene.freq.luad = percent(df.sum.luad$Gene.freq.luad/df.sum.luad$Gene.num.luad) -> Percent.Gene.freq.luad
df.sum.luad$Percent.Mut.freq.luad = percent(df.sum.luad$Mut.freq.luad/df.sum.luad$Gene.num.luad) -> Percent.Mut.freq.luad
df.sum.luad$Percent.Low.freq.luad = percent(df.sum.luad$Low.freq.luad/df.sum.luad$Gene.num.luad) -> Percent.Low.freq.luad
df.sum.luad$Percent.Moderate.freq.luad = percent(df.sum.luad$Moderate.freq.luad/df.sum.luad$Gene.num.luad) -> Percent.Moderate.freq.luad
df.sum.luad$Percent.High.freq.luad = percent(df.sum.luad$High.freq.luad/df.sum.luad$Gene.num.luad) -> Percent.High.freq.luad
df.sum.luad$Percent.Modifier.freq.luad = percent(df.sum.luad$Modifier.freq.luad/df.sum.luad$Gene.num.luad) -> Percent.Modifier.freq.luad
df.sum.luad$Percent.Miss.freq.luad = percent(df.sum.luad$Miss.freq.luad/df.sum.luad$Gene.num.luad) -> Percent.Miss.freq.luad
df.sum.luad$Percent.Prob.dam.luad = percent(df.sum.luad$Prob.dam.luad/df.sum.luad$Gene.num.luad) -> Percent.Prob.dam.luad
df.sum.luad$Percent.Poss.dam.luad = percent(df.sum.luad$Poss.dam.luad/df.sum.luad$Gene.num.luad) -> Percent.Poss.dam.luad
df.sum.luad$Percent.sift.del.luad = percent(df.sum.luad$sift.del.luad/df.sum.luad$Gene.num.luad) -> Percent.sift.del.luad
df.sum.luad$Percent.sift.del.lc.luad = percent(df.sum.luad$sift.del.lc.luad/df.sum.luad$Gene.num.luad) -> Percent.sift.del.lc.luad

df.sum.luad2 <- df.sum.luad %>% dplyr::select(cancer.type.luad,
                                              Number.Patients.luad,
                                              Gene.num.luad,
                                              Percent.Mut.freq.luad,
                                              Percent.Miss.freq.luad,
                                              Percent.Prob.dam.luad,
                                              Percent.Poss.dam.luad,
                                              Percent.sift.del.luad,
                                              Percent.sift.del.lc.luad
)
View(df.sum.luad2)

colnames(df.sum.luad2)[1] <- "Cancer Type"
colnames(df.sum.luad2)[2] <- "Number of Patients"
colnames(df.sum.luad2)[3] <- "Number of Genes"
colnames(df.sum.luad2)[4] <- "Percent of Genes that are Mutated in at least 1 Patient"
colnames(df.sum.luad2)[5] <- "Percent of Genes that are Missense Mutated in at least 1 Patient"
colnames(df.sum.luad2)[6] <- "Percent of Genes that are Probably Damaging"
colnames(df.sum.luad2)[7] <- "Percent of Genes that are Possibly Damaging"
colnames(df.sum.luad2)[8] <- "Percent of Genes that are Deleterious"
colnames(df.sum.luad2)[9] <- "Percent of Genes that are Deleterious (Low Confidence)"

View(df.sum.luad2)
saveRDS(df.sum.luad2, file = "~/tcga_biolinks1/RDS/df.sum.luad2")

#### STAD ######################################################################

readRDS(file = "~/tcga_biolinks1/RDS/STAD.final2") -> STAD.final2
readRDS(file = "~/tcga_biolinks1/RDS/stad.data") -> stad.data

"STAD" -> cancer.type.stad

length(STAD.final2$Gene) -> Gene.num.stad
sum(!is.na(STAD.final2$`Frequency Gene is Mutated`)) -> Gene.freq.stad
sum(!is.na(STAD.final2$`Number of Patients With Mutation`)) -> Mut.freq.stad
sum(!is.na(STAD.final2$`Number of High Impact Mutations`)) -> High.freq.stad
sum(!is.na(STAD.final2$`Number of Low Impact Mutations`)) -> Low.freq.stad
sum(!is.na(STAD.final2$`Number of Moderate Impact Mutations`)) -> Moderate.freq.stad
sum(!is.na(STAD.final2$`Number of Modifier Mutations`)) -> Modifier.freq.stad
sum(!is.na(STAD.final2$`Number of Patients with a Missense Mutation`)) -> Miss.freq.stad
length(which(STAD.final2$`PolyPhen Probably Damaging` != 0)) -> Prob.dam.stad
length(which(STAD.final2$`PolyPhen Possibly Damaging` != 0)) -> Poss.dam.stad
length(which(STAD.final2$`SIFT Deleterious` != 0)) -> sift.del.stad
length(which(STAD.final2$`SIFT Deleterious Low Confidence` !=0)) -> sift.del.lc.stad

stad.data$mutational %>% dplyr::count(Gene,case_id) -> case.id.gene.stad
data.frame(case.id.gene.stad)
df_uniq.stad <- unique(case.id.gene.stad$case_id)
length(df_uniq.stad) -> Number.Patients.stad



df.sum.stad <- data.frame(cancer.type.stad,
                          Number.Patients.stad,
                          Gene.num.stad,
                          Gene.freq.stad,
                          Mut.freq.stad,
                          Low.freq.stad,
                          Moderate.freq.stad,
                          High.freq.stad,
                          Modifier.freq.stad,
                          Miss.freq.stad,
                          Prob.dam.stad,
                          Poss.dam.stad,
                          sift.del.stad,
                          sift.del.lc.stad)
View(df.sum.stad)
# Calculate %
df.sum.stad$Percent.Gene.freq.stad = percent(df.sum.stad$Gene.freq.stad/df.sum.stad$Gene.num.stad) -> Percent.Gene.freq.stad
df.sum.stad$Percent.Mut.freq.stad = percent(df.sum.stad$Mut.freq.stad/df.sum.stad$Gene.num.stad) -> Percent.Mut.freq.stad
df.sum.stad$Percent.Low.freq.stad = percent(df.sum.stad$Low.freq.stad/df.sum.stad$Gene.num.stad) -> Percent.Low.freq.stad
df.sum.stad$Percent.Moderate.freq.stad = percent(df.sum.stad$Moderate.freq.stad/df.sum.stad$Gene.num.stad) -> Percent.Moderate.freq.stad
df.sum.stad$Percent.High.freq.stad = percent(df.sum.stad$High.freq.stad/df.sum.stad$Gene.num.stad) -> Percent.High.freq.stad
df.sum.stad$Percent.Modifier.freq.stad = percent(df.sum.stad$Modifier.freq.stad/df.sum.stad$Gene.num.stad) -> Percent.Modifier.freq.stad
df.sum.stad$Percent.Miss.freq.stad = percent(df.sum.stad$Miss.freq.stad/df.sum.stad$Gene.num.stad) -> Percent.Miss.freq.stad
df.sum.stad$Percent.Prob.dam.stad = percent(df.sum.stad$Prob.dam.stad/df.sum.stad$Gene.num.stad) -> Percent.Prob.dam.stad
df.sum.stad$Percent.Poss.dam.stad = percent(df.sum.stad$Poss.dam.stad/df.sum.stad$Gene.num.stad) -> Percent.Poss.dam.stad
df.sum.stad$Percent.sift.del.stad = percent(df.sum.stad$sift.del.stad/df.sum.stad$Gene.num.stad) -> Percent.sift.del.stad
df.sum.stad$Percent.sift.del.lc.stad = percent(df.sum.stad$sift.del.lc.stad/df.sum.stad$Gene.num.stad) -> Percent.sift.del.lc.stad

df.sum.stad2 <- df.sum.stad %>% dplyr::select(cancer.type.stad,
                                              Number.Patients.stad,
                                              Gene.num.stad,
                                              Percent.Mut.freq.stad,
                                              Percent.Miss.freq.stad,
                                              Percent.Prob.dam.stad,
                                              Percent.Poss.dam.stad,
                                              Percent.sift.del.stad,
                                              Percent.sift.del.lc.stad
)
View(df.sum.stad2)

colnames(df.sum.stad2)[1] <- "Cancer Type"
colnames(df.sum.stad2)[2] <- "Number of Patients"
colnames(df.sum.stad2)[3] <- "Number of Genes"
colnames(df.sum.stad2)[4] <- "Percent of Genes that are Mutated in at least 1 Patient"
colnames(df.sum.stad2)[5] <- "Percent of Genes that are Missense Mutated in at least 1 Patient"
colnames(df.sum.stad2)[6] <- "Percent of Genes that are Probably Damaging"
colnames(df.sum.stad2)[7] <- "Percent of Genes that are Possibly Damaging"
colnames(df.sum.stad2)[8] <- "Percent of Genes that are Deleterious"
colnames(df.sum.stad2)[9] <- "Percent of Genes that are Deleterious (Low Confidence)"

View(df.sum.stad2)
saveRDS(df.sum.stad2, file = "~/tcga_biolinks1/RDS/df.sum.stad2")

#### KIRP ######################################################################
readRDS(file = "~/tcga_biolinks1/RDS/KIRP.final2") -> KIRP.final2
readRDS(file = "~/tcga_biolinks1/RDS/kirp.data") -> kirp.data

"KIRP" -> cancer.type.kirp

length(KIRP.final2$Gene) -> Gene.num.kirp
sum(!is.na(KIRP.final2$`Frequency Gene is Mutated`)) -> Gene.freq.kirp
sum(!is.na(KIRP.final2$`Number of Patients With Mutation`)) -> Mut.freq.kirp
sum(!is.na(KIRP.final2$`Number of High Impact Mutations`)) -> High.freq.kirp
sum(!is.na(KIRP.final2$`Number of Low Impact Mutations`)) -> Low.freq.kirp
sum(!is.na(KIRP.final2$`Number of Moderate Impact Mutations`)) -> Moderate.freq.kirp
sum(!is.na(KIRP.final2$`Number of Modifier Mutations`)) -> Modifier.freq.kirp
sum(!is.na(KIRP.final2$`Number of Patients with a Missense Mutation`)) -> Miss.freq.kirp
length(which(KIRP.final2$`PolyPhen Probably Damaging` != 0)) -> Prob.dam.kirp
length(which(KIRP.final2$`PolyPhen Possibly Damaging` != 0)) -> Poss.dam.kirp
length(which(KIRP.final2$`SIFT Deleterious` != 0)) -> sift.del.kirp
length(which(KIRP.final2$`SIFT Deleterious Low Confidence` !=0)) -> sift.del.lc.kirp

kirp.data$mutational %>% dplyr::count(Gene,case_id) -> case.id.gene.kirp
data.frame(case.id.gene.kirp)
df_uniq.kirp <- unique(case.id.gene.kirp$case_id)
length(df_uniq.kirp) -> Number.Patients.kirp



df.sum.kirp <- data.frame(cancer.type.kirp,
                          Number.Patients.kirp,
                          Gene.num.kirp,
                          Gene.freq.kirp,
                          Mut.freq.kirp,
                          Low.freq.kirp,
                          Moderate.freq.kirp,
                          High.freq.kirp,
                          Modifier.freq.kirp,
                          Miss.freq.kirp,
                          Prob.dam.kirp,
                          Poss.dam.kirp,
                          sift.del.kirp,
                          sift.del.lc.kirp)
View(df.sum.kirp)
# Calculate %
df.sum.kirp$Percent.Gene.freq.kirp = percent(df.sum.kirp$Gene.freq.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.Gene.freq.kirp
df.sum.kirp$Percent.Mut.freq.kirp = percent(df.sum.kirp$Mut.freq.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.Mut.freq.kirp
df.sum.kirp$Percent.Low.freq.kirp = percent(df.sum.kirp$Low.freq.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.Low.freq.kirp
df.sum.kirp$Percent.Moderate.freq.kirp = percent(df.sum.kirp$Moderate.freq.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.Moderate.freq.kirp
df.sum.kirp$Percent.High.freq.kirp = percent(df.sum.kirp$High.freq.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.High.freq.kirp
df.sum.kirp$Percent.Modifier.freq.kirp = percent(df.sum.kirp$Modifier.freq.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.Modifier.freq.kirp
df.sum.kirp$Percent.Miss.freq.kirp = percent(df.sum.kirp$Miss.freq.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.Miss.freq.kirp
df.sum.kirp$Percent.Prob.dam.kirp = percent(df.sum.kirp$Prob.dam.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.Prob.dam.kirp
df.sum.kirp$Percent.Poss.dam.kirp = percent(df.sum.kirp$Poss.dam.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.Poss.dam.kirp
df.sum.kirp$Percent.sift.del.kirp = percent(df.sum.kirp$sift.del.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.sift.del.kirp
df.sum.kirp$Percent.sift.del.lc.kirp = percent(df.sum.kirp$sift.del.lc.kirp/df.sum.kirp$Gene.num.kirp) -> Percent.sift.del.lc.kirp

df.sum.kirp2 <- df.sum.kirp %>% dplyr::select(cancer.type.kirp,
                                              Number.Patients.kirp,
                                              Gene.num.kirp,
                                              Percent.Mut.freq.kirp,
                                              Percent.Miss.freq.kirp,
                                              Percent.Prob.dam.kirp,
                                              Percent.Poss.dam.kirp,
                                              Percent.sift.del.kirp,
                                              Percent.sift.del.lc.kirp
)
View(df.sum.kirp2)

colnames(df.sum.kirp2)[1] <- "Cancer Type"
colnames(df.sum.kirp2)[2] <- "Number of Patients"
colnames(df.sum.kirp2)[3] <- "Number of Genes"
colnames(df.sum.kirp2)[4] <- "Percent of Genes that are Mutated in at least 1 Patient"
colnames(df.sum.kirp2)[5] <- "Percent of Genes that are Missense Mutated in at least 1 Patient"
colnames(df.sum.kirp2)[6] <- "Percent of Genes that are Probably Damaging"
colnames(df.sum.kirp2)[7] <- "Percent of Genes that are Possibly Damaging"
colnames(df.sum.kirp2)[8] <- "Percent of Genes that are Deleterious"
colnames(df.sum.kirp2)[9] <- "Percent of Genes that are Deleterious (Low Confidence)"

View(df.sum.kirp2)
saveRDS(df.sum.kirp2, file = "~/tcga_biolinks1/RDS/df.sum.kirp2")

#### KIRC #######################################################################


readRDS(file = "~/tcga_biolinks1/RDS/KIRC.final2") -> KIRC.final2
readRDS(file = "~/tcga_biolinks1/RDS/kirc.data") -> kirc.data

"KIRC" -> cancer.type.kirc

length(KIRC.final2$Gene) -> Gene.num.kirc
sum(!is.na(KIRC.final2$`Frequency Gene is Mutated`)) -> Gene.freq.kirc
sum(!is.na(KIRC.final2$`Number of Patients With Mutation`)) -> Mut.freq.kirc
sum(!is.na(KIRC.final2$`Number of High Impact Mutations`)) -> High.freq.kirc
sum(!is.na(KIRC.final2$`Number of Low Impact Mutations`)) -> Low.freq.kirc
sum(!is.na(KIRC.final2$`Number of Moderate Impact Mutations`)) -> Moderate.freq.kirc
sum(!is.na(KIRC.final2$`Number of Modifier Mutations`)) -> Modifier.freq.kirc
sum(!is.na(KIRC.final2$`Number of Patients with a Missense Mutation`)) -> Miss.freq.kirc
length(which(KIRC.final2$`PolyPhen Probably Damaging` != 0)) -> Prob.dam.kirc
length(which(KIRC.final2$`PolyPhen Possibly Damaging` != 0)) -> Poss.dam.kirc
length(which(KIRC.final2$`SIFT Deleterious` != 0)) -> sift.del.kirc
length(which(KIRC.final2$`SIFT Deleterious Low Confidence` !=0)) -> sift.del.lc.kirc

kirc.data$mutational %>% dplyr::count(Gene,case_id) -> case.id.gene.kirc
data.frame(case.id.gene.kirc)
df_uniq.kirc <- unique(case.id.gene.kirc$case_id)
length(df_uniq.kirc) -> Number.Patients.kirc



df.sum.kirc <- data.frame(cancer.type.kirc,
                          Number.Patients.kirc,
                          Gene.num.kirc,
                          Gene.freq.kirc,
                          Mut.freq.kirc,
                          Low.freq.kirc,
                          Moderate.freq.kirc,
                          High.freq.kirc,
                          Modifier.freq.kirc,
                          Miss.freq.kirc,
                          Prob.dam.kirc,
                          Poss.dam.kirc,
                          sift.del.kirc,
                          sift.del.lc.kirc)
View(df.sum.kirc)
# Calculate %
df.sum.kirc$Percent.Gene.freq.kirc = percent(df.sum.kirc$Gene.freq.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.Gene.freq.kirc
df.sum.kirc$Percent.Mut.freq.kirc = percent(df.sum.kirc$Mut.freq.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.Mut.freq.kirc
df.sum.kirc$Percent.Low.freq.kirc = percent(df.sum.kirc$Low.freq.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.Low.freq.kirc
df.sum.kirc$Percent.Moderate.freq.kirc = percent(df.sum.kirc$Moderate.freq.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.Moderate.freq.kirc
df.sum.kirc$Percent.High.freq.kirc = percent(df.sum.kirc$High.freq.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.High.freq.kirc
df.sum.kirc$Percent.Modifier.freq.kirc = percent(df.sum.kirc$Modifier.freq.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.Modifier.freq.kirc
df.sum.kirc$Percent.Miss.freq.kirc = percent(df.sum.kirc$Miss.freq.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.Miss.freq.kirc
df.sum.kirc$Percent.Prob.dam.kirc = percent(df.sum.kirc$Prob.dam.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.Prob.dam.kirc
df.sum.kirc$Percent.Poss.dam.kirc = percent(df.sum.kirc$Poss.dam.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.Poss.dam.kirc
df.sum.kirc$Percent.sift.del.kirc = percent(df.sum.kirc$sift.del.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.sift.del.kirc
df.sum.kirc$Percent.sift.del.lc.kirc = percent(df.sum.kirc$sift.del.lc.kirc/df.sum.kirc$Gene.num.kirc) -> Percent.sift.del.lc.kirc

df.sum.kirc2 <- df.sum.kirc %>% dplyr::select(cancer.type.kirc,
                                              Number.Patients.kirc,
                                              Gene.num.kirc,
                                              Percent.Mut.freq.kirc,
                                              Percent.Miss.freq.kirc,
                                              Percent.Prob.dam.kirc,
                                              Percent.Poss.dam.kirc,
                                              Percent.sift.del.kirc,
                                              Percent.sift.del.lc.kirc
)
View(df.sum.kirc2)

colnames(df.sum.kirc2)[1] <- "Cancer Type"
colnames(df.sum.kirc2)[2] <- "Number of Patients"
colnames(df.sum.kirc2)[3] <- "Number of Genes"
colnames(df.sum.kirc2)[4] <- "Percent of Genes that are Mutated in at least 1 Patient"
colnames(df.sum.kirc2)[5] <- "Percent of Genes that are Missense Mutated in at least 1 Patient"
colnames(df.sum.kirc2)[6] <- "Percent of Genes that are Probably Damaging"
colnames(df.sum.kirc2)[7] <- "Percent of Genes that are Possibly Damaging"
colnames(df.sum.kirc2)[8] <- "Percent of Genes that are Deleterious"
colnames(df.sum.kirc2)[9] <- "Percent of Genes that are Deleterious (Low Confidence)"

View(df.sum.kirc2)
saveRDS(df.sum.kirc2, file = "~/tcga_biolinks1/RDS/df.sum.kirc2")


rbind(df.sum.brca2, df.sum.kirc2, df.sum.kirp2, df.sum.luad2, df.sum.prad2, df.sum.stad2) -> summary.table.mut
View(summary.table.mut)

colnames(summary.table.mut)[6] <- "Percent of Missense Mutated Genes that are Probably Damaging"
colnames(summary.table.mut)[7] <- "Percent of Missense Mutated Genes that are Possibly Damaging"
colnames(summary.table.mut)[8] <- "Percent of Missense Mutated Genes that are Deleterious"
colnames(summary.table.mut)[9] <- "Percent of Missense Mutated Genes that are Deleterious (Low Confidence)"
saveRDS(summary.table.mut, file = "~/tcga_biolinks1/RDS/summary.table.mut")
View(summary.table.mut)

write.csv(summary.table.mut, file = "~/tcga_biolinks1/Mutations/summary.table.mut.csv")
