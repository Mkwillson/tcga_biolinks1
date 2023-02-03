# show the names of the assays
names(brca.data$expression@assays)
# extract the BRCA expression data specifically TPM (Transcripts Per Million) into a standard dataframe
assay(brca.data$expression, "tpm_unstrand") -> exp.df.tpm.brca
dim(exp.df.tpm.brca)
head(exp.df.tpm.brca[,1:20])

sd.limit = 10
dim(exp.df.tpm.brca[rowSds(exp.df.tpm.brca)>sd.limit,])
exp.df.tpm.brca[rowSds(exp.df.tpm.brca)>sd.limit,] -> filt.exp.df.tpm.brca


data.frame(gene = rownames(filt.exp.df.tpm.brca),filt.exp.df.tpm.brca) -> filt.exp.df.tpm.brca.first.col
head(filt.exp.df.tpm.brca.first.col[1:20,1:20])

write.table(filt.exp.df.tpm.brca.first.col, file = "./filt.exp.df.tpm.brca.txt", sep = "\t", quote = F, row.names = F)

# extract the luad expression data specifically TPM (Transcripts Per Million) into a standard dataframe
assay(prad.data$expression, "tpm_unstrand") -> exp.df.tpm.prad
dim(exp.df.tpm.prad)
head(exp.df.tpm.prad[,1:20])

sd.limit = 10
dim(exp.df.tpm.prad[rowSds(exp.df.tpm.prad)>sd.limit,])
exp.df.tpm.prad[rowSds(exp.df.tpm.prad)>sd.limit,] -> filt.exp.df.tpm.prad


data.frame(gene = rownames(filt.exp.df.tpm.prad),filt.exp.df.tpm.prad) -> filt.exp.df.tpm.prad.first.col
head(filt.exp.df.tpm.prad.first.col[1:20,1:20])

write.table(filt.exp.df.tpm.prad.first.col, file = "./filt.exp.df.tpm.prad10.txt", sep = "\t", quote = F, row.names = F)

# extract the LUAD expression data specifically TPM (Transcripts Per Million) into a standard dataframe
assay(luad.data$expression, "tpm_unstrand") -> exp.df.tpm.luad
dim(exp.df.tpm.luad)
head(exp.df.tpm.luad[,1:20])

sd.limit = 50
dim(exp.df.tpm.luad[rowSds(exp.df.tpm.luad)>sd.limit,])
exp.df.tpm.luad[rowSds(exp.df.tpm.luad)>sd.limit,] -> filt.exp.df.tpm.luad


data.frame(gene = rownames(filt.exp.df.tpm.luad),filt.exp.df.tpm.luad) -> filt.exp.df.tpm.luad.first.col
head(filt.exp.df.tpm.luad.first.col[1:20,1:20])

write.table(filt.exp.df.tpm.luad.first.col, file = "./filt.exp.df.tpm.luad.txt", sep = "\t", quote = F, row.names = F)

# extract the STAD expression data specifically TPM (Transcripts Per Million) into a standard dataframe
assay(stad.data$expression, "tpm_unstrand") -> exp.df.tpm.stad
dim(exp.df.tpm.stad)
head(exp.df.tpm.stad[,1:20])

sd.limit = 50
dim(exp.df.tpm.stad[rowSds(exp.df.tpm.stad)>sd.limit,])
exp.df.tpm.stad[rowSds(exp.df.tpm.stad)>sd.limit,] -> filt.exp.df.tpm.stad


data.frame(gene = rownames(filt.exp.df.tpm.stad),filt.exp.df.tpm.stad) -> filt.exp.df.tpm.stad.first.col
head(filt.exp.df.tpm.stad.first.col[1:20,1:20])

write.table(filt.exp.df.tpm.stad.first.col, file = "./filt.exp.df.tpm.stad.txt", sep = "\t", quote = F, row.names = F)

# extract the kirp expression data specifically TPM (Transcripts Per Million) into a standard dataframe
assay(kirp.data$expression, "tpm_unstrand") -> exp.df.tpm.kirp
dim(exp.df.tpm.kirp)
head(exp.df.tpm.kirp[,1:20])

sd.limit = 50
dim(exp.df.tpm.kirp[rowSds(exp.df.tpm.kirp)>sd.limit,])
exp.df.tpm.kirp[rowSds(exp.df.tpm.kirp)>sd.limit,] -> filt.exp.df.tpm.kirp


data.frame(gene = rownames(filt.exp.df.tpm.kirp),filt.exp.df.tpm.kirp) -> filt.exp.df.tpm.kirp.first.col
head(filt.exp.df.tpm.kirp.first.col[1:20,1:20])

write.table(filt.exp.df.tpm.kirp.first.col, file = "./filt.exp.df.tpm.kirp.txt", sep = "\t", quote = F, row.names = F)

# extract the kirc expression data specifically TPM (Transcripts Per Million) into a standard dataframe
assay(kirc.data$expression, "tpm_unstrand") -> exp.df.tpm.kirc
dim(exp.df.tpm.kirc)
head(exp.df.tpm.kirc[,1:20])

sd.limit = 50
dim(exp.df.tpm.kirc[rowSds(exp.df.tpm.kirc)>sd.limit,])
exp.df.tpm.kirc[rowSds(exp.df.tpm.kirc)>sd.limit,] -> filt.exp.df.tpm.kirc


data.frame(gene = rownames(filt.exp.df.tpm.kirc),filt.exp.df.tpm.kirc) -> filt.exp.df.tpm.kirc.first.col
head(filt.exp.df.tpm.kirc.first.col[1:20,1:20])

write.table(filt.exp.df.tpm.kirc.first.col, file = "./filt.exp.df.tpm.kirc.txt", sep = "\t", quote = F, row.names = F)
