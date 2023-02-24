# show the names of the assays
names(brca.data$expression@assays)
# extract the expression data specifically TPM (Transcripts Per Million) into a standard dataframe
assay(brca.data$expression, "tpm_unstrand") -> exp.df.tpm.brca
dim(exp.df.tpm.brca)
head(exp.df.tpm.brca[,1:20])

sd.limit = 50
dim(exp.df.tpm.brca[rowSds(exp.df.tpm.brca)>sd.limit,])
exp.df.tpm.brca[rowSds(exp.df.tpm.brca)>sd.limit,] -> filt.exp.df.tpm.brca


data.frame(gene = rownames(filt.exp.df.tpm.brca),filt.exp.df.tpm.brca) -> filt.exp.df.tpm.brca.first.col
head(filt.exp.df.tpm.brca.first.col[1:20,1:20])

write.table(filt.exp.df.tpm.brca.first.col, file = "./filt.exp.df.tpm.brca.txt", sep = "\t", quote = F, row.names = F)
