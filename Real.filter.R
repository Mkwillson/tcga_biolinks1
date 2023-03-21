#### gp.style.filter function ####
# Construct a function for filtering data
gp.style.filter <- function (x,fold.change, delta, prop, base, prop.base, na.rm = TRUE, neg.rm = TRUE) 
{
  if (na.rm == TRUE){x <- x[!is.na(x)]}
  if (neg.rm == TRUE){x <- x[x > 0]}
  lenx <- length(x)
  if (lenx < 4 || lenx < prop + 1){return(FALSE)}
  srtd <- sort(x)
  if (prop < 1) {
    bot <- lenx * prop
  }else {bot <- prop}
  top <- lenx - bot
  fold.filter<-FALSE
  delta.filter<-FALSE
  base.filter<-FALSE
  if (max(srtd[bot:top]) / min(srtd[bot:top]) > fold.change){fold.filter<-TRUE}
  if (max(srtd[bot:top]) - min(srtd[bot:top]) > delta){delta.filter<-TRUE}
  if (length(which(srtd > base))>(prop.base*lenx)){base.filter<-TRUE}
  if (fold.filter ==TRUE & delta.filter == TRUE & base.filter == TRUE){
    return(TRUE)
  }else{
    return(FALSE)
  }
}


#### BRCA ####
# Filter at various levels until a subset of sufficient size is created (5000-1300 approx)
# Select data to filter
brca.data$expression
assay(brca.data$expression, "tpm_unstrand") -> exp.df.tpm.brca

# Look at data distribution
hist(exp.df.tpm.brca[3,])
hist(log10(rowMeans(exp.df.tpm.brca)))
hist(log10(rowMeans(exp.df.tpm.brca)))
hist(log10(rowSds(exp.df.tpm.brca)))


#### Filter attempt 1 - fold.change 1, delta 100, base 3 ####
brca.filt.exp1 <- exp.df.tpm.brca[apply(exp.df.tpm.brca,1,gp.style.filter,fold.change=1, delta=5, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(brca.filt.exp1)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(brca.filt.exp1),brca.filt.exp1) -> brca.filt.exp1.first.col
head(brca.filt.exp1.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(brca.filt.exp1.first.col, file = "./brca.filt.exp1.txt", sep = "\t", quote = F, row.names = F)
write.table(brca.filt.exp1.first.col[,1], file = "./brca.filt.exp1.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2 - fold.change 2, delta 10, base 3 ####
brca.filt.exp2 <- exp.df.tpm.brca[apply(exp.df.tpm.brca,1,gp.style.filter,fold.change=2, delta=10, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(brca.filt.exp2)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(brca.filt.exp2),brca.filt.exp2) -> brca.filt.exp2.first.col
head(brca.filt.exp2.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(brca.filt.exp2.first.col, file = "./brca.filt.exp2.txt", sep = "\t", quote = F, row.names = F)
write.table(brca.filt.exp2.first.col[,1], file = "./brca.filt.exp2.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2.5 - fold.change 2.5, delta 12.5, base 3 ####
brca.filt.exp2.5 <- exp.df.tpm.brca[apply(exp.df.tpm.brca,1,gp.style.filter,fold.change=2.5, delta=12.5, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(brca.filt.exp2.5)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(brca.filt.exp2.5),brca.filt.exp2.5) -> brca.filt.exp2.5.first.col
head(brca.filt.exp2.5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(brca.filt.exp2.5.first.col, file = "./brca.filt.exp2.5.txt", sep = "\t", quote = F, row.names = F)
write.table(brca.filt.exp2.5.first.col[,1], file = "./brca.filt.exp2.5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 3 - fold.change 3, delta 15, base 3 ####
brca.filt.exp3 <- exp.df.tpm.brca[apply(exp.df.tpm.brca,1,gp.style.filter,fold.change=3, delta=15, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(brca.filt.exp3)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(brca.filt.exp3),brca.filt.exp3) -> brca.filt.exp3.first.col
head(brca.filt.exp3.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(brca.filt.exp3.first.col, file = "./brca.filt.exp3.txt", sep = "\t", quote = F, row.names = F)
write.table(brca.filt.exp3.first.col[,1], file = "./brca.filt.exp3.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 4 - fold.change 4, delta 20, base 3 ####
brca.filt.exp4 <- exp.df.tpm.brca[apply(exp.df.tpm.brca,1,gp.style.filter,fold.change=4, delta=20, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(brca.filt.exp4)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(brca.filt.exp4),brca.filt.exp4) -> brca.filt.exp4.first.col
head(brca.filt.exp4.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(brca.filt.exp4.first.col, file = "./brca.filt.exp4.txt", sep = "\t", quote = F, row.names = F)
write.table(brca.filt.exp4.first.col[,1], file = "./brca.filt.exp4.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
 
#### Filter attempt 5 - fold.change 5, delta 25, base 3 ####
brca.filt.exp5 <- exp.df.tpm.brca[apply(exp.df.tpm.brca,1,gp.style.filter,fold.change=4, delta=20, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(brca.filt.exp5)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(brca.filt.exp5),brca.filt.exp5) -> brca.filt.exp5.first.col
head(brca.filt.exp5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(brca.filt.exp5.first.col, file = "./brca.filt.exp5.txt", sep = "\t", quote = F, row.names = F)
write.table(brca.filt.exp5.first.col[,1], file = "./brca.filt.exp5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### PRAD ###################################################################################################################
# Select data to filter
prad.data$expression
assay(prad.data$expression, "tpm_unstrand") -> exp.df.tpm.prad

# Look at data distribution
hist(exp.df.tpm.prad)

#### Filter attempt 1 - fold.change 1, delta 100, base 20 ####
prad.filt.exp1 <- exp.df.tpm.prad[apply(exp.df.tpm.prad,1,gp.style.filter,fold.change=1, delta=100, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(prad.filt.exp1)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(prad.filt.exp1),prad.filt.exp1) -> prad.filt.exp1.first.col
head(prad.filt.exp1.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(prad.filt.exp1.first.col, file = "./prad.filt.exp1.txt", sep = "\t", quote = F, row.names = F)
write.table(prad.filt.exp1.first.col[,1], file = "./prad.filt.exp1.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2 - fold.change 2, delta 200, base 20 ####
prad.filt.exp2 <- exp.df.tpm.prad[apply(exp.df.tpm.prad,1,gp.style.filter,fold.change=2, delta=200, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(prad.filt.exp2)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(prad.filt.exp2),prad.filt.exp2) -> prad.filt.exp2.first.col
head(prad.filt.exp2.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(prad.filt.exp2.first.col, file = "./prad.filt.exp2.txt", sep = "\t", quote = F, row.names = F)
write.table(prad.filt.exp2.first.col[,1], file = "./prad.filt.exp2.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2.5 - fold.change 2.5, delta 20, base 3 ####
prad.filt.exp2.5 <- exp.df.tpm.prad[apply(exp.df.tpm.prad,1,gp.style.filter,fold.change=2.5, delta=250, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(prad.filt.exp2.5)


# Look at data distribution
hist(prad.filt.exp2.5)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(prad.filt.exp2.5),prad.filt.exp2.5) -> prad.filt.exp2.5.first.col
head(prad.filt.exp2.5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(prad.filt.exp2.5.first.col, file = "./prad.filt.exp2.5.txt", sep = "\t", quote = F, row.names = F)
write.table(prad.filt.exp2.5.first.col[,1], file = "./prad.filt.exp2.5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 3 -  fold.change 3, delta 30, base 3 ####
prad.filt.exp3 <- exp.df.tpm.prad[apply(exp.df.tpm.prad,1,gp.style.filter,fold.change=3, delta=30, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(prad.filt.exp3)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(prad.filt.exp3),prad.filt.exp3) -> prad.filt.exp3.first.col
head(prad.filt.exp3.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(prad.filt.exp3.first.col, file = "./prad.filt.exp3.txt", sep = "\t", quote = F, row.names = F)
write.table(prad.filt.exp3.first.col[,1], file = "./prad.filt.exp3.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 4 - fold.change 3, delta 20, base 3 ####
prad.filt.exp4 <- exp.df.tpm.prad[apply(exp.df.tpm.prad,1,gp.style.filter,fold.change=3, delta=20, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(prad.filt.exp4)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(prad.filt.exp4),prad.filt.exp4) -> prad.filt.exp4.first.col
head(prad.filt.exp4.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(prad.filt.exp4.first.col, file = "./prad.filt.exp4.txt", sep = "\t", quote = F, row.names = F)
write.table(prad.filt.exp4.first.col[,1], file = "./prad.filt.exp4.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 5 - fold.change 4, delta 10, base 3 ####
prad.filt.exp5 <- exp.df.tpm.prad[apply(exp.df.tpm.prad,1,gp.style.filter,fold.change=4, delta=20, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(prad.filt.exp5)
# 504, 554

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(prad.filt.exp5),prad.filt.exp5) -> prad.filt.exp5.first.col
head(prad.filt.exp5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(prad.filt.exp5.first.col, file = "./prad.filt.exp5.txt", sep = "\t", quote = F, row.names = F)
write.table(prad.filt.exp5.first.col[,1], file = "./prad.filt.exp5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 6 - fold.change 3.5, delta 10, base 3 ####
prad.filt.exp6 <- exp.df.tpm.prad[apply(exp.df.tpm.prad,1,gp.style.filter,fold.change=3.5, delta=10, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(prad.filt.exp6)


# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(prad.filt.exp6),prad.filt.exp6) -> prad.filt.exp6.first.col
head(prad.filt.exp5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(prad.filt.exp6.first.col, file = "./prad.filt.exp6.txt", sep = "\t", quote = F, row.names = F)
write.table(prad.filt.exp6.first.col[,1], file = "./prad.filt.exp6.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)


#### LUAD ###################################################################################################################
# Select data to filter
luad.data$expression
assay(luad.data$expression, "tpm_unstrand") -> exp.df.tpm.luad

# Look at data distribution
hist(exp.df.tpm.luad)

#### Filter attempt 1 - fold.change 1, delta 100, base 20 ####
luad.filt.exp1 <- exp.df.tpm.luad[apply(exp.df.tpm.luad,1,gp.style.filter,fold.change=1, delta=100, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(luad.filt.exp1)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(luad.filt.exp1),luad.filt.exp1) -> luad.filt.exp1.first.col
head(luad.filt.exp1.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(luad.filt.exp1.first.col, file = "./luad.filt.exp1.txt", sep = "\t", quote = F, row.names = F)
write.table(luad.filt.exp1.first.col[,1], file = "./luad.filt.exp1.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2 - fold.change 2, delta 200, base 20 ####
luad.filt.exp2 <- exp.df.tpm.luad[apply(exp.df.tpm.luad,1,gp.style.filter,fold.change=2, delta=200, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(luad.filt.exp2)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(luad.filt.exp2),luad.filt.exp2) -> luad.filt.exp2.first.col
head(luad.filt.exp2.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(luad.filt.exp2.first.col, file = "./luad.filt.exp2.txt", sep = "\t", quote = F, row.names = F)
write.table(luad.filt.exp2.first.col[,1], file = "./luad.filt.exp2.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2.5 - fold.change 2.5, delta 250, base 20 ####
luad.filt.exp2.5 <- exp.df.tpm.luad[apply(exp.df.tpm.luad,1,gp.style.filter,fold.change=2.5, delta=250, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(luad.filt.exp2.5)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(luad.filt.exp2.5),luad.filt.exp2.5) -> luad.filt.exp2.5.first.col
head(luad.filt.exp2.5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(luad.filt.exp2.5.first.col, file = "./luad.filt.exp2.5.txt", sep = "\t", quote = F, row.names = F)
write.table(luad.filt.exp2.5.first.col[,1], file = "./luad.filt.exp2.5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 3 - fold.change 3, delta 30, base 3 ####
luad.filt.exp3 <- exp.df.tpm.luad[apply(exp.df.tpm.luad,1,gp.style.filter,fold.change=3, delta=30, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(luad.filt.exp3)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(luad.filt.exp3),luad.filt.exp3) -> luad.filt.exp3.first.col
head(luad.filt.exp3.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(luad.filt.exp3.first.col, file = "./luad.filt.exp3.txt", sep = "\t", quote = F, row.names = F)
write.table(luad.filt.exp3.first.col[,1], file = "./luad.filt.exp3.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 4 - fold.change 3, delta 20, base 3 ####
luad.filt.exp4 <- exp.df.tpm.luad[apply(exp.df.tpm.luad,1,gp.style.filter,fold.change=3, delta=20, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(luad.filt.exp4)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(luad.filt.exp4),luad.filt.exp4) -> luad.filt.exp4.first.col
head(luad.filt.exp4.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(luad.filt.exp4.first.col, file = "./luad.filt.exp4.txt", sep = "\t", quote = F, row.names = F)
write.table(luad.filt.exp4.first.col[,1], file = "./luad.filt.exp4.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)


#### KIRP ###################################################################################################################
# Select data to filter
kirp.data$expression
assay(kirp.data$expression, "tpm_unstrand") -> exp.df.tpm.kirp

# Look at data distribution
hist(exp.df.tpm.kirp)

#### Filter attempt 1 - fold.change 1, delta 100, base 20 ####
kirp.filt.exp1 <- exp.df.tpm.kirp[apply(exp.df.tpm.kirp,1,gp.style.filter,fold.change=1, delta=100, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirp.filt.exp1)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirp.filt.exp1),kirp.filt.exp1) -> kirp.filt.exp1.first.col
head(kirp.filt.exp1.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirp.filt.exp1.first.col, file = "./kirp.filt.exp1.txt", sep = "\t", quote = F, row.names = F)
write.table(kirp.filt.exp1.first.col[,1], file = "./kirp.filt.exp1.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2 - fold.change 2, delta 200, base 20 ####
kirp.filt.exp2 <- exp.df.tpm.kirp[apply(exp.df.tpm.kirp,1,gp.style.filter,fold.change=2, delta=200, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirp.filt.exp2)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirp.filt.exp2),kirp.filt.exp2) -> kirp.filt.exp2.first.col
head(kirp.filt.exp2.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirp.filt.exp2.first.col, file = "./kirp.filt.exp2.txt", sep = "\t", quote = F, row.names = F)
write.table(kirp.filt.exp2.first.col[,1], file = "./kirp.filt.exp2.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2.5 - fold.change 2.5, delta 250, base 20 ####
kirp.filt.exp2.5 <- exp.df.tpm.kirp[apply(exp.df.tpm.kirp,1,gp.style.filter,fold.change=2.5, delta=250, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirp.filt.exp2.5)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirp.filt.exp2.5),kirp.filt.exp2.5) -> kirp.filt.exp2.5.first.col
head(kirp.filt.exp2.5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirp.filt.exp2.5.first.col, file = "./kirp.filt.exp2.5.txt", sep = "\t", quote = F, row.names = F)
write.table(kirp.filt.exp2.5.first.col[,1], file = "./kirp.filt.exp2.5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 3 - fold.change 3, delta 30, base 20 ####
kirp.filt.exp3 <- exp.df.tpm.kirp[apply(exp.df.tpm.kirp,1,gp.style.filter,fold.change=3, delta=30, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirp.filt.exp3)

# Look at data distribution
hist(kirp.filt.exp3)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirp.filt.exp3),kirp.filt.exp3) -> kirp.filt.exp3.first.col
head(kirp.filt.exp3.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirp.filt.exp3.first.col, file = "./kirp.filt.exp3.txt", sep = "\t", quote = F, row.names = F)
write.table(kirp.filt.exp3.first.col[,1], file = "./kirp.filt.exp3.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 4 - fold.change 3, delta 20, base 5 ####
kirp.filt.exp4 <- exp.df.tpm.kirp[apply(exp.df.tpm.kirp,1,gp.style.filter,fold.change=3, delta=30, prop=0.05, base=5, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirp.filt.exp4)

# Look at data distribution
hist(kirp.filt.exp4)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirp.filt.exp4),kirp.filt.exp4) -> kirp.filt.exp4.first.col
head(kirp.filt.exp4.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirp.filt.exp4.first.col, file = "./kirp.filt.exp4.txt", sep = "\t", quote = F, row.names = F)
write.table(kirp.filt.exp4.first.col[,1], file = "./kirp.filt.exp4.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
#### Filter attempt 5 - fold.change 3, delta 20, base 5 ####
kirp.filt.exp5 <- exp.df.tpm.kirp[apply(exp.df.tpm.kirp,1,gp.style.filter,fold.change=3, delta=20, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirp.filt.exp5)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirp.filt.exp5),kirp.filt.exp5) -> kirp.filt.exp5.first.col
head(kirp.filt.exp5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirp.filt.exp5.first.col, file = "./kirp.filt.exp5.txt", sep = "\t", quote = F, row.names = F)
write.table(kirp.filt.exp5.first.col[,1], file = "./kirp.filt.exp5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)
#### KIRC ###################################################################################################################
# Select data to filter
kirc.data$expression
assay(kirc.data$expression, "tpm_unstrand") -> exp.df.tpm.kirc

# Look at data distribution
hist(exp.df.tpm.kirc)

#### Filter attempt 1 -  fold.change 1, delta 100, base 20 ####
kirc.filt.exp1 <- exp.df.tpm.kirc[apply(exp.df.tpm.kirc,1,gp.style.filter,fold.change=1, delta=100, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]

dim(kirc.filt.exp1)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirc.filt.exp1),kirc.filt.exp1) -> kirc.filt.exp1.first.col
head(kirc.filt.exp1.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirc.filt.exp1.first.col, file = "./kirc.filt.exp1.txt", sep = "\t", quote = F, row.names = F)
write.table(kirc.filt.exp1.first.col[,1], file = "./kirc.filt.exp1.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2 - fold.change 2, delta 200, base 20 ####
kirc.filt.exp2 <- exp.df.tpm.kirc[apply(exp.df.tpm.kirc,1,gp.style.filter,fold.change=2, delta=200, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirc.filt.exp2)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirc.filt.exp2),kirc.filt.exp2) -> kirc.filt.exp2.first.col
head(kirc.filt.exp2.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirc.filt.exp2.first.col, file = "./kirc.filt.exp2.txt", sep = "\t", quote = F, row.names = F)
write.table(kirc.filt.exp2.first.col[,1], file = "./kirc.filt.exp2.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2.5 - fold.change 2.5, delta 250, base 20 ####
kirc.filt.exp2.5 <- exp.df.tpm.kirc[apply(exp.df.tpm.kirc,1,gp.style.filter,fold.change=2.5, delta=250, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirc.filt.exp2.5)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirc.filt.exp2.5),kirc.filt.exp2.5) -> kirc.filt.exp2.5.first.col
head(kirc.filt.exp2.5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirc.filt.exp2.5.first.col, file = "./kirc.filt.exp2.5.txt", sep = "\t", quote = F, row.names = F)
write.table(kirc.filt.exp2.5.first.col[,1], file = "./kirc.filt.exp2.5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 3 - fold.change 3, delta 100, base 20 ####
kirc.filt.exp3 <- exp.df.tpm.kirc[apply(exp.df.tpm.kirc,1,gp.style.filter,fold.change=3, delta=100, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirc.filt.exp3)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirc.filt.exp3),kirc.filt.exp3) -> kirc.filt.exp3.first.col
head(kirc.filt.exp3.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirc.filt.exp3.first.col, file = "./kirc.filt.exp3.txt", sep = "\t", quote = F, row.names = F)
write.table(kirc.filt.exp3.first.col[,1], file = "./kirc.filt.exp3.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 4 - fold.change 3, delta 20, base 5 ####
kirc.filt.exp4 <- exp.df.tpm.kirc[apply(exp.df.tpm.kirc,1,gp.style.filter,fold.change=3, delta=20, prop=0.05, base=5, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirc.filt.exp4)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirc.filt.exp4),kirc.filt.exp4) -> kirc.filt.exp4.first.col
head(kirc.filt.exp4.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirc.filt.exp4.first.col, file = "./kirc.filt.exp4.txt", sep = "\t", quote = F, row.names = F)
write.table(kirc.filt.exp4.first.col[,1], file = "./kirc.filt.exp4.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### STAD ###################################################################################################################
# Select data to filter
stad.data$expression
assay(stad.data$expression, "tpm_unstrand") -> exp.df.tpm.stad

# Look at data distribution
hist(exp.df.tpm.stad)

#### Filter attempt 1 - fold.change 1, delta 100, base 20 ####
stad.filt.exp1 <- exp.df.tpm.stad[apply(exp.df.tpm.stad,1,gp.style.filter,fold.change=1, delta=100, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(stad.filt.exp1)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(stad.filt.exp1),stad.filt.exp1) -> stad.filt.exp1.first.col
head(stad.filt.exp1.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(stad.filt.exp1.first.col, file = "./stad.filt.exp1.txt", sep = "\t", quote = F, row.names = F)
write.table(stad.filt.exp1.first.col[,1], file = "./stad.filt.exp1.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2 - fold.change 2, delta 200, base 20 ####
stad.filt.exp2 <- exp.df.tpm.stad[apply(exp.df.tpm.stad,1,gp.style.filter,fold.change=2, delta=200, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(stad.filt.exp2)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(stad.filt.exp2),stad.filt.exp2) -> stad.filt.exp2.first.col
head(stad.filt.exp2.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(stad.filt.exp2.first.col, file = "./stad.filt.exp2.txt", sep = "\t", quote = F, row.names = F)
write.table(stad.filt.exp2.first.col[,1], file = "./stad.filt.exp2.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 2.5 - fold.change 2.5, delta 250, base 20 ####
stad.filt.exp2.5 <- exp.df.tpm.stad[apply(exp.df.tpm.stad,1,gp.style.filter,fold.change=2.5, delta=250, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(stad.filt.exp2.5)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(stad.filt.exp2.5),stad.filt.exp2.5) -> stad.filt.exp2.5.first.col
head(stad.filt.exp2.5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(stad.filt.exp2.5.first.col, file = "./stad.filt.exp2.5.txt", sep = "\t", quote = F, row.names = F)
write.table(stad.filt.exp2.5.first.col[,1], file = "./stad.filt.exp2.5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 3 - fold.change 3 , delta 30, base 5 ####
stad.filt.exp3 <- exp.df.tpm.stad[apply(exp.df.tpm.stad,1,gp.style.filter,fold.change=3, delta=30, prop=0.05, base=5, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(stad.filt.exp3)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(stad.filt.exp3),stad.filt.exp3) -> stad.filt.exp3.first.col
head(stad.filt.exp3.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(stad.filt.exp3.first.col, file = "./stad.filt.exp3.txt", sep = "\t", quote = F, row.names = F)
write.table(stad.filt.exp3.first.col[,1], file = "./stad.filt.exp3.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 4 - fold.change 4 , delta 20, base 3 ####
stad.filt.exp4 <- exp.df.tpm.stad[apply(exp.df.tpm.stad,1,gp.style.filter,fold.change=4, delta=20, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(stad.filt.exp4)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(stad.filt.exp4),stad.filt.exp4) -> stad.filt.exp4.first.col
head(stad.filt.exp4.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(stad.filt.exp4.first.col, file = "./stad.filt.exp4.txt", sep = "\t", quote = F, row.names = F)
write.table(stad.filt.exp4.first.col[,1], file = "./stad.filt.exp4.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### Filter attempt 5 - fold.change 3, delta 20, base 3 ####
stad.filt.exp5 <- exp.df.tpm.stad[apply(exp.df.tpm.stad,1,gp.style.filter,fold.change=3, delta=20, prop=0.05, base=3, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(stad.filt.exp5)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(stad.filt.exp5),stad.filt.exp5) -> stad.filt.exp5.first.col
head(stad.filt.exp5.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(stad.filt.exp5.first.col, file = "./stad.filt.exp5.txt", sep = "\t", quote = F, row.names = F)
write.table(stad.filt.exp5.first.col[,1], file = "./stad.filt.exp5.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)



#### NOT R CODE THIS IS JAVA FOR USE IN TERMINAL #####################################

#### BRCA5 ####

# Calculate threshold for MI
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/BRCA5/brca.filt.exp5.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/BRCA5 \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/BRCA5/brca.filt.exp5.tfs.txt \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20 

# 100 bootstrap
for i in {1..100}
do 
time java -Xmx25G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/BRCA5/brca.filt.exp5.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/BRCA5 --tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/BRCA5/brca.filt.exp5.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate
java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -o /home/ug_megan/tcga_biolinks1/ARACNe.data/BRCA5 --consolidate

#### PRAD6 #####
# Calculate threshold for MI - PRAD6
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD6/prad.filt.exp6.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD6 \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD6/prad.filt.exp6.tfs.txt \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20

# 100 bootstrap - PRAD2
for i in {1..100}
do 
time java -Xmx25G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD6/prad.filt.exp6.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD6 --tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD6/prad.filt.exp6.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate
java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -o /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD6 --consolidate

#### LUAD4 ############################################################################

# Calculate threshold for MI
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD4/luad.filt.exp4.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD4 \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD4/luad.filt.exp4.tfs.txt \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20

# 100 bootstrap - LUAD2
for i in {1..100}
do 
time java -Xmx25G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD4/luad.filt.exp4.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD4 --tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD4/luad.filt.exp4.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate
java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -o /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD4 --consolidate

#### KIRP5 ############################################################################
# Calculate threshold for MI
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP5/kirp.filt.exp5.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP5 \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP5/kirp.filt.exp5.tfs.txt \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20

# 100 bootstrap
for i in {1..100}
do 
time java -Xmx25G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP5/kirp.filt.exp5.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP5 --tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP5/kirp.filt.exp5.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate
java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP5 --consolidate

#### KIRC4 ############################################################################
# Calculate threshold for MI
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4/kirc.filt.exp4.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4 \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4/kirc.filt.exp4.tfs.txt \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20

# 100 bootstrap
for i in {1..100}
do 
time java -Xmx25G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4/kirc.filt.exp4.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4 --tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4/kirc.filt.exp4.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate
java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4 --consolidate

#### STAD5 #####
# Calculate threshold for MI
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5/stad.filt.exp5.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5 \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5/stad.filt.exp5.tfs.txt \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20

# 100 bootstrap
for i in {1..100}
do 
time java -Xmx25G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5/stad.filt.exp5.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5 --tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5/stad.filt.exp5.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate
java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -o /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5 --consolidate


