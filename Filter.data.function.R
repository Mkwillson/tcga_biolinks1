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
kirp.data$expression
assay(kirp.data$expression, "tpm_unstrand") -> exp.df.tpm.kirp
kirp.filt.exp <- exp.df.tpm.kirp[apply(exp.df.tpm.kirp,1,gp.style.filter,fold.change=3, delta=300, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]
dim(kirp.filt.exp)

# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(kirp.filt.exp),kirp.filt.exp) -> kirp.filt.exp.first.col
head(kirp.filt.exp.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(kirp.filt.exp.first.col, file = "./kirp.filt.exp.txt", sep = "\t", quote = F, row.names = F)
write.table(kirp.filt.exp.first.col[,1], file = "./kirp.filt.exp.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)