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

# Select data to filter
stad.data$expression
assay(stad.data$expression, "tpm_unstrand") -> exp.df.tpm.stad

# Look at data distribution
hist(exp.df.tpm.prad)

# Filter appropriatley
stad.filt.exp <- exp.df.tpm.stad[apply(exp.df.tpm.stad,1,gp.style.filter,fold.change=3, delta=300, prop=0.05, base=20, prop.base=0.05, na.rm = TRUE, neg.rm = TRUE),]

# Check size of data
dim(stad.filt.exp)
# Look at data distribution
hist(stad.filt.exp, breaks=8)


# Take just the names of genes to act as TFs to test against
data.frame(gene = rownames(stad.filt.exp),stad.filt.exp) -> stad.filt.exp.first.col
head(stad.filt.exp.first.col[1:20,1:20])

# Save expression data and TF data as txt files
write.table(stad.filt.exp.first.col, file = "./stad.filt.exp.txt", sep = "\t", quote = F, row.names = F)
write.table(stad.filt.exp.first.col[,1], file = "./stad.filt.exp.tfs.txt", sep = "\t", quote = F, row.names = F, col.names = F)

#### NOT R CODE THIS IS JAVA FOR USE IN TERMINAL ####

# Calculate threshold for MI - BRCA5
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD/stad.filt.exp.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD/stad.filt.exp.tfs.txt \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20

# 1 bootstrap
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD/stad.filt.exp.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD/stad.filt.exp.tfs.txt \
--pvalue 1E-8 \
--seed 1 \
--threads 20

# 100 bootstrap
for i in {1..100}
do 
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD/stad.filt.exp.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD/stad.filt.exp.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate
java -Xmx5G -jar Aracne.jar -o outputFolder --consolidate