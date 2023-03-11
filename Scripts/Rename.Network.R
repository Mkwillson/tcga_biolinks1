# Rename ensembl ids to match mutational data
# Note that the network files from consolidation have been manually renamed, using rename button in file section to avoid confusion.
# Construct function to rename ensemble ids without .version
rename.network.file <- function(x){
  temp.network <- read.delim(file = x)
  gsub("\\..*", '', temp.network$Regulator) -> temp.network$Regulator
  gsub("\\..*", '', temp.network$Target) -> temp.network$Target
  write.table(temp.network, file = gsub("network.txt", "network.rename.txt", x), sep = "\t")
}


# BRCA
x <- "~/tcga_biolinks1/ARACNe.data/BRCA5/brca5.network.txt"
rename.network.file(x)

# PRAD
x <- "~/tcga_biolinks1/ARACNe.data/PRAD5/prad5.network.txt"
rename.network.file(x)

# LUAD
x <- "~/tcga_biolinks1/ARACNe.data/LUAD4/luad4.network.txt"
rename.network.file(x)

# STAD - to do
x <- "~/tcga_biolinks1/ARACNe.data/STAD5/stad5.network.txt"
rename.network.file(x)

# KIRP
x <- "~/tcga_biolinks1/ARACNe.data/KIRP5/kirp5.network.txt"
rename.network.file(x)

# KIRC - to do
x <- "~/tcga_biolinks1/ARACNe.data/KIRC4/kirc4.network.txt"
rename.network.file(x)

# Export to cytoscape

################################

# Read back into R
read.graph("~/tcga_biolinks1/GRAPHML/brca5.1.graphml", format = "graphml") -> brca.graph
read.graph("~/tcga_biolinks1/GRAPHML/prad5.1.graphml", format = "graphml") -> prad.graph

read.csv(file = "~/tcga_biolinks1/GRAPHML/brca5.1.csv") -> brca.graph.df

#################################################
brca.graph.df %>% 
  filter(SIFT.Tolerated_low_confidence.brca==1) -> brca.graph.df.filt

unique(brca.graph.df.filt$name) -> all.genes

library(foreach)
library(doParallel)
registerDoParallel(cores = 30)

#i = 1
#out <- foreach(i=1:10, .combine = cbind)%dopar%{
out <- foreach(i=1:length(all.genes), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = all.genes[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes,])}

#sample(attr(V(GBMCOSMIC),"name"), length(all.genes)-1) -> all.genes.rand
sample(attr(V(brca.graph),"name"), length(all.genes)) -> all.genes.rand

out.rand.brca.graph <- foreach(i=1:length(all.genes.rand), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = all.genes.rand[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand,])}

saveRDS(out.rand.brca.graph, file = "./wherever?")
saveRDS(out, file = "./wherever?")

## draw a picture
hist(out)
## KS test
ks.test(as.numeric(out), as.numeric(out.rand.brca.graph)) -> ks.result

ks.result$statistic
ks.result$p.value

#############

hist(colMeans(out))
hist(colMeans(out.rand.brca.graph))
### turn inifnite into NAs
out[which(is.infinite(out))] <- NA

quantile(colMeans(out.rand.brca.graph, na.rm = T),0.01) -> random.cutoff

dat <- data.frame(Distance = c(colMeans(out), colMeans(out.rand.brca.graphc, na.rm = T)),
                  Distribution = factor(c(rep("Real",nrow(out)),rep("Random",nrow(out.rand.brca.graph)))))


ks.test(colMeans(out), colMeans(out.rand.brca.graph, na.rm = T)) -> ks.result

ks.result$statistic
ks.result$p.value


paste0("D = ",round(ks.result$statistic,2), ifelse(ks.result$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result$p.value,3)))) -> text.annot
paste0("D = ",round(ks.result$statistic,2)) -> text.annot1
paste0(ifelse(ks.result$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result$p.value,4)))) -> text.annot2
dim(out)
dim(out.rand.brca.graph)

ggplot(dat, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat, aes(xintercept=random.cutoff,  colour=Distribution),
             linetype="dashed", size=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = "N = 970") +
  ggtitle("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA Network")


##################
### create an input to run BRCA genes on a prad network

filtered.brca.genes.prad <- filter(brca.graph.df.filt, name %in% prad.graph.df.filt$name)
unique(filtered.brca.genes.prad$name) -> all.genes.brca.prad


