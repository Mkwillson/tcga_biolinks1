# Rename ensembl ids to match mutational data
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
