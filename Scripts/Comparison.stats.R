
### create an input to run BRCA genes on a prad network

filtered.brca.genes.prad <- filter(brca.graph.df.filt, name %in% prad.graph.df.filt$name)
unique(filtered.brca.genes.prad$name) -> all.genes.brca.prad