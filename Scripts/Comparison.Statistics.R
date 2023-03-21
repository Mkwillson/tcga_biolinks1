
## Filter BRCA
brca.graph <- readRDS(file = "~/tcga_biolinks1/RDS/brca.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/brca5.node.csv") -> brca.graph.df 

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> brca.graph.df.filt0.5

dim(brca.graph.df.filt0.5)
## Filter LUAD

luad.graph <- readRDS(file = "~/tcga_biolinks1/RDS/luad.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/luad4.node.csv") -> luad.graph.df

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 0.5) -> luad.graph.df.filt0.5

dim(luad.graph.df.filt0.5)

## Select genes that appear in both
filtered.brca.genes.luad <- filter(brca.graph.df.filt0.5, name %in% luad.graph.df$name)
unique(filtered.brca.genes.luad$name) -> unique.genes.brca.luad

dim(filtered.brca.genes.luad)

library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

## BRCA ON BRCA
# Calculate distances between all genes
#i = 1
brca.on.brca <- foreach(i=1:length(unique.genes.brca.luad), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.brca2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = unique.genes.brca.luad[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca.luad,])}

colnames(brca.on.brca) <- unique.genes.brca.luad

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(brca.graph),"name"), length(unique.genes.brca.luad)) -> rand.genes.brca.on.brca

random.brca.on.brca <- foreach(i=1:length(rand.genes.brca.on.brca), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = rand.genes.brca.on.brca[i],
    mode = c("all")
  ))
  return(temp.distance.table[rand.genes.brca.on.brca,])}

saveRDS(brca.on.brca, file = "~/tcga_biolinks1/stats/brca.on.brca2")
saveRDS(random.brca.on.brca, file = "~/tcga_biolinks1/stats/random.brca.on.brca")

## BRCA ON LUAD
# Calculate distances between all genes
#i = 1
brca.on.luad <- foreach(i=1:length(unique.genes.brca.luad), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.brca.luad[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca.luad,])}

 colnames(brca.on.luad) <- unique.genes.brca.luad

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.brca.luad)) -> rand.genes.brca.on.luad

random.brca.on.luad <- foreach(i=1:length(rand.genes.brca.on.luad), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = rand.genes.brca.on.luad[i],
    mode = c("all")
  ))
  return(temp.distance.table[rand.genes.brca.on.luad,])}

saveRDS(brca.on.luad, file = "~/tcga_biolinks1/stats/brca.on.luad")
saveRDS(random.brca.on.luad, file = "~/tcga_biolinks1/stats/random.brca.on.luad")

# readRDS(file = "~/tcga_biolinks1/stats/brca.on.luad") -> brca.on.luad
# readRDS(file = "~/tcga_biolinks1/stats/random.brca.on.luad") -> random.brca.on.luad
# readRDS(file = "~/tcga_biolinks1/stats/brca.on.brca") -> brca.on.brca
# readRDS(file = "~/tcga_biolinks1/stats/random.brca.on.brca") -> random.brca.on.brca

### STATS

brca.on.brca[which(is.infinite(brca.on.brca))] <- NA

# Create a random cutoff
quantile(colMeans(random.brca.on.brca, na.rm = T),0.01) -> random.cutoff.brca.on.brca
quantile(colMeans(random.brca.on.luad, na.rm = T),0.01) -> random.cutoff.brca.on.luad

dat.brca.luad <- data.frame(Distance = c(colMeans(brca.on.brca), colMeans(random.brca.on.brca, na.rm = T), colMeans(brca.on.luad), colMeans(random.brca.on.luad, na.rm = T)),
                            Distribution = factor(c(rep("BRCA Genes on BRCA Network",nrow(brca.on.brca)),rep("Random BRCA Genes on BRCA Network",nrow(random.brca.on.brca)),rep("BRCA Genes on LUAD Network",nrow(brca.on.luad)),rep("Random BRCA Genes on LUAD Network",nrow(random.brca.on.luad)))))

# Perform KS BRCA on BRCA v Random on BRCA
ks.test(colMeans(brca.on.brca), colMeans(random.brca.on.brca, na.rm = T)) -> ks.result.brca.on.brca

ks.result.brca.on.brca

# Print
ks.result.brca.on.brca$statistic
ks.result.brca.on.brca$p.value

# Perform KS BRCA on LUAD v Random on LUAD
ks.test(colMeans(brca.on.luad), colMeans(random.brca.on.luad, na.rm = T)) -> ks.result.brca.on.luad

ks.result.brca.on.luad

# Print
ks.result.brca.on.luad$statistic
ks.result.brca.on.luad$p.value

# Perform KS BRCA on BRCA v BRCA on LUAD
ks.test(colMeans(brca.on.brca), colMeans(brca.on.luad, na.rm = T)) -> ks.result.brca.v.luad

ks.result.brca.v.luad

# Print
ks.result.brca.v.luad$statistic
ks.result.brca.v.luad$p.value



# Create a text annotation of results
# paste0("D = ",round(ks.result0.5.brca$statistic,2), ifelse(ks.result0.5.brca$p.value==0, " p < 0.0001", paste0("p = >0.0001",round(ks.result0.5.brca$p.value,2)))) -> text.annot
# paste0("D = ",round(ks.result0.5.brca$statistic,2)) -> text.annot1
# paste0(ifelse(ks.result0.5.brca$p.value==0, " p < 0.0001", paste0("p = >0.001",round(ks.result0.5.brca$p.value,4)))) -> text.annot2
# length(ks.result.brca.on.brca$data$x) -> text.annot3

dim(brca.on.brca)
dim(brca.on.luad)
dim(random.brca.on.brca)
dim(random.brca.on.luad)

# Create gg plot
ggplot(dat.brca.luad, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca.luad, aes(xintercept=random.cutoff.brca.on.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  geom_vline(data=dat.brca.luad, aes(xintercept=random.cutoff.brca.on.luad,  colour=Distribution),
             linetype="dotted", linewidth=1) +
  # annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  # annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  theme(text = element_text(size = 20)) +
  ggtitle("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA and LUAD Networks") -> gg.brca.luad

gg.brca.luad
ggsave("~/tcga_biolinks1/Plots/gg.brca.luad.png", plot = gg.brca.luad, width = 18, height = 10)

######################################################################################################
## Filter BRCA
brca.graph <- readRDS(file = "~/tcga_biolinks1/RDS/brca.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/brca5.node.csv") -> brca.graph.df 

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> brca.graph.df.filt0.5

dim(brca.graph.df.filt0.5)
## Filter LUAD

kirc.graph <- readRDS(file = "~/tcga_biolinks1/RDS/kirc.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/kirc4.node.csv") -> kirc.graph.df

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> kirc.graph.df.filt0.5

dim(kirc.graph.df.filt0.5)

## Select genes that appear in both
filtered.brca.genes.kirc <- filter(brca.graph.df.filt0.5, name %in% kirc.graph.df$name)
unique(filtered.brca.genes.kirc$name) -> unique.genes.brca.kirc

dim(filtered.brca.genes.kirc)

library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

## BRCA ON BRCA
# Calculate distances between all genes
#i = 1
brca.on.brca <- foreach(i=1:length(unique.genes.brca.kirc), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.brca2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = unique.genes.brca.kirc[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca.kirc,])}

colnames(brca.on.brca) <- unique.genes.brca.kirc

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(brca.graph),"name"), length(unique.genes.brca.kirc)) -> rand.genes.brca.on.brca

random.brca.on.brca <- foreach(i=1:length(rand.genes.brca.on.brca), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = rand.genes.brca.on.brca[i],
    mode = c("all")
  ))
  return(temp.distance.table[rand.genes.brca.on.brca,])}

saveRDS(brca.on.brca, file = "~/tcga_biolinks1/stats/brca.on.brca")
saveRDS(random.brca.on.brca, file = "~/tcga_biolinks1/stats/random.brca.on.brca")

## BRCA ON LUAD
# Calculate distances between all genes
#i = 1
brca.on.kirc <- foreach(i=1:length(unique.genes.brca.kirc), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = unique.genes.brca.kirc[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca.kirc,])}

colnames(brca.on.kirc) <- unique.genes.brca.kirc

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirc.graph),"name"), length(unique.genes.brca.kirc)) -> rand.genes.brca.on.kirc

random.brca.on.kirc <- foreach(i=1:length(rand.genes.brca.on.kirc), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = rand.genes.brca.on.kirc[i],
    mode = c("all")
  ))
  return(temp.distance.table[rand.genes.brca.on.kirc,])}

saveRDS(brca.on.kirc, file = "~/tcga_biolinks1/stats/brca.on.kirc")
saveRDS(random.brca.on.kirc, file = "~/tcga_biolinks1/stats/random.brca.on.kirc")

# readRDS(file = "~/tcga_biolinks1/stats/brca.on.brca") -> brca.on.brca
# readRDS(file = "~/tcga_biolinks1/stats/random.brca.on.brca") -> random.brca.on.brca
# readRDS(file = "~/tcga_biolinks1/stats/brca.on.kirc") -> brca.on.kirc
# readRDS(file = "~/tcga_biolinks1/stats/random.brca.on.kirc") -> random.brca.on.kirc


### STATS

brca.on.brca[which(is.infinite(brca.on.brca))] <- NA

# Create a random cutoff
quantile(colMeans(random.brca.on.brca, na.rm = T),0.01) -> random.cutoff.brca.on.brca
quantile(colMeans(random.brca.on.kirc, na.rm = T),0.01) -> random.cutoff.brca.on.kirc

dat.brca.kirc <- data.frame(Distance = c(colMeans(brca.on.brca), colMeans(random.brca.on.brca, na.rm = T), colMeans(brca.on.kirc), colMeans(random.brca.on.kirc, na.rm = T)),
                            Distribution = factor(c(rep("BRCA Genes on BRCA Network",nrow(brca.on.brca)),rep("Random BRCA Genes on BRCA Network",nrow(random.brca.on.brca)),rep("BRCA Genes on KIRC Network",nrow(brca.on.kirc)),rep("Random BRCA Genes on KIRC Network",nrow(random.brca.on.kirc)))))

# Perform KS BRCA on BRCA v Random on BRCA
ks.test(colMeans(brca.on.brca), colMeans(random.brca.on.brca, na.rm = T)) -> ks.result.brca.on.brca

ks.result.brca.on.brca

# Print
ks.result.brca.on.brca$statistic
ks.result.brca.on.brca$p.value

# Perform KS BRCA on LUAD v Random on KIRC
ks.test(colMeans(brca.on.kirc), colMeans(random.brca.on.kirc, na.rm = T)) -> ks.result.brca.on.kirc

ks.result.brca.on.kirc

# Print
ks.result.brca.on.kirc$statistic
ks.result.brca.on.kirc$p.value

# Perform KS BRCA on BRCA v BRCA on KIRC
ks.test(colMeans(brca.on.brca), colMeans(brca.on.kirc, na.rm = T)) -> ks.result.brca.v.kirc

ks.result.brca.v.kirc

# Print
ks.result.brca.v.kirc$statistic
ks.result.brca.v.kirc$p.value



# Create a text annotation of results
# paste0("D = ",round(ks.result0.5.brca$statistic,2), ifelse(ks.result0.5.brca$p.value==0, " p < 0.0001", paste0("p = >0.0001",round(ks.result0.5.brca$p.value,2)))) -> text.annot
# paste0("D = ",round(ks.result0.5.brca$statistic,2)) -> text.annot1
# paste0(ifelse(ks.result0.5.brca$p.value==0, " p < 0.0001", paste0("p = >0.001",round(ks.result0.5.brca$p.value,4)))) -> text.annot2
# length(ks.result.brca.on.brca$data$x) -> text.annot3

dim(brca.on.brca)
dim(brca.on.kirc)
dim(random.brca.on.brca)
dim(random.brca.on.kirc)

# Create gg plot
ggplot(dat.brca.kirc, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca.kirc, aes(xintercept=random.cutoff.brca.on.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  geom_vline(data=dat.brca.kirc, aes(xintercept=random.cutoff.brca.on.kirc,  colour=Distribution),
             linetype="dotted", linewidth=1) +
  # annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  # annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  theme(text = element_text(size = 20)) +
  ggtitle("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA and KIRC Networks") -> gg.brca.kirc

gg.brca.kirc
ggsave("~/tcga_biolinks1/Plots/gg.brca.kirc.png", plot = gg.brca.kirc, width = 18, height = 10)

