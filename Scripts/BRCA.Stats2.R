
## Read back into R
# library(igraph)
# read.graph("~/tcga_biolinks1/GRAPHML/brca5.final.graphml", format = "graphml") -> brca.graph
# saveRDS(brca.graph, file = "~/tcga_biolinks1/RDS/brca.graph")

brca.graph <- readRDS(file = "~/tcga_biolinks1/RDS/brca.graph")

# Read table into R
read.csv(file = "~/tcga_biolinks1/GRAPHML/brca5.node.csv") -> brca.graph.df
brca.graph.df$Percent.of.Patients.with.a.Missense.Mutation

##### 1.5% PATIENTS ##############################################################

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1.5) -> brca.graph.df.filt1.5

unique(brca.graph.df.filt1.5$name) -> unique.genes.brca1.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.brca1.5 <- foreach(i=1:length(unique.genes.brca1.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.brca2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = unique.genes.brca1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca1.5,])}

colnames(out.brca1.5) <- unique.genes.brca1.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(brca.graph),"name"), length(unique.genes.brca1.5)) -> all.genes.rand.brca1.5

out.rand.brca.graph1.5 <- foreach(i=1:length(all.genes.rand.brca1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = all.genes.rand.brca1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.brca1.5,])}

# # Save results of distances
saveRDS(out.rand.brca.graph1.5, file = "~/tcga_biolinks1/stats/out.rand.brca.graph1.5")
saveRDS(out.brca1.5, file = "~/tcga_biolinks1/stats/out.brca1.5")


# out.rand.brca.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph2")
# out.brca2 <- readRDS( file = "~/tcga_biolinks1/stats/out.brca2")

# Create histogram of column means
# hist(colMeans(out.brca2))
# hist(colMeans(out.rand.brca.graph2))

# Turn inifnite into NAs
out.brca1.5[which(is.infinite(out.brca1.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.brca.graph1.5, na.rm = T),0.01) -> random.cutoff1.5.brca

dat.brca1.5 <- data.frame(Distance = c(colMeans(out.brca1.5), colMeans(out.rand.brca.graph1.5, na.rm = T)),
                        Distribution = factor(c(rep("Real",nrow(out.brca1.5)),rep("Random",nrow(out.rand.brca.graph1.5)))))

# Perform KS test on column means
ks.test(colMeans(out.brca1.5), colMeans(out.rand.brca.graph1.5, na.rm = T)) -> ks.result1.5.brca

# Print
ks.result1.5.brca$statistic
ks.result1.5.brca$p.value
# length(ks.result2.brca$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.5.brca$statistic,2), ifelse(ks.result1.5.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.5.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.5.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.brca$p.value,4)))) -> text.annot2
length(ks.result1.5.brca$data$x) -> text.annot3

dim(out.brca1.5)
dim(out.rand.brca.graph1.5)

# Create gg plot
ggplot(dat.brca1.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca1.5, aes(xintercept=random.cutoff1.5.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA Network") -> gg.brca1.5

gg.brca1.5
ggsave("./gg.brca1.5.png", plot = gg.brca1.5, width = 10, height = 10)

#### 1% PATIENTS #############################################################################

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1) -> brca.graph.df.filt1

unique(brca.graph.df.filt1$name) -> unique.genes.brca1

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.brca1 <- foreach(i=1:length(unique.genes.brca1), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.brca1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = unique.genes.brca1[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca1,])}

colnames(out.brca1) <- unique.genes.brca1

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(brca.graph),"name"), length(unique.genes.brca1)) -> all.genes.rand.brca1

out.rand.brca.graph1 <- foreach(i=1:length(all.genes.rand.brca1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = all.genes.rand.brca1[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.brca1,])}

# # Save results of distances
saveRDS(out.rand.brca.graph1, file = "~/tcga_biolinks1/stats/out.rand.brca.graph1")
saveRDS(out.brca1, file = "~/tcga_biolinks1/stats/out.brca1")


# out.rand.brca.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph2")
# out.brca2 <- readRDS( file = "~/tcga_biolinks1/stats/out.brca2")

# Create histogram of column means
# hist(colMeans(out.brca2))
# hist(colMeans(out.rand.brca.graph2))

# Turn inifnite into NAs
out.brca1[which(is.infinite(out.brca1))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.brca.graph1, na.rm = T),0.01) -> random.cutoff1.brca

dat.brca1 <- data.frame(Distance = c(colMeans(out.brca1), colMeans(out.rand.brca.graph1, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.brca1)),rep("Random",nrow(out.rand.brca.graph1)))))

# Perform KS test on column means
ks.test(colMeans(out.brca1), colMeans(out.rand.brca.graph1, na.rm = T)) -> ks.result1.brca

# Print
ks.result1.brca$statistic
ks.result1.brca$p.value
# length(ks.result2.brca$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.brca$statistic,2), ifelse(ks.result1.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.brca$p.value,4)))) -> text.annot2
length(ks.result1.brca$data$x) -> text.annot3

dim(out.brca1)
dim(out.rand.brca.graph1)

# Create gg plot
ggplot(dat.brca1, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca1, aes(xintercept=random.cutoff1.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA Network") -> gg.brca1

gg.brca1
ggsave("./gg.brca1.png", plot = gg.brca1, width = 10, height = 10)

###### 0.5% PATIENTS #######################################################################

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> brca.graph.df.filt0.5

unique(brca.graph.df.filt0.5$name) -> unique.genes.brca0.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.brca0.5 <- foreach(i=1:length(unique.genes.brca0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.brca2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = unique.genes.brca0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca0.5,])}

colnames(out.brca0.5) <- unique.genes.brca0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(brca.graph),"name"), length(unique.genes.brca0.5)) -> all.genes.rand.brca0.5

out.rand.brca.graph0.5 <- foreach(i=1:length(all.genes.rand.brca0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = all.genes.rand.brca0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.brca0.5,])}

# # Save results of distances
saveRDS(out.rand.brca.graph0.5, file = "~/tcga_biolinks1/stats/out.rand.brca.graph0.5")
saveRDS(out.brca0.5, file = "~/tcga_biolinks1/stats/out.brca0.5")


# out.rand.brca.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph2")
# out.brca2 <- readRDS( file = "~/tcga_biolinks1/stats/out.brca2")

# Create histogram of column means
# hist(colMeans(out.brca2))
# hist(colMeans(out.rand.brca.graph2))

# Turn inifnite into NAs
out.brca0.5[which(is.infinite(out.brca0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.brca.graph0.5, na.rm = T),0.01) -> random.cutoff0.5.brca

dat.brca0.5 <- data.frame(Distance = c(colMeans(out.brca0.5), colMeans(out.rand.brca.graph0.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.brca0.5)),rep("Random",nrow(out.rand.brca.graph0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.brca0.5), colMeans(out.rand.brca.graph0.5, na.rm = T)) -> ks.result0.5.brca

# Print
ks.result0.5.brca$statistic
ks.result0.5.brca$p.value
# length(ks.result2.brca$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result0.5.brca$statistic,2), ifelse(ks.result0.5.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result0.5.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result0.5.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.brca$p.value,4)))) -> text.annot2
length(ks.result0.5.brca$data$x) -> text.annot3

dim(out.brca0.5)
dim(out.rand.brca.graph0.5)

# Create gg plot
ggplot(dat.brca0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca0.5, aes(xintercept=random.cutoff0.5.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA Network") -> gg.brca0.5

gg.brca0.5
ggsave("./gg.brca0.5.png", plot = gg.brca0.5, width = 10, height = 10)

##### PROBABLY DAMAGING >0.908 ##############################

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> brca.graph.df.filt.prob

unique(brca.graph.df.filt.prob$name) -> unique.genes.brca.prob

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.brca.prob <- foreach(i=1:length(unique.genes.brca.prob), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.brca.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = unique.genes.brca.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca.prob,])}

colnames(out.brca.prob) <- unique.genes.brca.prob

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(brca.graph),"name"), length(unique.genes.brca.prob)) -> all.genes.rand.brca.prob

out.rand.brca.graph.prob <- foreach(i=1:length(all.genes.rand.brca.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = all.genes.rand.brca.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.brca.prob,])}

# # Save results of distances
saveRDS(out.rand.brca.graph.prob, file = "~/tcga_biolinks1/stats/out.rand.brca.graph.prob")
saveRDS(out.brca.prob, file = "~/tcga_biolinks1/stats/out.brca.prob")

# Create histogram of column means
# hist(colMeans(out.brca.prob))
# hist(colMeans(out.rand.brca.graph.prob))

# Turn inifnite into NAs
out.brca.prob[which(is.infinite(out.brca.prob))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.brca.graph.prob, na.rm = T),0.01) -> random.cutoff.prob.brca

dat.brca.prob <- data.frame(Distance = c(colMeans(out.brca.prob), colMeans(out.rand.brca.graph.prob, na.rm = T)),
                        Distribution = factor(c(rep("Real",nrow(out.brca.prob)),rep("Random",nrow(out.rand.brca.graph.prob)))))

# Perform KS test on column means
ks.test(colMeans(out.brca.prob), colMeans(out.rand.brca.graph.prob, na.rm = T)) -> ks.result.prob.brca

# Print
ks.result.prob.brca$statistic
ks.result.prob.brca$p.value
# length(ks.result2.brca$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.prob.brca$statistic,2), ifelse(ks.result.prob.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.prob.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.prob.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.brca$p.value,4)))) -> text.annot2
length(ks.result.prob.brca$data$x) -> text.annot3

dim(out.brca.prob)
dim(out.rand.brca.graph.prob)

# Create gg plot
ggplot(dat.brca.prob, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca.prob, aes(xintercept=random.cutoff.prob.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA Network") -> gg.brca.prob

gg.brca.prob
ggsave("./gg.brca.prob.png", plot = gg.brca.prob, width = 10, height = 10)

##### POSSIBLY DAMAGING =<0.908 >0.466 ##############################
brca.graph.df$PolyPhen.Mean
brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean <=0.908) %>%
  filter(PolyPhen.Mean >0.466)-> brca.graph.df.filt.poss

unique(brca.graph.df.filt.poss$name) -> unique.genes.brca.poss

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.brca.poss <- foreach(i=1:length(unique.genes.brca.poss), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.brca.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = unique.genes.brca.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca.poss,])}

colnames(out.brca.poss) <- unique.genes.brca.poss

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(brca.graph),"name"), length(unique.genes.brca.poss)) -> all.genes.rand.brca.poss

out.rand.brca.graph.poss <- foreach(i=1:length(all.genes.rand.brca.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = all.genes.rand.brca.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.brca.poss,])}

# # Save results of distances
saveRDS(out.rand.brca.graph.poss, file = "~/tcga_biolinks1/stats/out.rand.brca.graph.poss")
saveRDS(out.brca.poss, file = "~/tcga_biolinks1/stats/out.brca.poss")

# Create histogram of column means
# hist(colMeans(out.brca.poss))
# hist(colMeans(out.rand.brca.graph.poss))

# Turn inifnite into NAs
out.brca.poss[which(is.infinite(out.brca.poss))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.brca.graph.poss, na.rm = T),0.01) -> random.cutoff.poss.brca

dat.brca.poss <- data.frame(Distance = c(colMeans(out.brca.poss), colMeans(out.rand.brca.graph.poss, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.brca.poss)),rep("Random",nrow(out.rand.brca.graph.poss)))))

# Perform KS test on column means
ks.test(colMeans(out.brca.poss), colMeans(out.rand.brca.graph.poss, na.rm = T)) -> ks.result.poss.brca

# Print
ks.result.poss.brca$statistic
ks.result.poss.brca$p.value
# length(ks.result2.brca$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.poss.brca$statistic,2), ifelse(ks.result.poss.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.poss.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.poss.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.brca$p.value,4)))) -> text.annot2
length(ks.result.poss.brca$data$x) -> text.annot3

dim(out.brca.poss)
dim(out.rand.brca.graph.poss)

# Create gg plot
ggplot(dat.brca.poss, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca.poss, aes(xintercept=random.cutoff.poss.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA Network") -> gg.brca.poss

gg.brca.poss
ggsave("./gg.brca.poss.png", plot = gg.brca.poss, width = 10, height = 10)

##### ALL DAMAGING >0.466 ##########################################################################

brca.graph.df$PolyPhen.Mean
brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> brca.graph.df.filt.dam

unique(brca.graph.df.filt.dam$name) -> unique.genes.brca.dam

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.brca.dam <- foreach(i=1:length(unique.genes.brca.dam), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.brca.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = unique.genes.brca.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.brca.dam,])}

colnames(out.brca.dam) <- unique.genes.brca.dam

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(brca.graph),"name"), length(unique.genes.brca.dam)) -> all.genes.rand.brca.dam

out.rand.brca.graph.dam <- foreach(i=1:length(all.genes.rand.brca.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    brca.graph,
    to = all.genes.rand.brca.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.brca.dam,])}

# # Save results of distances
saveRDS(out.rand.brca.graph.dam, file = "~/tcga_biolinks1/stats/out.rand.brca.graph.dam")
saveRDS(out.brca.dam, file = "~/tcga_biolinks1/stats/out.brca.dam")

# Create histogram of column means
# hist(colMeans(out.brca.dam))
# hist(colMeans(out.rand.brca.graph.dam))

# Turn inifnite into NAs
out.brca.dam[which(is.infinite(out.brca.dam))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.brca.graph.dam, na.rm = T),0.01) -> random.cutoff.dam.brca

dat.brca.dam <- data.frame(Distance = c(colMeans(out.brca.dam), colMeans(out.rand.brca.graph.dam, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.brca.dam)),rep("Random",nrow(out.rand.brca.graph.dam)))))

# Perform KS test on column means
ks.test(colMeans(out.brca.dam), colMeans(out.rand.brca.graph.dam, na.rm = T)) -> ks.result.dam.brca

# Print
ks.result.dam.brca$statistic
ks.result.dam.brca$p.value
# length(ks.result2.brca$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.brca$statistic,2), ifelse(ks.result.dam.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.brca$p.value,4)))) -> text.annot2
length(ks.result.dam.brca$data$x) -> text.annot3

dim(out.brca.dam)
dim(out.rand.brca.graph.dam)

# Create gg plot
ggplot(dat.brca.dam, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca.dam, aes(xintercept=random.cutoff.dam.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA Network") -> gg.brca.dam

gg.brca.dam
ggsave("./gg.brca.dam.png", plot = gg.brca.dam, width = 10, height = 10)



#######################################

## Other method for saving
# pdf(file = "./test.pdf")
# gg.brca2
# dev.off()
# png(file = "./test.png")
# gg.brca2
# dev.off()
#### close all the connections
# graphics.off()

