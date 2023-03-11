
## Read back into R
# library(igraph)
# read.graph("~/tcga_biolinks1/GRAPHML/prad5.final.graphml", format = "graphml") -> prad.graph
# saveRDS(prad.graph, file = "~/tcga_biolinks1/RDS/prad.graph")

prad.graph <- readRDS(file = "~/tcga_biolinks1/RDS/prad.graph")

# Read table into R
read.csv(file = "~/tcga_biolinks1/GRAPHML/prad5.node.csv") -> prad.graph.df

##### 1.5% PATIENTS ##############################################################

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1.5) -> prad.graph.df.filt1.5

unique(prad.graph.df.filt1.5$name) -> unique.genes.prad1.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad1.5 <- foreach(i=1:length(unique.genes.prad1.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad1.5,])}

colnames(out.prad1.5) <- unique.genes.prad1.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad1.5)) -> all.genes.rand.prad1.5

out.rand.prad.graph1.5 <- foreach(i=1:length(all.genes.rand.prad1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad1.5,])}

# # Save results of distances
saveRDS(out.rand.prad.graph1.5, file = "~/tcga_biolinks1/stats/out.rand.prad.graph1.5")
saveRDS(out.prad1.5, file = "~/tcga_biolinks1/stats/out.prad1.5")


# out.rand.prad.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph2")
# out.prad2 <- readRDS( file = "~/tcga_biolinks1/stats/out.prad2")

# Create histogram of column means
# hist(colMeans(out.prad2))
# hist(colMeans(out.rand.prad.graph2))

# Turn inifnite into NAs
out.prad1.5[which(is.infinite(out.prad1.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph1.5, na.rm = T),0.01) -> random.cutoff1.5.prad

dat.prad1.5 <- data.frame(Distance = c(colMeans(out.prad1.5), colMeans(out.rand.prad.graph1.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.prad1.5)),rep("Random",nrow(out.rand.prad.graph1.5)))))

# Perform KS test on column means
ks.test(colMeans(out.prad1.5), colMeans(out.rand.prad.graph1.5, na.rm = T)) -> ks.result1.5.prad

# Print
ks.result1.5.prad$statistic
ks.result1.5.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.5.prad$statistic,2), ifelse(ks.result1.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.5.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.prad$p.value,4)))) -> text.annot2
length(ks.result1.5.prad$data$x) -> text.annot3

dim(out.prad1.5)
dim(out.rand.prad.graph1.5)

# Create gg plot
ggplot(dat.prad1.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad1.5, aes(xintercept=random.cutoff1.5.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and prad Genes on the prad Network") -> gg.prad1.5

gg.prad1.5
ggsave("~/tcga_biolinks1/gg.prad1.5.png", plot = gg.prad1.5, width = 10, height = 10)

#### 1% PATIENTS #############################################################################

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1) -> prad.graph.df.filt1

unique(prad.graph.df.filt1$name) -> unique.genes.prad1

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad1 <- foreach(i=1:length(unique.genes.prad1), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad1[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad1,])}

colnames(out.prad1) <- unique.genes.prad1

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad1)) -> all.genes.rand.prad1

out.rand.prad.graph1 <- foreach(i=1:length(all.genes.rand.prad1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad1[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad1,])}

# # Save results of distances
saveRDS(out.rand.prad.graph1, file = "~/tcga_biolinks1/stats/out.rand.prad.graph1")
saveRDS(out.prad1, file = "~/tcga_biolinks1/stats/out.prad1")


# out.rand.prad.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph2")
# out.prad2 <- readRDS( file = "~/tcga_biolinks1/stats/out.prad2")

# Create histogram of column means
# hist(colMeans(out.prad2))
# hist(colMeans(out.rand.prad.graph2))

# Turn inifnite into NAs
out.prad1[which(is.infinite(out.prad1))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph1, na.rm = T),0.01) -> random.cutoff1.prad

dat.prad1 <- data.frame(Distance = c(colMeans(out.prad1), colMeans(out.rand.prad.graph1, na.rm = T)),
                        Distribution = factor(c(rep("Real",nrow(out.prad1)),rep("Random",nrow(out.rand.prad.graph1)))))

# Perform KS test on column means
ks.test(colMeans(out.prad1), colMeans(out.rand.prad.graph1, na.rm = T)) -> ks.result1.prad

# Print
ks.result1.prad$statistic
ks.result1.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.prad$statistic,2), ifelse(ks.result1.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.prad$p.value,4)))) -> text.annot2
length(ks.result1.prad$data$x) -> text.annot3

dim(out.prad1)
dim(out.rand.prad.graph1)

# Create gg plot
ggplot(dat.prad1, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad1, aes(xintercept=random.cutoff1.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and prad Genes on the prad Network") -> gg.prad1

gg.prad1
ggsave("~/tcga_biolinks1/gg.prad1.png", plot = gg.prad1, width = 10, height = 10)

###### 0.5% PATIENTS #######################################################################

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> prad.graph.df.filt0.5

unique(prad.graph.df.filt0.5$name) -> unique.genes.prad0.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad0.5 <- foreach(i=1:length(unique.genes.prad0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad0.5,])}

colnames(out.prad0.5) <- unique.genes.prad0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad0.5)) -> all.genes.rand.prad0.5

out.rand.prad.graph0.5 <- foreach(i=1:length(all.genes.rand.prad0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad0.5,])}

# # Save results of distances
saveRDS(out.rand.prad.graph0.5, file = "~/tcga_biolinks1/stats/out.rand.prad.graph0.5")
saveRDS(out.prad0.5, file = "~/tcga_biolinks1/stats/out.prad0.5")


# out.rand.prad.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph2")
# out.prad2 <- readRDS( file = "~/tcga_biolinks1/stats/out.prad2")

# Create histogram of column means
# hist(colMeans(out.prad2))
# hist(colMeans(out.rand.prad.graph2))

# Turn inifnite into NAs
out.prad0.5[which(is.infinite(out.prad0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph0.5, na.rm = T),0.01) -> random.cutoff0.5.prad

dat.prad0.5 <- data.frame(Distance = c(colMeans(out.prad0.5), colMeans(out.rand.prad.graph0.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.prad0.5)),rep("Random",nrow(out.rand.prad.graph0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.prad0.5), colMeans(out.rand.prad.graph0.5, na.rm = T)) -> ks.result0.5.prad

# Print
ks.result0.5.prad$statistic
ks.result0.5.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result0.5.prad$statistic,2), ifelse(ks.result0.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result0.5.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result0.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.prad$p.value,4)))) -> text.annot2
length(ks.result0.5.prad$data$x) -> text.annot3

dim(out.prad0.5)
dim(out.rand.prad.graph0.5)

# Create gg plot
ggplot(dat.prad0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad0.5, aes(xintercept=random.cutoff0.5.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and prad Genes on the prad Network") -> gg.prad0.5

gg.prad0.5
ggsave("~/tcga_biolinks1/gg.prad0.5.png", plot = gg.prad0.5, width = 10, height = 10)

##### PROBABLY DAMAGING >0.908 ##############################

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> prad.graph.df.filt.prob

unique(prad.graph.df.filt.prob$name) -> unique.genes.prad.prob

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad.prob <- foreach(i=1:length(unique.genes.prad.prob), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad.prob,])}

colnames(out.prad.prob) <- unique.genes.prad.prob

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad.prob)) -> all.genes.rand.prad.prob

out.rand.prad.graph.prob <- foreach(i=1:length(all.genes.rand.prad.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad.prob,])}

# # Save results of distances
saveRDS(out.rand.prad.graph.prob, file = "~/tcga_biolinks1/stats/out.rand.prad.graph.prob")
saveRDS(out.prad.prob, file = "~/tcga_biolinks1/stats/out.prad.prob")

# Create histogram of column means
# hist(colMeans(out.prad.prob))
# hist(colMeans(out.rand.prad.graph.prob))

# Turn inifnite into NAs
out.prad.prob[which(is.infinite(out.prad.prob))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph.prob, na.rm = T),0.01) -> random.cutoff.prob.prad

dat.prad.prob <- data.frame(Distance = c(colMeans(out.prad.prob), colMeans(out.rand.prad.graph.prob, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.prad.prob)),rep("Random",nrow(out.rand.prad.graph.prob)))))

# Perform KS test on column means
ks.test(colMeans(out.prad.prob), colMeans(out.rand.prad.graph.prob, na.rm = T)) -> ks.result.prob.prad

# Print
ks.result.prob.prad$statistic
ks.result.prob.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.prob.prad$statistic,2), ifelse(ks.result.prob.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.prob.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.prob.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.prad$p.value,4)))) -> text.annot2
length(ks.result.prob.prad$data$x) -> text.annot3

dim(out.prad.prob)
dim(out.rand.prad.graph.prob)

# Create gg plot
ggplot(dat.prad.prob, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad.prob, aes(xintercept=random.cutoff.prob.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and prad Genes on the prad Network") -> gg.prad.prob

gg.prad.prob
ggsave("~/tcga_biolinks1/gg.prad.prob.png", plot = gg.prad.prob, width = 10, height = 10)

##### POSSIBLY DAMAGING =<0.908 >0.466 ##############################
prad.graph.df$PolyPhen.Mean
prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean <=0.908) %>%
  filter(PolyPhen.Mean >0.466)-> prad.graph.df.filt.poss

unique(prad.graph.df.filt.poss$name) -> unique.genes.prad.poss

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad.poss <- foreach(i=1:length(unique.genes.prad.poss), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad.poss,])}

colnames(out.prad.poss) <- unique.genes.prad.poss

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad.poss)) -> all.genes.rand.prad.poss

out.rand.prad.graph.poss <- foreach(i=1:length(all.genes.rand.prad.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad.poss,])}

# # Save results of distances
saveRDS(out.rand.prad.graph.poss, file = "~/tcga_biolinks1/stats/out.rand.prad.graph.poss")
saveRDS(out.prad.poss, file = "~/tcga_biolinks1/stats/out.prad.poss")

# Create histogram of column means
# hist(colMeans(out.prad.poss))
# hist(colMeans(out.rand.prad.graph.poss))

# Turn inifnite into NAs
out.prad.poss[which(is.infinite(out.prad.poss))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph.poss, na.rm = T),0.01) -> random.cutoff.poss.prad

dat.prad.poss <- data.frame(Distance = c(colMeans(out.prad.poss), colMeans(out.rand.prad.graph.poss, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.prad.poss)),rep("Random",nrow(out.rand.prad.graph.poss)))))

# Perform KS test on column means
ks.test(colMeans(out.prad.poss), colMeans(out.rand.prad.graph.poss, na.rm = T)) -> ks.result.poss.prad

# Print
ks.result.poss.prad$statistic
ks.result.poss.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.poss.prad$statistic,2), ifelse(ks.result.poss.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.poss.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.poss.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.prad$p.value,4)))) -> text.annot2
length(ks.result.poss.prad$data$x) -> text.annot3

dim(out.prad.poss)
dim(out.rand.prad.graph.poss)

# Create gg plot
ggplot(dat.prad.poss, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad.poss, aes(xintercept=random.cutoff.poss.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and prad Genes on the prad Network") -> gg.prad.poss

gg.prad.poss
ggsave("~/tcga_biolinks1/gg.prad.poss.png", plot = gg.prad.poss, width = 10, height = 10)

##### ALL DAMAGING >0.466 ##########################################################################

prad.graph.df$PolyPhen.Mean
prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> prad.graph.df.filt.dam

unique(prad.graph.df.filt.dam$name) -> unique.genes.prad.dam

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad.dam <- foreach(i=1:length(unique.genes.prad.dam), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad.dam,])}

colnames(out.prad.dam) <- unique.genes.prad.dam

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad.dam)) -> all.genes.rand.prad.dam

out.rand.prad.graph.dam <- foreach(i=1:length(all.genes.rand.prad.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad.dam,])}

# # Save results of distances
saveRDS(out.rand.prad.graph.dam, file = "~/tcga_biolinks1/stats/out.rand.prad.graph.dam")
saveRDS(out.prad.dam, file = "~/tcga_biolinks1/stats/out.prad.dam")

# Create histogram of column means
# hist(colMeans(out.prad.dam))
# hist(colMeans(out.rand.prad.graph.dam))

# Turn inifnite into NAs
out.prad.dam[which(is.infinite(out.prad.dam))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph.dam, na.rm = T),0.01) -> random.cutoff.dam.prad

dat.prad.dam <- data.frame(Distance = c(colMeans(out.prad.dam), colMeans(out.rand.prad.graph.dam, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.prad.dam)),rep("Random",nrow(out.rand.prad.graph.dam)))))

# Perform KS test on column means
ks.test(colMeans(out.prad.dam), colMeans(out.rand.prad.graph.dam, na.rm = T)) -> ks.result.dam.prad

# Print
ks.result.dam.prad$statistic
ks.result.dam.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.prad$statistic,2), ifelse(ks.result.dam.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.prad$p.value,4)))) -> text.annot2
length(ks.result.dam.prad$data$x) -> text.annot3

dim(out.prad.dam)
dim(out.rand.prad.graph.dam)

# Create gg plot
ggplot(dat.prad.dam, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad.dam, aes(xintercept=random.cutoff.dam.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and prad Genes on the prad Network") -> gg.prad.dam

gg.prad.dam
ggsave("~/tcga_biolinks1/gg.prad.dam.png", plot = gg.prad.dam, width = 10, height = 10)

##### ALL DAMAGING >0.466 and %>0.5 ##########################################################################

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(PolyPhen.Mean >0.466)-> prad.graph.df.filt.dam.0.5

unique(prad.graph.df.filt.dam.0.5$name) -> unique.genes.prad.dam.0.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad.dam.0.5 <- foreach(i=1:length(unique.genes.prad.dam.0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad.dam.0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad.dam.0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad.dam.0.5,])}

colnames(out.prad.dam.0.5) <- unique.genes.prad.dam.0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad.dam.0.5)) -> all.genes.rand.prad.dam.0.5

out.rand.prad.graph.dam.0.5 <- foreach(i=1:length(all.genes.rand.prad.dam.0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad.dam.0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad.dam.0.5,])}

# # Save results of distances
saveRDS(out.rand.prad.graph.dam.0.5, file = "~/tcga_biolinks1/stats/out.rand.prad.graph.dam.0.5")
saveRDS(out.prad.dam.0.5, file = "~/tcga_biolinks1/stats/out.prad.dam.0.5")
# readRDS(file = "~/tcga_biolinks1/stats/out.prad.dam.0.5") -> out.prad.dam.0.5
# readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph.dam.0.5") -> out.rand.prad.graph.dam.0.5

# Create histogram of column means
# hist(colMeans(out.prad.dam.0.5))
# hist(colMeans(out.rand.prad.graph.dam.0.5))

# Turn inifnite into NAs
out.prad.dam.0.5[which(is.infinite(out.prad.dam.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph.dam.0.5, na.rm = T),0.01) -> random.cutoff.dam.0.5.prad

dat.prad.dam.0.5 <- data.frame(Distance = c(colMeans(out.prad.dam.0.5), colMeans(out.rand.prad.graph.dam.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.prad.dam.0.5)),rep("Random",nrow(out.rand.prad.graph.dam.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.prad.dam.0.5), colMeans(out.rand.prad.graph.dam.0.5, na.rm = T)) -> ks.result.dam.0.5.prad

# Print
ks.result.dam.0.5.prad$statistic
ks.result.dam.0.5.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.0.5.prad$statistic,2), ifelse(ks.result.dam.0.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.0.5.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.0.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.prad$p.value,4)))) -> text.annot2
length(ks.result.dam.0.5.prad$data$x) -> text.annot3

dim(out.prad.dam.0.5)
dim(out.rand.prad.graph.dam.0.5)

# Create gg plot
ggplot(dat.prad.dam.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad.dam.0.5, aes(xintercept=random.cutoff.dam.0.5.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and PRAD Genes on the PRAD Network") -> gg.prad.dam.0.5

gg.prad.dam.0.5
ggsave("~/tcga_biolinks1/Plots/gg.prad.dam.0.5.png", plot = gg.prad.dam.0.5, width = 10, height = 10)

##### ALL DAMAGING >0.466 and %>1.5 ##########################################################################

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1.5) %>%
  filter(PolyPhen.Mean >0.466)-> prad.graph.df.filt.dam.1.5

unique(prad.graph.df.filt.dam.1.5$name) -> unique.genes.prad.dam.1.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad.dam.1.5 <- foreach(i=1:length(unique.genes.prad.dam.1.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad.dam.1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad.dam.1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad.dam.1.5,])}

colnames(out.prad.dam.1.5) <- unique.genes.prad.dam.1.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad.dam.1.5)) -> all.genes.rand.prad.dam.1.5

out.rand.prad.graph.dam.1.5 <- foreach(i=1:length(all.genes.rand.prad.dam.1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad.dam.1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad.dam.1.5,])}

# # Save results of distances
saveRDS(out.rand.prad.graph.dam.1.5, file = "~/tcga_biolinks1/stats/out.rand.prad.graph.dam.1.5")
saveRDS(out.prad.dam.1.5, file = "~/tcga_biolinks1/stats/out.prad.dam.1.5")
# readRDS(file = "~/tcga_biolinks1/stats/out.prad.dam.1.5") -> out.prad.dam.1.5
# readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph.dam.1.5") -> out.rand.prad.graph.dam.1.5

# Create histogram of column means
# hist(colMeans(out.prad.dam.1.5))
# hist(colMeans(out.rand.prad.graph.dam.1.5))

# Turn inifnite into NAs
out.prad.dam.1.5[which(is.infinite(out.prad.dam.1.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph.dam.1.5, na.rm = T),0.01) -> random.cutoff.dam.1.5.prad

dat.prad.dam.1.5 <- data.frame(Distance = c(colMeans(out.prad.dam.1.5), colMeans(out.rand.prad.graph.dam.1.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.prad.dam.1.5)),rep("Random",nrow(out.rand.prad.graph.dam.1.5)))))

# Perform KS test on column means
ks.test(colMeans(out.prad.dam.1.5), colMeans(out.rand.prad.graph.dam.1.5, na.rm = T)) -> ks.result.dam.1.5.prad

# Print
ks.result.dam.1.5.prad$statistic
ks.result.dam.1.5.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.1.5.prad$statistic,2), ifelse(ks.result.dam.1.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.1.5.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.1.5.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.1.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.1.5.prad$p.value,4)))) -> text.annot2
length(ks.result.dam.1.5.prad$data$x) -> text.annot3

dim(out.prad.dam.1.5)
dim(out.rand.prad.graph.dam.1.5)

# Create gg plot
ggplot(dat.prad.dam.1.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad.dam.1.5, aes(xintercept=random.cutoff.dam.1.5.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and PRAD Genes on the PRAD Network") -> gg.prad.dam.1.5

gg.prad.dam.1.5
ggsave("~/tcga_biolinks1/Plots/gg.prad.dam.1.5.png", plot = gg.prad.dam.1.5, width = 10, height = 10)

#### DELETERIOUS <0.05 ##########################################################################

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(SIFT.Mean <0.05)-> prad.graph.df.filt.del

unique(prad.graph.df.filt.del$name) -> unique.genes.prad.del

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad.del <- foreach(i=1:length(unique.genes.prad.del), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad.del), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad.del[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad.del,])}

colnames(out.prad.del) <- unique.genes.prad.del

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad.del)) -> all.genes.rand.prad.del

out.rand.prad.graph.del <- foreach(i=1:length(all.genes.rand.prad.del), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad.del[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad.del,])}

# # Save results of distances
saveRDS(out.rand.prad.graph.del, file = "~/tcga_biolinks1/stats/out.rand.prad.graph.del")
saveRDS(out.prad.del, file = "~/tcga_biolinks1/stats/out.prad.del")

# Create histogram of column means
# hist(colMeans(out.prad.del))
# hist(colMeans(out.rand.prad.graph.del))

# Turn inifnite into NAs
out.prad.del[which(is.infinite(out.prad.del))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph.del, na.rm = T),0.01) -> random.cutoff.del.prad

dat.prad.del <- data.frame(Distance = c(colMeans(out.prad.del), colMeans(out.rand.prad.graph.del, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.prad.del)),rep("Random",nrow(out.rand.prad.graph.del)))))

# Perform KS test on column means
ks.test(colMeans(out.prad.del), colMeans(out.rand.prad.graph.del, na.rm = T)) -> ks.result.del.prad

# Print
ks.result.del.prad$statistic
ks.result.del.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.prad$statistic,2), ifelse(ks.result.del.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.prad$p.value,4)))) -> text.annot2
length(ks.result.del.prad$data$x) -> text.annot3

dim(out.prad.del)
dim(out.rand.prad.graph.del)

# Create gg plot
ggplot(dat.prad.del, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad.del, aes(xintercept=random.cutoff.del.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and PRAD Genes on the PRAD Network") -> gg.prad.del

gg.prad.del
ggsave("~/tcga_biolinks1/Plots/gg.prad.del.png", plot = gg.prad.del, width = 10, height = 10)

##### ALL DELETERIOUS <0.05 and %>0.5 ##########################################################################

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(SIFT.Mean <0.05)-> prad.graph.df.filt.del.0.5

unique(prad.graph.df.filt.del.0.5$name) -> unique.genes.prad.del.0.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.prad.del.0.5 <- foreach(i=1:length(unique.genes.prad.del.0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.prad.del.0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = unique.genes.prad.del.0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.prad.del.0.5,])}

colnames(out.prad.del.0.5) <- unique.genes.prad.del.0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(prad.graph),"name"), length(unique.genes.prad.del.0.5)) -> all.genes.rand.prad.del.0.5

out.rand.prad.graph.del.0.5 <- foreach(i=1:length(all.genes.rand.prad.del.0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    prad.graph,
    to = all.genes.rand.prad.del.0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.prad.del.0.5,])}

# # Save results of distances
saveRDS(out.rand.prad.graph.del.0.5, file = "~/tcga_biolinks1/stats/out.rand.prad.graph.del.0.5")
saveRDS(out.prad.del.0.5, file = "~/tcga_biolinks1/stats/out.prad.del.0.5")
# readRDS(file = "~/tcga_biolinks1/stats/out.prad.del.0.5") -> out.prad.del.0.5
# readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph.del.0.5") -> out.rand.prad.graph.del.0.5

# Create histogram of column means
# hist(colMeans(out.prad.del.0.5))
# hist(colMeans(out.rand.prad.graph.del.0.5))

# Turn inifnite into NAs
out.prad.del.0.5[which(is.infinite(out.prad.del.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.prad.graph.del.0.5, na.rm = T),0.01) -> random.cutoff.del.0.5.prad

dat.prad.del.0.5 <- data.frame(Distance = c(colMeans(out.prad.del.0.5), colMeans(out.rand.prad.graph.del.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.prad.del.0.5)),rep("Random",nrow(out.rand.prad.graph.del.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.prad.del.0.5), colMeans(out.rand.prad.graph.del.0.5, na.rm = T)) -> ks.result.del.0.5.prad

# Print
ks.result.del.0.5.prad$statistic
ks.result.del.0.5.prad$p.value
# length(ks.result2.prad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.0.5.prad$statistic,2), ifelse(ks.result.del.0.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.prad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.0.5.prad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.0.5.prad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.prad$p.value,4)))) -> text.annot2
length(ks.result.del.0.5.prad$data$x) -> text.annot3

dim(out.prad.del.0.5)
dim(out.rand.prad.graph.del.0.5)

# Create gg plot
ggplot(dat.prad.del.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.prad.del.0.5, aes(xintercept=random.cutoff.del.0.5.prad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and PRAD Genes on the PRAD Network") -> gg.prad.del.0.5

gg.prad.del.0.5
ggsave("~/tcga_biolinks1/Plots/gg.prad.del.0.5.png", plot = gg.prad.del.0.5, width = 10, height = 10)



#######################################

## Other method for saving
# pdf(file = "./test.pdf")
# gg.prad2
# dev.off()
# png(file = "./test.png")
# gg.prad2
# dev.off()
#### close all the connections
# graphics.off()

