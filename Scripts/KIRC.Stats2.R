## Read back into R
# library(igraph)
read.graph("~/tcga_biolinks1/GRAPHML/kirc4.final.graphml", format = "graphml") -> kirc.graph
# saveRDS(kirc.graph, file = "~/tcga_biolinks1/RDS/kirc.graph")

kirc.graph <- readRDS(file = "~/tcga_biolinks1/RDS/kirc.graph")

# Read table into R
read.csv(file = "~/tcga_biolinks1/GRAPHML/kirc4.node.csv") -> kirc.graph.df
kirc.graph.df$Percent.of.Patients.with.a.Missense.Mutation

##### 1.5% PATIENTS ##############################################################

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1.5) -> kirc.graph.df.filt1.5

unique(kirc.graph.df.filt1.5$name) -> unique.genes.kirc1.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirc1.5 <- foreach(i=1:length(unique.genes.kirc1.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirc2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = unique.genes.kirc1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirc1.5,])}

colnames(out.kirc1.5) <- unique.genes.kirc1.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirc.graph),"name"), length(unique.genes.kirc1.5)) -> all.genes.rand.kirc1.5

out.rand.kirc.graph1.5 <- foreach(i=1:length(all.genes.rand.kirc1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = all.genes.rand.kirc1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirc1.5,])}

# # Save results of distances
saveRDS(out.rand.kirc.graph1.5, file = "~/tcga_biolinks1/stats/out.rand.kirc.graph1.5")
saveRDS(out.kirc1.5, file = "~/tcga_biolinks1/stats/out.kirc1.5")


# out.rand.kirc.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph2")
# out.kirc2 <- readRDS( file = "~/tcga_biolinks1/stats/out.kirc2")

# Create histogram of column means
# hist(colMeans(out.kirc2))
# hist(colMeans(out.rand.kirc.graph2))

# Turn inifnite into NAs
out.kirc1.5[which(is.infinite(out.kirc1.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirc.graph1.5, na.rm = T),0.01) -> random.cutoff1.5.kirc

dat.kirc1.5 <- data.frame(Distance = c(colMeans(out.kirc1.5), colMeans(out.rand.kirc.graph1.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.kirc1.5)),rep("Random",nrow(out.rand.kirc.graph1.5)))))

# Perform KS test on column means
ks.test(colMeans(out.kirc1.5), colMeans(out.rand.kirc.graph1.5, na.rm = T)) -> ks.result1.5.kirc

# Print
ks.result1.5.kirc$statistic
ks.result1.5.kirc$p.value
# length(ks.result2.kirc$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.5.kirc$statistic,2), ifelse(ks.result1.5.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.kirc$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.5.kirc$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.5.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.kirc$p.value,4)))) -> text.annot2
length(ks.result1.5.kirc$data$x) -> text.annot3

dim(out.kirc1.5)
dim(out.rand.kirc.graph1.5)

# Create gg plot
ggplot(dat.kirc1.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirc1.5, aes(xintercept=random.cutoff1.5.kirc,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirc Genes on the kirc Network") -> gg.kirc1.5

gg.kirc1.5
ggsave("~/tcga_biolinks1/gg.kirc1.5.png", plot = gg.kirc1.5, width = 10, height = 10)

#### 1% PATIENTS #############################################################################

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1) -> kirc.graph.df.filt1

unique(kirc.graph.df.filt1$name) -> unique.genes.kirc1

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirc1 <- foreach(i=1:length(unique.genes.kirc1), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirc1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = unique.genes.kirc1[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirc1,])}

colnames(out.kirc1) <- unique.genes.kirc1

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirc.graph),"name"), length(unique.genes.kirc1)) -> all.genes.rand.kirc1

out.rand.kirc.graph1 <- foreach(i=1:length(all.genes.rand.kirc1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = all.genes.rand.kirc1[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirc1,])}

# # Save results of distances
saveRDS(out.rand.kirc.graph1, file = "~/tcga_biolinks1/stats/out.rand.kirc.graph1")
saveRDS(out.kirc1, file = "~/tcga_biolinks1/stats/out.kirc1")


# out.rand.kirc.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph2")
# out.kirc2 <- readRDS( file = "~/tcga_biolinks1/stats/out.kirc2")

# Create histogram of column means
# hist(colMeans(out.kirc2))
# hist(colMeans(out.rand.kirc.graph2))

# Turn inifnite into NAs
out.kirc1[which(is.infinite(out.kirc1))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirc.graph1, na.rm = T),0.01) -> random.cutoff1.kirc

dat.kirc1 <- data.frame(Distance = c(colMeans(out.kirc1), colMeans(out.rand.kirc.graph1, na.rm = T)),
                        Distribution = factor(c(rep("Real",nrow(out.kirc1)),rep("Random",nrow(out.rand.kirc.graph1)))))

# Perform KS test on column means
ks.test(colMeans(out.kirc1), colMeans(out.rand.kirc.graph1, na.rm = T)) -> ks.result1.kirc

# Print
ks.result1.kirc$statistic
ks.result1.kirc$p.value
# length(ks.result2.kirc$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.kirc$statistic,2), ifelse(ks.result1.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.kirc$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.kirc$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.kirc$p.value,4)))) -> text.annot2
length(ks.result1.kirc$data$x) -> text.annot3

dim(out.kirc1)
dim(out.rand.kirc.graph1)

# Create gg plot
ggplot(dat.kirc1, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirc1, aes(xintercept=random.cutoff1.kirc,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirc Genes on the kirc Network") -> gg.kirc1

gg.kirc1
ggsave("~/tcga_biolinks1/gg.kirc1.png", plot = gg.kirc1, width = 10, height = 10)

###### 0.5% PATIENTS #######################################################################

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> kirc.graph.df.filt0.5

unique(kirc.graph.df.filt0.5$name) -> unique.genes.kirc0.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirc0.5 <- foreach(i=1:length(unique.genes.kirc0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirc2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = unique.genes.kirc0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirc0.5,])}

colnames(out.kirc0.5) <- unique.genes.kirc0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirc.graph),"name"), length(unique.genes.kirc0.5)) -> all.genes.rand.kirc0.5

out.rand.kirc.graph0.5 <- foreach(i=1:length(all.genes.rand.kirc0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = all.genes.rand.kirc0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirc0.5,])}

# # Save results of distances
saveRDS(out.rand.kirc.graph0.5, file = "~/tcga_biolinks1/stats/out.rand.kirc.graph0.5")
saveRDS(out.kirc0.5, file = "~/tcga_biolinks1/stats/out.kirc0.5")


# out.rand.kirc.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph2")
# out.kirc2 <- readRDS( file = "~/tcga_biolinks1/stats/out.kirc2")

# Create histogram of column means
# hist(colMeans(out.kirc2))
# hist(colMeans(out.rand.kirc.graph2))

# Turn inifnite into NAs
out.kirc0.5[which(is.infinite(out.kirc0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirc.graph0.5, na.rm = T),0.01) -> random.cutoff0.5.kirc

dat.kirc0.5 <- data.frame(Distance = c(colMeans(out.kirc0.5), colMeans(out.rand.kirc.graph0.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.kirc0.5)),rep("Random",nrow(out.rand.kirc.graph0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.kirc0.5), colMeans(out.rand.kirc.graph0.5, na.rm = T)) -> ks.result0.5.kirc

# Print
ks.result0.5.kirc$statistic
ks.result0.5.kirc$p.value
# length(ks.result2.kirc$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result0.5.kirc$statistic,2), ifelse(ks.result0.5.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.kirc$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result0.5.kirc$statistic,2)) -> text.annot1
paste0(ifelse(ks.result0.5.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.kirc$p.value,4)))) -> text.annot2
length(ks.result0.5.kirc$data$x) -> text.annot3

dim(out.kirc0.5)
dim(out.rand.kirc.graph0.5)

# Create gg plot
ggplot(dat.kirc0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirc0.5, aes(xintercept=random.cutoff0.5.kirc,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirc Genes on the kirc Network") -> gg.kirc0.5

gg.kirc0.5
ggsave("~/tcga_biolinks1/gg.kirc0.5.png", plot = gg.kirc0.5, width = 10, height = 10)

##### PROBABLY DAMAGING >0.908 ##############################

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> kirc.graph.df.filt.prob

unique(kirc.graph.df.filt.prob$name) -> unique.genes.kirc.prob

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirc.prob <- foreach(i=1:length(unique.genes.kirc.prob), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirc.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = unique.genes.kirc.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirc.prob,])}

colnames(out.kirc.prob) <- unique.genes.kirc.prob

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirc.graph),"name"), length(unique.genes.kirc.prob)) -> all.genes.rand.kirc.prob

out.rand.kirc.graph.prob <- foreach(i=1:length(all.genes.rand.kirc.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = all.genes.rand.kirc.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirc.prob,])}

# # Save results of distances
saveRDS(out.rand.kirc.graph.prob, file = "~/tcga_biolinks1/stats/out.rand.kirc.graph.prob")
saveRDS(out.kirc.prob, file = "~/tcga_biolinks1/stats/out.kirc.prob")

# Create histogram of column means
# hist(colMeans(out.kirc.prob))
# hist(colMeans(out.rand.kirc.graph.prob))

# Turn inifnite into NAs
out.kirc.prob[which(is.infinite(out.kirc.prob))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirc.graph.prob, na.rm = T),0.01) -> random.cutoff.prob.kirc

dat.kirc.prob <- data.frame(Distance = c(colMeans(out.kirc.prob), colMeans(out.rand.kirc.graph.prob, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.kirc.prob)),rep("Random",nrow(out.rand.kirc.graph.prob)))))

# Perform KS test on column means
ks.test(colMeans(out.kirc.prob), colMeans(out.rand.kirc.graph.prob, na.rm = T)) -> ks.result.prob.kirc

# Print
ks.result.prob.kirc$statistic
ks.result.prob.kirc$p.value
# length(ks.result2.kirc$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.prob.kirc$statistic,2), ifelse(ks.result.prob.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.kirc$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.prob.kirc$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.prob.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.kirc$p.value,4)))) -> text.annot2
length(ks.result.prob.kirc$data$x) -> text.annot3

dim(out.kirc.prob)
dim(out.rand.kirc.graph.prob)

# Create gg plot
ggplot(dat.kirc.prob, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirc.prob, aes(xintercept=random.cutoff.prob.kirc,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirc Genes on the kirc Network") -> gg.kirc.prob

gg.kirc.prob
ggsave("~/tcga_biolinks1/gg.kirc.prob.png", plot = gg.kirc.prob, width = 10, height = 10)

##### POSSIBLY DAMAGING =<0.908 >0.466 ##############################

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean <=0.908) %>%
  filter(PolyPhen.Mean >0.466)-> kirc.graph.df.filt.poss

unique(kirc.graph.df.filt.poss$name) -> unique.genes.kirc.poss

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirc.poss <- foreach(i=1:length(unique.genes.kirc.poss), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirc.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = unique.genes.kirc.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirc.poss,])}

colnames(out.kirc.poss) <- unique.genes.kirc.poss

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirc.graph),"name"), length(unique.genes.kirc.poss)) -> all.genes.rand.kirc.poss

out.rand.kirc.graph.poss <- foreach(i=1:length(all.genes.rand.kirc.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = all.genes.rand.kirc.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirc.poss,])}

# # Save results of distances
saveRDS(out.rand.kirc.graph.poss, file = "~/tcga_biolinks1/stats/out.rand.kirc.graph.poss")
saveRDS(out.kirc.poss, file = "~/tcga_biolinks1/stats/out.kirc.poss")

# Create histogram of column means
# hist(colMeans(out.kirc.poss))
# hist(colMeans(out.rand.kirc.graph.poss))

# Turn inifnite into NAs
out.kirc.poss[which(is.infinite(out.kirc.poss))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirc.graph.poss, na.rm = T),0.01) -> random.cutoff.poss.kirc

dat.kirc.poss <- data.frame(Distance = c(colMeans(out.kirc.poss), colMeans(out.rand.kirc.graph.poss, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.kirc.poss)),rep("Random",nrow(out.rand.kirc.graph.poss)))))

# Perform KS test on column means
ks.test(colMeans(out.kirc.poss), colMeans(out.rand.kirc.graph.poss, na.rm = T)) -> ks.result.poss.kirc

# Print
ks.result.poss.kirc$statistic
ks.result.poss.kirc$p.value
# length(ks.result2.kirc$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.poss.kirc$statistic,2), ifelse(ks.result.poss.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.kirc$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.poss.kirc$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.poss.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.kirc$p.value,4)))) -> text.annot2
length(ks.result.poss.kirc$data$x) -> text.annot3

dim(out.kirc.poss)
dim(out.rand.kirc.graph.poss)

# Create gg plot
ggplot(dat.kirc.poss, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirc.poss, aes(xintercept=random.cutoff.poss.kirc,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirc Genes on the kirc Network") -> gg.kirc.poss

gg.kirc.poss
ggsave("~/tcga_biolinks1/gg.kirc.poss.png", plot = gg.kirc.poss, width = 10, height = 10)

##### ALL DAMAGING >0.466 ##########################################################################

kirc.graph.df$PolyPhen.Mean
kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> kirc.graph.df.filt.dam

unique(kirc.graph.df.filt.dam$name) -> unique.genes.kirc.dam

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirc.dam <- foreach(i=1:length(unique.genes.kirc.dam), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirc.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = unique.genes.kirc.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirc.dam,])}

colnames(out.kirc.dam) <- unique.genes.kirc.dam

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirc.graph),"name"), length(unique.genes.kirc.dam)) -> all.genes.rand.kirc.dam

out.rand.kirc.graph.dam <- foreach(i=1:length(all.genes.rand.kirc.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirc.graph,
    to = all.genes.rand.kirc.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirc.dam,])}

# # Save results of distances
saveRDS(out.rand.kirc.graph.dam, file = "~/tcga_biolinks1/stats/out.rand.kirc.graph.dam")
saveRDS(out.kirc.dam, file = "~/tcga_biolinks1/stats/out.kirc.dam")

# Create histogram of column means
# hist(colMeans(out.kirc.dam))
# hist(colMeans(out.rand.kirc.graph.dam))

# Turn inifnite into NAs
out.kirc.dam[which(is.infinite(out.kirc.dam))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirc.graph.dam, na.rm = T),0.01) -> random.cutoff.dam.kirc

dat.kirc.dam <- data.frame(Distance = c(colMeans(out.kirc.dam), colMeans(out.rand.kirc.graph.dam, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.kirc.dam)),rep("Random",nrow(out.rand.kirc.graph.dam)))))

# Perform KS test on column means
ks.test(colMeans(out.kirc.dam), colMeans(out.rand.kirc.graph.dam, na.rm = T)) -> ks.result.dam.kirc

# Print
ks.result.dam.kirc$statistic
ks.result.dam.kirc$p.value
# length(ks.result2.kirc$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.kirc$statistic,2), ifelse(ks.result.dam.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.kirc$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.kirc$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.kirc$p.value,4)))) -> text.annot2
length(ks.result.dam.kirc$data$x) -> text.annot3

dim(out.kirc.dam)
dim(out.rand.kirc.graph.dam)

# Create gg plot
ggplot(dat.kirc.dam, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirc.dam, aes(xintercept=random.cutoff.dam.kirc,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirc Genes on the kirc Network") -> gg.kirc.dam

gg.kirc.dam
ggsave("~/tcga_biolinks1/Plots/gg.kirc.dam.png", plot = gg.kirc.dam, width = 10, height = 10)

