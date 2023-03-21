
## Read back into R
# library(igraph)
# read.graph("~/tcga_biolinks1/GRAPHML/luad4.final.graphml", format = "graphml") -> luad.graph
# saveRDS(luad.graph, file = "~/tcga_biolinks1/RDS/luad.graph")

luad.graph <- readRDS(file = "~/tcga_biolinks1/RDS/luad.graph")

# Read table into R
read.csv(file = "~/tcga_biolinks1/GRAPHML/luad4.node.csv") -> luad.graph.df
luad.graph.df$Percent.of.Patients.with.a.Missense.Mutation

##### 1.5% PATIENTS ##############################################################

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 1.5) -> luad.graph.df.filt1.5

unique(luad.graph.df.filt1.5$name) -> unique.genes.luad1.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad1.5 <- foreach(i=1:length(unique.genes.luad1.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad1.5,])}

colnames(out.luad1.5) <- unique.genes.luad1.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad1.5)) -> all.genes.rand.luad1.5

out.rand.luad.graph1.5 <- foreach(i=1:length(all.genes.rand.luad1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad1.5,])}

# # Save results of distances
saveRDS(out.rand.luad.graph1.5, file = "~/tcga_biolinks1/stats/out.rand.luad.graph1.5")
saveRDS(out.luad1.5, file = "~/tcga_biolinks1/stats/out.luad1.5")


# out.rand.luad.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph2")
# out.luad2 <- readRDS( file = "~/tcga_biolinks1/stats/out.luad2")

# Create histogram of column means
# hist(colMeans(out.luad2))
# hist(colMeans(out.rand.luad.graph2))

# Turn inifnite into NAs
out.luad1.5[which(is.infinite(out.luad1.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph1.5, na.rm = T),0.01) -> random.cutoff1.5.luad

dat.luad1.5 <- data.frame(Distance = c(colMeans(out.luad1.5), colMeans(out.rand.luad.graph1.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.luad1.5)),rep("Random",nrow(out.rand.luad.graph1.5)))))

# Perform KS test on column means
ks.test(colMeans(out.luad1.5), colMeans(out.rand.luad.graph1.5, na.rm = T)) -> ks.result1.5.luad

# Print
ks.result1.5.luad$statistic
ks.result1.5.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.5.luad$statistic,2), ifelse(ks.result1.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.5.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.luad$p.value,4)))) -> text.annot2
length(ks.result1.5.luad$data$x) -> text.annot3

dim(out.luad1.5)
dim(out.rand.luad.graph1.5)

# Create gg plot
ggplot(dat.luad1.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad1.5, aes(xintercept=random.cutoff1.5.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad1.5

gg.luad1.5
ggsave("./gg.luad1.5.png", plot = gg.luad1.5, width = 10, height = 10)

#### 1% PATIENTS #############################################################################

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 1) -> luad.graph.df.filt1

unique(luad.graph.df.filt1$name) -> unique.genes.luad1

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad1 <- foreach(i=1:length(unique.genes.luad1), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad1[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad1,])}

colnames(out.luad1) <- unique.genes.luad1

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad1)) -> all.genes.rand.luad1

out.rand.luad.graph1 <- foreach(i=1:length(all.genes.rand.luad1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad1[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad1,])}

# # Save results of distances
saveRDS(out.rand.luad.graph1, file = "~/tcga_biolinks1/stats/out.rand.luad.graph1")
saveRDS(out.luad1, file = "~/tcga_biolinks1/stats/out.luad1")


# out.rand.luad.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph2")
# out.luad2 <- readRDS( file = "~/tcga_biolinks1/stats/out.luad2")

# Create histogram of column means
# hist(colMeans(out.luad2))
# hist(colMeans(out.rand.luad.graph2))

# Turn inifnite into NAs
out.luad1[which(is.infinite(out.luad1))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph1, na.rm = T),0.01) -> random.cutoff1.luad

dat.luad1 <- data.frame(Distance = c(colMeans(out.luad1), colMeans(out.rand.luad.graph1, na.rm = T)),
                        Distribution = factor(c(rep("Real",nrow(out.luad1)),rep("Random",nrow(out.rand.luad.graph1)))))

# Perform KS test on column means
ks.test(colMeans(out.luad1), colMeans(out.rand.luad.graph1, na.rm = T)) -> ks.result1.luad

# Print
ks.result1.luad$statistic
ks.result1.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.luad$statistic,2), ifelse(ks.result1.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.luad$p.value,4)))) -> text.annot2
length(ks.result1.luad$data$x) -> text.annot3

dim(out.luad1)
dim(out.rand.luad.graph1)

# Create gg plot
ggplot(dat.luad1, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad1, aes(xintercept=random.cutoff1.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad1

gg.luad1
ggsave("./gg.luad1.png", plot = gg.luad1, width = 10, height = 10)

###### 0.5% PATIENTS #######################################################################


luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 0.5) -> luad.graph.df.filt0.5

unique(luad.graph.df.filt0.5$name) -> unique.genes.luad0.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad0.5 <- foreach(i=1:length(unique.genes.luad0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad0.5,])}

colnames(out.luad0.5) <- unique.genes.luad0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad0.5)) -> all.genes.rand.luad0.5

out.rand.luad.graph0.5 <- foreach(i=1:length(all.genes.rand.luad0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad0.5,])}

# # Save results of distances
saveRDS(out.rand.luad.graph0.5, file = "~/tcga_biolinks1/stats/out.rand.luad.graph0.5")
saveRDS(out.luad0.5, file = "~/tcga_biolinks1/stats/out.luad0.5")


# out.rand.luad.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph2")
# out.luad2 <- readRDS( file = "~/tcga_biolinks1/stats/out.luad2")

# Create histogram of column means
# hist(colMeans(out.luad2))
# hist(colMeans(out.rand.luad.graph2))

# Turn inifnite into NAs
out.luad0.5[which(is.infinite(out.luad0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph0.5, na.rm = T),0.01) -> random.cutoff0.5.luad

dat.luad0.5 <- data.frame(Distance = c(colMeans(out.luad0.5), colMeans(out.rand.luad.graph0.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.luad0.5)),rep("Random",nrow(out.rand.luad.graph0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.luad0.5), colMeans(out.rand.luad.graph0.5, na.rm = T)) -> ks.result0.5.luad

# Print
ks.result0.5.luad$statistic
ks.result0.5.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result0.5.luad$statistic,2), ifelse(ks.result0.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result0.5.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result0.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.luad$p.value,4)))) -> text.annot2
length(ks.result0.5.luad$data$x) -> text.annot3

dim(out.luad0.5)
dim(out.rand.luad.graph0.5)

# Create gg plot
ggplot(dat.luad0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad0.5, aes(xintercept=random.cutoff0.5.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad0.5

gg.luad0.5
ggsave("./gg.luad0.5.png", plot = gg.luad0.5, width = 10, height = 10)

##### PROBABLY DAMAGING >0.908 ##############################

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> luad.graph.df.filt.prob

unique(luad.graph.df.filt.prob$name) -> unique.genes.luad.prob

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad.prob <- foreach(i=1:length(unique.genes.luad.prob), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad.prob,])}

colnames(out.luad.prob) <- unique.genes.luad.prob

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad.prob)) -> all.genes.rand.luad.prob

out.rand.luad.graph.prob <- foreach(i=1:length(all.genes.rand.luad.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad.prob,])}

# # Save results of distances
saveRDS(out.rand.luad.graph.prob, file = "~/tcga_biolinks1/stats/out.rand.luad.graph.prob")
saveRDS(out.luad.prob, file = "~/tcga_biolinks1/stats/out.luad.prob")

# Create histogram of column means
# hist(colMeans(out.luad.prob))
# hist(colMeans(out.rand.luad.graph.prob))

# Turn inifnite into NAs
out.luad.prob[which(is.infinite(out.luad.prob))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph.prob, na.rm = T),0.01) -> random.cutoff.prob.luad

dat.luad.prob <- data.frame(Distance = c(colMeans(out.luad.prob), colMeans(out.rand.luad.graph.prob, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.luad.prob)),rep("Random",nrow(out.rand.luad.graph.prob)))))

# Perform KS test on column means
ks.test(colMeans(out.luad.prob), colMeans(out.rand.luad.graph.prob, na.rm = T)) -> ks.result.prob.luad

# Print
ks.result.prob.luad$statistic
ks.result.prob.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.prob.luad$statistic,2), ifelse(ks.result.prob.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.prob.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.prob.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.luad$p.value,4)))) -> text.annot2
length(ks.result.prob.luad$data$x) -> text.annot3

dim(out.luad.prob)
dim(out.rand.luad.graph.prob)

# Create gg plot
ggplot(dat.luad.prob, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad.prob, aes(xintercept=random.cutoff.prob.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad.prob

gg.luad.prob
ggsave("./gg.luad.prob.png", plot = gg.luad.prob, width = 10, height = 10)

##### POSSIBLY DAMAGING =<0.908 >0.466 ##############################
luad.graph.df$PolyPhen.Mean
luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean <=0.908) %>%
  filter(PolyPhen.Mean >0.466)-> luad.graph.df.filt.poss

unique(luad.graph.df.filt.poss$name) -> unique.genes.luad.poss

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad.poss <- foreach(i=1:length(unique.genes.luad.poss), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad.poss,])}

colnames(out.luad.poss) <- unique.genes.luad.poss

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad.poss)) -> all.genes.rand.luad.poss

out.rand.luad.graph.poss <- foreach(i=1:length(all.genes.rand.luad.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad.poss,])}

# # Save results of distances
saveRDS(out.rand.luad.graph.poss, file = "~/tcga_biolinks1/stats/out.rand.luad.graph.poss")
saveRDS(out.luad.poss, file = "~/tcga_biolinks1/stats/out.luad.poss")

# Create histogram of column means
# hist(colMeans(out.luad.poss))
# hist(colMeans(out.rand.luad.graph.poss))

# Turn inifnite into NAs
out.luad.poss[which(is.infinite(out.luad.poss))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph.poss, na.rm = T),0.01) -> random.cutoff.poss.luad

dat.luad.poss <- data.frame(Distance = c(colMeans(out.luad.poss), colMeans(out.rand.luad.graph.poss, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.luad.poss)),rep("Random",nrow(out.rand.luad.graph.poss)))))

# Perform KS test on column means
ks.test(colMeans(out.luad.poss), colMeans(out.rand.luad.graph.poss, na.rm = T)) -> ks.result.poss.luad

# Print
ks.result.poss.luad$statistic
ks.result.poss.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.poss.luad$statistic,2), ifelse(ks.result.poss.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.poss.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.poss.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.luad$p.value,4)))) -> text.annot2
length(ks.result.poss.luad$data$x) -> text.annot3

dim(out.luad.poss)
dim(out.rand.luad.graph.poss)

# Create gg plot
ggplot(dat.luad.poss, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad.poss, aes(xintercept=random.cutoff.poss.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad.poss

gg.luad.poss
ggsave("./gg.luad.poss.png", plot = gg.luad.poss, width = 10, height = 10)

##### ALL DAMAGING >0.466 ##########################################################################

luad.graph.df$PolyPhen.Mean
luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> luad.graph.df.filt.dam

unique(luad.graph.df.filt.dam$name) -> unique.genes.luad.dam

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad.dam <- foreach(i=1:length(unique.genes.luad.dam), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad.dam,])}

colnames(out.luad.dam) <- unique.genes.luad.dam

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad.dam)) -> all.genes.rand.luad.dam

out.rand.luad.graph.dam <- foreach(i=1:length(all.genes.rand.luad.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad.dam,])}

# # Save results of distances
saveRDS(out.rand.luad.graph.dam, file = "~/tcga_biolinks1/stats/out.rand.luad.graph.dam")
saveRDS(out.luad.dam, file = "~/tcga_biolinks1/stats/out.luad.dam")

# Create histogram of column means
# hist(colMeans(out.luad.dam))
# hist(colMeans(out.rand.luad.graph.dam))

# Turn inifnite into NAs
out.luad.dam[which(is.infinite(out.luad.dam))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph.dam, na.rm = T),0.01) -> random.cutoff.dam.luad

dat.luad.dam <- data.frame(Distance = c(colMeans(out.luad.dam), colMeans(out.rand.luad.graph.dam, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.luad.dam)),rep("Random",nrow(out.rand.luad.graph.dam)))))

# Perform KS test on column means
ks.test(colMeans(out.luad.dam), colMeans(out.rand.luad.graph.dam, na.rm = T)) -> ks.result.dam.luad

# Print
ks.result.dam.luad$statistic
ks.result.dam.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.luad$statistic,2), ifelse(ks.result.dam.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.luad$p.value,4)))) -> text.annot2
length(ks.result.dam.luad$data$x) -> text.annot3

dim(out.luad.dam)
dim(out.rand.luad.graph.dam)

# Create gg plot
ggplot(dat.luad.dam, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad.dam, aes(xintercept=random.cutoff.dam.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad.dam

gg.luad.dam
ggsave("./gg.luad.dam.png", plot = gg.luad.dam, width = 10, height = 10)

##### ALL DAMAGING >0.466 and %>0.5 ##########################################################################

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 0.5) %>%
  filter(PolyPhen.Mean >0.466)-> luad.graph.df.filt.dam.0.5

unique(luad.graph.df.filt.dam.0.5$name) -> unique.genes.luad.dam.0.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad.dam.0.5 <- foreach(i=1:length(unique.genes.luad.dam.0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad.dam.0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad.dam.0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad.dam.0.5,])}

colnames(out.luad.dam.0.5) <- unique.genes.luad.dam.0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad.dam.0.5)) -> all.genes.rand.luad.dam.0.5

out.rand.luad.graph.dam.0.5 <- foreach(i=1:length(all.genes.rand.luad.dam.0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad.dam.0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad.dam.0.5,])}

# # Save results of distances
saveRDS(out.rand.luad.graph.dam.0.5, file = "~/tcga_biolinks1/stats/out.rand.luad.graph.dam.0.5")
saveRDS(out.luad.dam.0.5, file = "~/tcga_biolinks1/stats/out.luad.dam.0.5")
# readRDS(file = "~/tcga_biolinks1/stats/out.luad.dam.0.5") -> out.luad.dam.0.5
# readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph.dam.0.5") -> out.rand.luad.graph.dam.0.5

# Create histogram of column means
# hist(colMeans(out.luad.dam.0.5))
# hist(colMeans(out.rand.luad.graph.dam.0.5))

# Turn inifnite into NAs
out.luad.dam.0.5[which(is.infinite(out.luad.dam.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph.dam.0.5, na.rm = T),0.01) -> random.cutoff.dam.0.5.luad

dat.luad.dam.0.5 <- data.frame(Distance = c(colMeans(out.luad.dam.0.5), colMeans(out.rand.luad.graph.dam.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.luad.dam.0.5)),rep("Random",nrow(out.rand.luad.graph.dam.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.luad.dam.0.5), colMeans(out.rand.luad.graph.dam.0.5, na.rm = T)) -> ks.result.dam.0.5.luad

# Print
ks.result.dam.0.5.luad$statistic
ks.result.dam.0.5.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.0.5.luad$statistic,2), ifelse(ks.result.dam.0.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.0.5.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.0.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.luad$p.value,4)))) -> text.annot2
length(ks.result.dam.0.5.luad$data$x) -> text.annot3

dim(out.luad.dam.0.5)
dim(out.rand.luad.graph.dam.0.5)

# Create gg plot
ggplot(dat.luad.dam.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad.dam.0.5, aes(xintercept=random.cutoff.dam.0.5.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad.dam.0.5

gg.luad.dam.0.5
ggsave("~/tcga_biolinks1/Plots/gg.luad.dam.0.5.png", plot = gg.luad.dam.0.5, width = 10, height = 10)

##### ALL DAMAGING >0.466 and %>1.5 ##########################################################################

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 1.5) %>%
  filter(PolyPhen.Mean >0.466)-> luad.graph.df.filt.dam.1.5

unique(luad.graph.df.filt.dam.1.5$name) -> unique.genes.luad.dam.1.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad.dam.1.5 <- foreach(i=1:length(unique.genes.luad.dam.1.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad.dam.1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad.dam.1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad.dam.1.5,])}

colnames(out.luad.dam.1.5) <- unique.genes.luad.dam.1.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad.dam.1.5)) -> all.genes.rand.luad.dam.1.5

out.rand.luad.graph.dam.1.5 <- foreach(i=1:length(all.genes.rand.luad.dam.1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad.dam.1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad.dam.1.5,])}

# # Save results of distances
saveRDS(out.rand.luad.graph.dam.1.5, file = "~/tcga_biolinks1/stats/out.rand.luad.graph.dam.1.5")
saveRDS(out.luad.dam.1.5, file = "~/tcga_biolinks1/stats/out.luad.dam.1.5")
# readRDS(file = "~/tcga_biolinks1/stats/out.luad.dam.1.5") -> out.luad.dam.1.5
# readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph.dam.1.5") -> out.rand.luad.graph.dam.1.5

# Create histogram of column means
# hist(colMeans(out.luad.dam.1.5))
# hist(colMeans(out.rand.luad.graph.dam.1.5))

# Turn inifnite into NAs
out.luad.dam.1.5[which(is.infinite(out.luad.dam.1.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph.dam.1.5, na.rm = T),0.01) -> random.cutoff.dam.1.5.luad

dat.luad.dam.1.5 <- data.frame(Distance = c(colMeans(out.luad.dam.1.5), colMeans(out.rand.luad.graph.dam.1.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.luad.dam.1.5)),rep("Random",nrow(out.rand.luad.graph.dam.1.5)))))

# Perform KS test on column means
ks.test(colMeans(out.luad.dam.1.5), colMeans(out.rand.luad.graph.dam.1.5, na.rm = T)) -> ks.result.dam.1.5.luad

# Print
ks.result.dam.1.5.luad$statistic
ks.result.dam.1.5.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.1.5.luad$statistic,2), ifelse(ks.result.dam.1.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.1.5.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.1.5.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.1.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.1.5.luad$p.value,4)))) -> text.annot2
length(ks.result.dam.1.5.luad$data$x) -> text.annot3

dim(out.luad.dam.1.5)
dim(out.rand.luad.graph.dam.1.5)

# Create gg plot
ggplot(dat.luad.dam.1.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad.dam.1.5, aes(xintercept=random.cutoff.dam.1.5.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad.dam.1.5

gg.luad.dam.1.5
ggsave("~/tcga_biolinks1/Plots/gg.luad.dam.1.5.png", plot = gg.luad.dam.1.5, width = 10, height = 10)

######

#### DELETERIOUS <0.05 ##########################################################################

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(SIFT.Mean <0.05)-> luad.graph.df.filt.del

unique(luad.graph.df.filt.del$name) -> unique.genes.luad.del

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad.del <- foreach(i=1:length(unique.genes.luad.del), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad.del), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad.del[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad.del,])}

colnames(out.luad.del) <- unique.genes.luad.del

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad.del)) -> all.genes.rand.luad.del

out.rand.luad.graph.del <- foreach(i=1:length(all.genes.rand.luad.del), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad.del[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad.del,])}

# # Save results of distances
saveRDS(out.rand.luad.graph.del, file = "~/tcga_biolinks1/stats/out.rand.luad.graph.del")
saveRDS(out.luad.del, file = "~/tcga_biolinks1/stats/out.luad.del")

# Create histogram of column means
# hist(colMeans(out.luad.del))
# hist(colMeans(out.rand.luad.graph.del))

# Turn inifnite into NAs
out.luad.del[which(is.infinite(out.luad.del))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph.del, na.rm = T),0.01) -> random.cutoff.del.luad

dat.luad.del <- data.frame(Distance = c(colMeans(out.luad.del), colMeans(out.rand.luad.graph.del, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.luad.del)),rep("Random",nrow(out.rand.luad.graph.del)))))

# Perform KS test on column means
ks.test(colMeans(out.luad.del), colMeans(out.rand.luad.graph.del, na.rm = T)) -> ks.result.del.luad

# Print
ks.result.del.luad$statistic
ks.result.del.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.luad$statistic,2), ifelse(ks.result.del.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.luad$p.value,4)))) -> text.annot2
length(ks.result.del.luad$data$x) -> text.annot3

dim(out.luad.del)
dim(out.rand.luad.graph.del)

# Create gg plot
ggplot(dat.luad.del, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad.del, aes(xintercept=random.cutoff.del.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad.del

gg.luad.del
ggsave("./gg.luad.del.png", plot = gg.luad.del, width = 10, height = 10)

##### ALL DELETERIOUS <0.05 and %>0.5 ##########################################################################

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 0.5) %>%
  filter(SIFT.Mean <0.05)-> luad.graph.df.filt.del.0.5

unique(luad.graph.df.filt.del.0.5$name) -> unique.genes.luad.del.0.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.luad.del.0.5 <- foreach(i=1:length(unique.genes.luad.del.0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.luad.del.0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = unique.genes.luad.del.0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.luad.del.0.5,])}

colnames(out.luad.del.0.5) <- unique.genes.luad.del.0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(luad.graph),"name"), length(unique.genes.luad.del.0.5)) -> all.genes.rand.luad.del.0.5

out.rand.luad.graph.del.0.5 <- foreach(i=1:length(all.genes.rand.luad.del.0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    luad.graph,
    to = all.genes.rand.luad.del.0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.luad.del.0.5,])}

# # Save results of distances
saveRDS(out.rand.luad.graph.del.0.5, file = "~/tcga_biolinks1/stats/out.rand.luad.graph.del.0.5")
saveRDS(out.luad.del.0.5, file = "~/tcga_biolinks1/stats/out.luad.del.0.5")
# readRDS(file = "~/tcga_biolinks1/stats/out.luad.del.0.5") -> out.luad.del.0.5
# readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph.del.0.5") -> out.rand.luad.graph.del.0.5

# Create histogram of column means
# hist(colMeans(out.luad.del.0.5))
# hist(colMeans(out.rand.luad.graph.del.0.5))

# Turn inifnite into NAs
out.luad.del.0.5[which(is.infinite(out.luad.del.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.luad.graph.del.0.5, na.rm = T),0.01) -> random.cutoff.del.0.5.luad

dat.luad.del.0.5 <- data.frame(Distance = c(colMeans(out.luad.del.0.5), colMeans(out.rand.luad.graph.del.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.luad.del.0.5)),rep("Random",nrow(out.rand.luad.graph.del.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.luad.del.0.5), colMeans(out.rand.luad.graph.del.0.5, na.rm = T)) -> ks.result.del.0.5.luad

# Print
ks.result.del.0.5.luad$statistic
ks.result.del.0.5.luad$p.value
# length(ks.result2.luad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.0.5.luad$statistic,2), ifelse(ks.result.del.0.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.luad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.0.5.luad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.0.5.luad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.luad$p.value,4)))) -> text.annot2
length(ks.result.del.0.5.luad$data$x) -> text.annot3

dim(out.luad.del.0.5)
dim(out.rand.luad.graph.del.0.5)

# Create gg plot
ggplot(dat.luad.del.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.luad.del.0.5, aes(xintercept=random.cutoff.del.0.5.luad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.luad.del.0.5

gg.luad.del.0.5
ggsave("~/tcga_biolinks1/Plots/gg.luad.del.0.5.png", plot = gg.luad.del.0.5, width = 10, height = 10)

#######################################

## Other method for saving
# pdf(file = "./test.pdf")
# gg.luad2
# dev.off()
# png(file = "./test.png")
# gg.luad2
# dev.off()
#### close all the connections
# graphics.off()

