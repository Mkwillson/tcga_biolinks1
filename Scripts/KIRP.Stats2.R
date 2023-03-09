
## Read back into R
# library(igraph)
# read.graph("~/tcga_biolinks1/GRAPHML/kirp5.final.graphml", format = "graphml") -> kirp.graph
# saveRDS(kirp.graph, file = "~/tcga_biolinks1/RDS/kirp.graph")

kirp.graph <- readRDS(file = "~/tcga_biolinks1/RDS/kirp.graph")

# Read table into R
read.csv(file = "~/tcga_biolinks1/GRAPHML/kirp5.node.csv") -> kirp.graph.df
kirp.graph.df$Percent.of.Patients.with.a.Missense.Mutation

##### 1.5% PATIENTS ##############################################################

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutations >= 1.5) -> kirp.graph.df.filt1.5

unique(kirp.graph.df.filt1.5$name) -> unique.genes.kirp1.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirp1.5 <- foreach(i=1:length(unique.genes.kirp1.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirp2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = unique.genes.kirp1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirp1.5,])}

colnames(out.kirp1.5) <- unique.genes.kirp1.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirp.graph),"name"), length(unique.genes.kirp1.5)) -> all.genes.rand.kirp1.5

out.rand.kirp.graph1.5 <- foreach(i=1:length(all.genes.rand.kirp1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = all.genes.rand.kirp1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirp1.5,])}

# # Save results of distances
saveRDS(out.rand.kirp.graph1.5, file = "~/tcga_biolinks1/stats/out.rand.kirp.graph1.5")
saveRDS(out.kirp1.5, file = "~/tcga_biolinks1/stats/out.kirp1.5")


# out.rand.kirp.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph2")
# out.kirp2 <- readRDS( file = "~/tcga_biolinks1/stats/out.kirp2")

# Create histogram of column means
# hist(colMeans(out.kirp2))
# hist(colMeans(out.rand.kirp.graph2))

# Turn inifnite into NAs
out.kirp1.5[which(is.infinite(out.kirp1.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirp.graph1.5, na.rm = T),0.01) -> random.cutoff1.5.kirp

dat.kirp1.5 <- data.frame(Distance = c(colMeans(out.kirp1.5), colMeans(out.rand.kirp.graph1.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.kirp1.5)),rep("Random",nrow(out.rand.kirp.graph1.5)))))

# Perform KS test on column means
ks.test(colMeans(out.kirp1.5), colMeans(out.rand.kirp.graph1.5, na.rm = T)) -> ks.result1.5.kirp

# Print
ks.result1.5.kirp$statistic
ks.result1.5.kirp$p.value
# length(ks.result2.kirp$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.5.kirp$statistic,2), ifelse(ks.result1.5.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.kirp$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.5.kirp$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.5.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.kirp$p.value,4)))) -> text.annot2
length(ks.result1.5.kirp$data$x) -> text.annot3

dim(out.kirp1.5)
dim(out.rand.kirp.graph1.5)

# Create gg plot
ggplot(dat.kirp1.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirp1.5, aes(xintercept=random.cutoff1.5.kirp,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirp Genes on the kirp Network") -> gg.kirp1.5

gg.kirp1.5
ggsave("~/tcga_biolinks1/gg.kirp1.5.png", plot = gg.kirp1.5, width = 10, height = 10)

#### 1% PATIENTS #############################################################################

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutations >= 1) -> kirp.graph.df.filt1

unique(kirp.graph.df.filt1$name) -> unique.genes.kirp1

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirp1 <- foreach(i=1:length(unique.genes.kirp1), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirp1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = unique.genes.kirp1[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirp1,])}

colnames(out.kirp1) <- unique.genes.kirp1

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirp.graph),"name"), length(unique.genes.kirp1)) -> all.genes.rand.kirp1

out.rand.kirp.graph1 <- foreach(i=1:length(all.genes.rand.kirp1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = all.genes.rand.kirp1[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirp1,])}

# # Save results of distances
saveRDS(out.rand.kirp.graph1, file = "~/tcga_biolinks1/stats/out.rand.kirp.graph1")
saveRDS(out.kirp1, file = "~/tcga_biolinks1/stats/out.kirp1")


# out.rand.kirp.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph2")
# out.kirp2 <- readRDS( file = "~/tcga_biolinks1/stats/out.kirp2")

# Create histogram of column means
# hist(colMeans(out.kirp2))
# hist(colMeans(out.rand.kirp.graph2))

# Turn inifnite into NAs
out.kirp1[which(is.infinite(out.kirp1))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirp.graph1, na.rm = T),0.01) -> random.cutoff1.kirp

dat.kirp1 <- data.frame(Distance = c(colMeans(out.kirp1), colMeans(out.rand.kirp.graph1, na.rm = T)),
                        Distribution = factor(c(rep("Real",nrow(out.kirp1)),rep("Random",nrow(out.rand.kirp.graph1)))))

# Perform KS test on column means
ks.test(colMeans(out.kirp1), colMeans(out.rand.kirp.graph1, na.rm = T)) -> ks.result1.kirp

# Print
ks.result1.kirp$statistic
ks.result1.kirp$p.value
# length(ks.result2.kirp$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.kirp$statistic,2), ifelse(ks.result1.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.kirp$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.kirp$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.kirp$p.value,4)))) -> text.annot2
length(ks.result1.kirp$data$x) -> text.annot3

dim(out.kirp1)
dim(out.rand.kirp.graph1)

# Create gg plot
ggplot(dat.kirp1, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirp1, aes(xintercept=random.cutoff1.kirp,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirp Genes on the kirp Network") -> gg.kirp1

gg.kirp1
ggsave("./gg.kirp1.png", plot = gg.kirp1, width = 10, height = 10)

###### 0.5% PATIENTS #######################################################################

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutations >= 0.5) -> kirp.graph.df.filt0.5

unique(kirp.graph.df.filt0.5$name) -> unique.genes.kirp0.5
# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirp0.5 <- foreach(i=1:length(unique.genes.kirp0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirp2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = unique.genes.kirp0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirp0.5,])}

colnames(out.kirp0.5) <- unique.genes.kirp0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirp.graph),"name"), length(unique.genes.kirp0.5)) -> all.genes.rand.kirp0.5

out.rand.kirp.graph0.5 <- foreach(i=1:length(all.genes.rand.kirp0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = all.genes.rand.kirp0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirp0.5,])}

# # Save results of distances
saveRDS(out.rand.kirp.graph0.5, file = "~/tcga_biolinks1/stats/out.rand.kirp.graph0.5")
saveRDS(out.kirp0.5, file = "~/tcga_biolinks1/stats/out.kirp0.5")


# out.rand.kirp.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph2")
# out.kirp2 <- readRDS( file = "~/tcga_biolinks1/stats/out.kirp2")

# Create histogram of column means
# hist(colMeans(out.kirp2))
# hist(colMeans(out.rand.kirp.graph2))

# Turn inifnite into NAs
out.kirp0.5[which(is.infinite(out.kirp0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirp.graph0.5, na.rm = T),0.01) -> random.cutoff0.5.kirp

dat.kirp0.5 <- data.frame(Distance = c(colMeans(out.kirp0.5), colMeans(out.rand.kirp.graph0.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.kirp0.5)),rep("Random",nrow(out.rand.kirp.graph0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.kirp0.5), colMeans(out.rand.kirp.graph0.5, na.rm = T)) -> ks.result0.5.kirp

# Print
ks.result0.5.kirp$statistic
ks.result0.5.kirp$p.value
# length(ks.result2.kirp$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result0.5.kirp$statistic,2), ifelse(ks.result0.5.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.kirp$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result0.5.kirp$statistic,2)) -> text.annot1
paste0(ifelse(ks.result0.5.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.kirp$p.value,4)))) -> text.annot2
length(ks.result0.5.kirp$data$x) -> text.annot3

dim(out.kirp0.5)
dim(out.rand.kirp.graph0.5)

# Create gg plot
ggplot(dat.kirp0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirp0.5, aes(xintercept=random.cutoff0.5.kirp,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirp Genes on the kirp Network") -> gg.kirp0.5

gg.kirp0.5
ggsave("~/tcga_biolinks1/gg.kirp0.5.png", plot = gg.kirp0.5, width = 10, height = 10)

##### PROBABLY DAMAGING >0.908 ##############################

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> kirp.graph.df.filt.prob

unique(kirp.graph.df.filt.prob$name) -> unique.genes.kirp.prob

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirp.prob <- foreach(i=1:length(unique.genes.kirp.prob), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirp.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = unique.genes.kirp.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirp.prob,])}

colnames(out.kirp.prob) <- unique.genes.kirp.prob

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirp.graph),"name"), length(unique.genes.kirp.prob)) -> all.genes.rand.kirp.prob

out.rand.kirp.graph.prob <- foreach(i=1:length(all.genes.rand.kirp.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = all.genes.rand.kirp.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirp.prob,])}

# # Save results of distances
saveRDS(out.rand.kirp.graph.prob, file = "~/tcga_biolinks1/stats/out.rand.kirp.graph.prob")
saveRDS(out.kirp.prob, file = "~/tcga_biolinks1/stats/out.kirp.prob")

# Create histogram of column means
# hist(colMeans(out.kirp.prob))
# hist(colMeans(out.rand.kirp.graph.prob))

# Turn inifnite into NAs
out.kirp.prob[which(is.infinite(out.kirp.prob))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirp.graph.prob, na.rm = T),0.01) -> random.cutoff.prob.kirp

dat.kirp.prob <- data.frame(Distance = c(colMeans(out.kirp.prob), colMeans(out.rand.kirp.graph.prob, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.kirp.prob)),rep("Random",nrow(out.rand.kirp.graph.prob)))))

# Perform KS test on column means
ks.test(colMeans(out.kirp.prob), colMeans(out.rand.kirp.graph.prob, na.rm = T)) -> ks.result.prob.kirp

# Print
ks.result.prob.kirp$statistic
ks.result.prob.kirp$p.value
# length(ks.result2.kirp$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.prob.kirp$statistic,2), ifelse(ks.result.prob.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.kirp$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.prob.kirp$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.prob.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.kirp$p.value,4)))) -> text.annot2
length(ks.result.prob.kirp$data$x) -> text.annot3

dim(out.kirp.prob)
dim(out.rand.kirp.graph.prob)

# Create gg plot
ggplot(dat.kirp.prob, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirp.prob, aes(xintercept=random.cutoff.prob.kirp,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirp Genes on the kirp Network") -> gg.kirp.prob

gg.kirp.prob
ggsave("~/tcga_biolinks1/gg.kirp.prob.png", plot = gg.kirp.prob, width = 10, height = 10)
ggsave("~/tcga_biolinks1/gg.kirp0.5.png", plot = gg.kirp0.5, width = 10, height = 10)

##### POSSIBLY DAMAGING =<0.908 >0.466 ##############################

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean <=0.908) %>%
  filter(PolyPhen.Mean >0.466)-> kirp.graph.df.filt.poss

unique(kirp.graph.df.filt.poss$name) -> unique.genes.kirp.poss

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirp.poss <- foreach(i=1:length(unique.genes.kirp.poss), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirp.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = unique.genes.kirp.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirp.poss,])}

colnames(out.kirp.poss) <- unique.genes.kirp.poss

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirp.graph),"name"), length(unique.genes.kirp.poss)) -> all.genes.rand.kirp.poss

out.rand.kirp.graph.poss <- foreach(i=1:length(all.genes.rand.kirp.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = all.genes.rand.kirp.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirp.poss,])}

# # Save results of distances
saveRDS(out.rand.kirp.graph.poss, file = "~/tcga_biolinks1/stats/out.rand.kirp.graph.poss")
saveRDS(out.kirp.poss, file = "~/tcga_biolinks1/stats/out.kirp.poss")

# Create histogram of column means
# hist(colMeans(out.kirp.poss))
# hist(colMeans(out.rand.kirp.graph.poss))

# Turn inifnite into NAs
out.kirp.poss[which(is.infinite(out.kirp.poss))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirp.graph.poss, na.rm = T),0.01) -> random.cutoff.poss.kirp

dat.kirp.poss <- data.frame(Distance = c(colMeans(out.kirp.poss), colMeans(out.rand.kirp.graph.poss, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.kirp.poss)),rep("Random",nrow(out.rand.kirp.graph.poss)))))

# Perform KS test on column means
ks.test(colMeans(out.kirp.poss), colMeans(out.rand.kirp.graph.poss, na.rm = T)) -> ks.result.poss.kirp

# Print
ks.result.poss.kirp$statistic
ks.result.poss.kirp$p.value
# length(ks.result2.kirp$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.poss.kirp$statistic,2), ifelse(ks.result.poss.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.kirp$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.poss.kirp$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.poss.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.kirp$p.value,4)))) -> text.annot2
length(ks.result.poss.kirp$data$x) -> text.annot3

dim(out.kirp.poss)
dim(out.rand.kirp.graph.poss)

# Create gg plot
ggplot(dat.kirp.poss, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirp.poss, aes(xintercept=random.cutoff.poss.kirp,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirp Genes on the kirp Network") -> gg.kirp.poss

gg.kirp.poss
ggsave("~/tcga_biolinks1/gg.kirp.poss.png", plot = gg.kirp.poss, width = 10, height = 10)

##### ALL DAMAGING >0.466 ##########################################################################

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> kirp.graph.df.filt.dam

unique(kirp.graph.df.filt.dam$name) -> unique.genes.kirp.dam

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.kirp.dam <- foreach(i=1:length(unique.genes.kirp.dam), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.kirp.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = unique.genes.kirp.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.kirp.dam,])}

colnames(out.kirp.dam) <- unique.genes.kirp.dam

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(kirp.graph),"name"), length(unique.genes.kirp.dam)) -> all.genes.rand.kirp.dam

out.rand.kirp.graph.dam <- foreach(i=1:length(all.genes.rand.kirp.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    kirp.graph,
    to = all.genes.rand.kirp.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.kirp.dam,])}

# # Save results of distances
saveRDS(out.rand.kirp.graph.dam, file = "~/tcga_biolinks1/stats/out.rand.kirp.graph.dam")
saveRDS(out.kirp.dam, file = "~/tcga_biolinks1/stats/out.kirp.dam")

# Create histogram of column means
# hist(colMeans(out.kirp.dam))
# hist(colMeans(out.rand.kirp.graph.dam))

# Turn inifnite into NAs
out.kirp.dam[which(is.infinite(out.kirp.dam))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirp.graph.dam, na.rm = T),0.01) -> random.cutoff.dam.kirp

dat.kirp.dam <- data.frame(Distance = c(colMeans(out.kirp.dam), colMeans(out.rand.kirp.graph.dam, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.kirp.dam)),rep("Random",nrow(out.rand.kirp.graph.dam)))))

# Perform KS test on column means
ks.test(colMeans(out.kirp.dam), colMeans(out.rand.kirp.graph.dam, na.rm = T)) -> ks.result.dam.kirp

# Print
ks.result.dam.kirp$statistic
ks.result.dam.kirp$p.value
# length(ks.result2.kirp$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.kirp$statistic,2), ifelse(ks.result.dam.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.kirp$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.kirp$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.kirp$p.value,4)))) -> text.annot2
length(ks.result.dam.kirp$data$x) -> text.annot3

dim(out.kirp.dam)
dim(out.rand.kirp.graph.dam)

# Create gg plot
ggplot(dat.kirp.dam, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirp.dam, aes(xintercept=random.cutoff.dam.kirp,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and kirp Genes on the kirp Network") -> gg.kirp.dam

gg.kirp.dam
ggsave("~/tcga_biolinks1/gg.kirp.dam.png", plot = gg.kirp.dam, width = 10, height = 10)

