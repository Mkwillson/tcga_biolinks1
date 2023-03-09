## Read back into R
# library(igraph)
# read.graph("~/tcga_biolinks1/GRAPHML/stad5.final.graphml", format = "graphml") -> stad.graph
# saveRDS(stad.graph, file = "~/tcga_biolinks1/RDS/stad.graph")

stad.graph <- readRDS(file = "~/tcga_biolinks1/RDS/stad.graph")

# Read table into R
read.csv(file = "~/tcga_biolinks1/GRAPHML/stad5.node.csv") -> stad.graph.df
stad.graph.df$Percent.of.Patients.with.a.Missense.Mutation

##### 1.5% PATIENTS ##############################################################

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1.5) -> stad.graph.df.filt1.5

unique(stad.graph.df.filt1.5$name) -> unique.genes.stad1.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.stad1.5 <- foreach(i=1:length(unique.genes.stad1.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.stad2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = unique.genes.stad1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.stad1.5,])}

colnames(out.stad1.5) <- unique.genes.stad1.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(stad.graph),"name"), length(unique.genes.stad1.5)) -> all.genes.rand.stad1.5

out.rand.stad.graph1.5 <- foreach(i=1:length(all.genes.rand.stad1.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = all.genes.rand.stad1.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.stad1.5,])}

# # Save results of distances
saveRDS(out.rand.stad.graph1.5, file = "~/tcga_biolinks1/stats/out.rand.stad.graph1.5")
saveRDS(out.stad1.5, file = "~/tcga_biolinks1/stats/out.stad1.5")


# out.rand.stad.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph2")
# out.stad2 <- readRDS( file = "~/tcga_biolinks1/stats/out.stad2")

# Create histogram of column means
# hist(colMeans(out.stad2))
# hist(colMeans(out.rand.stad.graph2))

# Turn inifnite into NAs
out.stad1.5[which(is.infinite(out.stad1.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph1.5, na.rm = T),0.01) -> random.cutoff1.5.stad

dat.stad1.5 <- data.frame(Distance = c(colMeans(out.stad1.5), colMeans(out.rand.stad.graph1.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.stad1.5)),rep("Random",nrow(out.rand.stad.graph1.5)))))

# Perform KS test on column means
ks.test(colMeans(out.stad1.5), colMeans(out.rand.stad.graph1.5, na.rm = T)) -> ks.result1.5.stad

# Print
ks.result1.5.stad$statistic
ks.result1.5.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.5.stad$statistic,2), ifelse(ks.result1.5.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.5.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.5.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.5.stad$p.value,4)))) -> text.annot2
length(ks.result1.5.stad$data$x) -> text.annot3

dim(out.stad1.5)
dim(out.rand.stad.graph1.5)

# Create gg plot
ggplot(dat.stad1.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad1.5, aes(xintercept=random.cutoff1.5.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and stad Genes on the stad Network") -> gg.stad1.5

gg.stad1.5
ggsave("~/tcga_biolinks1/Plots/gg.stad1.5.png", plot = gg.stad1.5, width = 10, height = 10)

#### 1% PATIENTS #############################################################################

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1) -> stad.graph.df.filt1

unique(stad.graph.df.filt1$name) -> unique.genes.stad1

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.stad1 <- foreach(i=1:length(unique.genes.stad1), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.stad1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = unique.genes.stad1[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.stad1,])}

colnames(out.stad1) <- unique.genes.stad1

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(stad.graph),"name"), length(unique.genes.stad1)) -> all.genes.rand.stad1

out.rand.stad.graph1 <- foreach(i=1:length(all.genes.rand.stad1), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = all.genes.rand.stad1[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.stad1,])}

# # Save results of distances
saveRDS(out.rand.stad.graph1, file = "~/tcga_biolinks1/stats/out.rand.stad.graph1")
saveRDS(out.stad1, file = "~/tcga_biolinks1/stats/out.stad1")


# out.rand.stad.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph2")
# out.stad2 <- readRDS( file = "~/tcga_biolinks1/stats/out.stad2")

# Create histogram of column means
# hist(colMeans(out.stad2))
# hist(colMeans(out.rand.stad.graph2))

# Turn inifnite into NAs
out.stad1[which(is.infinite(out.stad1))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph1, na.rm = T),0.01) -> random.cutoff1.stad

dat.stad1 <- data.frame(Distance = c(colMeans(out.stad1), colMeans(out.rand.stad.graph1, na.rm = T)),
                        Distribution = factor(c(rep("Real",nrow(out.stad1)),rep("Random",nrow(out.rand.stad.graph1)))))

# Perform KS test on column means
ks.test(colMeans(out.stad1), colMeans(out.rand.stad.graph1, na.rm = T)) -> ks.result1.stad

# Print
ks.result1.stad$statistic
ks.result1.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result1.stad$statistic,2), ifelse(ks.result1.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result1.stad$p.value,4)))) -> text.annot2
length(ks.result1.stad$data$x) -> text.annot3

dim(out.stad1)
dim(out.rand.stad.graph1)

# Create gg plot
ggplot(dat.stad1, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad1, aes(xintercept=random.cutoff1.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and stad Genes on the stad Network") -> gg.stad1

gg.stad1
ggsave("~/tcga_biolinks1/Plots/gg.stad1.png", plot = gg.stad1, width = 10, height = 10)

###### 0.5% PATIENTS #######################################################################

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> stad.graph.df.filt0.5

unique(stad.graph.df.filt0.5$name) -> unique.genes.stad0.5

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.stad0.5 <- foreach(i=1:length(unique.genes.stad0.5), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.stad2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = unique.genes.stad0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.stad0.5,])}

colnames(out.stad0.5) <- unique.genes.stad0.5

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(stad.graph),"name"), length(unique.genes.stad0.5)) -> all.genes.rand.stad0.5

out.rand.stad.graph0.5 <- foreach(i=1:length(all.genes.rand.stad0.5), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = all.genes.rand.stad0.5[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.stad0.5,])}

# # Save results of distances
saveRDS(out.rand.stad.graph0.5, file = "~/tcga_biolinks1/stats/out.rand.stad.graph0.5")
saveRDS(out.stad0.5, file = "~/tcga_biolinks1/stats/out.stad0.5")


# out.rand.stad.graph2 <- readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph2")
# out.stad2 <- readRDS( file = "~/tcga_biolinks1/stats/out.stad2")

# Create histogram of column means
# hist(colMeans(out.stad2))
# hist(colMeans(out.rand.stad.graph2))

# Turn inifnite into NAs
out.stad0.5[which(is.infinite(out.stad0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph0.5, na.rm = T),0.01) -> random.cutoff0.5.stad

dat.stad0.5 <- data.frame(Distance = c(colMeans(out.stad0.5), colMeans(out.rand.stad.graph0.5, na.rm = T)),
                          Distribution = factor(c(rep("Real",nrow(out.stad0.5)),rep("Random",nrow(out.rand.stad.graph0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.stad0.5), colMeans(out.rand.stad.graph0.5, na.rm = T)) -> ks.result0.5.stad

# Print
ks.result0.5.stad$statistic
ks.result0.5.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result0.5.stad$statistic,2), ifelse(ks.result0.5.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result0.5.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result0.5.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result0.5.stad$p.value,4)))) -> text.annot2
length(ks.result0.5.stad$data$x) -> text.annot3

dim(out.stad0.5)
dim(out.rand.stad.graph0.5)

# Create gg plot
ggplot(dat.stad0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad0.5, aes(xintercept=random.cutoff0.5.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and stad Genes on the stad Network") -> gg.stad0.5

gg.stad0.5
ggsave("~/tcga_biolinks1/gg.stad0.5.png", plot = gg.stad0.5, width = 10, height = 10)

########### %2 #################################################################
stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 2) -> stad.graph.df.filt2

unique(stad.graph.df.filt2$name) -> unique.genes.stad2

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.stad2 <- foreach(i=1:length(unique.genes.stad2), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.stad2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = unique.genes.stad2[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.stad2,])}

colnames(out.stad2) <- unique.genes.stad2

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(stad.graph),"name"), length(unique.genes.stad2)) -> all.genes.rand.stad2

out.rand.stad.graph2 <- foreach(i=1:length(all.genes.rand.stad2), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = all.genes.rand.stad2[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.stad2,])}

# # Save results of distances
saveRDS(out.rand.stad.graph2, file = "~/tcga_biolinks1/stats/out.rand.stad.graph2")
saveRDS(out.stad2, file = "~/tcga_biolinks1/stats/out.stad2")

# Turn inifnite into NAs
out.stad2[which(is.infinite(out.stad2))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph2, na.rm = T),0.01) -> random.cutoff2.stad

dat.stad2 <- data.frame(Distance = c(colMeans(out.stad2), colMeans(out.rand.stad.graph2, na.rm = T)),
                        Distribution = factor(c(rep("Real",nrow(out.stad2)),rep("Random",nrow(out.rand.stad.graph2)))))

# Perform KS test on column means
ks.test(colMeans(out.stad2), colMeans(out.rand.stad.graph2, na.rm = T)) -> ks.result2.stad

# Print
ks.result2.stad$statistic
ks.result2.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result2.stad$statistic,2), ifelse(ks.result2.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result2.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result2.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result2.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result2.stad$p.value,4)))) -> text.annot2
length(ks.result2.stad$data$x) -> text.annot3

dim(out.stad2)
dim(out.rand.stad.graph2)

# Create gg plot
ggplot(dat.stad2, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad2, aes(xintercept=random.cutoff2.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and stad Genes on the stad Network") -> gg.stad2

gg.stad2
ggsave("~/tcga_biolinks1/gg.stad2.png", plot = gg.stad2, width = 10, height = 10)


##### PROBABLY DAMAGING >0.908 ##############################

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> stad.graph.df.filt.prob

unique(stad.graph.df.filt.prob$name) -> unique.genes.stad.prob

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.stad.prob <- foreach(i=1:length(unique.genes.stad.prob), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.stad.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = unique.genes.stad.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.stad.prob,])}

colnames(out.stad.prob) <- unique.genes.stad.prob

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(stad.graph),"name"), length(unique.genes.stad.prob)) -> all.genes.rand.stad.prob

out.rand.stad.graph.prob <- foreach(i=1:length(all.genes.rand.stad.prob), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = all.genes.rand.stad.prob[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.stad.prob,])}

# # Save results of distances
saveRDS(out.rand.stad.graph.prob, file = "~/tcga_biolinks1/stats/out.rand.stad.graph.prob")
saveRDS(out.stad.prob, file = "~/tcga_biolinks1/stats/out.stad.prob")

# Create histogram of column means
# hist(colMeans(out.stad.prob))
# hist(colMeans(out.rand.stad.graph.prob))

# Turn inifnite into NAs
out.stad.prob[which(is.infinite(out.stad.prob))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph.prob, na.rm = T),0.01) -> random.cutoff.prob.stad

dat.stad.prob <- data.frame(Distance = c(colMeans(out.stad.prob), colMeans(out.rand.stad.graph.prob, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.stad.prob)),rep("Random",nrow(out.rand.stad.graph.prob)))))

# Perform KS test on column means
ks.test(colMeans(out.stad.prob), colMeans(out.rand.stad.graph.prob, na.rm = T)) -> ks.result.prob.stad

# Print
ks.result.prob.stad$statistic
ks.result.prob.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.prob.stad$statistic,2), ifelse(ks.result.prob.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.prob.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.prob.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.prob.stad$p.value,4)))) -> text.annot2
length(ks.result.prob.stad$data$x) -> text.annot3

dim(out.stad.prob)
dim(out.rand.stad.graph.prob)

# Create gg plot
ggplot(dat.stad.prob, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad.prob, aes(xintercept=random.cutoff.prob.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and stad Genes on the stad Network") -> gg.stad.prob

gg.stad.prob
ggsave("./gg.stad.prob.png", plot = gg.stad.prob, width = 10, height = 10)

##### POSSIBLY DAMAGING =<0.908 >0.466 ##############################

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean <=0.908) %>%
  filter(PolyPhen.Mean >0.466)-> stad.graph.df.filt.poss

unique(stad.graph.df.filt.poss$name) -> unique.genes.stad.poss

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.stad.poss <- foreach(i=1:length(unique.genes.stad.poss), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.stad.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = unique.genes.stad.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.stad.poss,])}

colnames(out.stad.poss) <- unique.genes.stad.poss

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(stad.graph),"name"), length(unique.genes.stad.poss)) -> all.genes.rand.stad.poss

out.rand.stad.graph.poss <- foreach(i=1:length(all.genes.rand.stad.poss), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = all.genes.rand.stad.poss[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.stad.poss,])}

# # Save results of distances
saveRDS(out.rand.stad.graph.poss, file = "~/tcga_biolinks1/stats/out.rand.stad.graph.poss")
saveRDS(out.stad.poss, file = "~/tcga_biolinks1/stats/out.stad.poss")

# Create histogram of column means
# hist(colMeans(out.stad.poss))
# hist(colMeans(out.rand.stad.graph.poss))

# Turn inifnite into NAs
out.stad.poss[which(is.infinite(out.stad.poss))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph.poss, na.rm = T),0.01) -> random.cutoff.poss.stad

dat.stad.poss <- data.frame(Distance = c(colMeans(out.stad.poss), colMeans(out.rand.stad.graph.poss, na.rm = T)),
                            Distribution = factor(c(rep("Real",nrow(out.stad.poss)),rep("Random",nrow(out.rand.stad.graph.poss)))))

# Perform KS test on column means
ks.test(colMeans(out.stad.poss), colMeans(out.rand.stad.graph.poss, na.rm = T)) -> ks.result.poss.stad

# Print
ks.result.poss.stad$statistic
ks.result.poss.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.poss.stad$statistic,2), ifelse(ks.result.poss.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.poss.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.poss.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.poss.stad$p.value,4)))) -> text.annot2
length(ks.result.poss.stad$data$x) -> text.annot3

dim(out.stad.poss)
dim(out.rand.stad.graph.poss)

# Create gg plot
ggplot(dat.stad.poss, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad.poss, aes(xintercept=random.cutoff.poss.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and stad Genes on the stad Network") -> gg.stad.poss

gg.stad.poss
ggsave("./gg.stad.poss.png", plot = gg.stad.poss, width = 10, height = 10)

##### ALL DAMAGING >0.466 ##########################################################################

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> stad.graph.df.filt.dam

unique(stad.graph.df.filt.dam$name) -> unique.genes.stad.dam

# Install packages
library(foreach)
library(doParallel)

# Run in parallel
registerDoParallel(cores = 30)

# Calculate distances between all genes
#i = 1
out.stad.dam <- foreach(i=1:length(unique.genes.stad.dam), .combine = cbind)%dopar%{
  # out <- foreach(i=1:length(unique.genes.stad.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = unique.genes.stad.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[unique.genes.stad.dam,])}

colnames(out.stad.dam) <- unique.genes.stad.dam

# Calculate distances between random genes of same length as filtered genes

sample(attr(V(stad.graph),"name"), length(unique.genes.stad.dam)) -> all.genes.rand.stad.dam

out.rand.stad.graph.dam <- foreach(i=1:length(all.genes.rand.stad.dam), .combine = cbind)%dopar%{
  temp.distance.table <- as.data.frame(distances(
    stad.graph,
    to = all.genes.rand.stad.dam[i],
    mode = c("all")
  ))
  return(temp.distance.table[all.genes.rand.stad.dam,])}

# # Save results of distances
saveRDS(out.rand.stad.graph.dam, file = "~/tcga_biolinks1/stats/out.rand.stad.graph.dam")
saveRDS(out.stad.dam, file = "~/tcga_biolinks1/stats/out.stad.dam")

# Create histogram of column means
# hist(colMeans(out.stad.dam))
# hist(colMeans(out.rand.stad.graph.dam))

# Turn inifnite into NAs
out.stad.dam[which(is.infinite(out.stad.dam))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph.dam, na.rm = T),0.01) -> random.cutoff.dam.stad

dat.stad.dam <- data.frame(Distance = c(colMeans(out.stad.dam), colMeans(out.rand.stad.graph.dam, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.stad.dam)),rep("Random",nrow(out.rand.stad.graph.dam)))))

# Perform KS test on column means
ks.test(colMeans(out.stad.dam), colMeans(out.rand.stad.graph.dam, na.rm = T)) -> ks.result.dam.stad

# Print
ks.result.dam.stad$statistic
ks.result.dam.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.stad$statistic,2), ifelse(ks.result.dam.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.stad$p.value,4)))) -> text.annot2
length(ks.result.dam.stad$data$x) -> text.annot3

dim(out.stad.dam)
dim(out.rand.stad.graph.dam)

# Create gg plot
ggplot(dat.stad.dam, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad.dam, aes(xintercept=random.cutoff.dam.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle("Mean Gene Distance Distributions of Both Random and stad Genes on the stad Network") -> gg.stad.dam

gg.stad.dam
ggsave("~/tcga_biolinks1/Plots/gg.stad.dam.png", plot = gg.stad.dam, width = 10, height = 10)

