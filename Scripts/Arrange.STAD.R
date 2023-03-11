#### ARRANGE GG PLOTS INTO FIGURE ##################################
# Create gg plots
#### 0.5% ######################################################################
readRDS(file = "~/tcga_biolinks1/stats/out.stad0.5") -> out.stad0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph0.5") -> out.rand.stad.graph0.5

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
  ggtitle(NULL) -> gg.stad0.5

gg.stad0.5

#### 1% ########################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.stad1") -> out.stad1
readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph1") -> out.rand.stad.graph1

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
  ggtitle(NULL) -> gg.stad1

gg.stad1

#### 1.5% ######################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.stad1.5") -> out.stad1.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph1.5") -> out.rand.stad.graph1.5

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
  ggtitle(NULL) -> gg.stad1.5

gg.stad1.5

#### PROB DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.stad.prob") -> out.stad.prob
readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph.prob") -> out.rand.stad.graph.prob
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
  ggtitle(NULL) -> gg.stad.prob

gg.stad.prob
#### POSS DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.stad.poss") -> out.stad.poss
readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph.poss") -> out.rand.stad.graph.poss

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
  ggtitle(NULL) -> gg.stad.poss

gg.stad.poss

#### ALL DAMAGING ##############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.stad.dam") -> out.stad.dam
readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph.dam") -> out.rand.stad.graph.dam

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
  ggtitle(NULL) -> gg.stad.dam

gg.stad.dam

##### ALL DAMAGING AND >=0.5% ##################################################

readRDS(file = "~/tcga_biolinks1/stats/out.stad.dam.0.5") -> out.stad.dam.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph.dam.0.5") -> out.rand.stad.graph.dam.0.5

out.stad.dam.0.5[which(is.infinite(out.stad.dam.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph.dam.0.5, na.rm = T),0.01) -> random.cutoff.dam.0.5.stad

dat.stad.dam.0.5 <- data.frame(Distance = c(colMeans(out.stad.dam.0.5), colMeans(out.rand.stad.graph.dam.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.stad.dam.0.5)),rep("Random",nrow(out.rand.stad.graph.dam.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.stad.dam.0.5), colMeans(out.rand.stad.graph.dam.0.5, na.rm = T)) -> ks.result.dam.0.5.stad

# Print
ks.result.dam.0.5.stad$statistic
ks.result.dam.0.5.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.0.5.stad$statistic,2), ifelse(ks.result.dam.0.5.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.0.5.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.0.5.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.stad$p.value,4)))) -> text.annot2
length(ks.result.dam.0.5.stad$data$x) -> text.annot3

dim(out.stad.dam.0.5)
dim(out.rand.stad.graph.dam.0.5)

# Create gg plot
ggplot(dat.stad.dam.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad.dam.0.5, aes(xintercept=random.cutoff.dam.0.5.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.stad.dam.0.5

gg.stad.dam.0.5

##### SIFT DELETERIOUS #########################################################

readRDS(file = "~/tcga_biolinks1/stats/out.stad.del") -> out.stad.del
readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph.del") -> out.rand.stad.graph.del

out.stad.del[which(is.infinite(out.stad.del))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph.del, na.rm = T),0.01) -> random.cutoff.del.stad

dat.stad.del <- data.frame(Distance = c(colMeans(out.stad.del), colMeans(out.rand.stad.graph.del, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.stad.del)),rep("Random",nrow(out.rand.stad.graph.del)))))

# Perform KS test on column means
ks.test(colMeans(out.stad.del), colMeans(out.rand.stad.graph.del, na.rm = T)) -> ks.result.del.stad

# Print
ks.result.del.stad$statistic
ks.result.del.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.stad$statistic,2), ifelse(ks.result.del.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.stad$p.value,4)))) -> text.annot2
length(ks.result.del.stad$data$x) -> text.annot3

dim(out.stad.del)
dim(out.rand.stad.graph.del)

# Create gg plot
ggplot(dat.stad.del, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad.del, aes(xintercept=random.cutoff.del.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.stad.del

gg.stad.del
##### SIFT DELETERIOUS AND >= 0.5% ##############################################

readRDS(file = "~/tcga_biolinks1/stats/out.stad.del.0.5") -> out.stad.del.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.stad.graph.del.0.5") -> out.rand.stad.graph.del.0.5

out.stad.del.0.5[which(is.infinite(out.stad.del.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.stad.graph.del.0.5, na.rm = T),0.01) -> random.cutoff.del.0.5.stad

dat.stad.del.0.5 <- data.frame(Distance = c(colMeans(out.stad.del.0.5), colMeans(out.rand.stad.graph.del.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.stad.del.0.5)),rep("Random",nrow(out.rand.stad.graph.del.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.stad.del.0.5), colMeans(out.rand.stad.graph.del.0.5, na.rm = T)) -> ks.result.del.0.5.stad

# Print
ks.result.del.0.5.stad$statistic
ks.result.del.0.5.stad$p.value
# length(ks.result2.stad$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.0.5.stad$statistic,2), ifelse(ks.result.del.0.5.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.stad$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.0.5.stad$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.0.5.stad$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.stad$p.value,4)))) -> text.annot2
length(ks.result.del.0.5.stad$data$x) -> text.annot3

dim(out.stad.del.0.5)
dim(out.rand.stad.graph.del.0.5)

# Create gg plot
ggplot(dat.stad.del.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.stad.del.0.5, aes(xintercept=random.cutoff.del.0.5.stad,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.stad.del.0.5

gg.stad.del.0.5
##### GGARRANGE ################################################################
library(ggpubr)
gg.stad0.5 + labs(subtitle = "Only Genes Missense Mutated in >= 0.5% of All Patients") -> gg.stad0.5t
gg.stad.0.5t + theme(text = element_text(size = 10)) -> p1
gg.stad1 + labs(subtitle = "Only Genes Missense Mutated in >= 1% of All Patients") -> gg.stad1t
gg.stad1t + theme(text = element_text(size = 10)) -> p2
gg.stad1.5 + labs(subtitle = "Only Genes Missense Mutated in >= 1.5% of All Patients") -> gg.stad1.5t
gg.stad1.5t + theme(text = element_text(size = 10)) -> p3
gg.stad.prob + labs(subtitle = "Only Genes with Probably Damaging Missense Mutations") -> gg.stad.probt
gg.stad.probt + theme(text = element_text(size = 10)) ->p4
gg.stad.poss + labs(subtitle = "Only Genes with Possibly Damaging Missense Mutations") -> gg.stad.posst
gg.stad.posst + theme(text = element_text(size = 10)) -> p5
gg.stad.dam + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations") -> gg.stad.damt
gg.stad.damt + theme(text = element_text(size = 10)) -> p6
gg.stad.dam.0.5 + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations,
in >=0.5% of All Patients") -> gg.stad.dam.0.5t
gg.stad.dam.0.5t + theme(text = element_text(size = 10)) -> p7
gg.stad.del + labs(subtitle = "Only Genes with Deleterious Missense Mutations") -> gg.stad.delt
gg.stad.delt + theme(text = element_text(size = 10)) -> p8
gg.stad.del.0.5 + labs(subtitle = "Only Genes with Deleterious Missense Mutations, in >0.5% of All Patients") -> gg.stad.del.0.5t
gg.stad.del.0.5t + theme(text = element_text(size = 10)) -> p9

stad.figure <- ggarrange(p1, p2, p3, 
                         p4, p5, p6,
                         p7, p8, p9,
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                         common.legend = TRUE, legend = "bottom",
                         ncol = 3, nrow = 3)
annotate_figure(stad.figure, top = text_grob("Mean Gene Distance Distributions of Both Random and STAD Genes on the STAD Network", size = 20, face = "bold")) -> stad.plot

ggsave("~/tcga_biolinks1/Plots/stad.plot.png", plot = stad.plot, width = 18, height = 10)

