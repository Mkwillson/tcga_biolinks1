#### ARRANGE GG PLOTS INTO FIGURE ##################################
# Create gg plots
#### 0.5% ######################################################################
readRDS(file = "~/tcga_biolinks1/stats/out.luad0.5") -> out.luad0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph0.5") -> out.rand.luad.graph0.5

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
  ggtitle(NULL) -> gg.luad0.5

gg.luad0.5

#### 1% ########################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.luad1") -> out.luad1
readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph1") -> out.rand.luad.graph1

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
  ggtitle(NULL) -> gg.luad1

gg.luad1

#### 1.5% ######################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.luad1.5") -> out.luad1.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph1.5") -> out.rand.luad.graph1.5

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
  ggtitle(NULL) -> gg.luad1.5

gg.luad1.5

#### PROB DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.luad.prob") -> out.luad.prob
readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph.prob") -> out.rand.luad.graph.prob
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
  ggtitle(NULL) -> gg.luad.prob

gg.luad.prob
#### POSS DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.luad.poss") -> out.luad.poss
readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph.poss") -> out.rand.luad.graph.poss

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
  ggtitle(NULL) -> gg.luad.poss

gg.luad.poss

#### ALL DAMAGING ##############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.luad.dam") -> out.luad.dam
readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph.dam") -> out.rand.luad.graph.dam

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
  ggtitle(NULL) -> gg.luad.dam

gg.luad.dam

##### ALL DAMAGING AND >=0.5% ##################################################

readRDS(file = "~/tcga_biolinks1/stats/out.luad.dam.0.5") -> out.luad.dam.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph.dam.0.5") -> out.rand.luad.graph.dam.0.5

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
  ggtitle(NULL) -> gg.luad.dam.0.5

gg.luad.dam.0.5

##### SIFT DELETERIOUS #########################################################

readRDS(file = "~/tcga_biolinks1/stats/out.luad.del") -> out.luad.del
readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph.del") -> out.rand.luad.graph.del

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
  ggtitle(NULL) -> gg.luad.del

gg.luad.del
##### SIFT DELETERIOUS AND >= 0.5% ##############################################

readRDS(file = "~/tcga_biolinks1/stats/out.luad.del.0.5") -> out.luad.del.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.luad.graph.del.0.5") -> out.rand.luad.graph.del.0.5

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
  ggtitle(NULL) -> gg.luad.del.0.5

gg.luad.del.0.5
##### GGARRANGE ################################################################
library(ggpubr)
gg.luad0.5 + labs(subtitle = "Only Genes Missense Mutated in >= 0.5% of All Patients") -> gg.luad0.5t
gg.luad0.5t + theme(text = element_text(size = 10)) -> p1
gg.luad1 + labs(subtitle = "Only Genes Missense Mutated in >= 1% of All Patients") -> gg.luad1t
gg.luad1t + theme(text = element_text(size = 10)) -> p2
gg.luad1.5 + labs(subtitle = "Only Genes Missense Mutated in >= 1.5% of All Patients") -> gg.luad1.5t
gg.luad1.5t + theme(text = element_text(size = 10)) -> p3
gg.luad.prob + labs(subtitle = "Only Genes with Probably Damaging Missense Mutations") -> gg.luad.probt
gg.luad.probt + theme(text = element_text(size = 10)) ->p4
gg.luad.poss + labs(subtitle = "Only Genes with Possibly Damaging Missense Mutations") -> gg.luad.posst
gg.luad.posst + theme(text = element_text(size = 10)) -> p5
gg.luad.dam + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations") -> gg.luad.damt
gg.luad.damt + theme(text = element_text(size = 10)) -> p6
gg.luad.dam.0.5 + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations,
in >=0.5% of All Patients") -> gg.luad.dam.0.5t
gg.luad.dam.0.5t + theme(text = element_text(size = 10)) -> p7
gg.luad.del + labs(subtitle = "Only Genes with Deleterious Missense Mutations") -> gg.luad.delt
gg.luad.delt + theme(text = element_text(size = 10)) -> p8
gg.luad.del.0.5 + labs(subtitle = "Only Genes with Deleterious Missense Mutations, in >0.5% of All Patients") -> gg.luad.del.0.5t
gg.luad.del.0.5t + theme(text = element_text(size = 10)) -> p9

luad.figure <- ggarrange(p1, p2, p3, 
                         p4, p5, p6,
                         p7, p8, p9,
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                         common.legend = TRUE, legend = "bottom",
                         ncol = 3, nrow = 3)
annotate_figure(luad.figure, top = text_grob("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network", size = 20, face = "bold")) -> luad.plot

ggsave("~/tcga_biolinks1/Plots/luad.plot.png", plot = luad.plot, width = 18, height = 10)

