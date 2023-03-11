#### ARRANGE GG PLOTS INTO FIGURE ##################################
# Create gg plots
#### 0.5% ######################################################################
readRDS(file = "~/tcga_biolinks1/stats/out.prad0.5") -> out.prad0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph0.5") -> out.rand.prad.graph0.5

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
  ggtitle(NULL) -> gg.prad0.5

gg.prad0.5

#### 1% ########################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.prad1") -> out.prad1
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph1") -> out.rand.prad.graph1

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
  ggtitle(NULL) -> gg.prad1

gg.prad1

#### 1.5% ######################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.prad1.5") -> out.prad1.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph1.5") -> out.rand.prad.graph1.5

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
  annotate("text", x = 3.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 3.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.prad1.5

gg.prad1.5

#### PROB DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.prad.prob") -> out.prad.prob
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph.prob") -> out.rand.prad.graph.prob
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
  ggtitle(NULL) -> gg.prad.prob

gg.prad.prob
#### POSS DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.prad.poss") -> out.prad.poss
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph.poss") -> out.rand.prad.graph.poss

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
  ggtitle(NULL) -> gg.prad.poss

gg.prad.poss

#### ALL DAMAGING ##############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.prad.dam") -> out.prad.dam
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph.dam") -> out.rand.prad.graph.dam

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
  ggtitle(NULL) -> gg.prad.dam

gg.prad.dam

##### ALL DAMAGING AND >=0.5% ##################################################

readRDS(file = "~/tcga_biolinks1/stats/out.prad.dam.0.5") -> out.prad.dam.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph.dam.0.5") -> out.rand.prad.graph.dam.0.5

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
  ggtitle(NULL) -> gg.prad.dam.0.5

gg.prad.dam.0.5

##### SIFT DELETERIOUS #########################################################

readRDS(file = "~/tcga_biolinks1/stats/out.prad.del") -> out.prad.del
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph.del") -> out.rand.prad.graph.del

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
  ggtitle(NULL) -> gg.prad.del

gg.prad.del
##### SIFT DELETERIOUS AND >= 0.5% ##############################################

readRDS(file = "~/tcga_biolinks1/stats/out.prad.del.0.5") -> out.prad.del.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph.del.0.5") -> out.rand.prad.graph.del.0.5

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
  ggtitle(NULL) -> gg.prad.del.0.5

gg.prad.del.0.5
##### GGARRANGE ################################################################

library(ggpubr)
gg.prad0.5 + labs(subtitle = "Only Genes Missense Mutated in >= 0.5% of All Patients") -> gg.prad0.5t
gg.prad0.5t + theme(text = element_text(size = 10)) -> p1
gg.prad1 + labs(subtitle = "Only Genes Missense Mutated in >= 1% of All Patients") -> gg.prad1t
gg.prad1t + theme(text = element_text(size = 10)) -> p2
gg.prad1.5 + labs(subtitle = "Only Genes Missense Mutated in >= 1.5% of All Patients") -> gg.prad1.5t
gg.prad1.5t + theme(text = element_text(size = 10)) -> p3
gg.prad.prob + labs(subtitle = "Only Genes with Probably Damaging Missense Mutations") -> gg.prad.probt
gg.prad.probt + theme(text = element_text(size = 10)) ->p4
gg.prad.poss + labs(subtitle = "Only Genes with Possibly Damaging Missense Mutations") -> gg.prad.posst
gg.prad.posst + theme(text = element_text(size = 10)) -> p5
gg.prad.dam + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations") -> gg.prad.damt
gg.prad.damt + theme(text = element_text(size = 10)) -> p6
gg.prad.dam.0.5 + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations,
in >=0.5% of All Patients") -> gg.prad.dam.0.5t
gg.prad.dam.0.5t + theme(text = element_text(size = 10)) -> p7
gg.prad.del + labs(subtitle = "Only Genes with Deleterious Missense Mutations") -> gg.prad.delt
gg.prad.delt + theme(text = element_text(size = 10)) -> p8
gg.prad.del.0.5 + labs(subtitle = "Only Genes with Deleterious Missense Mutations, in >0.5% of All Patients") -> gg.prad.del.0.5t
gg.prad.del.0.5t + theme(text = element_text(size = 10)) -> p9

prad.figure <- ggarrange(p1, p2, p3, 
                         p4, p5, p6,
                         p7, p8, p9,
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                         common.legend = TRUE, legend = "bottom",
                         ncol = 3, nrow = 3)
annotate_figure(prad.figure, top = text_grob("Mean Gene Distance Distributions of Both Random and PRAD Genes on the PRAD Network", size = 20, face = "bold")) -> prad.plot

ggsave("~/tcga_biolinks1/Plots/prad.plot.png", plot = prad.plot, width = 18, height = 10)

