#### ARRANGE GG PLOTS INTO FIGURE ##################################
# Create gg plots
#### 0.5% ######################################################################
readRDS(file = "~/tcga_biolinks1/stats/out.kirc0.5") -> out.kirc0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph0.5") -> out.rand.kirc.graph0.5

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
  ggtitle(NULL) -> gg.kirc0.5

gg.kirc0.5

#### 1% ########################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirc1") -> out.kirc1
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph1") -> out.rand.kirc.graph1

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
  ggtitle(NULL) -> gg.kirc1

gg.kirc1

#### 1.5% ######################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirc1.5") -> out.kirc1.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph1.5") -> out.rand.kirc.graph1.5

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
  ggtitle(NULL) -> gg.kirc1.5

gg.kirc1.5

#### PROB DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirc.prob") -> out.kirc.prob
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph.prob") -> out.rand.kirc.graph.prob
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
  ggtitle(NULL) -> gg.kirc.prob

gg.kirc.prob
#### POSS DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirc.poss") -> out.kirc.poss
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph.poss") -> out.rand.kirc.graph.poss

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
  ggtitle(NULL) -> gg.kirc.poss

gg.kirc.poss

#### ALL DAMAGING ##############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirc.dam") -> out.kirc.dam
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph.dam") -> out.rand.kirc.graph.dam

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
  ggtitle(NULL) -> gg.kirc.dam

gg.kirc.dam

##### ALL DAMAGING AND >=0.5% ##################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirc.dam.0.5") -> out.kirc.dam.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph.dam.0.5") -> out.rand.kirc.graph.dam.0.5

out.kirc.dam.0.5[which(is.infinite(out.kirc.dam.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirc.graph.dam.0.5, na.rm = T),0.01) -> random.cutoff.dam.0.5.kirc

dat.kirc.dam.0.5 <- data.frame(Distance = c(colMeans(out.kirc.dam.0.5), colMeans(out.rand.kirc.graph.dam.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.kirc.dam.0.5)),rep("Random",nrow(out.rand.kirc.graph.dam.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.kirc.dam.0.5), colMeans(out.rand.kirc.graph.dam.0.5, na.rm = T)) -> ks.result.dam.0.5.kirc

# Print
ks.result.dam.0.5.kirc$statistic
ks.result.dam.0.5.kirc$p.value
# length(ks.result2.kirc$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.0.5.kirc$statistic,2), ifelse(ks.result.dam.0.5.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.kirc$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.0.5.kirc$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.0.5.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.kirc$p.value,4)))) -> text.annot2
length(ks.result.dam.0.5.kirc$data$x) -> text.annot3

dim(out.kirc.dam.0.5)
dim(out.rand.kirc.graph.dam.0.5)

# Create gg plot
ggplot(dat.kirc.dam.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirc.dam.0.5, aes(xintercept=random.cutoff.dam.0.5.kirc,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.kirc.dam.0.5

gg.kirc.dam.0.5

##### SIFT DELETERIOUS #########################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirc.del") -> out.kirc.del
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph.del") -> out.rand.kirc.graph.del

out.kirc.del[which(is.infinite(out.kirc.del))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirc.graph.del, na.rm = T),0.01) -> random.cutoff.del.kirc

dat.kirc.del <- data.frame(Distance = c(colMeans(out.kirc.del), colMeans(out.rand.kirc.graph.del, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.kirc.del)),rep("Random",nrow(out.rand.kirc.graph.del)))))

# Perform KS test on column means
ks.test(colMeans(out.kirc.del), colMeans(out.rand.kirc.graph.del, na.rm = T)) -> ks.result.del.kirc

# Print
ks.result.del.kirc$statistic
ks.result.del.kirc$p.value
# length(ks.result2.kirc$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.kirc$statistic,2), ifelse(ks.result.del.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.kirc$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.kirc$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.kirc$p.value,4)))) -> text.annot2
length(ks.result.del.kirc$data$x) -> text.annot3

dim(out.kirc.del)
dim(out.rand.kirc.graph.del)

# Create gg plot
ggplot(dat.kirc.del, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirc.del, aes(xintercept=random.cutoff.del.kirc,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.kirc.del

gg.kirc.del
##### SIFT DELETERIOUS AND >= 0.5% ##############################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirc.del.0.5") -> out.kirc.del.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirc.graph.del.0.5") -> out.rand.kirc.graph.del.0.5

out.kirc.del.0.5[which(is.infinite(out.kirc.del.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirc.graph.del.0.5, na.rm = T),0.01) -> random.cutoff.del.0.5.kirc

dat.kirc.del.0.5 <- data.frame(Distance = c(colMeans(out.kirc.del.0.5), colMeans(out.rand.kirc.graph.del.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.kirc.del.0.5)),rep("Random",nrow(out.rand.kirc.graph.del.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.kirc.del.0.5), colMeans(out.rand.kirc.graph.del.0.5, na.rm = T)) -> ks.result.del.0.5.kirc

# Print
ks.result.del.0.5.kirc$statistic
ks.result.del.0.5.kirc$p.value
# length(ks.result2.kirc$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.0.5.kirc$statistic,2), ifelse(ks.result.del.0.5.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.kirc$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.0.5.kirc$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.0.5.kirc$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.kirc$p.value,4)))) -> text.annot2
length(ks.result.del.0.5.kirc$data$x) -> text.annot3

dim(out.kirc.del.0.5)
dim(out.rand.kirc.graph.del.0.5)

# Create gg plot
ggplot(dat.kirc.del.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirc.del.0.5, aes(xintercept=random.cutoff.del.0.5.kirc,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.kirc.del.0.5

gg.kirc.del.0.5
##### GGARRANGE ################################################################
library(ggpubr)
gg.kirc0.5 + labs(subtitle = "Only Genes Missense Mutated in >= 0.5% of All Patients") -> gg.kirc0.5t
gg.kirc0.5t + theme(text = element_text(size = 10)) -> p1
gg.kirc1 + labs(subtitle = "Only Genes Missense Mutated in >= 1% of All Patients") -> gg.kirc1t
gg.kirc1t + theme(text = element_text(size = 10)) -> p2
gg.kirc1.5 + labs(subtitle = "Only Genes Missense Mutated in >= 1.5% of All Patients") -> gg.kirc1.5t
gg.kirc1.5t + theme(text = element_text(size = 10)) -> p3
gg.kirc.prob + labs(subtitle = "Only Genes with Probably Damaging Missense Mutations") -> gg.kirc.probt
gg.kirc.probt + theme(text = element_text(size = 10)) ->p4
gg.kirc.poss + labs(subtitle = "Only Genes with Possibly Damaging Missense Mutations") -> gg.kirc.posst
gg.kirc.posst + theme(text = element_text(size = 10)) -> p5
gg.kirc.dam + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations") -> gg.kirc.damt
gg.kirc.damt + theme(text = element_text(size = 10)) -> p6
gg.kirc.dam.0.5 + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations,
in >=0.5% of All Patients") -> gg.kirc.dam.0.5t
gg.kirc.dam.0.5t + theme(text = element_text(size = 10)) -> p7
gg.kirc.del + labs(subtitle = "Only Genes with Deleterious Missense Mutations") -> gg.kirc.delt
gg.kirc.delt + theme(text = element_text(size = 10)) -> p8
gg.kirc.del.0.5 + labs(subtitle = "Only Genes with Deleterious Missense Mutations, in >0.5% of All Patients") -> gg.kirc.del.0.5t
gg.kirc.del.0.5t + theme(text = element_text(size = 10)) -> p9

kirc.figure <- ggarrange(p1, p2, p3, 
                         p4, p5, p6,
                         p7, p8, p9,
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                         common.legend = TRUE, legend = "bottom",
                         ncol = 3, nrow = 3)
annotate_figure(kirc.figure, top = text_grob("Mean Gene Distance Distributions of Both Random and KIRC Genes on the KIRC Network", size = 20, face = "bold")) -> kirc.plot

ggsave("~/tcga_biolinks1/Plots/kirc.plot.png", plot = kirc.plot, width = 18, height = 10)

