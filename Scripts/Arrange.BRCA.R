#### ARRANGE GG PLOTS INTO FIGURE ##################################
# Create gg plots
#### 0.5% ######################################################################
readRDS(file = "~/tcga_biolinks1/stats/out.brca0.5") -> out.brca0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph0.5") -> out.rand.brca.graph0.5

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
paste0("D = ",round(ks.result0.5.brca$statistic,2), ifelse(ks.result0.5.brca$p.value <0.0001, " p < 0.0001", paste0(" p = ",signif(ks.result0.5.brca$p.value,2)))) -> text.annot
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
  ggtitle(NULL) -> gg.brca0.5

gg.brca0.5

#### 1% ########################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.brca1") -> out.brca1
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph1") -> out.rand.brca.graph1

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
paste0("D = ",round(ks.result1.brca$statistic,2), ifelse(ks.result1.brca$p.value <0.0001, " p < 0.0001", paste0(" p = ",signif(ks.result1.brca$p.value,2)))) -> text.annot
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
  ggtitle(NULL) -> gg.brca1

gg.brca1

#### 1.5% ######################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.brca1.5") -> out.brca1.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph1.5") -> out.rand.brca.graph1.5

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

paste0("D = ",round(ks.result1.5.brca$statistic,2), ifelse(ks.result1.5.brca$p.value <0.0001, " p < 0.0001", paste0(" p = ",signif(ks.result1.5.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result1.5.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result1.5.brca$p.value==0, " p < 0.0001", paste0(" p = ",round(ks.result1.5.brca$p.value,8)))) -> text.annot2
length(ks.result1.5.brca$data$x) -> text.annot3

dim(out.brca1.5)
dim(out.rand.brca.graph1.5)

# Create gg plot
ggplot(dat.brca1.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca1.5, aes(xintercept=random.cutoff1.5.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  # annotate("text", x = 4.45, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.brca1.5

gg.brca1.5

#### PROB DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.brca.prob") -> out.brca.prob
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph.prob") -> out.rand.brca.graph.prob
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
paste0("D = ",round(ks.result.prob.brca$statistic,2), ifelse(ks.result.prob.brca$p.value <0.0001, " p < 0.0001", paste0(" p = ",signif(ks.result.prob.brca$p.value,2)))) -> text.annot
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
  ggtitle(NULL) -> gg.brca.prob

gg.brca.prob

#### POSS DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.brca.poss") -> out.brca.poss
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph.poss") -> out.rand.brca.graph.poss

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
paste0("D = ",round(ks.result.poss.brca$statistic,2), ifelse(ks.result.poss.brca$p.value <0.0001, " p < 0.0001", paste0(" p = ",signif(ks.result.poss.brca$p.value,2)))) -> text.annot
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
  ggtitle(NULL) -> gg.brca.poss

gg.brca.poss

#### ALL DAMAGING ##############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.brca.dam") -> out.brca.dam
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph.dam") -> out.rand.brca.graph.dam

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
paste0("D = ",round(ks.result.dam.brca$statistic,2), ifelse(ks.result.dam.brca$p.value <0.0001, " p < 0.0001", paste0(" p = ",signif(ks.result.dam.brca$p.value,2)))) -> text.annot
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
  ggtitle(NULL) -> gg.brca.dam

gg.brca.dam

##### ALL DAMAGING AND >=0.5% ##################################################

readRDS(file = "~/tcga_biolinks1/stats/out.brca.dam.0.5") -> out.brca.dam.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph.dam.0.5") -> out.rand.brca.graph.dam.0.5

out.brca.dam.0.5[which(is.infinite(out.brca.dam.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.brca.graph.dam.0.5, na.rm = T),0.01) -> random.cutoff.dam.0.5.brca

dat.brca.dam.0.5 <- data.frame(Distance = c(colMeans(out.brca.dam.0.5), colMeans(out.rand.brca.graph.dam.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.brca.dam.0.5)),rep("Random",nrow(out.rand.brca.graph.dam.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.brca.dam.0.5), colMeans(out.rand.brca.graph.dam.0.5, na.rm = T)) -> ks.result.dam.0.5.brca

# Print
ks.result.dam.0.5.brca$statistic
ks.result.dam.0.5.brca$p.value
# length(ks.result2.brca$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.0.5.brca$statistic,2), ifelse(ks.result.dam.0.5.brca$p.value <0.0001, " p < 0.0001", paste0(" p = ",signif(ks.result.dam.0.5.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.0.5.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.0.5.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.brca$p.value,4)))) -> text.annot2
length(ks.result.dam.0.5.brca$data$x) -> text.annot3

dim(out.brca.dam.0.5)
dim(out.rand.brca.graph.dam.0.5)

# Create gg plot
ggplot(dat.brca.dam.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca.dam.0.5, aes(xintercept=random.cutoff.dam.0.5.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.brca.dam.0.5

gg.brca.dam.0.5

##### SIFT DELETERIOUS #########################################################

readRDS(file = "~/tcga_biolinks1/stats/out.brca.del") -> out.brca.del
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph.del") -> out.rand.brca.graph.del

out.brca.del[which(is.infinite(out.brca.del))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.brca.graph.del, na.rm = T),0.01) -> random.cutoff.del.brca

dat.brca.del <- data.frame(Distance = c(colMeans(out.brca.del), colMeans(out.rand.brca.graph.del, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.brca.del)),rep("Random",nrow(out.rand.brca.graph.del)))))

# Perform KS test on column means
ks.test(colMeans(out.brca.del), colMeans(out.rand.brca.graph.del, na.rm = T)) -> ks.result.del.brca

# Print
ks.result.del.brca$statistic
ks.result.del.brca$p.value
# length(ks.result2.brca$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.brca$statistic,2), ifelse(ks.result.del.brca$p.value <0.0001, " p < 0.0001", paste0(" p = ",signif(ks.result.del.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.brca$p.value,4)))) -> text.annot2
length(ks.result.del.brca$data$x) -> text.annot3

dim(out.brca.del)
dim(out.rand.brca.graph.del)

# Create gg plot
ggplot(dat.brca.del, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca.del, aes(xintercept=random.cutoff.del.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.brca.del

gg.brca.del
##### SIFT DELETERIOUS AND >= 0.5% ##############################################

readRDS(file = "~/tcga_biolinks1/stats/out.brca.del.0.5") -> out.brca.del.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph.del.0.5") -> out.rand.brca.graph.del.0.5

out.brca.del.0.5[which(is.infinite(out.brca.del.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.brca.graph.del.0.5, na.rm = T),0.01) -> random.cutoff.del.0.5.brca

dat.brca.del.0.5 <- data.frame(Distance = c(colMeans(out.brca.del.0.5), colMeans(out.rand.brca.graph.del.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.brca.del.0.5)),rep("Random",nrow(out.rand.brca.graph.del.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.brca.del.0.5), colMeans(out.rand.brca.graph.del.0.5, na.rm = T)) -> ks.result.del.0.5.brca

# Print
ks.result.del.0.5.brca$statistic
ks.result.del.0.5.brca$p.value
# length(ks.result2.brca$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.0.5.brca$statistic,2), ifelse(ks.result.del.0.5.brca$p.value <0.0001, " p < 0.0001", paste0(" p = ",signif(ks.result.del.0.5.brca$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.0.5.brca$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.0.5.brca$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.brca$p.value,4)))) -> text.annot2
length(ks.result.del.0.5.brca$data$x) -> text.annot3

dim(out.brca.del.0.5)
dim(out.rand.brca.graph.del.0.5)

# Create gg plot
ggplot(dat.brca.del.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.brca.del.0.5, aes(xintercept=random.cutoff.del.0.5.brca,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.brca.del.0.5

gg.brca.del.0.5
##### GGARRANGE ################################################################
library(ggpubr)
gg.brca0.5 + labs(subtitle = "Only Genes Missense Mutated in >= 0.5% of All Patients") -> gg.brca0.5t
gg.brca0.5t + theme(text = element_text(size = 10)) -> p1
gg.brca1 + labs(subtitle = "Only Genes Missense Mutated in >= 1% of All Patients") -> gg.brca1t
gg.brca1t + theme(text = element_text(size = 10)) -> p2
gg.brca1.5 + labs(subtitle = "Only Genes Missense Mutated in >= 1.5% of All Patients") -> gg.brca1.5t
gg.brca1.5t + theme(text = element_text(size = 10)) -> p3
gg.brca.prob + labs(subtitle = "Only Genes with Probably Damaging Missense Mutations") -> gg.brca.probt
gg.brca.probt + theme(text = element_text(size = 10)) ->p4
gg.brca.poss + labs(subtitle = "Only Genes with Possibly Damaging Missense Mutations") -> gg.brca.posst
gg.brca.posst + theme(text = element_text(size = 10)) -> p5
gg.brca.dam + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations") -> gg.brca.damt
gg.brca.damt + theme(text = element_text(size = 10)) -> p6
gg.brca.dam.0.5 + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations,
in >=0.5% of All Patients") -> gg.brca.dam.0.5t
gg.brca.dam.0.5t + theme(text = element_text(size = 10)) -> p7
gg.brca.del + labs(subtitle = "Only Genes with Deleterious Missense Mutations") -> gg.brca.delt
gg.brca.delt + theme(text = element_text(size = 10)) -> p8
gg.brca.del.0.5 + labs(subtitle = "Only Genes with Deleterious Missense Mutations, in >= 0.5% of All Patients") -> gg.brca.del.0.5t
gg.brca.del.0.5t + theme(text = element_text(size = 10)) -> p9

brca.figure <- ggarrange(p1, p2, p3, 
                         p4, p5, p6,
                         p7, p8, p9,
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                         common.legend = TRUE, legend = "bottom",
                         ncol = 3, nrow = 3)


annotate_figure(brca.figure, top = text_grob("Mean Gene Distance Distributions of Both Random and BRCA Genes on the BRCA Network", size = 20, face = "bold")) -> brca.plot

ggsave("~/tcga_biolinks1/Plots/brca.plot.png", plot = brca.plot, width = 18, height = 10)

