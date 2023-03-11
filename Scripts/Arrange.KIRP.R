#### ARRANGE GG PLOTS INTO FIGURE ##################################
# Create gg plots
#### 0.5% ######################################################################
readRDS(file = "~/tcga_biolinks1/stats/out.kirp0.5") -> out.kirp0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph0.5") -> out.rand.kirp.graph0.5

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
  ggtitle(NULL) -> gg.kirp0.5

gg.kirp0.5

#### 1% ########################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirp1") -> out.kirp1
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph1") -> out.rand.kirp.graph1

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
  ggtitle(NULL) -> gg.kirp1

gg.kirp1

#### 1.5% ######################################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirp1.5") -> out.kirp1.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph1.5") -> out.rand.kirp.graph1.5

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
  ggtitle(NULL) -> gg.kirp1.5

gg.kirp1.5

#### PROB DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirp.prob") -> out.kirp.prob
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph.prob") -> out.rand.kirp.graph.prob
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
  ggtitle(NULL) -> gg.kirp.prob

gg.kirp.prob
#### POSS DAMAGING #############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirp.poss") -> out.kirp.poss
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph.poss") -> out.rand.kirp.graph.poss

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
  ggtitle(NULL) -> gg.kirp.poss

gg.kirp.poss

#### ALL DAMAGING ##############################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirp.dam") -> out.kirp.dam
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph.dam") -> out.rand.kirp.graph.dam

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
  ggtitle(NULL) -> gg.kirp.dam

gg.kirp.dam

##### ALL DAMAGING AND >=0.5% ##################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirp.dam.0.5") -> out.kirp.dam.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph.dam.0.5") -> out.rand.kirp.graph.dam.0.5

out.kirp.dam.0.5[which(is.infinite(out.kirp.dam.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirp.graph.dam.0.5, na.rm = T),0.01) -> random.cutoff.dam.0.5.kirp

dat.kirp.dam.0.5 <- data.frame(Distance = c(colMeans(out.kirp.dam.0.5), colMeans(out.rand.kirp.graph.dam.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.kirp.dam.0.5)),rep("Random",nrow(out.rand.kirp.graph.dam.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.kirp.dam.0.5), colMeans(out.rand.kirp.graph.dam.0.5, na.rm = T)) -> ks.result.dam.0.5.kirp

# Print
ks.result.dam.0.5.kirp$statistic
ks.result.dam.0.5.kirp$p.value
# length(ks.result2.kirp$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.dam.0.5.kirp$statistic,2), ifelse(ks.result.dam.0.5.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.kirp$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.dam.0.5.kirp$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.dam.0.5.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.dam.0.5.kirp$p.value,4)))) -> text.annot2
length(ks.result.dam.0.5.kirp$data$x) -> text.annot3

dim(out.kirp.dam.0.5)
dim(out.rand.kirp.graph.dam.0.5)

# Create gg plot
ggplot(dat.kirp.dam.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirp.dam.0.5, aes(xintercept=random.cutoff.dam.0.5.kirp,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.kirp.dam.0.5

gg.kirp.dam.0.5

##### SIFT DELETERIOUS #########################################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirp.del") -> out.kirp.del
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph.del") -> out.rand.kirp.graph.del

out.kirp.del[which(is.infinite(out.kirp.del))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirp.graph.del, na.rm = T),0.01) -> random.cutoff.del.kirp

dat.kirp.del <- data.frame(Distance = c(colMeans(out.kirp.del), colMeans(out.rand.kirp.graph.del, na.rm = T)),
                           Distribution = factor(c(rep("Real",nrow(out.kirp.del)),rep("Random",nrow(out.rand.kirp.graph.del)))))

# Perform KS test on column means
ks.test(colMeans(out.kirp.del), colMeans(out.rand.kirp.graph.del, na.rm = T)) -> ks.result.del.kirp

# Print
ks.result.del.kirp$statistic
ks.result.del.kirp$p.value
# length(ks.result2.kirp$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.kirp$statistic,2), ifelse(ks.result.del.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.kirp$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.kirp$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.kirp$p.value,4)))) -> text.annot2
length(ks.result.del.kirp$data$x) -> text.annot3

dim(out.kirp.del)
dim(out.rand.kirp.graph.del)

# Create gg plot
ggplot(dat.kirp.del, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirp.del, aes(xintercept=random.cutoff.del.kirp,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.kirp.del

gg.kirp.del
##### SIFT DELETERIOUS AND >= 0.5% ##############################################

readRDS(file = "~/tcga_biolinks1/stats/out.kirp.del.0.5") -> out.kirp.del.0.5
readRDS(file = "~/tcga_biolinks1/stats/out.rand.kirp.graph.del.0.5") -> out.rand.kirp.graph.del.0.5

out.kirp.del.0.5[which(is.infinite(out.kirp.del.0.5))] <- NA

# Create a random cutoff
quantile(colMeans(out.rand.kirp.graph.del.0.5, na.rm = T),0.01) -> random.cutoff.del.0.5.kirp

dat.kirp.del.0.5 <- data.frame(Distance = c(colMeans(out.kirp.del.0.5), colMeans(out.rand.kirp.graph.del.0.5, na.rm = T)),
                               Distribution = factor(c(rep("Real",nrow(out.kirp.del.0.5)),rep("Random",nrow(out.rand.kirp.graph.del.0.5)))))

# Perform KS test on column means
ks.test(colMeans(out.kirp.del.0.5), colMeans(out.rand.kirp.graph.del.0.5, na.rm = T)) -> ks.result.del.0.5.kirp

# Print
ks.result.del.0.5.kirp$statistic
ks.result.del.0.5.kirp$p.value
# length(ks.result2.kirp$data$x)


# Create a text annotation of results
paste0("D = ",round(ks.result.del.0.5.kirp$statistic,2), ifelse(ks.result.del.0.5.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.kirp$p.value,2)))) -> text.annot
paste0("D = ",round(ks.result.del.0.5.kirp$statistic,2)) -> text.annot1
paste0(ifelse(ks.result.del.0.5.kirp$p.value==0, " p < 0.0001", paste0("p = ",round(ks.result.del.0.5.kirp$p.value,4)))) -> text.annot2
length(ks.result.del.0.5.kirp$data$x) -> text.annot3

dim(out.kirp.del.0.5)
dim(out.rand.kirp.graph.del.0.5)

# Create gg plot
ggplot(dat.kirp.del.0.5, aes(x=Distance, colour=Distribution)) +
  geom_density() +
  geom_vline(data=dat.kirp.del.0.5, aes(xintercept=random.cutoff.del.0.5.kirp,  colour=Distribution),
             linetype="dashed", linewidth=1) +
  annotate("text", x = 4.45, y = 1.75, label = text.annot) +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Mean Distance") +
  labs(y = "Empirical Density") +
  theme_minimal() +
  ggtitle(NULL) -> gg.kirp.del.0.5

gg.kirp.del.0.5
##### GGARRANGE ################################################################
library(ggpubr)
gg.kirp0.5 + labs(subtitle = "Only Genes Missense Mutated in >= 0.5% of All Patients") -> gg.kirp0.5t
gg.kirp.0.5t + theme(text = element_text(size = 10)) -> p1
gg.kirp1 + labs(subtitle = "Only Genes Missense Mutated in >= 1% of All Patients") -> gg.kirp1t
gg.kirp1t + theme(text = element_text(size = 10)) -> p2
gg.kirp1.5 + labs(subtitle = "Only Genes Missense Mutated in >= 1.5% of All Patients") -> gg.kirp1.5t
gg.kirp1.5t + theme(text = element_text(size = 10)) -> p3
gg.kirp.prob + labs(subtitle = "Only Genes with Probably Damaging Missense Mutations") -> gg.kirp.probt
gg.kirp.probt + theme(text = element_text(size = 10)) ->p4
gg.kirp.poss + labs(subtitle = "Only Genes with Possibly Damaging Missense Mutations") -> gg.kirp.posst
gg.kirp.posst + theme(text = element_text(size = 10)) -> p5
gg.kirp.dam + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations") -> gg.kirp.damt
gg.kirp.damt + theme(text = element_text(size = 10)) -> p6
gg.kirp.dam.0.5 + labs(subtitle = "Only Genes with Either Probably or Possibly Damaging Missense Mutations,
in >=0.5% of All Patients") -> gg.kirp.dam.0.5t
gg.kirp.dam.0.5t + theme(text = element_text(size = 10)) -> p7
gg.kirp.del + labs(subtitle = "Only Genes with Deleterious Missense Mutations") -> gg.kirp.delt
gg.kirp.delt + theme(text = element_text(size = 10)) -> p8
gg.kirp.del.0.5 + labs(subtitle = "Only Genes with Deleterious Missense Mutations, in >0.5% of All Patients") -> gg.kirp.del.0.5t
gg.kirp.del.0.5t + theme(text = element_text(size = 10)) -> p9

kirp.figure <- ggarrange(p1, p2, p3, 
                         p4, p5, p6,
                         p7, p8, p9,
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                         common.legend = TRUE, legend = "bottom",
                         ncol = 3, nrow = 3)
annotate_figure(kirp.figure, top = text_grob("Mean Gene Distance Distributions of Both Random and KIRP Genes on the KIRP Network", size = 20, face = "bold")) -> kirp.plot

ggsave("~/tcga_biolinks1/Plots/kirp.plot.png", plot = kirp.plot, width = 18, height = 10)

