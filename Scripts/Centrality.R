#### there are different measures of centrality that you could test if different between categories

# ## Stress
# wilcox.test(Stress~(PolyPhen.Mean), data = luad.graph.df.filt0.5)
# boxplot(Stress~(PolyPhen.Mean >0.908), data = luad.graph.df.filt0.5 )
# tapply(luad.graph.df.filt0.5$Stress, (luad.graph.df.filt0.5$PolyPhen.Mean >0.908), median)
# 
# wilcox.test(Stress~PolyPhen.Probably.Damaging, data = luad.graph.df.filt0.5)
# boxplot(Stress~PolyPhen.Probably.Damaging, data = luad.graph.df.filt0.5 )
# tapply(luad.graph.df.filt0.5$Stress, luad.graph.df.filt0.5$PolyPhen.Probably.Damaging, median)
# 
# ## Closeness Centrality
# wilcox.test(ClosenessCentrality~(PolyPhen.Mean >0.908), data = luad.graph.df.filt0.5)
# boxplot(ClosenessCentrality~(PolyPhen.Mean >0.908), data = luad.graph.df.filt0.5 )
# tapply(luad.graph.df.filt0.5$ClosenessCentrality, (luad.graph.df.filt0.5$PolyPhen.Mean >0.908), median)
# 
# wilcox.test(ClosenessCentrality~PolyPhen.Probably.Damaging, data = luad.graph.df.filt0.5)
# boxplot(ClosenessCentrality~PolyPhen.Probably.Damaging, data = luad.graph.df.filt0.5 )
# tapply(luad.graph.df.filt0.5$ClosenessCentrality, luad.graph.df.filt0.5$PolyPhen.Probably.Damaging, median)
# 
# ## Degree
# wilcox.test(Degree~(PolyPhen.Mean >0.908), data = luad.graph.df.filt0.5)
# boxplot(Degree~(PolyPhen.Mean >0.908), data = luad.graph.df.filt0.5 )
# tapply(luad.graph.df.filt0.5$Degree, (luad.graph.df.filt0.5$PolyPhen.Mean >0.908), median)
# 
# wilcox.test(ClosenessCentrality~PolyPhen.Probably.Damaging, data = luad.graph.df.filt0.5)
# boxplot(ClosenessCentrality~PolyPhen.Probably.Damaging, data = luad.graph.df.filt0.5 )
# tapply(luad.graph.df.filt0.5$ClosenessCentrality, luad.graph.df.filt0.5$PolyPhen.Probably.Damaging, median)

#### LUAD ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/luad4.node.csv") -> luad.graph.df
luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 0.5) -> luad.graph.df.filt0.5

dim(luad.graph.df.filt0.5)
## Wilcoxon Rank Sum Test with continuity correction/Mann-Whitney Test
wilcox.test(Stress~PolyPhen.Probably.Damaging, data = luad.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~PolyPhen.Probably.Damaging, data = luad.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~PolyPhen.Probably.Damaging, data = luad.graph.df.filt0.5)
wilcox.test(Stress~SIFT.Deleterious, data = luad.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~SIFT.Deleterious, data = luad.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~SIFT.Deleterious, data = luad.graph.df.filt0.5)

## Spearman Correlation Test
cor.test(luad.graph.df.filt0.5$PolyPhen.Mean, log(luad.graph.df.filt0.5$Stress), method = "spearman")
cor.test(luad.graph.df.filt0.5$PolyPhen.Mean, log(luad.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(luad.graph.df.filt0.5$PolyPhen.Mean, log(luad.graph.df.filt0.5$ClosenessCentrality), method = "spearman")
cor.test(luad.graph.df.filt0.5$SIFT.Mean, log(luad.graph.df.filt0.5$Stress), method = "spearman")
cor.test(luad.graph.df.filt0.5$SIFT.Mean, log(luad.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(luad.graph.df.filt0.5$SIFT.Mean, log(luad.graph.df.filt0.5$ClosenessCentrality), method = "spearman")

# lm(luad.graph.df$Stress~luad.graph.df$SIFT.Mean, na.rm = T)
# plot(luad.graph.df$SIFT.Mean, luad.graph.df$Stress)


####################


data.frame(luad.graph.df.filt0.5$PolyPhen.Mean, luad.graph.df.filt0.5$BetweennessCentrality) -> luad.c1
ggplot(luad.c1, aes(x=luad.graph.df.filt0.5.PolyPhen.Mean, y=log(luad.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 15)) +
  ggtitle("PolyPhen Mean Against Betweeness Centrality") -> luad.c1.plot
luad.c1.plot

data.frame(luad.graph.df.filt0.5$PolyPhen.Mean, luad.graph.df.filt0.5$Stress) -> luad.c2
ggplot(luad.c2, aes(x=luad.graph.df.filt0.5.PolyPhen.Mean, y=log(luad.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Stress Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 15)) +
    ggtitle("PolyPhen Mean Against Stress Centrality") -> luad.c2.plot
luad.c2.plot

data.frame(luad.graph.df.filt0.5$PolyPhen.Mean, luad.graph.df.filt0.5$ClosenessCentrality) -> luad.c3
ggplot(luad.c3, aes(x=luad.graph.df.filt0.5.PolyPhen.Mean, y=log(luad.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 15)) +
  ggtitle("PolyPhen Mean Against Closeness Centrality") -> luad.c3.plot
luad.c3.plot

data.frame(luad.graph.df.filt0.5$SIFT.Mean, luad.graph.df.filt0.5$BetweennessCentrality) -> luad.c4
ggplot(luad.c4, aes(x=luad.graph.df.filt0.5.SIFT.Mean, y=log(luad.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 15)) +
  ggtitle("SIFT Mean Against Betweeness Centrality") -> luad.c4.plot
luad.c4.plot

data.frame(luad.graph.df.filt0.5$SIFT.Mean, luad.graph.df.filt0.5$Stress) -> luad.c5
ggplot(luad.c5, aes(x=luad.graph.df.filt0.5.SIFT.Mean, y=log(luad.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Stress Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 15)) +
  ggtitle("SIFT Mean Against Stress Centrality") -> luad.c5.plot
luad.c5.plot

data.frame(luad.graph.df.filt0.5$SIFT.Mean, luad.graph.df.filt0.5$ClosenessCentrality) -> luad.c6
ggplot(luad.c6, aes(x=luad.graph.df.filt0.5.SIFT.Mean, y=log(luad.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 15)) +
  ggtitle("SIFT Mean Against Closeness Centrality") -> luad.c6.plot
luad.c6.plot


luad.c.figure <- ggarrange(luad.c2.plot, luad.c1.plot, luad.c3.plot,luad.c5.plot, luad.c4.plot, luad.c6.plot,
                         labels = c("A", "B", "C", "D", "E", "F"),
                         common.legend = TRUE, legend = "bottom",
                         ncol = 2, nrow = 3)
luad.c.figure
annotate_figure(luad.c.figure, top = text_grob("Lung Adenocarcinoma Network Parameters Against SIFT and PolyPhen Means", size = 20, face = "bold")) -> luad.centrality

ggsave("~/tcga_biolinks1/Plots/luad.centrality.png", plot = luad.centrality, width = 12, height = 16)

#### PRAD ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/prad5.node.csv") -> prad.graph.df
prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> prad.graph.df.filt0.5

dim(prad.graph.df.filt0.5)
## Wilcoxon Rank Sum Test with continuity correction/Mann-Whitney Test
wilcox.test(Stress~PolyPhen.Probably.Damaging, data = prad.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~PolyPhen.Probably.Damaging, data = prad.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~PolyPhen.Probably.Damaging, data = prad.graph.df.filt0.5)
wilcox.test(Stress~SIFT.Deleterious, data = prad.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~SIFT.Deleterious, data = prad.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~SIFT.Deleterious, data = prad.graph.df.filt0.5)

## Spearman Correlation Test
cor.test(prad.graph.df.filt0.5$PolyPhen.Mean, log(prad.graph.df.filt0.5$Stress), method = "spearman")
cor.test(prad.graph.df.filt0.5$PolyPhen.Mean, log(prad.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(prad.graph.df.filt0.5$PolyPhen.Mean, log(prad.graph.df.filt0.5$ClosenessCentrality), method = "spearman")
cor.test(prad.graph.df.filt0.5$SIFT.Mean, log(prad.graph.df.filt0.5$Stress), method = "spearman")
cor.test(prad.graph.df.filt0.5$SIFT.Mean, log(prad.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(prad.graph.df.filt0.5$SIFT.Mean, log(prad.graph.df.filt0.5$ClosenessCentrality), method = "spearman")

# lm(prad.graph.df$Stress~prad.graph.df$SIFT.Mean, na.rm = T)
# plot(prad.graph.df$SIFT.Mean, prad.graph.df$Stress)


# Create Plots
data.frame(prad.graph.df.filt0.5$PolyPhen.Mean, prad.graph.df.filt0.5$BetweennessCentrality) -> prad.c1
ggplot(prad.c1, aes(x=prad.graph.df.filt0.5.PolyPhen.Mean, y=log(prad.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Prostate Adenocarcinoma Genes Against Betweeness Centrality") -> prad.c1.plot
prad.c1.plot

data.frame(prad.graph.df.filt0.5$PolyPhen.Mean, prad.graph.df.filt0.5$Stress) -> prad.c2
ggplot(prad.c2, aes(x=prad.graph.df.filt0.5.PolyPhen.Mean, y=log(prad.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Prostate Adenocarcinoma Genes Against Stress") -> prad.c2.plot
prad.c2.plot

data.frame(prad.graph.df.filt0.5$PolyPhen.Mean, prad.graph.df.filt0.5$ClosenessCentrality) -> prad.c3
ggplot(prad.c3, aes(x=prad.graph.df.filt0.5.PolyPhen.Mean, y=log(prad.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Prostate Adenocarcinoma Genes Against Closeness Centrality") -> prad.c3.plot
prad.c3.plot

data.frame(prad.graph.df.filt0.5$SIFT.Mean, prad.graph.df.filt0.5$BetweennessCentrality) -> prad.c4
ggplot(prad.c4, aes(x=prad.graph.df.filt0.5.SIFT.Mean, y=log(prad.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Prostate Adenocarcinoma Genes Against Betweeness Centrality") -> prad.c4.plot
prad.c4.plot

data.frame(prad.graph.df.filt0.5$SIFT.Mean, prad.graph.df.filt0.5$Stress) -> prad.c5
ggplot(prad.c5, aes(x=prad.graph.df.filt0.5.SIFT.Mean, y=log(prad.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Prostate Adenocarcinoma Genes Against Stress") -> prad.c5.plot
prad.c5.plot

data.frame(prad.graph.df.filt0.5$SIFT.Mean, prad.graph.df.filt0.5$ClosenessCentrality) -> prad.c6
ggplot(prad.c6, aes(x=prad.graph.df.filt0.5.SIFT.Mean, y=log(prad.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Prostate Adenocarcinoma Genes Against Closeness Centrality") -> prad.c6.plot
prad.c6.plot

prad.c.figure <- ggarrange(prad.c1.plot, prad.c2.plot, prad.c3.plot, prad.c4.plot, prad.c5.plot, prad.c6.plot,
                           labels = c("A", "B", "C", "D", "E", "F"),
                           common.legend = TRUE, legend = "bottom",
                           ncol = 3, nrow = 2)
annotate_figure(prad.c.figure, top = text_grob("Prostate Adenocarcinoma Network Parameters Against SIFT and PolyPhen Means", size = 20, face = "bold")) -> prad.centrality

ggsave("~/tcga_biolinks1/Plots/prad.centrality.png", plot = prad.centrality, width = 18, height = 10)

#### BRCA ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/brca5.node.csv") -> brca.graph.df
brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> brca.graph.df.filt0.5

dim(brca.graph.df.filt0.5)
## Wilcoxon Rank Sum Test with continuity correction/Mann-Whitney Test
wilcox.test(Stress~PolyPhen.Probably.Damaging, data = brca.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~PolyPhen.Probably.Damaging, data = brca.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~PolyPhen.Probably.Damaging, data = brca.graph.df.filt0.5)
wilcox.test(Stress~SIFT.Deleterious, data = brca.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~SIFT.Deleterious, data = brca.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~SIFT.Deleterious, data = brca.graph.df.filt0.5)

## Spearman Correlation Test
cor.test(brca.graph.df.filt0.5$PolyPhen.Mean, log(brca.graph.df.filt0.5$Stress), method = "spearman")
cor.test(brca.graph.df.filt0.5$PolyPhen.Mean, log(brca.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(brca.graph.df.filt0.5$PolyPhen.Mean, log(brca.graph.df.filt0.5$ClosenessCentrality), method = "spearman")
cor.test(brca.graph.df.filt0.5$SIFT.Mean, log(brca.graph.df.filt0.5$Stress), method = "spearman")
cor.test(brca.graph.df.filt0.5$SIFT.Mean, log(brca.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(brca.graph.df.filt0.5$SIFT.Mean, log(brca.graph.df.filt0.5$ClosenessCentrality), method = "spearman")

# lm(brca.graph.df$Stress~brca.graph.df$SIFT.Mean, na.rm = T)
# plot(brca.graph.df$SIFT.Mean, brca.graph.df$Stress)


####################
data.frame(brca.graph.df.filt0.5$PolyPhen.Mean, brca.graph.df.filt0.5$BetweennessCentrality) -> brca.c1
ggplot(brca.c1, aes(x=brca.graph.df.filt0.5.PolyPhen.Mean, y=log(brca.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Breast Invasive Carcinoma Genes Against Betweeness Centrality") -> brca.c1.plot
brca.c1.plot

data.frame(brca.graph.df.filt0.5$PolyPhen.Mean, brca.graph.df.filt0.5$Stress) -> brca.c2
ggplot(brca.c2, aes(x=brca.graph.df.filt0.5.PolyPhen.Mean, y=log(brca.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Breast Invasive Carcinoma Genes Against Stress") -> brca.c2.plot
brca.c2.plot

data.frame(brca.graph.df.filt0.5$PolyPhen.Mean, brca.graph.df.filt0.5$ClosenessCentrality) -> brca.c3
ggplot(brca.c3, aes(x=brca.graph.df.filt0.5.PolyPhen.Mean, y=log(brca.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Breast Invasive Carcinoma Genes Against Closeness Centrality") -> brca.c3.plot
brca.c3.plot

data.frame(brca.graph.df.filt0.5$SIFT.Mean, brca.graph.df.filt0.5$BetweennessCentrality) -> brca.c4
ggplot(brca.c4, aes(x=brca.graph.df.filt0.5.SIFT.Mean, y=log(brca.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Breast Invasive Carcinoma Genes Against Betweeness Centrality") -> brca.c4.plot
brca.c4.plot

data.frame(brca.graph.df.filt0.5$SIFT.Mean, brca.graph.df.filt0.5$Stress) -> brca.c5
ggplot(brca.c5, aes(x=brca.graph.df.filt0.5.SIFT.Mean, y=log(brca.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Breast Invasive Carcinoma Genes Against Stress") -> brca.c5.plot
brca.c5.plot

data.frame(brca.graph.df.filt0.5$SIFT.Mean, brca.graph.df.filt0.5$ClosenessCentrality) -> brca.c6
ggplot(brca.c6, aes(x=brca.graph.df.filt0.5.SIFT.Mean, y=log(brca.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Breast Invasive Carcinoma Against Closeness Centrality") -> brca.c6.plot
brca.c6.plot

brca.c.figure <- ggarrange(brca.c1.plot, brca.c2.plot, brca.c3.plot, brca.c4.plot, brca.c5.plot, brca.c6.plot,
                           labels = c("A", "B", "C", "D", "E", "F"),
                           common.legend = TRUE, legend = "bottom",
                           ncol = 3, nrow = 2)
annotate_figure(brca.c.figure, top = text_grob("Breast Invasive Carcinoma Network Parameters Against SIFT and PolyPhen Means", size = 20, face = "bold")) -> brca.centrality

ggsave("~/tcga_biolinks1/Plots/brca.centrality.png", plot = brca.centrality, width = 18, height = 10)


#### STAD ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/stad5.node.csv") -> stad.graph.df
stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> stad.graph.df.filt0.5

dim(stad.graph.df.filt0.5)
## Wilcoxon Rank Sum Test with continuity correction/Mann-Whitney Test
wilcox.test(Stress~PolyPhen.Probably.Damaging, data = stad.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~PolyPhen.Probably.Damaging, data = stad.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~PolyPhen.Probably.Damaging, data = stad.graph.df.filt0.5)
wilcox.test(Stress~SIFT.Deleterious, data = stad.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~SIFT.Deleterious, data = stad.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~SIFT.Deleterious, data = stad.graph.df.filt0.5)

## Spearman Correlation Test
cor.test(stad.graph.df.filt0.5$PolyPhen.Mean, log(stad.graph.df.filt0.5$Stress), method = "spearman")
cor.test(stad.graph.df.filt0.5$PolyPhen.Mean, log(stad.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(stad.graph.df.filt0.5$PolyPhen.Mean, log(stad.graph.df.filt0.5$ClosenessCentrality), method = "spearman")
cor.test(stad.graph.df.filt0.5$SIFT.Mean, log(stad.graph.df.filt0.5$Stress), method = "spearman")
cor.test(stad.graph.df.filt0.5$SIFT.Mean, log(stad.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(stad.graph.df.filt0.5$SIFT.Mean, log(stad.graph.df.filt0.5$ClosenessCentrality), method = "spearman")

# lm(stad.graph.df$Stress~stad.graph.df$SIFT.Mean, na.rm = T)
# plot(stad.graph.df$SIFT.Mean, stad.graph.df$Stress)


####################
data.frame(stad.graph.df.filt0.5$PolyPhen.Mean, stad.graph.df.filt0.5$BetweennessCentrality) -> stad.c1
ggplot(stad.c1, aes(x=stad.graph.df.filt0.5.PolyPhen.Mean, y=log(stad.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Stomach Adenocarcinoma Genes Against Betweeness Centrality") -> stad.c1.plot
stad.c1.plot

data.frame(stad.graph.df.filt0.5$PolyPhen.Mean, stad.graph.df.filt0.5$Stress) -> stad.c2
ggplot(stad.c2, aes(x=stad.graph.df.filt0.5.PolyPhen.Mean, y=log(stad.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Stomach Adenocarcinoma Genes Against Stress") -> stad.c2.plot
stad.c2.plot

data.frame(stad.graph.df.filt0.5$PolyPhen.Mean, stad.graph.df.filt0.5$ClosenessCentrality) -> stad.c3
ggplot(stad.c3, aes(x=stad.graph.df.filt0.5.PolyPhen.Mean, y=log(stad.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Stomach Adenocarcinoma Genes Against Closeness Centrality") -> stad.c3.plot
stad.c3.plot

data.frame(stad.graph.df.filt0.5$SIFT.Mean, stad.graph.df.filt0.5$BetweennessCentrality) -> stad.c4
ggplot(stad.c4, aes(x=stad.graph.df.filt0.5.SIFT.Mean, y=log(stad.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Stomach Adenocarcinoma Genes Against Betweeness Centrality") -> stad.c4.plot
stad.c4.plot

data.frame(stad.graph.df.filt0.5$SIFT.Mean, stad.graph.df.filt0.5$Stress) -> stad.c5
ggplot(stad.c5, aes(x=stad.graph.df.filt0.5.SIFT.Mean, y=log(stad.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Stomach Adenocarcinoma Genes Against Stress") -> stad.c5.plot
stad.c5.plot

data.frame(stad.graph.df.filt0.5$SIFT.Mean, stad.graph.df.filt0.5$ClosenessCentrality) -> stad.c6
ggplot(stad.c6, aes(x=stad.graph.df.filt0.5.SIFT.Mean, y=log(stad.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Stomach Adenocarcinoma Genes Against Closeness Centrality") -> stad.c6.plot
stad.c6.plot

stad.c.figure <- ggarrange(stad.c1.plot, stad.c2.plot, stad.c3.plot, stad.c4.plot, stad.c5.plot, stad.c6.plot,
                           labels = c("A", "B", "C", "D", "E", "F"),
                           common.legend = TRUE, legend = "bottom",
                           ncol = 3, nrow = 2)
annotate_figure(stad.c.figure, top = text_grob("Stomach Adenocarcinoma Network Parameters Against SIFT and PolyPhen Means", size = 20, face = "bold")) -> stad.centrality

ggsave("~/tcga_biolinks1/Plots/stad.centrality.png", plot = stad.centrality, width = 18, height = 10)

#### KIRP ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/kirp5.node.csv") -> kirp.graph.df
kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutations >= 0.5) -> kirp.graph.df.filt0.5

dim(kirp.graph.df.filt0.5)
## Wilcoxon Rank Sum Test with continuity correction/Mann-Whitney Test
wilcox.test(Stress~PolyPhen.Probably.Damaging, data = kirp.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~PolyPhen.Probably.Damaging, data = kirp.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~PolyPhen.Probably.Damaging, data = kirp.graph.df.filt0.5)
wilcox.test(Stress~SIFT.Deleterious, data = kirp.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~SIFT.Deleterious, data = kirp.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~SIFT.Deleterious, data = kirp.graph.df.filt0.5)

## Spearman Correlation Test
cor.test(kirp.graph.df.filt0.5$PolyPhen.Mean, log(kirp.graph.df.filt0.5$Stress), method = "spearman")
cor.test(kirp.graph.df.filt0.5$PolyPhen.Mean, log(kirp.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(kirp.graph.df.filt0.5$PolyPhen.Mean, log(kirp.graph.df.filt0.5$ClosenessCentrality), method = "spearman")
cor.test(kirp.graph.df.filt0.5$SIFT.Mean, log(kirp.graph.df.filt0.5$Stress), method = "spearman")
cor.test(kirp.graph.df.filt0.5$SIFT.Mean, log(kirp.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(kirp.graph.df.filt0.5$SIFT.Mean, log(kirp.graph.df.filt0.5$ClosenessCentrality), method = "spearman")

# lm(kirp.graph.df$Stress~kirp.graph.df$SIFT.Mean, na.rm = T)
# plot(kirp.graph.df$SIFT.Mean, kirp.graph.df$Stress)


####################
data.frame(kirp.graph.df.filt0.5$PolyPhen.Mean, kirp.graph.df.filt0.5$BetweennessCentrality) -> kirp.c1
ggplot(kirp.c1, aes(x=kirp.graph.df.filt0.5.PolyPhen.Mean, y=log(kirp.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Kidney Renal Papillary Cell Carcinoma Genes Against Betweeness Centrality") -> kirp.c1.plot
kirp.c1.plot

data.frame(kirp.graph.df.filt0.5$PolyPhen.Mean, kirp.graph.df.filt0.5$Stress) -> kirp.c2
ggplot(kirp.c2, aes(x=kirp.graph.df.filt0.5.PolyPhen.Mean, y=log(kirp.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Kidney Renal Papillary Cell Carcinomaa Genes Against Stress") -> kirp.c2.plot
kirp.c2.plot

data.frame(kirp.graph.df.filt0.5$PolyPhen.Mean, kirp.graph.df.filt0.5$ClosenessCentrality) -> kirp.c3
ggplot(kirp.c3, aes(x=kirp.graph.df.filt0.5.PolyPhen.Mean, y=log(kirp.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Kidney Renal Papillary Cell Carcinoma Genes Against Closeness Centrality") -> kirp.c3.plot
kirp.c3.plot

data.frame(kirp.graph.df.filt0.5$SIFT.Mean, kirp.graph.df.filt0.5$BetweennessCentrality) -> kirp.c4
ggplot(kirp.c4, aes(x=kirp.graph.df.filt0.5.SIFT.Mean, y=log(kirp.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Kidney Renal Papillary Cell Carcinoma Genes Against Betweeness Centrality") -> kirp.c4.plot
kirp.c4.plot

data.frame(kirp.graph.df.filt0.5$SIFT.Mean, kirp.graph.df.filt0.5$Stress) -> kirp.c5
ggplot(kirp.c5, aes(x=kirp.graph.df.filt0.5.SIFT.Mean, y=log(kirp.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Kidney Renal Papillary Cell Carcinoma Genes Against Stress") -> kirp.c5.plot
kirp.c5.plot

data.frame(kirp.graph.df.filt0.5$SIFT.Mean, kirp.graph.df.filt0.5$ClosenessCentrality) -> kirp.c6
ggplot(kirp.c6, aes(x=kirp.graph.df.filt0.5.SIFT.Mean, y=log(kirp.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Kidney Renal Papillary Cell Carcinoma Genes Against Closeness Centrality") -> kirp.c6.plot
kirp.c6.plot

kirp.c.figure <- ggarrange(kirp.c1.plot, kirp.c2.plot, kirp.c3.plot, kirp.c4.plot, kirp.c5.plot, kirp.c6.plot,
                           labels = c("A", "B", "C", "D", "E", "F"),
                           common.legend = TRUE, legend = "bottom",
                           ncol = 3, nrow = 2)
annotate_figure(kirp.c.figure, top = text_grob("Kidney Renal Papillary Cell Carcinoma Network Parameters Against SIFT and PolyPhen Means", size = 20, face = "bold")) -> kirp.centrality

ggsave("~/tcga_biolinks1/Plots/kirp.centrality.png", plot = kirp.centrality, width = 18, height = 10)



#### KIRC ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/kirc4.node.csv") -> kirc.graph.df
kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> kirc.graph.df.filt0.5

dim(kirc.graph.df.filt0.5)
## Wilcoxon Rank Sum Test with continuity correction/Mann-Whitney Test
wilcox.test(Stress~PolyPhen.Probably.Damaging, data = kirc.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~PolyPhen.Probably.Damaging, data = kirc.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~PolyPhen.Probably.Damaging, data = kirc.graph.df.filt0.5)
wilcox.test(Stress~SIFT.Deleterious, data = kirc.graph.df.filt0.5)
wilcox.test(BetweennessCentrality~SIFT.Deleterious, data = kirc.graph.df.filt0.5)
wilcox.test(ClosenessCentrality~SIFT.Deleterious, data = kirc.graph.df.filt0.5)

## Spearman Correlation Test
cor.test(kirc.graph.df.filt0.5$PolyPhen.Mean, log(kirc.graph.df.filt0.5$Stress), method = "spearman")
cor.test(kirc.graph.df.filt0.5$PolyPhen.Mean, log(kirc.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(kirc.graph.df.filt0.5$PolyPhen.Mean, log(kirc.graph.df.filt0.5$ClosenessCentrality), method = "spearman")
cor.test(kirc.graph.df.filt0.5$SIFT.Mean, log(kirc.graph.df.filt0.5$Stress), method = "spearman")
cor.test(kirc.graph.df.filt0.5$SIFT.Mean, log(kirc.graph.df.filt0.5$BetweennessCentrality), method = "spearman")
cor.test(kirc.graph.df.filt0.5$SIFT.Mean, log(kirc.graph.df.filt0.5$ClosenessCentrality), method = "spearman")

# lm(kirc.graph.df$Stress~kirc.graph.df$SIFT.Mean, na.rm = T)
# plot(kirc.graph.df$SIFT.Mean, kirc.graph.df$Stress)


####################
data.frame(kirc.graph.df.filt0.5$PolyPhen.Mean, kirc.graph.df.filt0.5$BetweennessCentrality) -> kirc.c1
ggplot(kirc.c1, aes(x=kirc.graph.df.filt0.5.PolyPhen.Mean, y=log(kirc.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Kidney Renal Clear Cell Carcinoma Genes Against Betweeness Centrality") -> kirc.c1.plot
kirc.c1.plot

data.frame(kirc.graph.df.filt0.5$PolyPhen.Mean, kirc.graph.df.filt0.5$Stress) -> kirc.c2
ggplot(kirc.c2, aes(x=kirc.graph.df.filt0.5.PolyPhen.Mean, y=log(kirc.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Kidney Renal Clear Cell Carcinomaa Genes Against Stress") -> kirc.c2.plot
kirc.c2.plot

data.frame(kirc.graph.df.filt0.5$PolyPhen.Mean, kirc.graph.df.filt0.5$ClosenessCentrality) -> kirc.c3
ggplot(kirc.c3, aes(x=kirc.graph.df.filt0.5.PolyPhen.Mean, y=log(kirc.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "PolyPhen Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("PolyPhen Mean of Kidney Renal Clear Cell Carcinoma Genes Against Closeness Centrality") -> kirc.c3.plot
kirc.c3.plot

data.frame(kirc.graph.df.filt0.5$SIFT.Mean, kirc.graph.df.filt0.5$BetweennessCentrality) -> kirc.c4
ggplot(kirc.c4, aes(x=kirc.graph.df.filt0.5.SIFT.Mean, y=log(kirc.graph.df.filt0.5.BetweennessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Betweeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Kidney Renal Clear Cell Carcinoma Genes Against Betweeness Centrality") -> kirc.c4.plot
kirc.c4.plot

data.frame(kirc.graph.df.filt0.5$SIFT.Mean, kirc.graph.df.filt0.5$Stress) -> kirc.c5
ggplot(kirc.c5, aes(x=kirc.graph.df.filt0.5.SIFT.Mean, y=log(kirc.graph.df.filt0.5.Stress))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Stress)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Kidney Renal Clear Cell Carcinoma Genes Against Stress") -> kirc.c5.plot
kirc.c5.plot

data.frame(kirc.graph.df.filt0.5$SIFT.Mean, kirc.graph.df.filt0.5$ClosenessCentrality) -> kirc.c6
ggplot(kirc.c6, aes(x=kirc.graph.df.filt0.5.SIFT.Mean, y=log(kirc.graph.df.filt0.5.ClosenessCentrality))) +
  geom_point() +
  labs(x = "SIFT Mean") +
  labs(y = "Log(Closeness Centrality)") +
  theme_minimal() +
  geom_smooth(method = "loess", size = 1.5, color = "red") +
  theme(text = element_text(size = 8)) +
  ggtitle("SIFT Mean of Kidney Renal Clear Cell Carcinoma Genes Against Closeness Centrality") -> kirc.c6.plot
kirc.c6.plot

kirc.c.figure <- ggarrange(kirc.c1.plot, kirc.c2.plot, kirc.c3.plot, kirc.c4.plot, kirc.c5.plot, kirc.c6.plot,
                           labels = c("A", "B", "C", "D", "E", "F"),
                           common.legend = TRUE, legend = "bottom",
                           ncol = 3, nrow = 2)
annotate_figure(kirc.c.figure, top = text_grob("Kidney Renal Clear Cell Carcinoma Network Parameters Against SIFT and PolyPhen Means", size = 20, face = "bold")) -> kirc.centrality

ggsave("~/tcga_biolinks1/Plots/kirc.centrality.png", plot = kirc.centrality, width = 18, height = 10)
