# Rank % of patients with Gene for 0.5%

#### BRCA ######################################################################

brca.graph <- readRDS(file = "~/tcga_biolinks1/RDS/brca.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/brca5.node.csv") -> brca.graph.df
brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> brca.graph.df.filt0.5

rank(-brca.graph.df.filt0.5$Percent.of.Patients.With.Mutation) -> brca.graph.df.filt0.5$rank
brca.graph.df.filt0.5$rank

library(ggrepel)
# Genes of interest = genes identified buy literature and my results to be important
genes.of.interest <- c( "PTEN", "TP53", "BRCA1", "BRCA2", "CDH1", "GATA3", "PIK3CA", "MUC16", "MAP3K1")
brca.graph.df.filt0.5$Hugo.Symbol -> brca.graph.df.filt0.5$label
brca.graph.df.filt0.5$label[!(brca.graph.df.filt0.5$label%in%genes.of.interest)] <- NA

ggplot(brca.graph.df.filt0.5, aes(x=rank, y = Percent.of.Patients.With.Mutation, label = label)) +
  geom_line() +
  geom_point(color = "red") +
  labs(x = "Rank") +
  labs(y = "Percent of Patients With Mutation") +
  theme_minimal() +
  geom_text_repel(nudge_x = 100, nudge_y = 2, max.overlaps = Inf) +
  ggtitle("Breast Invasive Carcinoma") +
  theme(plot.title = element_text(hjust = 0.5)) -> rank.brca0.5

rank.brca0.5
ggsave("~/tcga_biolinks1/Rank/rank.brca0.5.png", plot = rank.brca0.5, width = 10, height = 5)

#### PRAD ######################################################################


prad.graph <- readRDS(file = "~/tcga_biolinks1/RDS/prad.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/prad5.node.csv") -> prad.graph.df
prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> prad.graph.df.filt0.5

rank(-prad.graph.df.filt0.5$Percent.of.Patients.With.Mutation) -> prad.graph.df.filt0.5$rank
prad.graph.df.filt0.5$rank

library(ggrepel)
# Genes of interest = genes identified buy literature and my results to be important
genes.of.interest <- c( "SPOP", "TP53", "KMT2D", "MUC16", "EGFR", "MYC", "VEGFA", "PTEN")
prad.graph.df.filt0.5$Hugo.Symbol -> prad.graph.df.filt0.5$label
prad.graph.df.filt0.5$label[!(prad.graph.df.filt0.5$label%in%genes.of.interest)] <- NA

ggplot(prad.graph.df.filt0.5, aes(x=rank, y = Percent.of.Patients.With.Mutation, label = label)) +
  geom_line() +
  geom_point(color = "red") +
  labs(x = "Rank") +
  labs(y = "Percent of Patients With Mutation") +
  theme_minimal() +
  geom_text_repel(nudge_x = 50, nudge_y = 2, max.overlaps = Inf) +
  ggtitle("Prostate Adenocarcinoma") +
  theme(plot.title = element_text(hjust = 0.5)) -> rank.prad0.5

rank.prad0.5
ggsave("~/tcga_biolinks1/Rank/rank.prad0.5.png", plot = rank.prad0.5, width = 10, height = 5)

####  LUAD ######################################################################
luad.graph <- readRDS(file = "~/tcga_biolinks1/RDS/luad.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/luad4.node.csv") -> luad.graph.df

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 0.5) -> luad.graph.df.filt0.5

rank(-luad.graph.df.filt0.5$Percent.of.Patients.With.Mutation) -> luad.graph.df.filt0.5$rank
luad.graph.df.filt0.5$rank

library(ggrepel)
# Genes of interest = genes identified buy literature and my results to be important
genes.of.interest <- c( "EGFR", "ALK", "TP53", "KRAS", "MUC16", "CSMD3")
luad.graph.df.filt0.5$Hugo.Symbol -> luad.graph.df.filt0.5$label
luad.graph.df.filt0.5$label[!(luad.graph.df.filt0.5$label%in%genes.of.interest)] <- NA

ggplot(luad.graph.df.filt0.5, aes(x=rank, y = Percent.of.Patients.With.Mutation, label = label)) +
  geom_line() +
  geom_point(color = "red") +
  labs(x = "Rank") +
  labs(y = "Percent of Patients With Mutation") +
  theme_minimal() +
  geom_text_repel(nudge_x = 500, nudge_y = 2, max.overlaps = Inf) +
  ggtitle("Lung Adenocarcinoma") +
  theme(plot.title = element_text(hjust = 0.5)) -> rank.luad0.5
rank.luad0.5
ggsave("~/tcga_biolinks1/Rank/rank.luad0.5.png", plot = rank.luad0.5, width = 10, height = 5)

#### STAD ######################################################################
stad.graph <- readRDS(file = "~/tcga_biolinks1/RDS/stad.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/stad5.node.csv") -> stad.graph.df

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> stad.graph.df.filt0.5

rank(-stad.graph.df.filt0.5$Percent.of.Patients.With.Mutation) -> stad.graph.df.filt0.5$rank
stad.graph.df.filt0.5$rank

library(ggrepel)
# Genes of interest = genes identified buy literature and my results to be important
genes.of.interest <- c( "FN1", "TP53", "TIMP1", "SPP1", "APOE", "VCAN", "CDH1", "MUC16", "LRP1B", "ARID1A", "CSMD3")
stad.graph.df.filt0.5$Hugo.Symbol -> stad.graph.df.filt0.5$label
stad.graph.df.filt0.5$label[!(stad.graph.df.filt0.5$label%in%genes.of.interest)] <- NA

ggplot(stad.graph.df.filt0.5, aes(x=rank, y = Percent.of.Patients.With.Mutation, label = label)) +
  geom_line() +
  geom_point(color = "red") +
  labs(x = "Rank") +
  labs(y = "Percent of Patients With Mutation") +
  theme_minimal() +
  geom_text_repel(nudge_x = 500, nudge_y = 2, max.overlaps = Inf) +
  ggtitle("Stomach Adenocarcinoma") +
  theme(plot.title = element_text(hjust = 0.5)) -> rank.stad0.5

rank.stad0.5
ggsave("~/tcga_biolinks1/Rank/rank.stad0.5.png", plot = rank.stad0.5, width = 10, height = 5)

#### KIRP ######################################################################
kirp.graph <- readRDS(file = "~/tcga_biolinks1/RDS/kirp.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/kirp5.node.csv") -> kirp.graph.df

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutations >= 0.5) -> kirp.graph.df.filt0.5


rank(-kirp.graph.df.filt0.5$Percent.of.Patients.With.Mutation) -> kirp.graph.df.filt0.5$rank
kirp.graph.df.filt0.5$rank

library(ggrepel)
# Genes of interest = genes identified buy literature and my results to be important
genes.of.interest <- c( "VHL", "MET", "MUC16", "KMT2C", "KMT2D", "FAT1", "SETD2", "BAP1", "PBRM1", "TP53", "ARID1A", "PKHD1", "FH")
kirp.graph.df.filt0.5$Hugo.Symbol -> kirp.graph.df.filt0.5$label
kirp.graph.df.filt0.5$label[!(kirp.graph.df.filt0.5$label%in%genes.of.interest)] <- NA

ggplot(kirp.graph.df.filt0.5, aes(x=rank, y = Percent.of.Patients.With.Mutation, label = label)) +
  geom_line() +
  geom_point(color = "red") +
  labs(x = "Rank") +
  labs(y = "Percent of Patients With Mutation") +
  theme_minimal() +
  geom_text_repel(nudge_x = 50, nudge_y = 0.5, max.overlaps = Inf) +
  ggtitle("Kidney Renal Papillary Cell Carcinoma") +
  theme(plot.title = element_text(hjust = 0.5)) -> rank.kirp0.5

rank.kirp0.5
ggsave("~/tcga_biolinks1/Rank/rank.kirp0.5.png", plot = rank.kirp0.5, width = 10, height = 5)

#### KIRC ######################################################################
kirc.graph <- readRDS(file = "~/tcga_biolinks1/RDS/kirc.graph")
read.csv(file = "~/tcga_biolinks1/GRAPHML/kirc4.node.csv") -> kirc.graph.df

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> kirc.graph.df.filt0.5


rank(-kirc.graph.df.filt0.5$Percent.of.Patients.With.Mutation) -> kirc.graph.df.filt0.5$rank
kirc.graph.df.filt0.5$rank

library(ggrepel)
# Genes of interest = genes identified buy literature and my results to be important
genes.of.interest <- c( "VHL", "MET", "MUC16", "KMT2C", "KMT2D", "FAT1", "SETD2", "BAP1", "PBRM1", "TP53", "ARID1A", "PKHD1", "FH")
kirc.graph.df.filt0.5$Hugo.Symbol -> kirc.graph.df.filt0.5$label
kirc.graph.df.filt0.5$label[!(kirc.graph.df.filt0.5$label%in%genes.of.interest)] <- NA

ggplot(kirc.graph.df.filt0.5, aes(x=rank, y = Percent.of.Patients.With.Mutation, label = label, colour=)) +
  geom_line() +
  geom_point(color = "red") +
  labs(x = "Rank") +
  labs(y = "Percent of Patients With Mutation") +
  theme_minimal() +
  geom_text_repel(nudge_x = 150, nudge_y = 3, max.overlaps = Inf) +
  ggtitle("Kidney Renal Clear Cell Carcinoma") +
  theme(plot.title = element_text(hjust = 0.5)) -> rank.kirc0.5


rank.kirc0.5
ggsave("~/tcga_biolinks1/Rank/rank.kirc0.5.png", plot = rank.kirc0.5, width = 10, height = 5)

################################################################################
rank.brca0.5 + theme(text = element_text(size = 10)) -> p1
rank.kirc0.5 + theme(text = element_text(size = 10)) -> p2
rank.kirp0.5 + theme(text = element_text(size = 10)) -> p3
rank.luad0.5 + theme(text = element_text(size = 10)) -> p4
rank.prad0.5 + theme(text = element_text(size = 10)) -> p5
rank.stad0.5 + theme(text = element_text(size = 10)) -> p6
rank.figure <- ggarrange(p1, p2, p3, p4, p5, p6,
                         labels = c("A", "B", "C", "D", "E", "F"),
                         common.legend = TRUE, legend = "bottom",
                         ncol = 3, nrow = 1)
annotate_figure(rank.figure, top = text_grob("Genes Ranked by Percent of Patients With Mutated Gene", size = 20, face = "bold")) -> rank.plot
rank.plot
ggsave("~/tcga_biolinks1/Rank/rank.plot.png", plot = brca.plot, width = 28, height = 10)


