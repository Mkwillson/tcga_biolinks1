
##########################################


#### this is the part where you filter your wordart to contain the bits you want to show
luad.graph.df %>%
  filter(Number.of.Patients.With.Mutation >= 1) %>%
  filter(Hugo.Symbol != "") -> luad.graph.df.filt1


luad.graph.df.filt1$Number.of.Patients.With.Mutation


# library(foreach)

# unique(luad.graph.df.filt1$Hugo.Symbol) -> mut.levels
# 
# wordcloud.out <- foreach(i = 1:length(mut.levels), .combine = c)%do%{
# name.of.mut <- mut.levels[i] 
# luad.graph.df.filt1[which(luad.graph.df.filt1$Hugo.Symbol == mut.levels[i]),"Number.of.Patients.With.Mutation"] -> freq.of.mut
# return(rep(name.of.mut, freq.of.mut))
# }
# 
# write.csv(wordcloud.out, file = "./wordcloud.csv", quote = F, row.names = F, col.names = F)
#### BRCA ####
### Create data.frame for brca 0.5%
data.frame(brca.graph.df.filt0.5$Hugo.Symbol, brca.graph.df.filt0.5$Number.of.Patients.With.Mutation, brca.graph.df.filt0.5$PolyPhen.Mean) -> wordcloud.df.brca.test
### output as a .csv
write.csv(wordcloud.df.brca.test, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.brca.test.csv", quote = F, row.names = F)

### Create data.frame for brca 1.5%
data.frame(brca.graph.df.filt1.5$Hugo.Symbol, brca.graph.df.filt1.5$Number.of.Patients.With.Mutation) -> wordcloud.df.brca.1.5
### output as a .csv
write.csv(wordcloud.df.brca.1.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.brca.1.5.csv", quote = F, row.names = F)

### Create data.frame for brca Probably damaging
data.frame(brca.graph.df.filt.prob$Hugo.Symbol, brca.graph.df.filt.prob$Number.of.Patients.With.Mutation) -> wordcloud.df.brca.prob
### output as a .csv
write.csv(wordcloud.df.brca.prob, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.brca.prob.csv", quote = F, row.names = F)

### Create data.frame for  brca all damaging
data.frame( brca.graph.df.filt.dam$Hugo.Symbol,  brca.graph.df.filt.dam$Number.of.Patients.With.Mutation) -> wordcloud.df.brca.dam
### output as a .csv
write.csv(wordcloud.df.brca.dam, file = "~/tcga_biolinks1/wordcloud/wordcloud.df. brca.dam.csv", quote = F, row.names = F)

### PRAD ####
### Create data.frame for prad 0.5%
data.frame(prad.graph.df.filt0.5$Hugo.Symbol, prad.graph.df.filt0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.0.5
### output as a .csv
write.csv(wordcloud.df.prad.0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.prad.0.5.csv", quote = F, row.names = F)

### Create data.frame for prad 1.5%
data.frame(prad.graph.df.filt1.5$Hugo.Symbol, prad.graph.df.filt1.5$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.1.5
### output as a .csv
write.csv(wordcloud.df.prad.1.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.prad.1.5.csv", quote = F, row.names = F)

### Create data.frame for prad Probably damaging
data.frame(prad.graph.df.filt.prob$Hugo.Symbol, prad.graph.df.filt.prob$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.prob
### output as a .csv
write.csv(wordcloud.df.prad.prob, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.prad.prob.csv", quote = F, row.names = F)

### Create data.frame for   prad all damaging
data.frame(  prad.graph.df.filt.dam$Hugo.Symbol,prad.graph.df.filt.dam$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.dam
### output as a .csv
write.csv(wordcloud.df.prad.dam, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.  prad.dam.csv", quote = F, row.names = F)


### LUAD ####
### Create data.frame for luad 0.5%
data.frame(luad.graph.df.filt0.5$Hugo.Symbol, luad.graph.df.filt0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.luad.0.5
### output as a .csv
write.csv(wordcloud.df.luad.0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.luad.0.5.csv", quote = F, row.names = F)

### Create data.frame for luad 1.5%
data.frame(luad.graph.df.filt1.5$Hugo.Symbol, luad.graph.df.filt1.5$Number.of.Patients.With.Mutation, luad) -> wordcloud.df.luad.1.5
### output as a .csv
write.csv(wordcloud.df.luad.1.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.luad.1.5.csv", quote = F, row.names = F)

### Create data.frame for luad Probably damaging
data.frame(luad.graph.df.filt.prob$Hugo.Symbol, luad.graph.df.filt.prob$Number.of.Patients.With.Mutation) -> wordcloud.df.luad.prob
### output as a .csv
write.csv(wordcloud.df.luad.prob, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.luad.prob.csv", quote = F, row.names = F)

### Create data.frame for  luad all damaging
data.frame( luad.graph.df.filt.dam$Hugo.Symbol,  luad.graph.df.filt.dam$Number.of.Patients.With.Mutation) -> wordcloud.df.luad.dam
### output as a .csv
write.csv(wordcloud.df.luad.dam, file = "~/tcga_biolinks1/wordcloud/wordcloud.df. luad.dam.csv", quote = F, row.names = F)

####kirp####
### Create data.frame for  luad 0.5%
data.frame( kirp.graph.df.filt0.5$Hugo.Symbol,  kirp.graph.df.filt0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.0.5
### output as a .csv
write.csv(wordcloud.df.kirp.0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.0.5.csv", quote = F, row.names = F)

### Create data.frame for kirp 1.5%
data.frame(kirp.graph.df.filt1.5$Hugo.Symbol, kirp.graph.df.filt1.5$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.1.5
### output as a .csv
write.csv(wordcloud.df.kirp.1.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.1.5.csv", quote = F, row.names = F)

### Create data.frame for kirp Probably damaging
data.frame(kirp.graph.df.filt.prob$Hugo.Symbol, kirp.graph.df.filt.prob$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.prob
### output as a .csv
write.csv(wordcloud.df.kirp.prob, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.prob.csv", quote = F, row.names = F)

### Create data.frame for kirp all damaging
data.frame(kirp.graph.df.filt.dam$Hugo.Symbol, kirp.graph.df.filt.dam$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.dam
### output as a .csv
write.csv(wordcloud.df.kirp.dam, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.dam.csv", quote = F, row.names = F)


############

rank(-luad.graph.df.filt1$Number.of.Patients.With.Mutation) -> luad.graph.df.filt1$rank
# luad.graph.df.filt1$Number.of.Patients.With.Mutation/total.no.of.patients -> luad.graph.df.filt1$perc


library(ggrepel)

genes.of.interest <- c("TP53", "KRAS")
luad.graph.df.filt1$Hugo.Symbol -> luad.graph.df.filt1$label
luad.graph.df.filt1$label[!(luad.graph.df.filt1$label%in%genes.of.interest)] <- NA


ggplot(luad.graph.df.filt1, aes(x=rank, y = Number.of.Patients.With.Mutation, label = label)) +
  geom_line() +
  geom_point(color = "red") +
  # geom_vline(data=dat.luad2, aes(xintercept=random.cutoff2.luad,  colour=Distribution),
  #            linetype="dashed", linewidth=1) +
  annotate("text", x = 2000, y = 200, label = "Lung Cancer spectrum") +
  #annotate("text", x = 4.5, y = 1.65, label = text.annot2) +
  # annotate("text", x = 4.45, y = 1.55, label = paste0("N = ",text.annot3)) +
  labs(x = "Rank") +
  labs(y = "Number of Patients With Mutation") +
  theme_minimal() +
  geom_text_repel(arrow = T, nudge_x = 2) + 
  ggtitle("Mean Gene Distance Distributions of Both Random and LUAD Genes on the LUAD Network") -> gg.spectrum


###################

#### filtered graph as a starting point...
luad.graph.df.filt2
read.csv(file = "~/tcga_biolinks1/GRAPHML/luad4.stats.node.csv") -> luad.graph.df
hist(luad.graph.df$BetweennessCentrality)
hist(luad.graph.df.filt2$BetweennessCentrality)
luad.graph.df %>%
  filter(Frequency.Gene.is.Mutated >= 1) -> luad.graph.df.filt2

luad.graph.df$PolyPhen.Probably.Damaging

sample(attr(V(luad.graph),"name"), length(unique.genes.luad2)) -> all.genes.rand.luad2


#### there are different measures of centrality that you could test if different between categories

wilcox.test(Stress~PolyPhen.Probably.Damaging, data = luad.graph.df)
boxplot(Stress~PolyPhen.Probably.Damaging, data = luad.graph.df )
tapply(luad.graph.df$Stress, luad.graph.df$PolyPhen.Probably.Damaging, median)


wilcox.test(ClosenessCentrality~PolyPhen.Probably.Damaging, data = luad.graph.df)
boxplot(ClosenessCentrality~PolyPhen.Probably.Damaging, data = luad.graph.df )
tapply(luad.graph.df$ClosenessCentrality, luad.graph.df$PolyPhen.Probably.Damaging, median)


plot(luad.graph.df$SIFT.Mean, luad.graph.df$Degree)
plot(luad.graph.df$SIFT.Mean, log(luad.graph.df$BetweennessCentrality))

cor.test(luad.graph.df$SIFT.Mean, log(luad.graph.df$BetweennessCentrality), method = "spearman")

lm(luad.graph.df$Stress~luad.graph.df$SIFT.Mean, na.rm = T)
plot(luad.graph.df$SIFT.Mean, luad.graph.df$Stress)


####################

quantile(colMeans(out.rand.luad.graph2, na.rm = T),0.01) -> random.cutoff2.luad

### make sure you set the out.luad2 gene names (see above)
colMeans(out.luad2) -> col.mean.out.luad2

### output this object as .csv to export back in to cytoscape 
write.csv(col.mean.out.luad2, file = "./col.mean.out.luad2.csv")

#### these are the genes relevant at the 1% level
col.mean.out.luad2[col.mean.out.luad2 < random.cutoff2.luad]




