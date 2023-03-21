
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

read.csv(file = "~/tcga_biolinks1/GRAPHML/brca5.node.csv") -> brca.graph.df

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1.5) -> brca.graph.df.filt1.5

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> brca.graph.df.filt0.5

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> brca.graph.df.filt.prob

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> brca.graph.df.filt.dam

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(PolyPhen.Mean >0.466)-> brca.graph.df.filt.dam0.5

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(SIFT.Mean <0.05)-> brca.graph.df.filt.del

brca.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(SIFT.Mean <0.05)-> brca.graph.df.filt.del0.5



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

### Create data frame for brca all damaging and >= 0.5%
data.frame( brca.graph.df.filt.dam0.5$Hugo.Symbol,  brca.graph.df.filt.dam0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.brca.dam0.5
### output as a .csv
write.csv(wordcloud.df.brca.dam0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.brca.dam0.5.csv", quote = F, row.names = F)

### Create data frame for brca deleterious
data.frame( brca.graph.df.filt.del$Hugo.Symbol,  brca.graph.df.filt.del$Number.of.Patients.With.Mutation) -> wordcloud.df.brca.del
### output as a .csv
write.csv(wordcloud.df.brca.del, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.brca.del.csv", quote = F, row.names = F)

### Create data frame for brca deleterious and >= 0.5%
data.frame( brca.graph.df.filt.del0.5$Hugo.Symbol,  brca.graph.df.filt.del0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.brca.del0.5
### output as a .csv
write.csv(wordcloud.df.brca.del0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.brca.del0.5.csv", quote = F, row.names = F)



#### PRAD ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/prad5.node.csv") -> prad.graph.df

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1.5) -> prad.graph.df.filt1.5

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> prad.graph.df.filt0.5

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> prad.graph.df.filt.prob

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> prad.graph.df.filt.dam

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(PolyPhen.Mean >0.466)-> prad.graph.df.filt.dam0.5

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(SIFT.Mean <0.05)-> prad.graph.df.filt.del

prad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(SIFT.Mean <0.05)-> prad.graph.df.filt.del0.5



### Create data.frame for prad 0.5%
data.frame(prad.graph.df.filt0.5$Hugo.Symbol, prad.graph.df.filt0.5$Number.of.Patients.With.Mutation, prad.graph.df.filt0.5$PolyPhen.Mean) -> wordcloud.df.prad.test
### output as a .csv
write.csv(wordcloud.df.prad.test, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.prad.test.csv", quote = F, row.names = F)

### Create data.frame for prad 1.5%
data.frame(prad.graph.df.filt1.5$Hugo.Symbol, prad.graph.df.filt1.5$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.1.5
### output as a .csv
write.csv(wordcloud.df.prad.1.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.prad.1.5.csv", quote = F, row.names = F)

### Create data.frame for prad Probably damaging
data.frame(prad.graph.df.filt.prob$Hugo.Symbol, prad.graph.df.filt.prob$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.prob
### output as a .csv
write.csv(wordcloud.df.prad.prob, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.prad.prob.csv", quote = F, row.names = F)

### Create data.frame for  prad all damaging
data.frame( prad.graph.df.filt.dam$Hugo.Symbol,  prad.graph.df.filt.dam$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.dam
### output as a .csv
write.csv(wordcloud.df.prad.dam, file = "~/tcga_biolinks1/wordcloud/wordcloud.df. prad.dam.csv", quote = F, row.names = F)

### Create data frame for prad all damaging and >= 0.5%
data.frame( prad.graph.df.filt.dam0.5$Hugo.Symbol,  prad.graph.df.filt.dam0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.dam0.5
### output as a .csv
write.csv(wordcloud.df.prad.dam0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.prad.dam0.5.csv", quote = F, row.names = F)

### Create data frame for prad deleterious
data.frame( prad.graph.df.filt.del$Hugo.Symbol,  prad.graph.df.filt.del$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.del
### output as a .csv
write.csv(wordcloud.df.prad.del, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.prad.del.csv", quote = F, row.names = F)

### Create data frame for prad deleterious and >= 0.5%
data.frame( prad.graph.df.filt.del0.5$Hugo.Symbol,  prad.graph.df.filt.del0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.prad.del0.5
### output as a .csv
write.csv(wordcloud.df.prad.del0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.prad.del0.5.csv", quote = F, row.names = F)


#### LUAD ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/luad4.node.csv") -> luad.graph.df

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 1.5) -> luad.graph.df.filt1.5

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 0.5) -> luad.graph.df.filt0.5

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> luad.graph.df.filt.prob

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> luad.graph.df.filt.dam

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 0.5) %>%
  filter(PolyPhen.Mean >0.466)-> luad.graph.df.filt.dam0.5

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(SIFT.Mean <0.05)-> luad.graph.df.filt.del

luad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.With.a.Missense.Mutation >= 0.5) %>%
  filter(SIFT.Mean <0.05)-> luad.graph.df.filt.del0.5



### Create data.frame for luad 0.5%
data.frame(luad.graph.df.filt0.5$Hugo.Symbol, luad.graph.df.filt0.5$Number.of.Patients.With.Mutation, luad.graph.df.filt0.5$PolyPhen.Mean) -> wordcloud.df.luad.test
### output as a .csv
write.csv(wordcloud.df.luad.test, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.luad.test.csv", quote = F, row.names = F)

### Create data.frame for luad 1.5%
data.frame(luad.graph.df.filt1.5$Hugo.Symbol, luad.graph.df.filt1.5$Number.of.Patients.With.Mutation) -> wordcloud.df.luad.1.5
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

### Create data frame for luad all damaging and >= 0.5%
data.frame( luad.graph.df.filt.dam0.5$Hugo.Symbol,  luad.graph.df.filt.dam0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.luad.dam0.5
### output as a .csv
write.csv(wordcloud.df.luad.dam0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.luad.dam0.5.csv", quote = F, row.names = F)

### Create data frame for luad deleterious
data.frame( luad.graph.df.filt.del$Hugo.Symbol,  luad.graph.df.filt.del$Number.of.Patients.With.Mutation) -> wordcloud.df.luad.del
### output as a .csv
write.csv(wordcloud.df.luad.del, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.luad.del.csv", quote = F, row.names = F)

### Create data frame for luad deleterious and >= 0.5%
data.frame( luad.graph.df.filt.del0.5$Hugo.Symbol,  luad.graph.df.filt.del0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.luad.del0.5
### output as a .csv
write.csv(wordcloud.df.luad.del0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.luad.del0.5.csv", quote = F, row.names = F)


#### STAD ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/stad5.node.csv") -> stad.graph.df

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1.5) -> stad.graph.df.filt1.5

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> stad.graph.df.filt0.5

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> stad.graph.df.filt.prob

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> stad.graph.df.filt.dam

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(PolyPhen.Mean >0.466)-> stad.graph.df.filt.dam0.5

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(SIFT.Mean <0.05)-> stad.graph.df.filt.del

stad.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(SIFT.Mean <0.05)-> stad.graph.df.filt.del0.5



### Create data.frame for stad 0.5%
data.frame(stad.graph.df.filt0.5$Hugo.Symbol, stad.graph.df.filt0.5$Number.of.Patients.With.Mutation, stad.graph.df.filt0.5$PolyPhen.Mean) -> wordcloud.df.stad.test
### output as a .csv
write.csv(wordcloud.df.stad.test, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.stad.test.csv", quote = F, row.names = F)

### Create data.frame for stad 1.5%
data.frame(stad.graph.df.filt1.5$Hugo.Symbol, stad.graph.df.filt1.5$Number.of.Patients.With.Mutation) -> wordcloud.df.stad.1.5
### output as a .csv
write.csv(wordcloud.df.stad.1.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.stad.1.5.csv", quote = F, row.names = F)

### Create data.frame for stad Probably damaging
data.frame(stad.graph.df.filt.prob$Hugo.Symbol, stad.graph.df.filt.prob$Number.of.Patients.With.Mutation) -> wordcloud.df.stad.prob
### output as a .csv
write.csv(wordcloud.df.stad.prob, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.stad.prob.csv", quote = F, row.names = F)

### Create data.frame for  stad all damaging
data.frame( stad.graph.df.filt.dam$Hugo.Symbol,  stad.graph.df.filt.dam$Number.of.Patients.With.Mutation) -> wordcloud.df.stad.dam
### output as a .csv
write.csv(wordcloud.df.stad.dam, file = "~/tcga_biolinks1/wordcloud/wordcloud.df. stad.dam.csv", quote = F, row.names = F)

### Create data frame for stad all damaging and >= 0.5%
data.frame( stad.graph.df.filt.dam0.5$Hugo.Symbol,  stad.graph.df.filt.dam0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.stad.dam0.5
### output as a .csv
write.csv(wordcloud.df.stad.dam0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.stad.dam0.5.csv", quote = F, row.names = F)

### Create data frame for stad deleterious
data.frame( stad.graph.df.filt.del$Hugo.Symbol,  stad.graph.df.filt.del$Number.of.Patients.With.Mutation) -> wordcloud.df.stad.del
### output as a .csv
write.csv(wordcloud.df.stad.del, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.stad.del.csv", quote = F, row.names = F)

### Create data frame for stad deleterious and >= 0.5%
data.frame( stad.graph.df.filt.del0.5$Hugo.Symbol,  stad.graph.df.filt.del0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.stad.del0.5
### output as a .csv
write.csv(wordcloud.df.stad.del0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.stad.del0.5.csv", quote = F, row.names = F)




#### KIRP ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/kirp5.node.csv") -> kirp.graph.df

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutations >= 1.5) -> kirp.graph.df.filt1.5

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutations >= 0.5) -> kirp.graph.df.filt0.5

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> kirp.graph.df.filt.prob

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> kirp.graph.df.filt.dam

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutations >= 0.5) %>%
  filter(PolyPhen.Mean >0.466)-> kirp.graph.df.filt.dam0.5

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(SIFT.Mean <0.05)-> kirp.graph.df.filt.del

kirp.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutations >= 0.5) %>%
  filter(SIFT.Mean <0.05)-> kirp.graph.df.filt.del0.5



### Create data.frame for kirp 0.5%
data.frame(kirp.graph.df.filt0.5$Hugo.Symbol, kirp.graph.df.filt0.5$Number.of.Patients.With.Mutation, kirp.graph.df.filt0.5$PolyPhen.Mean) -> wordcloud.df.kirp.test
### output as a .csv
write.csv(wordcloud.df.kirp.test, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.test.csv", quote = F, row.names = F)

### Create data.frame for kirp 1.5%
data.frame(kirp.graph.df.filt1.5$Hugo.Symbol, kirp.graph.df.filt1.5$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.1.5
### output as a .csv
write.csv(wordcloud.df.kirp.1.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.1.5.csv", quote = F, row.names = F)

### Create data.frame for kirp Probably damaging
data.frame(kirp.graph.df.filt.prob$Hugo.Symbol, kirp.graph.df.filt.prob$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.prob
### output as a .csv
write.csv(wordcloud.df.kirp.prob, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.prob.csv", quote = F, row.names = F)

### Create data.frame for  kirp all damaging
data.frame( kirp.graph.df.filt.dam$Hugo.Symbol,  kirp.graph.df.filt.dam$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.dam
### output as a .csv
write.csv(wordcloud.df.kirp.dam, file = "~/tcga_biolinks1/wordcloud/wordcloud.df. kirp.dam.csv", quote = F, row.names = F)

### Create data frame for kirp all damaging and >= 0.5%
data.frame( kirp.graph.df.filt.dam0.5$Hugo.Symbol,  kirp.graph.df.filt.dam0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.dam0.5
### output as a .csv
write.csv(wordcloud.df.kirp.dam0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.dam0.5.csv", quote = F, row.names = F)

### Create data frame for kirp deleterious
data.frame( kirp.graph.df.filt.del$Hugo.Symbol,  kirp.graph.df.filt.del$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.del
### output as a .csv
write.csv(wordcloud.df.kirp.del, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.del.csv", quote = F, row.names = F)

### Create data frame for kirp deleterious and >= 0.5%
data.frame( kirp.graph.df.filt.del0.5$Hugo.Symbol,  kirp.graph.df.filt.del0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.kirp.del0.5
### output as a .csv
write.csv(wordcloud.df.kirp.del0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirp.del0.5.csv", quote = F, row.names = F)

#### KIRC ####

read.csv(file = "~/tcga_biolinks1/GRAPHML/kirc4.node.csv") -> kirc.graph.df

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 1.5) -> kirc.graph.df.filt1.5

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) -> kirc.graph.df.filt0.5

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.908) -> kirc.graph.df.filt.prob

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(PolyPhen.Mean >0.466)-> kirc.graph.df.filt.dam

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(PolyPhen.Mean >0.466)-> kirc.graph.df.filt.dam0.5

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(SIFT.Mean <0.05)-> kirc.graph.df.filt.del

kirc.graph.df %>%
  filter(Number.of.Missense.Mutations >=1) %>%
  filter(Percent.of.Patients.with.a.Missense.Mutation >= 0.5) %>%
  filter(SIFT.Mean <0.05)-> kirc.graph.df.filt.del0.5



### Create data.frame for kirc 0.5%
data.frame(kirc.graph.df.filt0.5$Hugo.Symbol, kirc.graph.df.filt0.5$Number.of.Patients.With.Mutation, kirc.graph.df.filt0.5$PolyPhen.Mean) -> wordcloud.df.kirc.test
### output as a .csv
write.csv(wordcloud.df.kirc.test, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirc.test.csv", quote = F, row.names = F)

### Create data.frame for kirc 1.5%
data.frame(kirc.graph.df.filt1.5$Hugo.Symbol, kirc.graph.df.filt1.5$Number.of.Patients.With.Mutation) -> wordcloud.df.kirc.1.5
### output as a .csv
write.csv(wordcloud.df.kirc.1.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirc.1.5.csv", quote = F, row.names = F)

### Create data.frame for kirc Probably damaging
data.frame(kirc.graph.df.filt.prob$Hugo.Symbol, kirc.graph.df.filt.prob$Number.of.Patients.With.Mutation) -> wordcloud.df.kirc.prob
### output as a .csv
write.csv(wordcloud.df.kirc.prob, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirc.prob.csv", quote = F, row.names = F)

### Create data.frame for  kirc all damaging
data.frame( kirc.graph.df.filt.dam$Hugo.Symbol,  kirc.graph.df.filt.dam$Number.of.Patients.With.Mutation) -> wordcloud.df.kirc.dam
### output as a .csv
write.csv(wordcloud.df.kirc.dam, file = "~/tcga_biolinks1/wordcloud/wordcloud.df. kirc.dam.csv", quote = F, row.names = F)

### Create data frame for kirc all damaging and >= 0.5%
data.frame( kirc.graph.df.filt.dam0.5$Hugo.Symbol,  kirc.graph.df.filt.dam0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.kirc.dam0.5
### output as a .csv
write.csv(wordcloud.df.kirc.dam0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirc.dam0.5.csv", quote = F, row.names = F)

### Create data frame for kirc deleterious
data.frame( kirc.graph.df.filt.del$Hugo.Symbol,  kirc.graph.df.filt.del$Number.of.Patients.With.Mutation) -> wordcloud.df.kirc.del
### output as a .csv
write.csv(wordcloud.df.kirc.del, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirc.del.csv", quote = F, row.names = F)

### Create data frame for kirc deleterious and >= 0.5%
data.frame( kirc.graph.df.filt.del0.5$Hugo.Symbol,  kirc.graph.df.filt.del0.5$Number.of.Patients.With.Mutation) -> wordcloud.df.kirc.del0.5
### output as a .csv
write.csv(wordcloud.df.kirc.del0.5, file = "~/tcga_biolinks1/wordcloud/wordcloud.df.kirc.del0.5.csv", quote = F, row.names = F)

#######################################

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


####################

## Writing CSV of only genes less than random cut off (%>=0.5, missense)
## BRCA
readRDS(file = "~/tcga_biolinks1/stats/out.rand.brca.graph0.5") -> out.rand.brca.graph0.5
readRDS(file = "~/tcga_biolinks1/stats/out.brca0.5") -> out.brca0.5

### Define random cutoff
quantile(colMeans(out.rand.brca.graph0.5, na.rm = T),0.01) -> random.cutoff0.5.brca

### Set the gene names

colMeans(out.brca0.5) -> col.mean.out.brca0.5

### output this object as .csv to export back in to cytoscape 
write.csv(col.mean.out.brca0.5, file = "~/tcga_biolinks1/Cytoscape/col.mean.out.brca0.5.csv")

#### these are the genes relevant at the 1% level
col.mean.out.brca0.5[col.mean.out.brca0.5 < random.cutoff0.5.brca] -> genes.cutoff.brca
write.csv(genes.cutoff.brca, file = "~/tcga_biolinks1/Cytoscape/genes.cutoff.brca.csv")

#### PRAD ####
readRDS(file = "~/tcga_biolinks1/stats/out.rand.prad.graph0.5") -> out.rand.prad.graph0.5
readRDS(file = "~/tcga_biolinks1/stats/out.prad0.5") -> out.prad0.5

### Define random cutoff
quantile(colMeans(out.rand.prad.graph0.5, na.rm = T),0.01) -> random.cutoff0.5.prad

### Set the gene names

colMeans(out.prad0.5) -> col.mean.out.prad0.5

### output this object as .csv to export back in to cytoscape 
write.csv(col.mean.out.prad0.5, file = "~/tcga_biolinks1/Cytoscape/col.mean.out.prad0.5.csv")

#### these are the genes relevant at the 1% level
col.mean.out.prad0.5[col.mean.out.prad0.5 < random.cutoff0.5.prad] -> genes.cutoff.prad
write.csv(genes.cutoff.prad, file = "~/tcga_biolinks1/Cytoscape/genes.cutoff.prad.csv")



