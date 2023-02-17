library(dplyr)

head(prad.data$mutational)

colnames(prad.data$mutational)

prad.data$mutational %>%
  count(Gene) -> no.mutations

prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
count(Gene) -> no.missense

prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  filter(grepl("deleterious", SIFT))%>%
  count(Gene) -> no.sift.del

prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  filter(grepl("possibly", PolyPhen))%>%
  count(Gene) -> no.possible.polyphen

prad.data$mutational %>%
  #group_by(Gene) %>%
  filter(Variant_Classification=="Missense_Mutation") %>%
  filter(grepl("probably", PolyPhen))%>%
  count(Gene) -> no.probable.polyphen

prad.data$mutational %>% dplyr::count(SIFT)
gsub("deleterious\\(",'',prad.data$mutational$SIFT) -> new.variable
new.variable
gsub("\\)","",new.variable) -> new.variable2
new.variable2
gsub("tolerated\\(",'',new.variable2) -> new.variable3
new.variable3
gsub("deleterious_low_confidence\\(",'',new.variable3) -> new.variable4
new.variable4
gsub("tolerated_low_confidence\\(",'',new.variable4) -> SIFT.number
SIFT.number

data.frame(no.mutations,
           no.missense,
           no.possible.polyphen,
           no.probable.polyphen,
           SIFT.number)

length(no.mutations) <- length(no.missense)
length(no.mutations) <- length(no.possible.polyphen)
length(no.mutations) <- length(no.probable.polyphen)
