##### BRCA histological type
# an object you want to summarise
colData(brca.data$expression) -> object_of_int

table(object_of_int$gender)
table(object_of_int$gender, object_of_int$primary_diagnosis)

brca.data$expression$sample

# create an index
match(brca.data$clinical$bcr_patient_barcode, brca.data$expression$patient) -> idx

brca.data$clinical[idx[!is.na(idx)],] -> matched.exp.clinical
dim(matched.exp.clinical)

table(matched.exp.clinical$histological_type)

#### PRAD histological type ################################

colData(prad.data$expression) -> object_of_int2

table(object_of_int2$gender)
table(object_of_int2$gender, object_of_int2$primary_diagnosis)

prad.data$expression$sample

# create an index
match(prad.data$clinical$bcr_patient_barcode, prad.data$expression$patient) -> idx

brca.data$clinical[idx[!is.na(idx)],] -> matched.prad.exp.clinical
dim(matched.prad.exp.clinical)

table(matched.prad.exp.clinical$histological_type)

#### LUAD ##################################################
# An object you want to summarise
colData(luad.data$expression) -> object_of_int3
object_of_int3

table(object_of_int3$gender)
table(object_of_int3$gender, object_of_int3$primary_diagnosis)

luad.data$expression$sample

# create an index
match(luad.data$clinical$bcr_patient_barcode, luad.data$expression$patient) -> idx3
idx3

luad.data$clinical[idx3[!is.na(idx3)],] -> matched.luad.exp.clinical
matched.luad.exp.clinical

dim(matched.luad.exp.clinical)

table(matched.luad.exp.clinical$histological_type)

#### STAD ##############################################
# An object you want to summarise
colData(stad.data$expression) -> object_of_int4
object_of_int4

table(object_of_int4$gender)
table(object_of_int4$gender, object_of_int4$primary_diagnosis)

stad.data$expression$sample

# create an index
match(stad.data$clinical$bcr_patient_barcode, stad.data$expression$patient) -> idx4
idx4

stad.data$clinical[idx4[!is.na(idx4)],] -> matched.stad.exp.clinical
matched.stad.exp.clinical

dim(matched.stad.exp.clinical)

table(matched.stad.exp.clinical$histological_type)

#### KIRC #############################################

# An object you want to summarise
colData(kirc.data$expression) -> object_of_int5
object_of_int5

table(object_of_int5$gender)
table(object_of_int5$gender, object_of_int5$primary_diagnosis)

kirc.data$expression$sample

# create an index
match(kirc.data$clinical$bcr_patient_barcode, kirc.data$expression$patient) -> idx5
idx5

kirc.data$clinical[idx5[!is.na(idx5)],] -> matched.kirc.exp.clinical
matched.kirc.exp.clinical

dim(matched.kirc.exp.clinical)

table(matched.kirc.exp.clinical$histological_type)

#### KIRP ############################################
# An object you want to summarise
colData(kirp.data$expression) -> object_of_int6
object_of_int6

table(object_of_int6$gender)
table(object_of_int6$gender, object_of_int6$primary_diagnosis)

kirp.data$expression$sample

# create an index
match(kirp.data$clinical$bcr_patient_barcode, kirp.data$expression$patient) -> idx6
idx6

kirp.data$clinical[idx6[!is.na(idx6)],] -> matched.kirp.exp.clinical
matched.kirp.exp.clinical

dim(matched.kirp.exp.clinical)

table(matched.kirp.exp.clinical$histological_type)

