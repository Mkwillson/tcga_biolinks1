
########## stad ##############
# MI threshold
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP/kirp.filt.exp.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP/kirp.filt.exp.tfs.txt \
--pvalue 1E-8 \
--seed 1 \
--threads 20 \
--calculateThreshold

# 1 bootstrap 
time java -Xmx13G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD/luad.filt.exp.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD/luad.filt.exp.tfs.txt \
--pvalue 1E-8 \
--seed 1 \
--threads 20

# 100 bootstrap threshold
for i in {1..100}
do 
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP/kirp.filt.exp.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP/kirp.filt.exp.tfs.txt \
--pvalue 1E-8 \
--seed $i \
--threads 20 \
done

# 100 bootstrap
for i in {1..100}
do 
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD/luad.filt.exp.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD --tfs  /home/ug_megan/tcga_biolinks1/ARACNe.data/LUAD/luad.filt.exp.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate bootstraps
java -Xmx5G -jar  /home/anafox/ARACNe-AP/dist/aracne.jar \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRP \
--consolidate

############################################################

for i in {1..100}
do 
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD/prad.filt.exp.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD --tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/PRAD/prad.filt.exp.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done
