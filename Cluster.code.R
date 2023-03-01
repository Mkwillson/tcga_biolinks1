# Calculate threshold for MI
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5/stad.filt.exp5.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5 \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5/stad.filt.exp5.tfs.txt \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20

# 100 bootstrap
for i in {1..100}
do 
time java -Xmx25G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5/stad.filt.exp5.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5 --tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5/stad.filt.exp5.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate
java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -o /home/ug_megan/tcga_biolinks1/ARACNe.data/STAD5 --consolidate

# Calculate threshold for MI
time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4/kirc.filt.exp4.txt \
-o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4 \
--tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4/kirc.filt.exp4.tfs.txt \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20

# 100 bootstrap
for i in {1..100}
do 
time java -Xmx25G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -e /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4/kirc.filt.exp4.txt -o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4 --tfs /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4/kirc.filt.exp4.tfs.txt --pvalue 1E-8 --seed $i --threads 20
done

# Consolidate
java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar -o /home/ug_megan/tcga_biolinks1/ARACNe.data/KIRC4 --consolidate

