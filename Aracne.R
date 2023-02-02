time java -Xmx5G -jar /home/anafox/ARACNe-AP/dist/aracne.jar \
-e /home/ug_megan/tcga_biolinks1/Tables/filt.exp.df.tpm.prad10.txt \
-o /home/ug_megan/tcga_biolinks1 \
--pvalue 1E-8 --seed 1 \
--calculateThreshold \
--threads 20 