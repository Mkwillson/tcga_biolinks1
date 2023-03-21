# tcga_biolinks
CMB3000 Research Project Code

Step 1 - real.pull.data script:
- This contains a singular function for the extraction of mutation, expression and clinical data.
- Data pulled for 6 cancer types

Step 2 - Real.filter script:
- Contains a function for filtering of expression data
- Contains various different filtering attempts for all cancer types
- At the bottom of this script is the JAVA script for running ARACNe-AP in a terminal. First in calculate threshold mode, then on 100 bootstraps, then consolidation of all files into 1 network file.

Step 3 - Rename.Network script:
- The expression networks have genes with ensembl id with versions (e.g. .number on the end) and must be converted to ensembl ids with no version on the end in order to annotate mutation data.
- This code just removes the .number on the end of the ensembl id

Step 5 - All.mutation.data script:
- This script manipulates mutation data in to a 'per gene' format for annotation onto expression networks
- Select 4 categories (High, Moderate, Low and Modifier Impact)
- Select missense mutations
- Create flags (0 or 1) for the SIFT and PolyPhen categories
- Calculate mean, median, maximum and standard deviations for each mutation Mean and PolyPhen scores in each genes
- Create data frame
- Must remove % and NAs to be able to annotate into Cytoscape and use continuous mapping.

Step 6 - Summary.Table script:
- A table to summarise the basic categories within the mutation data frame as this is too big to go in write up.

Step 7 - Centrality script:
- This perfroms a Spearmans rank correlartion test of PolyPhen and SIFT means agaisnt various network parameters
- Plots PolyPhen and SIFT means against network parameters using ggplot2 with LOESS smoothing

Step 8 - BRCA.Stats2, PRAD.Stats2, KIRC.Stats etc
- These scripts contain the code to perform a Kolmogorov-Smirnov (KS) test of filtered mean gene distance of mutated genes against a list of random genes of the size, from the same cancer type.
- Plots mean gene distance distributions

Step 9 - Arrange.BRCA, Arrange.PRAD, Arrange.KIRC etc
- These scripts contain the code to arrange all distribution plots for that cancer type into 1 plot

Step 10 - Comparison.Statistic script:
- Same principle as above using KS test, but on different cancer networks

Step 11 - Wordclouds script:
- Contains the code to select only genes mutated in greater than or equal to 0.5% and with a mean gene distance of less than a randomly generated 1% cut-off
- This list of genes can then be input into only word art software.

Step 12 - Rank.Graph script:
- Code to rank mutations by percent of patients they are seen to be mutated in 
- Plot as a graph with 'genes of interest' highlighted.
