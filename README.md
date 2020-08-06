# Outlier Detection Tool for Gene expression and Methylation Data
## Abstract
Expression and methylation datasets are standard genomic techniques and an increasing number of
computational methods are implemented to aid in analyzing the huge and complex amount of generated data. Such
generated datasets often contain a sizeable fraction of outliers that cause misleading results in downstream analysis.
Here, we present a comprehensive approach to detect sample and gene outliers in expression or methylation
datasets. The core algorithms detected most outliers that were artificially introduced by us. Sample outliers detected
by hierarchical clustering are validated by the Silhouette coefficient. At the gene level, the GESD, Boxplot, and MAD
algorithms detected with f-measure of at least 83% the simulated outlier genes in non-intersected distributions.
This combined approach detected many outliers in publicly available datasets from the TCGA and GEO portals.
Frequently, some functionally similar genes marked as outliers turned out to have outlier observations in common
samples. As such cases may be of special interest, they are labeled for further investigations. Expression and DNA
methylation datasets should clearly be checked for outlier points before proceeding with any further analysis. We
suggest that already 2 outlier observations are enough to label an outlier gene as they are enough to ruin a perfect
co-expression. Besides, outliers might also carry useful information and thus functionally similar outliers should be
labeled for further investigation. The presented software is freely available via github.

You can find more in the following link.
https://www.longdom.org/open-access/robust-detection-of-outlier-samples-and-genes-in-expression-datasets-jpb-1000387.pdf
