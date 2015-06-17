class helpmenuOutlierDetection():
    def helpscripts(self):
        scripts=('OUTLIER DETECTION MENU \n'
        '--- Gene Level ---\n\nSuppose user chooses "Gene Level"'
        '\nUser can apply four algorithms which are the generalized ESD, modified Z-score, adjusted box plot, median rule for gene level outlier detection. The generalized ESD and modified Z-score perform better than adjusted box plot and median rule when data follows normal distribution.'
        'So, if user applies more than one algorithm, algorithms are applied based on their distribution affinity.For example, suppose user chooses modified Z-score and adjusted box plot,'
        'then modified Z-score is applied for genes that follow normal distribution and adjusted box plot is applied for genes that do not follow normal distribution. However, If user chooses just one algorithm or two algorithms that work better when data follows a normal distribution, algorithms are applied to all genes regardless of distribution.'
        '\n\n ### Threshold Value ###\n'
        'If the number of the outlier gene expression observations is greater than the decided thereshold value, the gene is labbeled as an outlier.'
        '\n\n-- Sample Level --\n User selects the "Average Hierarchical Clustering" algorithm and applies it. The dendrogram will be automatically saved into selected folder where gene expression files exist. The user can also check the Silhouette Test which tells the how well the clustering is performed.')
        return scripts

hlpMENU=helpmenuOutlierDetection()


