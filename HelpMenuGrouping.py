class helpmenuGrouping():
    def helpscripts(self):
        scripts=('## Grouping Outliers ## \n\n--- Functional Similarity (FunSimMat) ---\n'
        'We have used FunSimMat database to compute functional similarity between outlier genes. We connected the server \n"http://funsimmat.bioinf.mpi-inf.mpg.de/xmlrpc.php" in order to query outlier genes. Therefore user"s machine has to connect internet while analyzing the dataset. The asked threshold is to determine whether two genes are functional similar. In our tool user can decide the threshold and enter it to corresponding place.\n'
        '\n--- Co-expression ---\n'
        'Pearson correlation was used to measure co-expression. For all pair of genes, Pearson correlation that ranges between -1 and 1 was computed. Only pair of co-expressed genes that are above the threshold are stored.\n'
        '\n--- KEGG Pathway Participation ---\n'
        'We used human genes involved in metabolic pathways derived from KEGG Pathway database (Downloaded from http://www.broadinstitute.org/gsea/msigdb/collections.jsp and converted to tab delimited text file). We used these data as a reference to match outlier genes with corresponding pathways. If outlier genes include more than 20% (which can be adjust by user) of the genes participate in a pathway, we label them as non-outlier.')
        #print scripts
        return scripts

helpMG=helpmenuGrouping()
helpMG.helpscripts()
