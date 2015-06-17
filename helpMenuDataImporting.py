class helpmenuForImporting():
    def helpscripts(self):
        scripts=(' IMPORTING DATASET \n\n1. Importing Dataset from TCGA\n'
        'You just need to select the path where TCGA expression files exist. Please make sure that there should not be any other text file except expression files in selected folder !!!\n'
        '\n2. Importing Dataset from GEO\n'
        'Our tool takes particular pattern of dataset as an input. You should modify your dataset prior to analysis. Dataset must be in tab delimited text file.'
        '\n\n--Example of GEO Dataset--\n\n'
        'Empty  Tumor  Tumor  Tumor   Normal\n'
        'Gene1\nGene2\nGene3\nGene4\n.\n.\n.\nGenen\n\n'
        '\n### P-value threshold ###\n\n'
        'We use Shapiro-Wilk test to calculate the normaliy of expression values of each gene. P-values is used to test the hypothesis. We let you to use your own p-value. Generally 0.05 is widely accepted and used however we did not force user to use p-value which we decided as default value. So user can decide p-value by themselves for their own project')
        #print scripts
        return scripts

hlpmenu=helpmenuForImporting()


