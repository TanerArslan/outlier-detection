import csv

class Extracting_outlierGenes():
    ''' After detecting outliers genes, they can be extracted along
        with detected algorithms and gene row.'''
    
    def ExtractOutlierGene(self,list_outliers):
        tutorial_out=open('Outliers.csv','wb')
        mywriter = csv.writer(tutorial_out)
        firstrow=['Algorithm','Gene symbol','Gene row']
        mywriter.writerow(firstrow)
        for item in list_outliers:
            mywriter.writerow(item)
        tutorial_out.close()
        return

extractOutliers=Extracting_outlierGenes()
