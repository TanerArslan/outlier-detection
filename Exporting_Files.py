### Exporting and Extracting files from tool ####

class ExportingIncompleteGenes():
    ''' Exports incomplete Genes'''
    def ExportIncmopleteGenes(self,incompleteGenes):
        IncompleteGenes=open('Incomplete Genes List.txt','a') 
        IncompleteGenes.write('Row  Genes')
        for line in incompleteGenes:
            IncompleteGenes.write('\n' + line)
        
        IncompleteGenes.close()

class ExportingKEGGgenes():
    ''' Exporting KEGG genes and pathways'''
    def ExportKeggPathways(self,KEGGgenes):
        keggpaths=open('Kegg Pathway Genekeggpaths.writes.txt','a')
        keggpaths.write('Outlier Genes with Participated Pathways')
        for path in KEGGgenes:
            newstr=','.join(path)
            keggpaths.write('\n' + newstr)
            keggpaths.write('\n------------------------------------------------------------------------------------------')
        keggpaths.close()
 


class ExtractingExpressionFiles():
    ''' Extracting Expression files '''
    # Extracting TCGA files #
    def Extracting_files_TCGA(self,lst):
        sample=zip(*lst)
        for i in range(0,len(sample)):
            if i==0:
                pass
            else:
                fl=open('%s.txt'%sample[i][0],"a")
                for k in range(0,len(sample[i])):
                    new_str="%s \t %s "%(sample[0][k],sample[i][k])
                    if '\n' in new_str:
                        new_str=new_str.rstrip()
                        fl.write(new_str)
                    else:
                        new_str=new_str.rstrip()
                        fl.write('\n' + new_str)
                fl.close()

    # Extracting GEO Expression file#                
    def Extracting_files_GEO(self,lst):
        pdata=open("PureDataset.txt","a")
        for gene in lst:
            string=''
            for i in gene:
                string=string + str(i)        
            string=string.replace(",", "\t")
            pdata.write('\n' + string)
        pdata.close()
        return
        


exportingKeggGenes = ExportingKEGGgenes()
exportincompgenes = ExportingIncompleteGenes()
extractFile = ExtractingExpressionFiles()
