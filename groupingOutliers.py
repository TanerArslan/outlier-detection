from scipy.stats.stats import pearsonr
import xmlrpclib
import mygene
server = xmlrpclib.ServerProxy('http://funsimmat.bioinf.mpi-inf.mpg.de/xmlrpc.php')

class CoExpression_OutlierGene():
    ''' Co-expression is calculated for each pair
    of outlier genes. '''

    ### the need of full expression values is obtained ### 
    def get_data_coexpress(self,fulldata,outlierdata,threshold):
        self.lstfull=[]
        for i in range(0,len(outlierdata)):
            for k in range(0,len(fulldata)):
                if type(outlierdata[i])==list:
                    if outlierdata[i][0]==fulldata[k][0]:
                        self.lstfull.append(fulldata[k])
                    else:
                        pass
                elif type(outlierdata[i])==str:
                    if outlierdata[i]==fulldata[k][0]:
                        self.lstfull.append(fulldata[k])
                    else:
                        pass
                else:

                    pass
        return self.lstfull

    ### main co-expression of each pair of outlier genes ###
    def co_express_calculation(self,completedata,outliergenes,threshold):
        gene_data=self.get_data_coexpress(completedata,outliergenes,threshold)
        self.co_expressed=[]
        for i in range(0,len(gene_data)):
            for k in range(i+1,len(gene_data)):
                result=pearsonr(gene_data[i][1:],gene_data[k][1:])
                if result[0]>float(threshold) or result[0] <-float(threshold):
                    if gene_data[i][0] in self.co_expressed:
                        pass
                    else:
                        self.co_expressed.append(gene_data[i][0])
                    if gene_data[k][0] in self.co_expressed:
                        pass
                    else:
                        self.co_expressed.append(gene_data[k][0])
                else:
                    pass
    
        if len(self.co_expressed)> (len(outliergenes)*2.0/100):
            #print 'greater than 2'
            return self.co_expressed
        else:
            #print 'not greater'
            self.co_expressed=[]
        #print len(self.co_expressed)
        #return self.co_expressed


class KEGG_Pathway_Participation():
    ''' Outlier genes are checked if they belong to
    same set of pathways'''

    # Opening reference KEGG reference file and matching pathways with outlier genes #
    def keggPathway_merge_data(self,genelist,filePath):
        finalpathway=[] # includes the genes and pathways of KEGG file
        pathways=open(filePath)
        for line in pathways.readlines():
            temp=[]
            data=line
            data1=data.split("\t")
            for i in data1:
                if len(i)>2:
                    temp.append(i)
                else:
                    pass
            finalpathway.append(temp)
    
        outgenepath=[] # outlier gene are matched corresponding pathways
        for k in range(0,len(finalpathway)):
            for j in range(0,len(genelist)):
                genepath=[]
                if type(genelist[j])==list:
                    if genelist[j][0] in finalpathway[k]:
                        genepath.append(finalpathway[k][0])
                        genepath.append(genelist[j][0])
                        outgenepath.append(genepath)
                    else:
                        pass
                elif type(genelist[j])==str:
                    if genelist[j] in finalpathway[k]:
                        genepath.append(finalpathway[k][0])
                        genepath.append(genelist[j])
                        outgenepath.append(genepath)
                    else:
                        pass
                else:
                    pass
        return outgenepath,finalpathway

    # outlier genes merged if they paricipate same set of pathways
    def merge_genes_kegg(self,genelist,filepath):
        outgenepath,finlpathwy=self.keggPathway_merge_data(genelist,filepath)
        listpathwaygene= []
        for row in outgenepath:
            for i, resrow in enumerate(listpathwaygene):
                if row[0]==resrow[0]:
                    listpathwaygene[i] += row[1:]
                    break
            else:
                listpathwaygene.append(row)
        #print listpathwaygene
        return listpathwaygene

    # Comparing outier genes if they pass certain thresholf #
    def keggPathwaySim(self,genelist,filepath,threshold):
        listpathwaygene=self.merge_genes_kegg(genelist,filepath)
        outgenepath,originalpath=self.keggPathway_merge_data(genelist,filepath)
        finallist=[]
        for i in range(0,len(listpathwaygene)):
            for k in range(0,len(originalpath)):
                if listpathwaygene[i][0]==originalpath[k][0]:
                    if len(listpathwaygene[i])>= ((len(originalpath[k])*(float(threshold)/100))): 
                        finallist.append(listpathwaygene[i])
                    else:
                        pass
                else:
                    pass
        #removing the same genes to avoid redundancy #
        self.kegg_outlier=[]
        for lst in finallist:
            for gene in lst[1:]:
                if gene in self.kegg_outlier:
                    pass
                else:   
                    self.kegg_outlier.append(gene)
        self.finallist=finallist
        
        return self.kegg_outlier

class FunSimMat_calculation():
    ''' Functional Similarity is calculated
        for each pair of genes.'''

    #Gene symbol is converted to uniprot id using mygene module #
    def get_uniprot_id(self,genelist):
        mg=mygene.MyGeneInfo()
        genesymbol=[]
        gene_uniprot_id=[]
        gene_uniprot_id_single=[]
        genesymbol_uniprotid={}
        for i in range(0,len(genelist)):
            if type(genelist[i])==list:
                genesymbol.append(genelist[i][0])
            elif type(genelist[i])==str:
                genesymbol.append(genelist[i])
            else:
                pass

        out = mg.querymany(genesymbol, scopes='symbol', fields='uniprot', species='human')
        for k in range(0,len(out)):
            try:
                gene_uniprot_id.append(out[k]['uniprot']['Swiss-Prot'])
                genesymbol_uniprotid[genelist[k]]=out[k]['uniprot']['Swiss-Prot']
            except:
                pass
        for j in range(0,len(gene_uniprot_id)):
            if type(gene_uniprot_id[j])==list:
                gene_uniprot_id_single.append(gene_uniprot_id[j][0])
            else:
                gene_uniprot_id_single.append(gene_uniprot_id[j])

        #print genesymbol_uniprotid    
        #print len(gene_uniprot_id_single)
        return gene_uniprot_id_single,genesymbol_uniprotid

    def funsimmat(self,genelist,thresholdFSM):
        uniprotId,genesymbol_uniprotid_dict=self.get_uniprot_id(genelist)
        simlist=[]
        strgene=""
        for k in uniprotId:
            strgene=",".join((strgene,str(k)))
        
        #Computing FunSimMat #
        for i in uniprotId:
            query=str(i)
            list = strgene
            result = server.Functional.getScoresGp2List(query, list)
            header = result[0]
            for k in range(1,len(result)):
                acc1 = result[k][0]
                acc2 = result[k][1]
                try:
                    score1 = float(result[k][2])
                except:
                    continue
                if score1 >float(thresholdFSM):
                    if acc1!=acc2:
                        if acc1 in simlist:
                            pass
                        else:  
                            simlist.append(acc1)
                        if acc2 in simlist:
                            pass
                        else:
                            simlist.append(acc2)
                    else:
                        pass
                else:
                    pass
        #print simlist

        #uniprot ids are converted back to gene symbols #  
        self.genelist_fromuniprot=[]

        for k,v in genesymbol_uniprotid_dict.items():
            if type(v)==unicode:
                if v in simlist:
                    self.genelist_fromuniprot.append(k)
                else:
                    pass
            else:
                self.genelist_fromuniprot.append(k)

        #print len(self.genelist_fromuniprot)
            
CoExpression = CoExpression_OutlierGene()
KEGGpp = KEGG_Pathway_Participation()
FunSimMat=FunSimMat_calculation()
