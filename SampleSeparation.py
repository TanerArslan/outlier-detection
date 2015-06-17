class sampleSeparation():
    '''Samples are separated for further
    tests and analysis'''
    # cancer sample for GEO data set #
    def cancersampleGEO(self,labelname,gene):
        self.cancer=[]
        for i in range(0,len(labelname)):
            if labelname[i]=='Tumor':
                self.cancer.append(gene[i])
            else:
                pass
        return self.cancer
    # normal sample for GEO data set #            
    def normalsampleGEO(self,labelname,gene):
        self.normal=[]
        for i in range(0,len(labelname)):
            if labelname[i]=='Normal':
                self.normal.append(gene[i])
            else:
                pass
        return self.normal
    # cancer sample for TCGA dataset #
    def cancersample(self,liste,gene):
        self.liste=liste
        self.gene=gene
        self.cancer=[]
        sampleid=['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29']
        for i in range(0,len(self.liste)):
            code=self.liste[i]
            for j in sampleid[0:9]:
                if j in code[13:15]:
                    self.cancer.append(self.gene[i])
                else:
                    continue
        #print self.cancer,'cancer sample'
        return self.cancer
    # normal sample for GEO data set #
    def normalsample(self,liste,gene):
        self.sample=[]
        self.liste=liste
        self.gene=gene
        sampleid=['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29']
        for i in range(0,len(self.liste)):
            code=liste[i]
            for j in sampleid[9:18]:
                if j in code[13:15]:
                    self.sample.append(self.gene[i])
                else:
                    continue
        #print(sample)
        return self.sample

sampleClass=sampleSeparation()
    
