
class removingSameGene():
    ''' Removing same outlier genes genes if normal and tumor samples
        are applied both'''
    
    def removeGene(self,lst):
        print lst
        self.temp=[] # temprorary list
        for gene in lst:
            if gene in self.temp:
                pass
            else:
                self.temp.append(gene)
        lst=self.temp
        print lst
        return

removeOutlierGene=removingSameGene()


