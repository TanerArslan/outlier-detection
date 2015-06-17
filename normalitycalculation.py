import scipy.stats as stats

class normalityCalculation():
    '''For each gene the normality is calculated '''
    def distribution(self,gene,thresholdNorm):
        self.z,self.pval=stats.shapiro(gene[1:])
        if self.pval<thresholdNorm:
            #print 'not normal distribution'
            return self.pval
        else:
            #print'normal'
            return self.pval

normality=normalityCalculation()


