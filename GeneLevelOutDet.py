import math
import numpy as np
from scipy import stats
from math import sqrt
from math import exp
import scipy.stats as stats

class GeneLevelOutDet():
    '''Gene level outlier detection algorithms
    are applied'''
    
    ########### Modified Z-Score ###########
    
    #Calculating the MAD #
    def mad(self,gene):
        gene = np.ma.array(gene[1:]).compressed() # should be faster to not use masked arrays.  
        med = np.median(gene)
        #print("mad",np.median(np.abs(gene - med)))
        return np.median(np.abs(gene - med))

    # Calculating median #
    def median(self,gene): 
        sortedgene = sorted(gene[1:]) 
        if len(sortedgene) % 2 == 0: 
            n = len(sortedgene)//2
            #print ("median",sortedgene[n]+sortedgene[n-1]/2)
            return (sortedgene[n]+sortedgene[n-1])/2
        else:
            #print("median:",sortedgene[len(sortedgene)//2])
            return sortedgene[len(sortedgene)//2]

    #Main Modified Z-Score calculation for gene #
    def modified_z_score(self,gene,threshold=3.5):
        z_c=[]
        medianvalue=self.median(gene)
        madvalue=self.mad(gene)
        for i in gene[1:]:
                mzscore= abs(0.6745*(i-medianvalue)/madvalue)
                if mzscore>threshold :
                    z_c.append(i)
                    #print(len(z_c))
                else:
                    pass
                    #print("no outlier detected",i)
        #print(mzscore>threshold)
        return z_c

    #Checking whether gene is outlier #
    def outlier_detection_modified_z_score(self,gene,thresholdGeneLevel):
        z_c=self.modified_z_score(gene,threshold=3.5)
        #print(len(z_c))
        count=len(gene[1:])-1
        if len(z_c) > float(thresholdGeneLevel):
            #print('gene is outlier')
            #print [gene[0]] + z_c
            return [gene[0]] + z_c
        else:
            pass
            #print('gene is not outlier')
        return


    ######### GESD Algorithm ########
        
    def main_gesd(self,gene,thresholdGeneLevel):
        n=len(gene[1:])
        diff=[]
        final_outlier_gesd=[]
        for k in range(0, int(n/2)+1):
            i=k+1
            mean=float(sum(gene[1:])/len(gene[1:]))      # mean calculation
            sd=sqrt(sum((x-mean)**2 for x in gene[1:])/(len(gene[1:])-1))   # standard deviation
            alpha=0.05
            p=1-alpha/(2*(n-i+1)) # calculating p-value
            t=stats.t.ppf(p,n-i-1)
            cv=t*(n-i)/sqrt((n-i-1+t**2)*(n-i+1)) #calculating r critical calues 
            po=[]
            for j in gene[1:]:
                if sd==0:
                    a=0
                    po.append(a)
                else:
                    R=(abs(float(j-mean)/sd)) # calculating r test for each expression values
                    po.append(R)
            Rmax=max(po)
            gene_index=po.index(max(po))
            final_outlier_gesd.append(gene[gene_index+1])
            difference=Rmax-cv
            diff.append(difference)
            #print(i,'=>',Rmax)
            index=po.index(Rmax)
            del gene[int(index)+1] # removing the value that maximize the R test 
            #print(gene)
        final=max(diff)
        indx=diff.index(final)
        if final > 0:
            b=5
            #print('there are ',int(indx+1),'outliers in the given data sets') 
            if int(indx+1)>float(thresholdGeneLevel):
                #print('gene is outlier')
                #print gene[0]
                return [gene[0]] + final_outlier_gesd[0:indx+1]
            else:
                pass
                #print('gene is not outlier')
        else:
            pass
            #print('there are no outliers')

        return

    ########### Adjusted Box Plot ##############
        
    #Calculating median#
    def median_abp(self,gene): 
        sortedgene = sorted(gene[1:]) 
        if len(sortedgene) % 2 == 0: 
            n = len(sortedgene)//2
            #print ("median",(sortedgene[n]+sortedgene[n-1])/2.0)
            return (sortedgene[n]+sortedgene[n-1])/2
        else:
            #print("median:",sortedgene[len(sortedgene)//2])
            return sortedgene[len(sortedgene)//2]

    #Calculating the lower quartile of box plot #
    def lower_q1(self,gene):
        sortedgene=sorted(gene[1:])
        mdn=self.median_abp(gene)
        nr=(len(gene[1:])-1)*0.25
        i=int(nr)
        f=nr-int(nr)
        if f==0:
            #print ('q1:',sortedgene[i])
            return sortedgene[i]
        else:
            #print('q1:',sortedgene[i] + f*(sortedgene[i+1]-sortedgene[i]))
            return sortedgene[i] + f*(sortedgene[i+1]-sortedgene[i])

    # Calculating the upper quartile #
    def upper_q3(self,gene):
        sortedgene=sorted(gene[1:])
        mdn=self.median_abp(gene)
        nr=(len(gene[1:])-1)*0.75
        i=int(nr)
        f=nr-int(nr)
        if f==0:
            #print ('q1:',sortedgene[i])
            return sortedgene[i]
        else:
            #print('q1:',sortedgene[i] + f*(sortedgene[i+1]-sortedgene[i]))
            return sortedgene[i] + f*(sortedgene[i+1]-sortedgene[i])

    # Calculating the inter quartile range # 
    def iqr(self,gene):
        sortedgene=sorted(gene[1:])
        mdn=self.median_abp(gene)
        Q1=self.lower_q1(gene)
        Q3=self.upper_q3(gene)
        IQR=float(Q3-Q1)
        #print("IQR",IQR)
        return IQR


    #Calculating the medcouple for skewness #
    def medcouple(self,gene):
        sortedgene=sorted(gene[1:])
        mdn=self.median_abp(gene)
        mdcouple=[]
        if len(sortedgene) %2 ==1:
            for j in sortedgene[len(sortedgene)//2:]:
                for i in sortedgene[0:len(sortedgene)//2]:
                    if i==j:
                        pass
                    else:
                        mc=((j-mdn)-(mdn-i))/(j-i)
                        mdcouple.append(mc)
        else:
            for j in sortedgene[len(sortedgene)/2:]:
                for i in sortedgene[0:len(sortedgene)/2]:
                    if i==j:
                        pass
                    else:
                        mc=((j-mdn)-(mdn-i))/(j-i)
                        mdcouple.append(mc)
        #print('mdcouple',mdcouple)
        return(mdcouple)

    #Calculating the median of medcouple #
    def medcouplemedian(self,gene):
        mdcouple=self.medcouple(gene)
        sortedmdcouple = sorted(mdcouple)
        if len(sortedmdcouple) % 2 == 0: 
            n = len(sortedmdcouple)//2
            #print ("median",(sortedmdcouple[n]+sortedmdcouple[n-1])/2.0)
            return (sortedmdcouple[n]+sortedmdcouple[n-1])/2.0
        else:
            #print("median:",sortedmdcouple[len(sortedmdcouple)//2])
            return sortedmdcouple[len(sortedmdcouple)//2]

    #Calculating the lower and upper fences for adjusted box plot#
    def lower_upperFence(self,gene):
        Q1=self.lower_q1(gene)
        Q3=self.upper_q3(gene)
        IQR=self.iqr(gene)
        MC=self.medcouplemedian(gene)
        #print(MC)
        if MC > 0:
            lower=(Q1 - 1.5*exp((-3.5)*(MC))*IQR)
            upper=(Q3 + 1.5*exp((4*MC))*IQR)
            #print 'lower:',lower,'upper:',upper
            return lower,upper
        elif MC < 0:
            lower=(Q1 - 1.5*exp((-4.0)*(MC))*IQR)
            upper=(Q3 + 1.5*exp((3.5)*(MC))*IQR)
            #print 'lower:',lower,'upper:',upper
            return lower,upper
        else:
            lower=(Q1-1.5*IQR)
            upper=(Q3+1.5*IQR)
            return lower,upper

    #Outlier detection in adjusted box plot for each gene expression value #
    def outlier_detection_adjusted_box_plot(self,gene,thresholdGeneLevel):
        values=self.lower_upperFence(gene)
        lowervalue,uppervalue=values[0],values[1]
        outlierlist=[]
        for i in gene[1:]:
            if i<lowervalue or i>uppervalue:
                outlierlist.append(i)
            else:
                continue
        #print(outlierlist)
        if len(outlierlist)>float(thresholdGeneLevel):
            #print'gene is outlier'
            return [gene[0]] + outlierlist 
        else:
            pass
            #print'gene is not outlier'

        return


    ########## Median Rule ###########
        
    # Calculating the intervals for median rule #
    # interquartile range and median are taken from adjusted box plot # 
    def median_rule(self,gene):
        IqR=self.iqr(gene)
        mdn=self.median(gene)
        C1=mdn - 2.3 *IqR
        C2=mdn + 2.3 *IqR
        return C1,C2

    #Checking for outliers #
    def median_rule_outlier_det(self,gene,thresholdGeneLevel):
        lower_interval,upper_interval=self.median_rule(gene)
        outlierlist_median_rule=[]
        for i in gene[1:]:
            if i<lower_interval or i>upper_interval:
                outlierlist_median_rule.append(i)
            else:
                continue
        #print(outlierlist)
        if len(outlierlist_median_rule)>float(thresholdGeneLevel):
            #print 'outlier'
            #print [gene[0]] +outlierlist_median_rule  
            return [gene[0]] +outlierlist_median_rule  
        else:
            pass
            #print'gene is not outlier'
        return

    ######### Common Outliers for Adjusted Box plot and Median Rule #######


    def common_outlier(self,list1,list2,gene):
        common_list=[]
        for sublist in list1:
            for sublist2 in list2:
                if sublist[0]==sublist2[0]:
                    commonOutlierss=len(set(sublist) & set(sublist2))-1
                    if commonOutlierss > len(gene[1:])*2.0/100:
                        if sublist2[0] in common_list:
                            pass
                        else:
                            common_list.append(sublist2[0])
                    else:
                        pass
                else:
                    pass
        #print len(common_list)
        return common_list

    ########### Combining Outlier genes for Generalized ESD and Modified Z-Score

    def singleListforMzsGesd(self,lst1,lst2):
        commonMzsGesd_list=[]
        for i in range (0,len(lst1)):
            if lst1[i] in commonMzsGesd_list:
                pass
            else:
                commonMzsGesd_list.append(lst1[i])

            for k in range(0,len(lst2)):
                if lst1[i][0]==lst2[k][0]:
                    if lst2[k] in commonMzsGesd_list:
                       pass
                else:
                    if lst2[k] in commonMzsGesd_list:
                        pass
                    else:
                        commonMzsGesd_list.append(lst2[k])
        

        #print len(commonMzsGesd_list)
        return commonMzsGesd_list

    ################# Deleting gene level outliers #############
    
    def deleting_outliers(self,data,genes):
        deleting_index=[]
        for gene in genes:
            for k in range(0,len(data)):
                if gene==data[k][0]:
                    deleting_index.append(k)
                else:
                    pass
        for index in sorted(deleting_index,reverse=True):
            del data[index]
        #print len(data)
        return data        
GeneLevel=GeneLevelOutDet()    
   
