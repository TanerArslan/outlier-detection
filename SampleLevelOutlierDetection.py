from hcluster import pdist, linkage, dendrogram
import numpy
from numpy.random import rand
import matplotlib
import scipy.cluster.hierarchy as sch
from matplotlib import pyplot as plt

from numpy import *
import numpy as np
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from sklearn.cluster import KMeans

class SampleLevelOutDetTCGA():
    ''' Average Hierarchical clustering and
    Silhouette Test is applied for TCGA dataset'''
    #Samples are separated into tumor and normal #
    def classLabel(self,lst):
        self.label=[]
        c='cancer'
        n='normal'
        # If selected data set is TCGA #
        sampleid=['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29']
        for i in range(0,len(lst)):
            code=lst[i]
            for j in sampleid[0:9]:
                if j in code[13:15]:
                    self.label.append(c)
                else:
                    pass
        
            for j in sampleid[9:18]:
                if j in code[13:15]:
                    self.label.append(n)
                else:
                    pass
                                       
        #print len(self.label)
        return self.label

    #Each sample is indexed with respect to its sample type#
    def numbering(self,lst):
        self.newlist=[]
        x=0
        y=0
        for i in self.label:
            if i=='cancer':
                x=x+1
                i='cancer%s'%x
                self.newlist.append(i)
            elif i=='normal':
                y=y+1
                i='normal%s'%y
                self.newlist.append(i)
        #print len(self.newlist)
        return self.newlist
    
    #Main average hierarchical algorithm #
    def hierarchical(self,lst,fulldataset):
        #Samples are colored according to its sample type #
        label_color={}
        for i in self.numbering(self.classLabel(lst)):
            r=('r')
            b=('b')
            if i[0:6]=='cancer':
                label_color[i]=r
                #print label_colors
            elif i[0:6]=='normal' :
                label_color[i]=b
                #print label_colors
            else:
                continue
        tg=zip(*fulldataset)
        Y = pdist(tg)
        #average linkage is applied #
        Z = linkage(Y,method='average')
        sch.set_link_color_palette(['black'])
        a=sch.dendrogram(Z,leaf_font_size=6,labels=self.newlist)
            

        #dendrogram is plotted #
        ax = plt.gca()
        xlbls = ax.get_xmajorticklabels()
    
        for lbl in xlbls:
            lbl.set_color(label_color[lbl.get_text()])
        plt.title("Average Hierarchical Clustering Algorithm")
        plt.savefig('Average Hierarchical Clustering.pdf',dpi=500)
        #plt.show()
        plt.close()

        self.labels=array([])
        c=array([1])
        n=array([0])

        #Silhouette Test #
        #Samples are converted into '0' or '1' for validation #
        for i in self.classLabel(lst):
            if i=='cancer':
                self.labels=np.concatenate([self.labels,c])
            else:
                self.labels=np.concatenate([self.labels,n])

        self.labels=np.delete(self.labels,self.labels[-1])
        self.score=metrics.silhouette_score(Z, self.labels, metric='euclidean')
        #print self.score

class SampleLevelOutDetGEO():
    ''' Average Hierarchical clustering and
    Silhouette Test is applied for TCGA dataset'''
    #Samples are separated into tumor and normal #
    def classLabel(self,lst):
        self.label=[]
        c='cancer'
        n='normal'
        # If selected data set is GEO #
        if lst[0]== 'Tumor' or 'Normal':
            for i in range(0,len(lst)):
                if lst[i]=='Tumor':
                    self.label.append(c)
                else:
                    pass
                
            for i in range(0,len(lst)):
                if lst[i]=='Normal':
                    self.label.append(n)
                else:
                    pass
                                       
        #print len(self.label)
        return self.label

    #Each sample is indexed with respect to its sample type#
    def numbering(self,lst):
        self.newlist=[]
        x=0
        y=0
        for i in self.label:
            if i=='cancer':
                x=x+1
                i='cancer%s'%x
                self.newlist.append(i)
            elif i=='normal':
                y=y+1
                i='normal%s'%y
                self.newlist.append(i)
        #print len(self.newlist)
        return self.newlist
    
    #Main average hierarchical algorithm #
    def hierarchical(self,lst,fulldataset):
        #Samples are colored according to its sample type #
        label_color={}
        for i in self.numbering(self.classLabel(lst)):
            r=('r')
            b=('b')
            if i[0:6]=='cancer':
                label_color[i]=r
                #print label_colors
            elif i[0:6]=='normal' :
                label_color[i]=b
                #print label_colors
            else:
                continue
        tg=zip(*fulldataset)
        Y = pdist(tg)
        #average linkage is applied #
        Z = linkage(Y,method='average')
        sch.set_link_color_palette(['black'])
        a=sch.dendrogram(Z,leaf_font_size=6,labels=self.newlist)
            

        #dendrogram is plotted #
        ax = plt.gca()
        xlbls = ax.get_xmajorticklabels()
    
        for lbl in xlbls:
            lbl.set_color(label_color[lbl.get_text()])
        plt.title("Average Hierarchical Clustering Algorithm")
        plt.savefig('Average Hierarchical Clustering.pdf',dpi=500)
        #plt.show()
        plt.close()

        self.labels=array([])
        c=array([1])
        n=array([0])

        #Silhouette Test #
        #Samples are converted into '0' or '1' for validation #
        for i in self.classLabel(lst):
            if i=='cancer':
                self.labels=np.concatenate([self.labels,c])
            else:
                self.labels=np.concatenate([self.labels,n])

        self.labels=np.delete(self.labels,self.labels[-1])
        self.score=metrics.silhouette_score(Z, self.labels, metric='euclidean')
        #print self.score

SampleLevelTCGA=SampleLevelOutDetTCGA()
SampleLevelGEO=SampleLevelOutDetGEO()        
