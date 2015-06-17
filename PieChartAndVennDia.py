import matplotlib.pyplot as plt 
from matplotlib_venn import venn2_unweighted
from Tkinter import *
import tkMessageBox
class PieCharts():
    ''' Pie charts for the detected outlier genes are generated'''
    
    ### Pie chart for Median rule, Generalized ESD, Adjusted box plt and Modified Z-Score###\
    def pie_chart_medianRule_gesd_abplot_mzs(self,x,y,z,t,q,k):
        labels = 'GESD ','Median R.','Common Outliers','Adjusted BoxPlot ','Non-Outliers','M. Z-Score'
        sizes=[len(x),len(y)-len(z),len(z),len(t)-len(z),len(q)+len(z)-len(x)-len(t)-len(y)-len(k),len(k)]
        #sizes=[len(self.gesd_outlier_list),len(self.medianRule_outlier_list)-len(self.common),len(self.common),len(self.abplot_outlier_list)-len(self.common),len(self.cancersets)+len(self.common)-len(self.medianRule_outlier_list)-len(self.abplot_outlier_list)-len(self.gesd_outlier_list)]
        colors = ['gold', 'red','yellowgreen','lightcoral','lightskyblue','blue']
        explode = (0, 0.1,0,0,0,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=270) 
        plt.axis('equal')
        plt.title('Gene Level Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_mzs_gesd_abplot_mzs.pdf') 
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')

    ### Pie chart for Median rule, Generalized ESD and Adjusted box plt ###
    def pie_chart_medianRule_gesd_abplot(self,x,y,z,t,q):
        labels = 'GESD ','Median R.','Common Outliers','Adjusted B.P. ','Non-Outliers'
        sizes=[len(x),len(y)-len(z),len(z),len(t)-len(z),len(q)+len(z)-len(x)-len(t)-len(y)]
        #sizes=[len(self.gesd_outlier_list),len(self.medianRule_outlier_list)-len(self.common),len(self.common),len(self.abplot_outlier_list)-len(self.common),len(self.cancersets)+len(self.common)-len(self.medianRule_outlier_list)-len(self.abplot_outlier_list)-len(self.gesd_outlier_list)]
        colors = ['gold', 'red','yellowgreen','lightcoral','lightskyblue']
        explode = (0, 0.1,0,0,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Level Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_mzs_gesd_abplot.pdf') 
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    ### Pie chart for Median Rule, Modified Z-Score and Adjusted box plt ###
    def pie_chart_medianRule_mzs_abplot(self,x,y,z,t,q):
        labels = 'M. Z-Score','Median R.','Common Outlier','Adjusted B.P. ','Non-Outliers'
        sizes=[len(x),len(y)-len(z),len(y),len(t)-len(z),len(q)+len(z)-len(y)-len(t)-len(x)]
        #sizes=[len(self.mzs_outlier_list),len(self.medianRule_outlier_list)-len(self.common),len(self.common),len(self.abplot_outlier_list)-len(self.common),len(self.cancersets)+len(self.common)-len(self.medianRule_outlier_list)-len(self.abplot_outlier_list)-len(self.mzs_outlier_list)]
        colors = ['gold', 'red','yellowgreen','lightcoral','lightskyblue']
        explode = (0, 0.1,0,0,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_mzs_gesd_abplot.pdf') 
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    ### Pie chart for Median Rule, Modified Z-Score and Generalized ESD ###
    def pie_chart_medianRule_mzs_gesd(self,x,y,z,t,q):
        labels = 'M. Z-Score','Median R.','Common Outlier','GESD ','Non-Outliers'
        sizes=[len(x)-len(z),len(y),len(z),len(t)-len(z),len(q)+len(z)-len(y)-len(x)-len(t)]
        #sizes=[len(self.mzs_outlier_list)-len(self.common),len(self.medianRule_outlier_list),len(self.common),len(self.gesd_outlier_list)-len(self.common),len(self.cancersets)+len(self.common)-len(self.medianRule_outlier_list)-len(self.mzs_outlier_list)-len(self.gesd_outlier_list)]
        colors = ['gold', 'red','yellowgreen','lightcoral','lightskyblue']
        explode = (0, 0.1,0,0,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=0) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_mzs_gesd_abplot.pdf') 
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    
    ### Pie chart for Median Rule and Generalized ESD ####
    def pie_chart_medianRule_gesd(self,x,y,z):
        labels = 'Median R.','GESD ','Non-Outliers'
        sizes=[len(x),len(y),len(z)-len(x)-len(y)]
        #sizes=[len(self.medianRule_outlier_list),len(self.gesd_outlier_list),len(self.cancersets)-len(self.medianRule_outlier_list)-len(self.gesd_outlier_list)]
        colors = ['gold', 'lightskyblue','yellowgreen']
        explode = (0, 0.1,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_abplot_gesd.pdf')
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    ### Pie chart for Median Rule and Modified Z-Score ####
    def pie_chart_medianRule_mzs(self,x,y,z):
        labels = 'M. Z-Score','Median R.','Non-Outliers'
        sizes=[len(x),len(y),len(z)-len(x)-len(y)]
        #sizes=[len(self.mzs_outlier_list),len(self.medianRule_outlier_list),len(self.cancersets)-len(self.medianRule_outlier_list)-len(self.mzs_outlier_list)]
        colors = ['gold', 'lightskyblue','lightcoral']
        explode = (0, 0.1,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_medianRule_mzs.pdf')
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    ### Pie chart for Median Rule ####
    def pie_chart_medianRule(self,x,y):
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        labels = 'Median R.','Non-Outliers'
        sizes=[len(x),len(y)-len(x)]
        #sizes=[len(self.medianRule_outlier_list),len(self.cancersets)-len(self.medianRule_outlier_list)]
        colors = ['gold', 'lightskyblue']
        explode = (0, 0.1)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.03)
        plt.savefig('pie_chart_medianrule.pdf')
        plt.close()
        return

    ### Pie chart for Modified Z-Score ####
    def pie_chart_mzs(self,x,y):
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        labels = 'M. Z-Score','Non-Outliers'
        sizes=[len(x),len(y)-len(x)]
        #sizes=[len(self.mzs_outlier_list),len(self.cancersets)-len(self.mzs_outlier_list)]
        colors = ['gold', 'lightskyblue']
        explode = (0, 0.1)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.03)
        plt.savefig('pie_chart_mzs.pdf')
        plt.close()
        return

    ### Pie chart for Generalized ESD ####
    def pie_chart_gesd(self,x,y):
        labels = 'GESD','Non-Outliers'
        sizes=[len(x),len(y)-len(x)]
        #sizes=[len(self.gesd_outlier_list),len(self.cancersets)-len(self.gesd_outlier_list)]
        colors = ['gold', 'lightskyblue']
        explode = (0, 0.1)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.03)
        plt.savefig('pie_chart_gesd.pdf')
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    ### Pie chart for Adjusted Box plot ####
    def pie_chart_abplot(self,x,y):
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        labels = 'Adjusted B.P.','Non-Outliers'
        sizes=[len(x),len(y)-len(x)]
        #sizes=[len(self.abplot_outlier_list),len(self.cancersets)-len(self.abplot_outlier_list)]
        colors = ['gold', 'lightskyblue']
        explode = (0, 0.1)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=120) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.03)
        plt.savefig('pie_chart_abplot.pdf')
        plt.close()
        
        return

    ### Pie chart for Adjusted Box plot and Mofidied Z-Score ####
    def pie_chart_abplot_mzs(self,x,y,z):
        labels = 'M. Z-Score','Adjusted B.P.','Non-Outliers'
        sizes=[len(x),len(y),len(z)-len(y)-len(x)]
        #sizes=[len(self.mzs_outlier_list),len(self.abplot_outlier_list),len(self.cancersets)-len(self.abplot_outlier_list)-len(self.mzs_outlier_list)]
        colors = ['gold', 'lightskyblue','lightcoral']
        explode = (0, 0.1,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_abplot_mzs.pdf')
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    ### Pie chart for Adjusted Box plot and Generalized Z-Score ####
    def pie_chart_abplot_gesd(self,x,y,z):
        labels = 'Adjusted B.P.','GESD ','Non-Outliers'
        sizes=[len(x),len(y),len(z)-len(x)-len(y)]
        #sizes=[len(self.abplot_outlier_list),len(self.gesd_outlier_list),len(self.cancersets)-len(self.abplot_outlier_list)-len(self.gesd_outlier_list)]
        colors = ['gold', 'lightskyblue','yellowgreen']
        explode = (0, 0.1,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_abplot_gesd.pdf')
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    ### Pie chart for Modified Z-score plot and Generalized Z-Score ####
    def pie_chart_mzs_gesd(self,x,y,z,t):
        labels = 'M. Z-Score','GESD ','Common Outlier','Non-Outliers'
        sizes=[len(x)-len(z),len(y)-len(z),len(z),len(t) + len(z)-len(x)-len(y)]
        #sizes=[len(self.mzs_outlier_list),len(self.gesd_outlier_list),len(self.cancersets)-len(self.mzs_outlier_list)-len(self.gesd_outlier_list)]
        colors = ['gold', 'lightskyblue','yellowgreen','red']
        explode = (0, 0.1,0,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_mzs_gesd.pdf')
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return
    
    ### Pie chart for Adjusted Box plot, Modified Z-Score and Generalized ESD ####
    def pie_chart_abplot_mzs_gesd(self,x,y,z,t,q):
        
        labels = 'M. Z-Score','Adjusted B.P.','Common Outlier','GESD','Non-Outliers'
        sizes=[len(x)-len(z),len(y),len(z),len(t)-len(z),len(q)+len(z)-len(y)-len(x)]
        #sizes=[len(self.mzs_outlier_list)-len(self.common),len(self.abplot_outlier_list),len(self.common),len(self.gesd_outlier_list)-len(self.common),len(self.cancersets)+len(self.common)-len(self.abplot_outlier_list)-len(self.mzs_outlier_list)]
        colors = ['gold', 'red','yellowgreen','lightcoral','lightskyblue']
        explode = (0, 0.1,0,0,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=180) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_mzs_gesd_abplot.pdf') 
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    ### Pie chart for Adjusted Box plot and Median Rule ####
    def pie_chart_abplot_medianrule(self,x,y,z,t):
        labels = 'Median R.','Adjusted B.P.','Common Outlier','Non-Outliers'
        sizes=[len(x)-len(z),len(y)-len(z),len(z),len(t)+len(z)-len(y)-len(x)]
        #sizes=[len(self.medianRule_outlier_list)-len(self.common),len(self.abplot_outlier_list),len(self.common),len(self.cancersets)+len(self.common)-len(self.abplot_outlier_list)-len(self.medianRule_outlier_list)]
        colors = [ 'red','yellowgreen','gold','lightskyblue']
        explode = (0, 0.1,0,0)
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90) 
        plt.axis('equal')
        plt.title('Gene Based Outlier Detection Pie Chart ', bbox={'facecolor':'white'},y=1.05)
        plt.savefig('pie_chart_median ruleabplot.pdf') 
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return


class VennDiagrams():
    '''Pie charts for the detected outlier genes are generated'''

    #### Venn Diagram for Modified Z-Score and Modified Z-Score ###
    def venn_diagram_gesd_mzs(self,algo1,algo2,commons):
        commonOutlier=len(commons)
        first=len(algo1)-commonOutlier
        second=len(lgo2)-commonOutlier
        venn2_unweighted(subsets=(first,second,commonOutlier),set_labels=('modified Z-score','Generalized ESD',''))
        plt.title('Venn Diagram of Dataset')
        plt.savefig('Venn_Diagram_abplot_mzs.pdf')
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return

    #### Venn Diagram for Adjusted box plot and Median Rule ###
    def venn_diagram_abplt_medianrule(self,algo1,algo2,commons):
        commonOutlier=len(commons)
        first=len(algo1)-commonOutlier
        second=len(algo2)-commonOutlier
        venn2_unweighted(subsets=(first,second,commonOutlier),set_labels=('Median Rule','Adjusted Box-Plot',''))
        plt.title('Venn Diagram of Dataset')
        plt.savefig('Venn_Diagram_abplot_mzs.pdf')
        plt.close()
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        return


venndiagram= VennDiagrams()
piechart= PieCharts()


        
