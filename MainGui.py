from Tkinter import *
import tkColorChooser
import tkFileDialog
import tkMessageBox
from tkFileDialog import askopenfilename
import glob
import os
from numpy import *
from SampleSeparation import *
from normalitycalculation import *
from Exporting_Files import *
from GeneLevelOutDet import *
from PieChartAndVennDia import *
from groupingOutliers import *
from SampleLevelOutlierDetection import *
from GetSingleOutlierGene import *
from helpMenuDataImporting import *
from HelpMenuOutlierMenu import *
from HelpMenuGrouping import *
from ExtractionOutliers import *
class Tool(Frame):
    ''' Main GUI for outlier detection'''

    # Initializing the tool #
    def __init__(self,master):
        Frame.__init__(self,master)
        self.grid()
        self.radiobutton()

        self.helpmenudatabutton=Button(mGui,text='Help',command=self.helpmenu1)
        self.helpmenudatabutton.place(relx=0.42,rely=0.9)

    # Help menu for the data selection and p value for distribution #
    def helpmenu1(self):
        tkMessageBox.showwarning(title='Help',message=hlpmenu.helpscripts())
        return
    
    # Creating radio button for data set selection #
    def radiobutton(self):
        self.radiobuttonLabel=Label(mGui,text='Please select the dataset source',font='Times 18 bold italic')
        self.radiobuttonLabel.place(relx=0.25,rely=0.25)

        self.favorite=StringVar()
        
        self.TCGARadiobutton=Radiobutton(mGui,text='TCGA',value=1,variable=self.favorite,command=self.TCGARadiobutton)
        self.TCGARadiobutton.place(relx=0.3,rely=0.4)

        self.GEORadiobutton=Radiobutton(mGui,text='GEO',value=2,variable=self.favorite,command=self.GEORadiobutton)
        self.GEORadiobutton.place(relx=0.5,rely=0.4)

        self.startRadiobutton=Button(mGui,text='Start',command=self.radiobuttonstart)
        self.startRadiobutton.place(relx=0.3,rely=0.6)

        self.cancelRadioButton=Button(mGui,text="Quit",command=self.mQuit)
        self.cancelRadioButton.place(relx=0.5,rely=0.6)

    # assigning TCGA radiobutton for further steps #
    def TCGARadiobutton(self):
        self.selectedRbutton=1

    # assigning GEO radiobutton for further steps #
    def GEORadiobutton(self):
        self.selectedRbutton=2

    # next step after selecting data type #
    def radiobuttonstart(self):
        try:
            if self.selectedRbutton==1 or self.selectedRbutton==2: 
                # removing buttons # 
                self.radiobuttonLabel.destroy()
                self.TCGARadiobutton.destroy()
                self.GEORadiobutton.destroy()
                self.startRadiobutton.destroy()
                self.cancelRadioButton.destroy()                                                    
                # TCGA dataset #
                if self.selectedRbutton==1:
                    self.browse()
                # GEO data set#
                else:
                    self.GEOselectbutton=Button(mGui,text='Select File',command=self.selectFileGEO)
                    self.GEOselectbutton.place(relx=0.73,rely=0.25)
                    self.GEOEntry=Entry(mGui,width=35)
                    self.GEOEntry.place(relx=0.13,rely=0.25)
        # warning when start button first is clicked, before selecting data source
        except:
            tkMessageBox.showwarning(title='Attention',message='Please select data source first !!!')
    # Selecting GEO data set #
    def selectFileGEO(self):
        self.GEOfilename=askopenfilename(parent=mGui)
        self.GEOEntry.delete(0,END)
        self.GEOEntry.insert(0,self.GEOfilename)

        #Creating buttons #
        self.startbuttonGEO=Button(mGui,text='Start',command=self.readdataGEO)
        self.startbuttonGEO.place(relx=0.3,rely=0.4)
        self.cancelbuttonGEO=Button(mGui,text="Quit",command=self.mQuit)
        self.cancelbuttonGEO.place(relx=0.6,rely=0.4)    

    #Opening and Reading GEO dataset and merging #
    def readdataGEO(self):
        fileGEO=open(self.GEOEntry.get())
        self.data_set=[]
        for i,line in enumerate(fileGEO.readlines()):
            data=line
            data1=data.split("\t")
            data2=[i.strip() for i in data1]
            self.data_set.append(data2)

        #P value for shapiro test #            
        self.distthresholdlabel=Label(mGui,text='Input Shapiro test treshold')
        self.distthresholdlabel.place(relx=0.25,rely=0.55)
        self.distthresholdEntry=Entry(mGui,width=6)
        self.distthresholdEntry.place(relx=0.61,rely=0.55)
        self.flterButton=Button(mGui,text='Filter Data',command=self.filter_data)
        self.flterButton.place(relx=0.4,rely=0.7)

    # Creating browse button for TCGA data set  #
    def browse(self):
        self.entry=Entry(mGui,width=30)
        self.entry.place(relx=0.2,rely=0.2)
        self.pathButton=Button(mGui)
        self.pathButton.config(text="Browse",font=10,command=self.path)
        self.pathButton.place(relx=0.75,rely=0.2)
        
    #Selecting the path where expression files are exist #
    def path(self):
        mGui.fileName =  tkFileDialog.askdirectory(parent=mGui, title='Select Directory')
        pathName=mGui.fileName
        self.entry.delete(0,END)
        self.entry.insert(0,pathName)
        self.StartButton=Button(mGui,text="Start",font=10,command=self.combining_genes)
        self.StartButton.place(relx=0.3,rely=0.4)
        self.CancelButton=Button(mGui,text="Quit",font=10,command=self.mQuit)
        self.CancelButton.place(relx=0.6,rely=0.4)

    #Opening single text files, reading and merging #
    def combining_genes(self):
        directory=self.entry.get()
        os.chdir(directory)
        list_of_files=glob.glob('*.txt')
        self.data_set=[]
        for i,filename in enumerate(list_of_files):
            data_list=open(filename,'r')
            if i==0:
                for j,line in enumerate(data_list.readlines()):
                    data=line
                    data1=data.split("\t")
                    self.data_set.append(data1)
            else:
                for j,line in enumerate(data_list.readlines()):
                    temp_data=[]
                    data=line
                    data1=data.split("\t")
                    for k in self.data_set[j]:
                        temp_data.append(k)
                    temp_data.append(data1[1])
                    self.data_set[j]=temp_data
            data_list.close()
        ##P value for shapiro test #
        self.distthresholdlabel=Label(mGui,text='Input Shapiro test treshold')
        self.distthresholdlabel.place(relx=0.25,rely=0.55)
        self.distthresholdEntry=Entry(mGui,width=6)
        self.distthresholdEntry.place(relx=0.61,rely=0.55)
        #Creating button for further step #
        self.flterButton=Button(mGui,text='Filter Data',command=self.filter_data)
        self.flterButton.place(relx=0.4,rely=0.65)

    #Filterign data set #
    def filter_data(self):
        #checking if p value is input #
        try:
            self.thresholdforPvalue=float(self.distthresholdEntry.get())
        except:
            tkMessageBox.showwarning(title='Attention',message='Please enter numeric value for threshold')
        self.cancerandnormaldatafull=[]
        self.incomplete_genes=[]
        self.dist=[]
        #checking if correct data set is selected #
        try:
            self.list1=self.data_set[0][1:]
        except:
            tkMessageBox.showwarning(title='Attention',message='Please select correct directory or check txt files and Restart the program')
        self.cancersets=[]
        self.full_dataset=[]
        self.normalsets=[]
        #data set filtering for TCGA  #
        if self.selectedRbutton==1:
            for i in range(0,len(self.data_set)):
                if i>=2:
                    data=[]
                    temporary=self.data_set[i]
                    for k in temporary[1:]:
                        try:
                            tempfloat=float(k)
                            data.append(tempfloat)
                        except:
                            pass
                    data=[x for x in data if x]
                    temp=[]
                    for j in data:
                        if type(j)==float:
                            temp.append(j)
                        else:
                            pass
                    if len(temp)==(len(temporary)-1):
    
                        self.full_dataset.append(data)
                        #separating samples into tumor and normal #
                        data2=sampleClass.cancersample(self.list1,data)
                        self.data3=[temporary[0]]+data2 # tumor sample
                        self.cancersets.append(self.data3)
                        data4=sampleClass.normalsample(self.list1,data)
                        self.data5=[temporary[0]]+data4 #normal sample
                        self.normalsets.append(self.data5)
                        self.cancerNormal=[temporary[0]]+data2 +data4 
                        self.cancerandnormaldatafull.append(self.cancerNormal)
                        try:
                            pval=normality.distribution(self.data3,self.thresholdforPvalue)
                        except:
                            tkMessageBox.showwarning(title='Attention',message='Please select correct directory or check txt files and Restart the program')
                        #normality calculation #
                        if pval>self.thresholdforPvalue:
                            self.dist.append(self.data3[0])
                        else:
                            pass
                    #incomplete genes are detected #
                    else:
                        line=str(i+1)+str('   ')+str(temporary[0])
                        self.incomplete_genes.append(line)

        #Data set filtering fot GEO #
        else:
            for i in range(0,len(self.data_set)):
                if i>=1:
                    data=[]
                    temporary=self.data_set[i]
                    for k in temporary[1:]:
                        try:
                            tempfloat=float(k)
                            data.append(tempfloat)
                        except:
                            pass
                    data=[x for x in data if x]
                    #data=[float(a) for a in data[1:]]
                    temp=[]
                    for j in data:
                        if type(j)==float:
                            temp.append(j)
                        else:
                            pass
                    if len(temp)==(len(temporary)-1):
    
                        self.full_dataset.append(data)
                        #separating samples into tumor and normal #
                        data2=sampleClass.cancersampleGEO(self.list1,data)
                        self.data3=[temporary[0]]+data2 # tumor samples
                        self.cancersets.append(self.data3)
                        data4=sampleClass.normalsampleGEO(self.list1,data)
                        self.data5=[temporary[0]]+data4 #normal samples
                        self.normalsets.append(self.data5)
                        self.cancerNormal=[temporary[0]]+data2 +data4
                        self.cancerandnormaldatafull.append(self.cancerNormal)
                        try: 
                            pval=normality.distribution(self.data3,self.thresholdforPvalue)
                        except:
                            tkMessageBox.showwarning(title='Attention',message='Please select correct file and Restart the program')
                        #normality calculation #
                        if pval>self.thresholdforPvalue:
                            self.dist.append(self.data3[0])
                    #incomplete genes are detected #
                    else:
                        line=str(i+1)+str('   ')+str(temporary[0])
                        self.incomplete_genes.append(line)
        #normal samples are added data if number of normal samples are > 50 #
        if len(self.normalsets[0])>50:
            self.cancersets=self.cancersets + self.normalsets
        else:
            pass
        #Destroying buttons #
        self.helpmenudatabutton.destroy()
        #if TCGA is selected #
        try:
            self.distthresholdlabel.destroy()
            self.distthresholdEntry.destroy()
            self.StartButton.destroy()
            self.CancelButton.destroy()
            self.pathButton.destroy()
            self.entry.destroy()
            self.flterButton.destroy()
        #if GEO is selected #
        except:
            self.GEOEntry.destroy()
            self.GEOselectbutton.destroy()
            self.startbuttonGEO.destroy()
            self.cancelbuttonGEO.destroy()
            self.distthresholdlabel.destroy()
            self.distthresholdEntry.destroy()
            self.flterButton.destroy()

        #Reporting the characteristic of data set #
        self.datahead=Label(mGui,text ='Characteristics of Data Set',font='Times 18 bold italic')
        self.datahead.place(relx=0.1,rely=0.2)
        self.geneText=Text(mGui,height=1,width=50,wrap=None)
        self.geneText.place(relx=0.1,rely=0.3)
        gText="Number of Genes :                           "+str(len(self.data_set))
        self.geneText.insert(0.0,gText)
        self.sampleText=Text(mGui,height=1,width=50)
        self.sampleText.place(relx=0.1,rely=0.35)
        sText="Number of Samples:                          "+ str(len(self.data_set[0][1:]))
        self.sampleText.insert(0.0,sText)
        self.cancerText=Text(mGui,height=1,width=50)
        self.cancerText.place(relx=0.1,rely=0.4)
        cText="Number of Cancer Samples :                  "+ str(len(self.data3[1:]))
        self.cancerText.insert(0.0,cText)
        self.normalText=Text(mGui,height=1,width=50)
        self.normalText.place(relx=0.1,rely=0.45)
        nText="Number of Normal Samples:                   "+ str(len(self.data5[1:]))
        self.normalText.insert(0.0,nText)
        self.incompleteText=Text(mGui,height=1,width=50)
        self.incompleteText.place(relx=0.1,rely=0.5)
        iText="Number of incomplete Genes:                 "+ str(len(self.incomplete_genes))
        self.incompleteText.insert(0.0,iText)
        self.distText=Text(mGui,height=1,width=50)
        self.distText.place(relx=0.1,rely=0.55)
        dText="Number of genes follow normal distribution: " +str(len(self.dist))
        self.distText.insert(0.0,dText)

        #Creating the buttons for next steps #
        self.goForAnalysis=Button(mGui,text='Outlier Menu',command=self.goForAnalysis)
        self.goForAnalysis.place(relx=0.07,rely=0.7)
        self.CancelButton=Button(mGui,text="Quit",command=self.mQuit)
        self.CancelButton.place(relx=0.72,rely=0.7)
        self.exportincompletegenes=Button(mGui,text='Export Incomplete Genes',command=self.ExportIncmopleteGene)
        self.exportincompletegenes.place(relx=0.32,rely=0.7)

    #Exporting Incomplete Genes where expression file exists #
    def ExportIncmopleteGene(self):
        exportincompgenes.ExportIncmopleteGenes(self.incomplete_genes)
        tkMessageBox.showinfo(title='Info',message='File were successfully exported')

    #Propceeding sample/gene level outlier detection #
    def goForAnalysis(self):
        #Deleting buttons #
        self.datahead.destroy()
        self.geneText.destroy()
        self.sampleText.destroy()
        self.cancerText.destroy()
        self.normalText.destroy()
        self.incompleteText.destroy()
        self.distText.destroy()
        self.goForAnalysis.destroy()
        self.CancelButton.destroy()
        self.exportincompletegenes.destroy()
        #Creating buttons#
        self.selectLabel=Label(mGui,text='Please select outlier detection level type',font='Times 14 bold italic')
        self.selectLabel.place(relx=0.22,rely=0.2)
        self.genelLevelButton=Button(mGui,text='Gene Level  ',command=self.checkbutgene)
        self.genelLevelButton.place(relx=0.22,rely=0.27)
        self.sampleLevelButton=Button(mGui,text='Sample Level',command=self.checkbutsample)
        self.sampleLevelButton.place(relx=0.52,rely=0.27)
        self.helpmenualgorithmsbutton=Button(mGui,text='Help',command=self.helpmenu2)
        self.helpmenualgorithmsbutton.place(relx=0.44,rely=0.93)

    #Sample level outlier detection menu #
    def checkbutsample(self):
        # Removing labels and buttons if gene level is selected before #
        try:
            self.selectalgogeneLabel.destroy()
            self.checkbutton_Mzs.destroy()
            self.checkbutton_Abp.destroy()
            self.checkbutton_Gesd.destroy()
            self.applygeneButton.destroy()
            self.checkbutton_MedianRule.destroy()
            self.selectedgeneLabel.destroy()
            self.cancelButgen.destroy()
            self.genethresholdEntry.destroy()
            self.genethresholdlabel.destroy()
            self.helpmenualgorithmsbutton.destroy()
            self.backToStatsButton.destroy()
        except:
            pass
        #Creating buttons and labels for sample level outlier detection #
        self.selectLabel.destroy()
        self.selectedsampleLabel=Label(mGui,text='Detecting Outliers at sample level is selected',font='Times 14 bold italic')
        self.selectedsampleLabel.place(relx=0.22,rely=0.2)
        self.selectalgogeneLabel=Label(mGui,text='Please select algorithm',font='Times 14 bold ')
        self.selectalgogeneLabel.place(relx=0.3,rely=0.4)
        #Adding button for sample level #
        self.ahc=BooleanVar()
        self.checkbutton_ahc=Checkbutton(mGui,text='Average Hierarchical Clustering',variable=self.ahc,command=self.checkbutsamplealgo)
        self.checkbutton_ahc.place(relx=0.2,rely=0.5)

        self.applysampleButton=Button(mGui,text='Apply',command=self.applycheckbutsamplealgo)
        self.applysampleButton.place(relx=0.2,rely=0.7)
        self.cancelButsample=Button(mGui,text='Quit',command=self.mQuit)
        self.cancelButsample.place(relx=0.7,rely=0.7)
        self.goBackStatsButton=Button(mGui, text='Back to Stats',command=self.backtoStatistics)
        self.goBackStatsButton.place(relx=0.4, rely=0.7)

    #Checking if sample level outlier detection is selected #
    def checkbutsamplealgo(self):
        if self.ahc.get():
            #print 'just mzscore'
            pass

    #Applying Average hierarchical algorithm and displaying the results#
    def applycheckbutsamplealgo(self):
        #TCGA dataset#
        if self.selectedRbutton==1:
            SampleLevelTCGA.hierarchical(self.list1,self.full_dataset)
        #GEO dataset#
        else:
            SampleLevelGEO.hierarchical(self.list1,self.full_dataset)
        tkMessageBox.showinfo(title='Attention',message='Image was saved into selected diretory')
        #Removing buttons #
        self.genelLevelButton.destroy()
        self.sampleLevelButton.destroy()
        self.selectedsampleLabel.destroy()
        self.selectalgogeneLabel.destroy()
        self.checkbutton_ahc.destroy()
        self.applysampleButton.destroy()
        self.cancelButsample.destroy()
        self.helpmenualgorithmsbutton.destroy()
        self.goBackStatsButton.destroy()
        #Results #
        self.SilhouetteCofText=Text(mGui,height=1,width=50,wrap=None)
        self.SilhouetteCofText.place(relx=0.2,rely=0.22)
        if self.selectedRbutton==1:
            SilCofText="Silhouette Coefficient : "+str(SampleLevelTCGA.score)
        else:
            SilCofText="Silhouette Coefficient : "+str(SampleLevelGEO.score)
        self.SilhouetteCofText.insert(0.0,SilCofText)
        self.sampleLevelOutlierLabel=Label(mGui,text='Please enter the outlier label',font='Times 14 bold ')
        self.sampleLevelOutlierLabel.place(relx=0.27,rely=0.32)
        self.sampleLevelOutlierEntry=Entry(mGui,width=25)
        self.sampleLevelOutlierEntry.place(relx=0.27,rely=0.4)
        self.sampleNameLabel=Label(mGui,text='If there is more than one outlier,Please seperate them by comma',font='8')
        self.sampleNameLabel.place(relx=0.05,rely=0.5)
        #Creating buttons for futher studies #
        self.GoGenelevelOutlierButton=Button(mGui,text='Go Gene Level',command=self.GoGeneLevel)
        self.GoGenelevelOutlierButton.place(relx=0.15,rely=0.6)
        self.removeSampleOutlierButton=Button(mGui,text='Remove Outlier',command=self.deleteOutlierfromSampleLevel)
        self.removeSampleOutlierButton.place(relx=0.45,rely=0.6)
        self.cancelSampleOutlierButton=Button(mGui,text='Quit',command=self.mQuit)
        self.cancelSampleOutlierButton.place(relx=0.75,rely=0.6)
    
    #deleting outliers#
    def deleteOutlierfromSampleLevel(self):
        self.sampleoutlierstring=self.sampleLevelOutlierEntry.get() # getting outliers from input box        
        #removing labels and buttons #
        self.sampleLevelOutlierLabel.destroy()
        self.sampleLevelOutlierEntry.destroy()
        self.sampleNameLabel.destroy()
        self.GoGenelevelOutlierButton.destroy()
        self.removeSampleOutlierButton.destroy()
        self.cancelSampleOutlierButton.destroy()
        self.SilhouetteCofText.destroy()
        #creating labels for the removing outlier section section #
        self.dataCharacterwoutOutlierSampleLabel=Label(mGui,text='Characteristics of Pure Data Set ',font='Times 18 bold italic')
        self.dataCharacterwoutOutlierSampleLabel.place(relx=0.2,rely=0.2)
        #when gene level outlier algorithms are applied too#
        try:
            seperateSampleOutlier=self.sampleoutlierstring.split(',')
            self.sampleindex=[]
            # reading sample outliers from text widget and preparing for removal #
            for sample in seperateSampleOutlier:
                if self.selectedRbutton==1:
                    if sample in SampleLevelTCGA.newlist:
                        index=SampleLevelTCGA.newlist.index(sample)
                        self.sampleindex.append(index)
                    else:
                        pass
                else:
                    if sample in SampleLevelGEO.newlist:
                        index=SampleLevelGEO.newlist.index(sample)
                        self.sampleindex.append(index)
                    else:
                        pass
                
                reverseData=zip(*self.cancerandnormaldatafull)
                #deleting sample level outliers #
                for sampleindex in sorted(self.sampleindex,reverse=True):
                    del reverseData[int(sampleindex+1)]
    
                self.backcancerandnormaldatafull=zip(*reverseData)
                #deleting gene level outliers #
                cleaneddatagenes=GeneLevel.deleting_outliers(self.backcancerandnormaldatafull,self.RealOutlier)
                self.cleaneddatagenes=cleaneddatagenes
    
                #results#
                self.numberofdeletedSampleText=Text(mGui,height=1,width=60,wrap=None)
                self.numberofdeletedSampleText.place(relx=0.1,rely=0.3)
                numberofdeletedSample='Number of Deleted Sample:            '+str(len(self.sampleindex))
                self.numberofdeletedSampleText.insert(0.0,numberofdeletedSample)
        
                self.numberofremainSampleText=Text(mGui,height=1,width=60,wrap=None)
                self.numberofremainSampleText.place(relx=0.1,rely=0.35)
                numberofremainSample='Number of Sample of Data Set:        '+str(len(reverseData)-1)
                self.numberofremainSampleText.insert(0.0,numberofremainSample)
    
                self.numberofdeletedgenesText=Text(mGui,height=1,width=60,wrap=None)
                self.numberofdeletedgenesText.place(relx=0.1,rely=0.40)
                numberofdeletedgenes='Number of Deleted Outlier Genes:     '+str(len(self.RealOutlier))
                self.numberofdeletedgenesText.insert(0.0,numberofdeletedgenes)
        
                self.numberofcleanedgenesText=Text(mGui,height=1,width=60,wrap=None)
                self.numberofcleanedgenesText.place(relx=0.1,rely=0.45)
                numberofcleangenes='Number of Genes of Data Set:         '+str(len(self.cleaneddatagenes))
                self.numberofcleanedgenesText.insert(0.0,numberofcleangenes)
                #creating buttons for next steps #
                self.extractsampleButton=Button(mGui,text='Extract Files',command=self.browseforsamplelevel)
                self.extractsampleButton.place(relx=0.25,rely=0.5)
    
                self.CancelfromSampleButton=Button(mGui,text='Quit',command=self.mQuit)
                self.CancelfromSampleButton.place(relx=0.55,rely=0.5)
        #when just sample level outlier detection algorithm is applied #
        except:
            seperateSampleOutlier=self.sampleoutlierstring.split(',')
            self.sampleindex=[]
            # reading sample outliers from text widget and preparing for removal #
            for sample in seperateSampleOutlier:
                if self.selectedRbutton==1:
                    if sample in SampleLevelTCGA.newlist:
                        index=SampleLevelTCGA.newlist.index(sample)
                        self.sampleindex.append(index)
                    else:
                        pass
                else:
                    if sample in SampleLevelGEO.newlist:
                        index=SampleLevelGEO.newlist.index(sample)
                        self.sampleindex.append(index)
                    else:
                        pass
            
            reverseData=zip(*self.cancerandnormaldatafull)
            #deleting sample level outliers #
            for sampleindex in sorted(self.sampleindex,reverse=True):
                del reverseData[int(sampleindex+1)]

            self.backcancerandnormaldatafull=zip(*reverseData)
            self.cleaneddatagenes=self.backcancerandnormaldatafull
            #result#
            self.numberofdeletedSampleText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofdeletedSampleText.place(relx=0.1,rely=0.3)
            numberofdeletedSample='Number of Deleted Sample:            '+str(len(self.sampleindex))
            self.numberofdeletedSampleText.insert(0.0,numberofdeletedSample)
    
            self.numberofremainSampleText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofremainSampleText.place(relx=0.1,rely=0.35)
            numberofremainSample='Number of Sample of Data Set:        '+str(len(reverseData)-1)
            self.numberofremainSampleText.insert(0.0,numberofremainSample)
            ##creating buttons for next steps ##
            self.extractsampleButton=Button(mGui,text='Extract Files',command=self.browseforsamplelevel)
            self.extractsampleButton.place(relx=0.25,rely=0.5)

            self.CancelfromSampleButton=Button(mGui,text='Quit',command=self.mQuit)
            self.CancelfromSampleButton.place(relx=0.55,rely=0.5)

    #Browsing path button to extract files#
    def browseforsamplelevel(self):
        self.browseSampleEntry=Entry(mGui,width=26)
        self.browseSampleEntry.place(relx=0.2,rely=0.6)
        self.browseSampleButton=Button(mGui)
        self.browseSampleButton.config(text="Browse",font=10,command=self.getdirectoryforSampleExtract)
        self.browseSampleButton.place(relx=0.7,rely=0.6)

    #Selecting path using pop-up #
    def getdirectoryforSampleExtract(self):        
        mGui.fileNameSampleExtract =  tkFileDialog.askdirectory(parent=mGui, title='Select Directory')
        pathname=mGui.fileNameSampleExtract
        self.browseSampleEntry.delete(0,END)
        self.browseSampleEntry.insert(0,pathname)
        #removing buttons #
        self.extractsampleButton.destroy()
        self.CancelfromSampleButton.destroy()
        #creating buttons for next steps#
        self.ExtractSampleButton=Button(mGui,text="Extract",font=10,command=self.FinalSampleExtract)
        self.ExtractSampleButton.place(relx=0.3,rely=0.68)
        self.CancelExtractSampleButton=Button(mGui,text="Quit",font=10,command=self.mQuit)
        self.CancelExtractSampleButton.place(relx=0.6,rely=0.68)

    #Extracting files#
    def FinalSampleExtract(self):
        new_directorySample=self.browseSampleEntry.get()
        os.chdir(new_directorySample)
        
        for sampleindex in sorted(self.sampleindex,reverse=True):
            del self.data_set[0][int(sampleindex+1)]

        for sampleindex in sorted(self.sampleindex,reverse=True):
            del self.data_set[1][int(sampleindex+1)]        
        #TCGA dataset#
        if self.selectedRbutton==1:
            final_data=[self.data_set[0]] + [self.data_set[1]] + self.cleaneddatagenes
            extractFile.Extracting_files_TCGA(final_data)
        #GEO dataset#
        else:
            final_data=[self.data_set[0]] + self.cleaneddatagenes
            extractFile.Extracting_files_GEO(final_data)
        tkMessageBox.showinfo(title='Info',message='Files were successfully extracted')


    #Back to gene level #
    def GoGeneLevel(self):
        self.sampleoutlierstring=self.sampleLevelOutlierEntry.get()
        
        self.sampleLevelOutlierLabel.destroy()
        self.sampleLevelOutlierEntry.destroy()
        self.sampleNameLabel.destroy()
        self.GoGenelevelOutlierButton.destroy()
        self.removeSampleOutlierButton.destroy()
        self.cancelSampleOutlierButton.destroy()
        self.SilhouetteCofText.destroy()

        self.selectLabel=Label(mGui,text='Please select outlier detection level type',font='Times 14 bold italic')
        self.selectLabel.place(relx=0.22,rely=0.2)
        self.genelLevelButton=Button(mGui,text='Gene Level  ',command=self.checkbutgene)
        self.genelLevelButton.place(relx=0.22,rely=0.27)
        self.sampleLevelButton=Button(mGui,text='Sample Level',command=self.checkbutsample)
        self.sampleLevelButton.place(relx=0.52,rely=0.27)
    #Gene level outlier detection menu #
    def checkbutgene(self):
        #destroying buttons and labels if sample is selected before #
        try:
            self.applysampleButton.destroy()
            self.checkbutton_ahc.destroy()
            self.applysampleButton.destroy()
            self.checkbutton_ahc.destroy()
            self.selectedsampleLabel.destroy()
            self.cancelButsample.destroy()
            self.selectalgogeneLabel.destroy()
            self.goBackStatsButton.destroy()
        except:
            pass
        
        self.selectLabel.destroy()
        self.selectedgeneLabel=Label(mGui,text='Detecting Outliers at gene level is selected',font='Times 14 bold italic')
        self.selectedgeneLabel.place(relx=0.22,rely=0.2)
        self.selectalgogeneLabel=Label(mGui,text='Please select algorithms that you would like to apply',font='Times 14 bold ')
        self.selectalgogeneLabel.place(relx=0.1,rely=0.4)

        #Creating buttons for gene level algorithms #
        self.mzscore=BooleanVar()
        self.checkbutton_Mzs=Checkbutton(mGui,text='Modified Z-Score (MAD)',variable=self.mzscore,command=self.checkbutgenealgo)
        self.checkbutton_Mzs.place(relx=0.2,rely=0.5)

        self.abplot=BooleanVar()
        self.checkbutton_Abp=Checkbutton(mGui,text='Adjusted Box Plot',variable=self.abplot,command=self.checkbutgenealgo)
        self.checkbutton_Abp.place(relx=0.2,rely=0.55)

        self.gesd=BooleanVar()
        self.checkbutton_Gesd=Checkbutton(mGui,text='Generalized Extreme Studentized Deviate ',variable=self.gesd,command=self.checkbutgenealgo)
        self.checkbutton_Gesd.place(relx=0.2,rely=0.6)

        self.medianRule=BooleanVar()
        self.checkbutton_MedianRule=Checkbutton(mGui,text='Median Rule ',variable=self.medianRule,command=self.checkbutgenealgo)
        self.checkbutton_MedianRule.place(relx=0.2,rely=0.65)

        # Threshold for gene level outlier detection 
        self.genethresholdlabel=Label(mGui,text='Outlier observation threshold')
        self.genethresholdlabel.place(relx=0.2,rely=0.75)
        self.genethresholdEntry=Entry(mGui,width=4)
        self.genethresholdEntry.place(relx=0.61,rely=0.75)

        #Buttons for further steps #
        self.applygeneButton=Button(mGui,text='Apply',command=self.applycheckbutgenealgo)
        self.applygeneButton.place(relx=0.2,rely=0.85)

        self.cancelButgen=Button(mGui,text='Quit',command=self.mQuit)
        self.cancelButgen.place(relx=0.7,rely=0.85)

        self.backToStatsButton=Button(mGui,text='Back to Stats',command=self.backtoStatistics)
        self.backToStatsButton.place(relx=0.4, rely=0.85)
    #Back to stats #
    def backtoStatistics(self):
        try:
            #deleting buttons and labels if sample level is selected #
            self.applysampleButton.destroy()             
            self.checkbutton_ahc.destroy()
            self.applysampleButton.destroy()
            self.checkbutton_ahc.destroy()
            self.selectedsampleLabel.destroy()
            self.cancelButsample.destroy()
            self.selectalgogeneLabel.destroy()
            self.goBackStatsButton.destroy()
        except:
            #deleting buttons and labels if gene level is selected #
            self.selectalgogeneLabel.destroy()
            self.checkbutton_Mzs.destroy()
            self.checkbutton_Abp.destroy()
            self.checkbutton_Gesd.destroy()
            self.applygeneButton.destroy()
            self.checkbutton_MedianRule.destroy()
            self.selectedgeneLabel.destroy()
            self.cancelButgen.destroy()
            self.genethresholdEntry.destroy()
            self.genethresholdlabel.destroy()
            self.helpmenualgorithmsbutton.destroy()
            self.backToStatsButton.destroy()
        self.genelLevelButton.destroy()
        self.sampleLevelButton.destroy()
        
        #Reporting the characteristic of data set #
        self.datahead=Label(mGui,text ='Characteristics of Data Set',font='Times 18 bold italic')
        self.datahead.place(relx=0.1,rely=0.2)
        self.geneText=Text(mGui,height=1,width=50,wrap=None)
        self.geneText.place(relx=0.1,rely=0.3)
        gText="Number of Genes :                           "+str(len(self.data_set))
        self.geneText.insert(0.0,gText)
        self.sampleText=Text(mGui,height=1,width=50)
        self.sampleText.place(relx=0.1,rely=0.35)
        sText="Number of Samples:                          "+ str(len(self.data_set[0][1:]))
        self.sampleText.insert(0.0,sText)
        self.cancerText=Text(mGui,height=1,width=50)
        self.cancerText.place(relx=0.1,rely=0.4)
        cText="Number of Cancer Samples :                  "+ str(len(self.data3[1:]))
        self.cancerText.insert(0.0,cText)
        self.normalText=Text(mGui,height=1,width=50)
        self.normalText.place(relx=0.1,rely=0.45)
        nText="Number of Normal Samples:                   "+ str(len(self.data5[1:]))
        self.normalText.insert(0.0,nText)
        self.incompleteText=Text(mGui,height=1,width=50)
        self.incompleteText.place(relx=0.1,rely=0.5)
        iText="Number of incomplete Genes:                 "+ str(len(self.incomplete_genes))
        self.incompleteText.insert(0.0,iText)
        self.distText=Text(mGui,height=1,width=50)
        self.distText.place(relx=0.1,rely=0.55)
        dText="Number of genes follow normal distribution: " +str(len(self.dist))
        self.distText.insert(0.0,dText)

        #Creating the buttons for next steps #
        self.goForAnalysis=Button(mGui,text='Outlier Menu',command=self.GoOutlierMenu)
        self.goForAnalysis.place(relx=0.07,rely=0.7)
        self.CancelButton=Button(mGui,text="Quit",command=self.mQuit)
        self.CancelButton.place(relx=0.72,rely=0.7)
        self.exportincompletegenes=Button(mGui,text='Export Incomplete Genes',command=self.ExportIncmopleteGene)
        self.exportincompletegenes.place(relx=0.32,rely=0.7)
        
    #Go outlier menu after clicking back to stats#
    def GoOutlierMenu(self):
        #Deleting buttons #
        self.datahead.destroy()
        self.geneText.destroy()
        self.sampleText.destroy()
        self.cancerText.destroy()
        self.normalText.destroy()
        self.incompleteText.destroy()
        self.distText.destroy()
        self.goForAnalysis.destroy()
        self.CancelButton.destroy()
        self.exportincompletegenes.destroy()
        #Creating buttons#
        self.selectLabel=Label(mGui,text='Please select outlier detection level type',font='Times 14 bold italic')
        self.selectLabel.place(relx=0.22,rely=0.2)
        self.genelLevelButton=Button(mGui,text='Gene Level  ',command=self.checkbutgene)
        self.genelLevelButton.place(relx=0.22,rely=0.27)
        self.sampleLevelButton=Button(mGui,text='Sample Level',command=self.checkbutsample)
        self.sampleLevelButton.place(relx=0.52,rely=0.27)
        self.helpmenualgorithmsbutton=Button(mGui,text='Help',command=self.helpmenu2)
        self.helpmenualgorithmsbutton.place(relx=0.42,rely=0.93)
    
    #checking which algorithms is selected #
    def checkbutgenealgo(self):
        if self.medianRule.get():
            pass
        
        if self.mzscore.get():
            pass

        if self.abplot.get():
            pass

        if self.gesd.get():
            pass

        if self.mzscore.get() and self.medianRule.get():
            pass
        
        if self.medianRule.get() and self.abplot.get():
            pass
        
        if self.mzscore.get() and self.abplot.get():
            pass
        
        if self.gesd.get() and self.medianRule.get():
            pass
        
        if self.gesd.get() and self.abplot.get():
            pass

        if self.gesd.get() and self.mzscore.get():
            pass

        if self.gesd.get() and self.mzscore.get() and self.medianRule.get():
            pass
        
        if self.medianRule.get() and self.mzscore.get() and self.abplot.get():
            pass
        
        if self.gesd.get() and self.mzscore.get() and self.abplot.get():
            pass
        
        if self.gesd.get() and self.medianRule.get() and self.abplot.get():
            pass

        if self.gesd.get() and self.medianRule.get() and self.abplot.get() and self.medianRule.get():
            pass

    #Applying algortihms for gene level outlier detection #
    def applycheckbutgenealgo(self):
        #Empty list to extract outliers #
        self.extractOut=[]
        
        #Generalized ESD, Adjusted box plot, Modified Z-Score Algorithms and Median Rule #
        if self.gesd.get() and self.mzscore.get() and self.abplot.get() and self.medianRule.get():
            self.abplot_outlier_list=[]
            self.gesd_outlier_list=[]
            self.medianRule_outlier_list=[]
            self.mzs_outlier_list=[]
            self.total_outlier=[]
            for i in range(0,len(self.cancersets)):
                #Normality Calculation#
                pval=normality.distribution(self.cancersets[i],self.thresholdforPvalue)
                if pval>self.thresholdforPvalue:
                    #Modified Z-Score and Generalized ESD  for normal distribution gene  #
                    self.out_mz=GeneLevel.outlier_detection_modified_z_score(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    self.out_gesd=GeneLevel.main_gesd(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    
                    if self.out_mz==None:
                        pass
                    else:
                        self.mzs_outlier_list.append(self.out_mz)
                        temp=['Modified Z-score',self.out_mz[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)

                    if self.out_gesd==None:
                        pass
                    else:
                        self.gesd_outlier_list.append(self.out_gesd)
                        temp=['GESD',self.out_gesd[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                #Adjusted box plot and Median Rule for non-normal distribution gene #
                else:
                    self.out_medianrule=GeneLevel.median_rule_outlier_det(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    self.out_abp=GeneLevel.outlier_detection_adjusted_box_plot(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_abp==None:
                        pass
                    else:
                        self.abplot_outlier_list.append(self.out_abp)
                        temp=['Adjusted BP',self.out_abp[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)

                    if self.out_medianrule==None:
                        pass
                    else:
                        self.medianRule_outlier_list.append(self.out_medianrule)
                        temp=['Median R.',self.out_medianrule[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)

            #common outliers for median rule and adjusted box plot#                        
            self.commonAbp_Median=GeneLevel.common_outlier(self.medianRule_outlier_list,self.abplot_outlier_list,self.cancersets[0])                        
            #common outliers for Generalized ESD and Modified Z-score #                        
            self.common=GeneLevel.common_outlier(self.mzs_outlier_list,self.gesd_outlier_list,self.cancersets[0])
            #merging outliers for Generalized ESD and Modified Z-score #
            self.single=GeneLevel.singleListforMzsGesd(self.mzs_outlier_list,self.gesd_outlier_list)
            #total outlier#
            self.total_outlier_list=self.single + self.commonAbp_Median

            #Removing labels and buttons #
            self.RemoveButtonandLabel()
            #Displaying results#
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)

            self.numberofMzsText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMzsText.place(relx=0.1,rely=0.4)
            numberoutliemMzs="Modified Z-Score            "+str(len(self.mzs_outlier_list))
            self.numberofMzsText.insert(0.0,numberoutliemMzs)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.55)
            numberoutliemabplot="Adjusted Box Plot           "+str(len(self.abplot_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)

            self.numberofcommonText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcommonText.place(relx=0.1,rely=0.5)
            numberoutliemcommon="Common (GESD & M.Z-score)   "+str(len(self.common))
            self.numberofcommonText.insert(0.0,numberoutliemcommon)

            self.numberofcommonAbp_MedianText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcommonAbp_MedianText.place(relx=0.1,rely=0.65)
            numberoutliemcommonAbp_Median="Common (A.BoxPlot & M.Rule) "+str(len(self.commonAbp_Median))
            self.numberofcommonAbp_MedianText.insert(0.0,numberoutliemcommonAbp_Median)
            
            self.numberofGesdText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofGesdText.place(relx=0.1,rely=0.45)
            numberoutliemGesd="Generalized ESD             "+str(len(self.gesd_outlier_list))
            self.numberofGesdText.insert(0.0,numberoutliemGesd)

            self.numberofMedianText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMedianText.place(relx=0.1,rely=0.6)
            numberoutlierMedian="Median Rule                 "+str(len(self.medianRule_outlier_list))
            self.numberofMedianText.insert(0.0,numberoutlierMedian)
            #Creating buttons for next steps#
            self.CreateButtonandLabel()
            self.venndiagramButton=Button(mGui,text='Venn Diagram',command=self.Venn_Diagram)
            self.venndiagramButton.place(relx=0.1,rely=0.8)
            return
        
        #Generalized ESD, Adjusted box plot and Modified Z-Score Algorithms #
        if self.gesd.get() and self.mzscore.get() and self.abplot.get():
            self.abplot_outlier_list=[]
            self.gesd_outlier_list=[]
            self.mzs_outlier_list=[]
            self.total_outlier=[]
            for i in range(0,len(self.cancersets)):
                #Normality Calculation#
                pval=normality.distribution(self.cancersets[i],self.thresholdforPvalue)
                if pval>self.thresholdforPvalue:
                    #Modified Z-Score and Generalized ESD  for normal distribution gene  #
                    self.out_mz=GeneLevel.outlier_detection_modified_z_score(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    self.out_gesd=GeneLevel.main_gesd(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    
                    if self.out_mz==None:
                        pass
                    else:
                        self.mzs_outlier_list.append(self.out_mz)
                        temp=['Modified Z-score',self.out_mz[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)

                    if self.out_gesd==None:
                        pass
                    else:
                        self.gesd_outlier_list.append(self.out_gesd)
                        temp=['GESD',self.out_gesd[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                #Adjusted box plot for non-normal distribution gene #
                else:
                    
                    self.out_abp=GeneLevel.outlier_detection_adjusted_box_plot(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_abp==None:
                        pass
                    else:
                        self.abplot_outlier_list.append(self.out_abp)
                        temp=['Adjusted BP',self.out_abp[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
            #common outliers for Generalized ESD and Modified Z-score #                        
            self.common=GeneLevel.common_outlier(self.mzs_outlier_list,self.gesd_outlier_list,self.cancersets[0])
            #merging outliers for Generalized ESD and Modified Z-score #
            self.single=GeneLevel.singleListforMzsGesd(self.mzs_outlier_list,self.gesd_outlier_list)
            #total outlier#
            self.total_outlier_list=self.single + self.abplot_outlier_list

            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying results#
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)

            self.numberofMzsText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMzsText.place(relx=0.1,rely=0.4)
            numberoutliemMzs="Modified Z-Score            "+str(len(self.mzs_outlier_list))
            self.numberofMzsText.insert(0.0,numberoutliemMzs)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.45)
            numberoutliemabplot="Adjusted Box Plot           "+str(len(self.abplot_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)

            self.numberofcommonText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcommonText.place(relx=0.1,rely=0.5)
            numberoutliemcommon="Common (GESD & M. Z-score)  "+str(len(self.common))
            self.numberofcommonText.insert(0.0,numberoutliemcommon)

            self.numberofGesdText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofGesdText.place(relx=0.1,rely=0.55)
            numberoutliemGesd="Generalized ESD             "+str(len(self.gesd_outlier_list))
            self.numberofGesdText.insert(0.0,numberoutliemGesd)
            #Creating buttons for next steps#
            self.CreateButtonandLabel()
            self.venndiagramButton=Button(mGui,text='Venn Diagram',command=self.Venn_Diagram)
            self.venndiagramButton.place(relx=0.1,rely=0.8)                           
            return
            
        #Median Rule, Adjusted box plot and Modified Z-Score Algorithms #
        if self.medianRule.get() and self.mzscore.get() and self.abplot.get():
            self.abplot_outlier_list=[]
            self.medianRule_outlier_list=[]
            self.mzs_outlier_list=[]
            self.total_outlier=[]
            for i in range(0,len(self.cancersets)):
                #Normality Calculation#
                pval=normality.distribution(self.cancersets[i],self.thresholdforPvalue)
                #Modified Z-Score for normal distribution gene #
                if pval>self.thresholdforPvalue:
                    self.out_mz=GeneLevel.outlier_detection_modified_z_score(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    
                    if self.out_mz==None:
                        pass
                    else:
                        self.mzs_outlier_list.append(self.out_mz)
                        temp=['Modified Z-score',self.out_mz[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                #Median Rule and Adjusted box plot for non-normal distribution gene #
                else:
                    self.out_abp=GeneLevel.outlier_detection_adjusted_box_plot(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    self.out_medianrule=GeneLevel.median_rule_outlier_det(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    
                    if self.out_abp==None:
                        pass
                    else:
                        self.abplot_outlier_list.append(self.out_abp)
                        temp=['Adjusted BP',self.out_abp[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)

                    if self.out_medianrule==None:
                        pass
                    else:
                        self.medianRule_outlier_list.append(self.out_medianrule)
                        temp=['Median R.',self.out_medianrule[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
            #common outliers for median rule and adjusted box plot#                        
            self.common=GeneLevel.common_outlier(self.medianRule_outlier_list,self.abplot_outlier_list,self.cancersets[0])
            #total outliers#
            self.total_outlier_list=self.common + self.mzs_outlier_list

            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying results#
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)

            self.numberofMzsText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMzsText.place(relx=0.1,rely=0.4)
            numberoutliemMzs="Modified Z-Score                 "+str(len(self.mzs_outlier_list))
            self.numberofMzsText.insert(0.0,numberoutliemMzs)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.45)
            numberoutliemabplot="Adjusted Box Plot                "+str(len(self.abplot_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)

            self.numberofcommonText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcommonText.place(relx=0.1,rely=0.55)
            numberoutliemcommon="Common (Median R. & A. Boxplot)  "+str(len(self.common))
            self.numberofcommonText.insert(0.0,numberoutliemcommon)

            self.numberofGesdText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofGesdText.place(relx=0.1,rely=0.5)
            numberoutliemGesd="Median Rule                      "+str(len(self.medianRule_outlier_list))
            self.numberofGesdText.insert(0.0,numberoutliemGesd)
            #Creating buttons for next steps#
            self.CreateButtonandLabel()
            self.venndiagramButton=Button(mGui,text='Venn Diagram',command=self.Venn_Diagram)
            self.venndiagramButton.place(relx=0.1,rely=0.8)

            return

        #Median Rule, Adjusted box plot and Generalized ESD Algorithms#
        if self.medianRule.get() and self.gesd.get() and self.abplot.get():
            self.abplot_outlier_list=[]
            self.medianRule_outlier_list=[]
            self.gesd_outlier_list=[]
            self.total_outlier=[]
            for i in range(0,len(self.cancersets)):
                #Normality Calculation #
                pval=normality.distribution(self.cancersets[i],self.thresholdforPvalue)
                if pval>self.thresholdforPvalue:
                    #Generalized ESD for normal distribution gene #
                    self.out_gesd=GeneLevel.main_gesd(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_gesd==None:
                        pass
                    else:
                        self.gesd_outlier_list.append(self.out_gesd)
                        temp=['GESD',self.out_gesd[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                #Median Rule and Adjusted box plot for non-normal distribution gene #
                else:
                    self.out_abp=GeneLevel.outlier_detection_adjusted_box_plot(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_abp==None:
                        pass
                    else:
                        self.abplot_outlier_list.append(self.out_abp)
                        temp=['Adjusted BP',self.out_abp[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)

                    self.out_medianrule=GeneLevel.median_rule_outlier_det(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())

                    if self.out_medianrule==None:
                        pass
                    else:
                        self.medianRule_outlier_list.append(self.out_medianrule)
                        temp=['Median R.',self.out_medianrule[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                        
            #common outliers for median rule and adjusted box plot #
            self.common=GeneLevel.common_outlier(self.medianRule_outlier_list,self.abplot_outlier_list,self.cancersets[0])
            #total outliers #
            self.total_outlier_list=self.common + self.gesd_outlier_list
            #Removing labels and buttons #
            self.RemoveButtonandLabel()
            #Dispalying Results #
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)

            self.numberofMzsText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMzsText.place(relx=0.1,rely=0.4)
            numberoutliemMzs="Generalized ESD                  "+str(len(self.gesd_outlier_list))
            self.numberofMzsText.insert(0.0,numberoutliemMzs)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.45)
            numberoutliemabplot="Adjusted Box PLot                "+str(len(self.abplot_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)

            self.numberofcommonText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcommonText.place(relx=0.1,rely=0.55)
            numberoutliemcommon="Common (Median R. & A. Boxplot)  "+str(len(self.common))
            self.numberofcommonText.insert(0.0,numberoutliemcommon)

            self.numberofGesdText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofGesdText.place(relx=0.1,rely=0.5)
            numberoutliemGesd="Median Rule                      "+str(len(self.medianRule_outlier_list))
            self.numberofGesdText.insert(0.0,numberoutliemGesd)
            #Creating buttons for next steps#
            self.CreateButtonandLabel()
            self.venndiagramButton=Button(mGui,text='Venn Diagram',command=self.Venn_Diagram)
            self.venndiagramButton.place(relx=0.1,rely=0.8)

            return 

        #Median Rule, Mofified Z-Score and Generalized ESD Algorithms#
        if self.gesd.get() and self.mzscore.get() and self.medianRule.get():
            self.medianRule_outlier_list=[]
            self.gesd_outlier_list=[]
            self.mzs_outlier_list=[]
            self.total_outlier=[]
            for i in range(0,len(self.cancersets)):
                #Normality Calculation#
                pval=normality.distribution(self.cancersets[i],self.thresholdforPvalue)
                if pval>self.thresholdforPvalue:
                    #Generalized ESD and Modified Z-score for normal distribution genes #
                    self.out_mz=GeneLevel.outlier_detection_modified_z_score(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    self.out_gesd=GeneLevel.main_gesd(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    
                    if self.out_mz==None:
                        pass
                    else:
                        self.mzs_outlier_list.append(self.out_mz)
                        temp=['Modified Z-Score',self.out_mz[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)

                    if self.out_gesd==None:
                        pass
                    else:
                        self.gesd_outlier_list.append(self.out_gesd)
                        temp=['GESD',self.out_gesd[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                        
                #Median Rule for non-normal distribution gene #
                else:
                    self.out_medianrule=GeneLevel.median_rule_outlier_det(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_medianrule==None:
                        pass
                    else:
                        self.medianRule_outlier_list.append(self.out_medianrule)
                        temp=['Median R.',self.out_medianrule[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
            #Common outliers #                        
            self.common=GeneLevel.common_outlier(self.mzs_outlier_list,self.gesd_outlier_list,self.cancersets[0])
            #Merging outliers for Generalized ESD and Modified Z-score #
            self.single=GeneLevel.singleListforMzsGesd(self.mzs_outlier_list,self.gesd_outlier_list)
            #Total outliers #
            self.total_outlier_list=self.single + self.medianRule_outlier_list

            #Removing labels and buttons #
            self.RemoveButtonandLabel()
            #Dispalying Results #
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)

            self.numberofMzsText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMzsText.place(relx=0.1,rely=0.4)
            numberoutliemMzs="Modified Z-Score      "+str(len(self.mzs_outlier_list))
            self.numberofMzsText.insert(0.0,numberoutliemMzs)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.45)
            numberoutliemabplot="Median Rule           "+str(len(self.medianRule_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)

            self.numberofcommonText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcommonText.place(relx=0.1,rely=0.5)
            numberoutliemcommon="Common (GESD & M.Zs)  "+str(len(self.common))
            self.numberofcommonText.insert(0.0,numberoutliemcommon)

            self.numberofGesdText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofGesdText.place(relx=0.1,rely=0.55)
            numberoutliemGesd="Generalized ESD       "+str(len(self.gesd_outlier_list))
            self.numberofGesdText.insert(0.0,numberoutliemGesd)
            
            #Creating buttons for next steps#
            self.CreateButtonandLabel()
            self.venndiagramButton=Button(mGui,text='Venn Diagram',command=self.self.Venn_Diagram)
            self.venndiagramButton.place(relx=0.1,rely=0.8)

            return
            

        #Mofified Z-Score and Generalized ESD Algorithms#
        if self.gesd.get() and self.mzscore.get():
            self.gesd_outlier_list=[]
            self.mzs_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                #Modified Z-Score#
                self.out_mz=GeneLevel.outlier_detection_modified_z_score(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                #Generalized ESD #
                self.out_gesd=GeneLevel.main_gesd(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    
                if self.out_mz==None:
                    pass
                else:
                    self.mzs_outlier_list.append(self.out_mz)
                    temp=['Modified Z-score',self.out_mz[0],i] # temprorary list to add main extraction list for outliers
                    self.extractOut.append(temp)

                if self.out_gesd==None:
                    pass
                else:
                    self.gesd_outlier_list.append(self.out_gesd)
                    temp=['GESD',self.out_gesd[0],i] # temprorary list to add main extraction list for outliers
                    self.extractOut.append(temp)
            #Common Outliers #
            self.common=GeneLevel.common_outlier(self.mzs_outlier_list,self.gesd_outlier_list,self.cancersets[0])            
            #Merging outliers#
            self.single=GeneLevel.singleListforMzsGesd(self.mzs_outlier_list,self.gesd_outlier_list)
            
            self.total_outlier_list=self.single
            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying Results#
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)

            self.numberofMzsText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMzsText.place(relx=0.1,rely=0.4)
            numberoutliemMzs="Modified Z-Score            "+str(len(self.mzs_outlier_list))
            self.numberofMzsText.insert(0.0,numberoutliemMzs)

            self.numberofGesdText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofGesdText.place(relx=0.1,rely=0.5)
            numberoutliemGesd="Generalized ESD             "+str(len(self.gesd_outlier_list))
            self.numberofGesdText.insert(0.0,numberoutliemGesd)

            self.numberofcommonText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcommonText.place(relx=0.1,rely=0.6)
            numberoutliemcommon="Common Outliers             "+str(len(self.common))
            self.numberofcommonText.insert(0.0,numberoutliemcommon)

            #Creating buttons for further steps#            
            self.CreateButtonandLabel()
            self.venndiagramButton=Button(mGui,text='Venn Diagram',command=self.Venn_Diagram)
            self.venndiagramButton.place(relx=0.1,rely=0.8)
            return
            
        #Adjusted box plot and Generalized ESD Algorithms#
        if self.abplot.get() and self.gesd.get():
            self.abplot_outlier_list=[]
            self.gesd_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                #normality calculation #
                pval=normality.distribution(self.cancersets[i],self.thresholdforPvalue)
                #GESD for normal gene #
                if pval>self.thresholdforPvalue:
                    self.out_gesd=GeneLevel.main_gesd(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_gesd==None:
                        pass
                    else:
                        self.gesd_outlier_list.append(self.out_gesd)
                        temp=['GESD',self.out_gesd[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                #Adjusted Box plot for normal gene #
                else:
                    self.out_abp=GeneLevel.outlier_detection_adjusted_box_plot(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_abp==None:
                        pass
                    else:
                        self.abplot_outlier_list.append(self.out_abp)
                        temp=['Adjusted BP',self.out_abp[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
            self.total_outlier_list=self.abplot_outlier_list + self.gesd_outlier_list

            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying Results#
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene L',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.45)
            numberoutliemabplot="Adjusted Box Plot          "+str(len(self.abplot_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)
            
            self.numberofGesdText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofGesdText.place(relx=0.1,rely=0.4)
            numberoutliemGesd="Generalized ESD            "+str(len(self.gesd_outlier_list))
            self.numberofGesdText.insert(0.0,numberoutliemGesd)

            #Creating buttons for further steps#            
            self.CreateButtonandLabel()

            return            
            

        #Median Rule and Generalized ESD Algorithms #
        if self.medianRule.get() and self.gesd.get():
            self.medianRule_outlier_list=[]
            self.gesd_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                #normality calculation #
                pval=normality.distribution(self.cancersets[i],self.thresholdforPvalue)
                #GESD for normal gene #
                if pval>self.thresholdforPvalue:
                    self.out_gesd=GeneLevel.main_gesd(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_gesd==None:
                        pass
                    else:
                        self.gesd_outlier_list.append(self.out_gesd)
                        temp=['GESD',self.out_gesd[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                #Median Rule for normal gene #
                else:
                    self.out_medianrule=GeneLevel.median_rule_outlier_det(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_medianrule==None:
                        pass
                    else:
                        self.medianRule_outlier_list.append(self.out_medianrule)
                        temp=['Median R.',self.out_medianrule[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
            self.total_outlier_list=self.medianRule_outlier_list + self.gesd_outlier_list
            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying Results#
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.45)
            numberoutliemabplot="Median Rule               "+str(len(self.medianRule_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)
            
            self.numberofGesdText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofGesdText.place(relx=0.1,rely=0.4)
            numberoutliemGesd="Generalized ESD           "+str(len(self.gesd_outlier_list))
            self.numberofGesdText.insert(0.0,numberoutliemGesd)

            #Creating buttons for further steps#            
            self.CreateButtonandLabel()

            return
            
        #Adjusted box plot  and Modified Z-Score algorithm #
        if self.mzscore.get() and self.abplot.get():
            self.mzs_outlier_list=[]
            self.abplot_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                #normality calculation #
                pval=normality.distribution(self.cancersets[i],self.thresholdforPvalue)
                #Modified Z-score for normal gene #
                if pval>self.thresholdforPvalue:
                    self.out_mz=GeneLevel.outlier_detection_modified_z_score(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    
                    if self.out_mz==None:
                        pass
                    else:
                        self.mzs_outlier_list.append(self.out_mz)
                        temp=['Modified Z.Score',self.out_mz[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                #Adjusted box plot for non-normal gene #
                else:
                    self.out_abp=GeneLevel.outlier_detection_adjusted_box_plot(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_abp==None:
                        pass
                    else:
                        self.abplot_outlier_list.append(self.out_abp)
                        temp=['Adjusted BP',self.out_abp[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
            self.total_outlier_list=self.mzs_outlier_list + self.abplot_outlier_list

            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying Results#
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)
            
            self.numberofMzsText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMzsText.place(relx=0.1,rely=0.4)
            numberoutliemMzs="Modified Z-Score              "+str(len(self.mzs_outlier_list))
            self.numberofMzsText.insert(0.0,numberoutliemMzs)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.5)
            numberoutliemabplot="Adjusted Box PLot             "+str(len(self.abplot_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)

            #Creating buttons for further steps#            
            self.CreateButtonandLabel()

            return
            

        #Median Rule and Modified Z-Score algorithm #
        if self.mzscore.get() and self.medianRule.get():
            self.mzs_outlier_list=[]
            self.medianRule_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                pval=normality.distribution(self.cancersets[i],self.thresholdforPvalue)
                #Modified Z-score for normal gene #
                if pval>self.thresholdforPvalue: 
                    self.out_mz=GeneLevel.outlier_detection_modified_z_score(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_mz==None:
                        pass
                    else:
                        self.mzs_outlier_list.append(self.out_mz)
                        temp=['Modified Z.Score',self.out_mz[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
                #Median Rule for non-normal gene #
                else:
                    self.out_medianrule=GeneLevel.median_rule_outlier_det(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                    if self.out_medianrule==None:
                        pass
                    else:
                        self.medianRule_outlier_list.append(self.out_medianrule)
                        temp=['Median R.',self.out_medianrule[0],i] # temprorary list to add main extraction list for outliers
                        self.extractOut.append(temp)
            self.total_outlier_list=self.mzs_outlier_list + self.medianRule_outlier_list

            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying Results#
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)
            
            self.numberofMzsText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMzsText.place(relx=0.1,rely=0.4)
            numberoutliemMzs="Modified Z-Score            "+str(len(self.mzs_outlier_list))
            self.numberofMzsText.insert(0.0,numberoutliemMzs)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.5)
            numberoutliemabplot="Median Rule                 "+str(len(self.medianRule_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)

            #Creating buttons for further steps#            
            self.CreateButtonandLabel()

            return

        #Median Rule and Adjusted box plot algorithm #
        if self.medianRule.get() and self.abplot.get():
            self.medianRule_outlier_list=[]
            self.abplot_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                self.out_medianrule=GeneLevel.median_rule_outlier_det(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                self.out_abp=GeneLevel.outlier_detection_adjusted_box_plot(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                if self.out_medianrule==None:
                    pass
                else:
                    self.medianRule_outlier_list.append(self.out_medianrule)
                    temp=['Median R.',self.out_medianrule[0],i] # temprorary list to add main extraction list for outliers
                    self.extractOut.append(temp)
                
                if self.out_abp==None:
                    pass
                else:
                    self.abplot_outlier_list.append(self.out_abp)
                    temp=['Adjusted BP',self.out_abp[0],i] # temprorary list to add main extraction list for outliers
                    self.extractOut.append(temp)
                    
            self.common=GeneLevel.common_outlier(self.abplot_outlier_list,self.medianRule_outlier_list,self.cancersets[0])
            self.total_outlier_list=self.common

            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying Results#
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.2)

            self.appliedalgoLabel=Label(mGui,text='Applied Algorithm',font='bold')
            self.appliedalgoLabel.place(relx=0.1,rely=0.3)

            self.numberofoutlierLabel=Label(mGui,text='Number of Potential Outliers',font='bold')
            self.numberofoutlierLabel.place(relx=0.45,rely=0.3)
            
            self.numberofMzsText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofMzsText.place(relx=0.1,rely=0.4)
            numberoutliemMzs="Median Rule                  "+str(len(self.medianRule_outlier_list))
            self.numberofMzsText.insert(0.0,numberoutliemMzs)

            self.numberofabplotText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofabplotText.place(relx=0.1,rely=0.5)
            numberoutliemabplot="Adjusted Box PLot            "+str(len(self.abplot_outlier_list))
            self.numberofabplotText.insert(0.0,numberoutliemabplot)

            self.numberofcommonText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcommonText.place(relx=0.1,rely=0.6)
            numberoutliemcommon="Common Outliers              "+str(len(self.common))
            self.numberofcommonText.insert(0.0,numberoutliemcommon)

            #Creating buttons for further steps#            
            self.CreateButtonandLabel()
            self.venndiagramButton=Button(mGui,text='Venn Diagram',command=self.Venn_Diagram)
            self.venndiagramButton.place(relx=0.1,rely=0.8)
            return

        #Modified Z-Score algorithm#
        if self.mzscore.get():
            self.mzs_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                self.out_mz=GeneLevel.outlier_detection_modified_z_score(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                if self.out_mz==None:
                    pass
                else:
                    self.mzs_outlier_list.append(self.out_mz)
                    temp=['Modified Z-score',self.out_mz[0],i] # temprorary list to add main extraction list for outliers
                    self.extractOut.append(temp)
            self.total_outlier_list=self.mzs_outlier_list

            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying results #
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.25)

            self.appliedalgoText=Text(mGui,height=1,width=50,wrap=None)
            self.appliedalgoText.place(relx=0.1,rely=0.4)
            appliedText="Applied Algorithm  : Modified Z-Score"
            self.appliedalgoText.insert(0.0,appliedText)
            
            self.numberofoutlierText=Text(mGui,height=1,width=50,wrap=None)
            self.numberofoutlierText.place(relx=0.1,rely=0.5)
            numberoutlierText="Number of Outliers : "+str(len(self.mzs_outlier_list))
            self.numberofoutlierText.insert(0.0,numberoutlierText)

            #Creating buttons for further steps#            
            self.CreateButtonandLabel()
            
            return self.mzs_outlier_list
        
        #Adjusted Box plot algorithm #
        if self.abplot.get():
            self.abplot_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                self.out_abp=GeneLevel.outlier_detection_adjusted_box_plot(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                if self.out_abp==None:
                    pass
                else:
                    self.abplot_outlier_list.append(self.out_abp)
                    temp=['Adjusted BP',self.out_abp[0],i] # temprorary list to add main extraction list for outliers
                    self.extractOut.append(temp)

            self.total_outlier_list=self.abplot_outlier_list

            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying results #
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.25)
            
            self.appliedalgoText=Text(mGui,height=1,width=50,wrap=None)
            self.appliedalgoText.place(relx=0.1,rely=0.4)
            appliedText="Applied Algorithm  : Adjusted Box Plot "
            self.appliedalgoText.insert(0.0,appliedText)

            self.numberofoutlierText=Text(mGui,height=1,width=50,wrap=None)
            self.numberofoutlierText.place(relx=0.1,rely=0.5)
            numberoutlierText="Number of Outliers : "+str(len(self.abplot_outlier_list))
            self.numberofoutlierText.insert(0.0,numberoutlierText)

            #Creating buttons for further steps#            
            self.CreateButtonandLabel()
            
            return self.abplot_outlier_list

        #Generalized ESD Algorithm #
        if self.gesd.get():
            self.gesd_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                self.out_gesd=GeneLevel.main_gesd(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                if self.out_gesd==None:
                    pass
                else:
                    self.gesd_outlier_list.append(self.out_gesd)
                    temp=['GESD',self.out_gesd[0],i] # temprorary list to add main extraction list for outliers
                    self.extractOut.append(temp)

            self.total_outlier_list=self.gesd_outlier_list

            #Removing labels and buttons #
            self.RemoveButtonandLabel()

            #Displaying results #
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.25)
            
            self.appliedalgoText=Text(mGui,height=1,width=50,wrap=None)
            self.appliedalgoText.place(relx=0.1,rely=0.4)
            appliedText="Applied Algorithm: Generalized Extreme Studentized Deviate "
            self.appliedalgoText.insert(0.0,appliedText)

            self.numberofoutlierText=Text(mGui,height=1,width=50,wrap=None)
            self.numberofoutlierText.place(relx=0.1,rely=0.5)
            numberoutlierText="Number of Outliers : "+str(len(self.gesd_outlier_list))
            self.numberofoutlierText.insert(0.0,numberoutlierText)

            #Creating buttons for further steps#            
            self.CreateButtonandLabel()

            return self.gesd_outlier_list
        
        # Median Rule Algorthim #
        if self.medianRule.get():
            self.medianRule_outlier_list=[]
            for i in range(0,len(self.cancersets)):
                self.out_medianrule=GeneLevel.median_rule_outlier_det(self.cancersets[i],thresholdGeneLevel=self.genethresholdEntry.get())
                if self.out_medianrule==None:
                    pass
                else:
                    self.medianRule_outlier_list.append(self.out_medianrule)
                    temp=['Median R',self.out_medianrule[0],i] # temprorary list to add main extraction list for outliers
                    self.extractOut.append(temp)
                    
            self.total_outlier_list=self.medianRule_outlier_list 

            #Removing labels and buttons #
            self.RemoveButtonandLabel()
            
            #Displaying results #
            self.resultgenelevelLabel=Label(mGui,text='Detected outliers at Gene Level',font='Times 18 bold')
            self.resultgenelevelLabel.place(relx=0.1,rely=0.25)
            
            self.appliedalgoText=Text(mGui,height=1,width=50,wrap=None)
            self.appliedalgoText.place(relx=0.1,rely=0.4)
            appliedText="Applied Algorithm : Median Rule "
            self.appliedalgoText.insert(0.0,appliedText)

            self.numberofoutlierText=Text(mGui,height=1,width=50,wrap=None)
            self.numberofoutlierText.place(relx=0.1,rely=0.5)
            numberoutlierText="Number of Outliers: "+str(len(self.medianRule_outlier_list))
            self.numberofoutlierText.insert(0.0,numberoutlierText)

            #Creating buttons for further steps#
            self.CreateButtonandLabel()
            
            return self.total_outlier_list

    # Help menu for outlier detection algorithms #
    def helpmenu2(self):
        tkMessageBox.showwarning(title='Help',message=hlpMENU.helpscripts())
        return
    
    # Removing buttons and labels for gene level outlier detection #
    def RemoveButtonandLabel(self):
        #Removing buttons and labels #
        self.sampleLevelButton.destroy()
        self.genelLevelButton.destroy()
        self.selectedgeneLabel.destroy()
        self.selectalgogeneLabel.destroy()
        self.checkbutton_Mzs.destroy()
        self.checkbutton_Abp.destroy()
        self.checkbutton_Gesd.destroy()
        self.checkbutton_MedianRule.destroy()
        self.applygeneButton.destroy()
        self.cancelButgen.destroy()
        self.genethresholdEntry.destroy()
        self.genethresholdlabel.destroy()
#        self.genethrsholdPercentagelabel.destroy()
        self.helpmenualgorithmsbutton.destroy()
        self.backToStatsButton.destroy()

    #Creating buttons for next steps after gene level outlier detection #        
    def CreateButtonandLabel(self):
        #Creating buttons for further steps#
        self.groupOutlierButton=Button(mGui,text='Group Outlier',command=self.groupOutlierMethod)
        self.groupOutlierButton.place(relx=0.1,rely=0.7)
            
        self.backmainOutlierMenuButton=Button(mGui,text='Back',command=self.backMainOutlierMenu)
        self.backmainOutlierMenuButton.place(relx=0.37,rely=0.7)
            
        self.piechartMzsButton=Button(mGui,text='Pie Chart',command=self.PieChart)
        self.piechartMzsButton.place(relx=0.51,rely=0.7)

        self.cancelgeneoutlierButton=Button(mGui,text='Quit',command=self.mQuit)
        self.cancelgeneoutlierButton.place(relx=0.72,rely=0.7)

        self.ExtractOutliersButton=Button (mGui,text='Extract Outliers',command=self.EXtractOUT)
        self.ExtractOutliersButton.place(relx=0.37,rely=0.8)

    #Extracting Outlier Files#
    def EXtractOUT(self):
        extractOutliers.ExtractOutlierGene(self.extractOut)
        tkMessageBox.showinfo(title='Attention',message='File was saved as Outliers.csv into same same directory where expression files exist')

    #Plotting Pie Chart  #
    def PieChart(self):
        if self.abplot.get() and self.medianRule.get() and self.gesd.get() and self.mzscore.get():
            piechart.pie_chart_medianRule_gesd_abplot_mzs(self.gesd_outlier_list,self.medianRule_outlier_list,self.common,self.abplot_outlier_list,self.cancersets,self.mzs_outlier_list)
            return
        
        if self.abplot.get() and self.medianRule.get() and self.gesd.get():
            piechart.pie_chart_medianRule_gesd_abplot(self.gesd_outlier_list,self.medianRule_outlier_list,self.common,self.abplot_outlier_list,self.cancersets)
            return
        
        if self.mzscore.get() and self.medianRule.get() and self.abplot.get():
            piechart.pie_chart_medianRule_mzs_abplot(self.mzs_outlier_list,self.medianRule_outlier_list,self.common,self.abplot_outlier_list,self.cancersets)
            return
        
        if self.mzscore.get() and self.medianRule.get() and self.gesd.get():
            piechart.pie_chart_medianRule_mzs_gesd(self.mzs_outlier_list,self.medianRule_outlier_list,self.common,self.gesd_outlier_list,self.cancersets)
            return
        
        if self.mzscore.get() and self.abplot.get() and self.gesd.get():
            piechart.pie_chart_abplot_mzs_gesd(self.mzs_outlier_list,self.abplot_outlier_list,self.common,self.gesd_outlier_list,self.cancersets)
            return

        if self.mzscore.get() and self.gesd.get():
            piechart.pie_chart_mzs_gesd(self.mzs_outlier_list,self.gesd_outlier_list,self.common,self.cancersets)
            return

        if self.medianRule.get() and self.gesd.get():
            piechart.pie_chart_medianRule_gesd(self.medianRule_outlier_list,self.gesd_outlier_list,self.cancersets)
            return
        
        if self.abplot.get() and self.gesd.get():
            piechart.pie_chart_abplot_gesd(self.abplot_outlier_list,self.gesd_outlier_list,self.cancersets)
            return
        
        if self.medianRule.get() and self.abplot.get():
            piechart.pie_chart_abplot_medianrule(self.medianRule_outlier_list,self.abplot_outlier_list, self.common,self.cancersets)
            return

        if self.mzscore.get() and self.medianRule.get():
            piechart.pie_chart_medianRule_mzs(self.mzs_outlier_list,self.medianRule_outlier_list,self.cancersets)
            return
        
        if self.mzscore.get() and self.abplot.get():
            piechart.pie_chart_abplot_mzs(self.mzs_outlier_list,self.abplot_outlier_list,self.cancersets)
            return
        
        if self.abplot.get():
            piechart.pie_chart_abplot(self.abplot_outlier_list,self.cancersets)
            return
        
        if self.medianRule.get():
            piechart.pie_chart_medianRule(self.medianRule_outlier_list,self.cancersets)
            return
        
        if self.mzscore.get():
            piechart.pie_chart_mzs(self.mzs_outlier_list,self.cancersets)
            return
        
        if self.gesd.get():
            piechart.pie_chart_gesd(self.gesd_outlier_list,self.cancersets)
            return
        
            
    #Plotting Venn diagram#
    def Venn_Diagram(self):
        if self.medianRule.get() and self.abplot.get():
            venndiagram.venn_diagram_abplt_medianrule(self.medianRule_outlier_list,self.abplot_outlier_list,self.common)
            return
        if self.mzscore.get() and self.gesd.get():
            venndiagram.venn_diagram_gesd_mzs(self.mzs_outlier_list,self.gesd_outlier_list,self.common)
            return
        if self.mzscore.get() and self.abplot.get() and self.gesd.get():
            venndiagram.venn_diagram_gesd_mzs(self.mzs_outlier_list,self.gesd_outlier_list,self.common)
            return
        if self.mzscore.get() and self.medianRule.get() and self.gesd.get():
            venndiagram.venn_diagram_gesd_mzs(self.mzs_outlier_list,self.gesd_outlier_list,self.common)
            return
        if self.mzscore.get() and self.medianRule.get() and self.abplot.get():
            venndiagram.venn_diagram_abplt_medianrule(self.medianRule_outlier_list,self.abplot_outlier_list,self.common)
            return
        if self.abplot.get() and self.medianRule.get() and self.gesd.get():
            venndiagram.venn_diagram_abplt_medianrule(self.medianRule_outlier_list,self.abplot_outlier_list,self.common)
            return
        if self.abplot.get() and self.medianRule.get() and self.gesd.get() and self.mzscore.get():
            venndiagram.venn_diagram_abplt_medianrule(self.medianRule_outlier_list,self.abplot_outlier_list,self.common)
            return
        
        
    #Grouping outlier main menu buttons and labels #
    def groupOutlierMethod(self):
        #Removing labels and buttons #
        self.resultgenelevelLabel.destroy()
        self.backmainOutlierMenuButton.destroy()
        self.cancelgeneoutlierButton.destroy()
        self.piechartMzsButton.destroy()
        self.groupOutlierButton.destroy()
        self.cancelgeneoutlierButton.destroy()
        self.ExtractOutliersButton.destroy()
        try:
            self.numberofoutlierText.destroy()
        except:
            pass
        try:
            self.numberofoutlierLabel.destroy()
        except:
            pass
        try:
            self.appliedalgoText.destroy()
        except:
            pass
        try:
            self.appliedalgoLabel.destroy()
        except:
            pass
        try:
            self.venndiagramButton.destroy()
        except:
            pass
        try:
            self.numberofMzsText.destroy()
        except:
            pass
        try:
            self.numberofcommonText.destroy()
        except:
            pass
        try:
            self.numberofGesdText.destroy()
        except:
            pass
        try:
            self.numberofabplotText.destroy()
        except:
            pass
        try:
            self.numberofcommonAbp_MedianText.destroy()
            self.numberofMedianText.destroy()
        except:
            pass
        #Creating labels and buttons for grouping outlier #
        self.groupOutliergeneLabel=Label(mGui,text='Finding Similarity of Outlier Genes',font='Times 18 bold')
        self.groupOutliergeneLabel.place(relx=0.2,rely=0.2)

        self.selectmethodgeneLabel=Label(mGui,text='Please select methods',font='Times 14 bold ')
        self.selectmethodgeneLabel.place(relx=0.35,rely=0.27)

        self.funsim=BooleanVar()
        self.checkbutton_funsim=Checkbutton(mGui,text='Functional Similarity',variable=self.funsim,command=self.checkbutsimMethods)
        self.checkbutton_funsim.place(relx=0.2,rely=0.35)

        self.keggPath=BooleanVar()
        self.checkbutton_keggPath=Checkbutton(mGui,text='Kegg Pathway Participation',variable=self.keggPath,command=self.checkbutsimMethods)
        self.checkbutton_keggPath.place(relx=0.2,rely=0.48)

        self.coexpress=BooleanVar()
        self.checkbutton_coexpress=Checkbutton(mGui,text='Co-expression (Pearson Correlation) ',variable=self.coexpress,command=self.checkbutsimMethods)
        self.checkbutton_coexpress.place(relx=0.2,rely=0.70)
        
        self.settingsimButton=Button(mGui,text='Settings',command=self.similaritysetting)
        self.settingsimButton.place(relx=0.25,rely=0.85)

        self.CancelsimButton=Button(mGui,text='Quit',command=self.mQuit)
        self.CancelsimButton.place(relx=0.6,rely=0.85)

        self.helpmenumethodbutton=Button(mGui,text='Help',command=self.helpmenu3)
        self.helpmenumethodbutton.place(relx=0.8,rely=0.90)

    #Going back to main outlier menu from gene level outlier result#
    def backMainOutlierMenu(self):
        self.resultgenelevelLabel.destroy()
        self.ExtractOutliersButton.destroy()
        self.groupOutlierButton.destroy()
        self.backmainOutlierMenuButton.destroy()
        self.cancelgeneoutlierButton.destroy()
        self.piechartMzsButton.destroy()
        try:
            self.appliedalgoText.destroy()
        except:
            pass
        try:
            self.numberofoutlierText.destroy()
        except:
            pass
        try:
            self.numberofMzsText.destroy()
        except:
            pass
        
        try:
            self.numberofabplotText.destroy()
        except:
            pass
        try:
            self.numberofcommonText.destroy()
        except:
            pass
        try:
            self.venndiagramButton.destroy()
        except:
            pass
        try:
            self.appliedalgoLabel.destroy()
        except:
            pass
        try:
            self.numberofoutlierLabel.destroy()
        except:
            pass
        try:
            self.numberofGesdText.destroy()
        except:
            pass
        self.selectLabel=Label(mGui,text='Please select outlier detection level type',font='Times 14 bold italic')
        self.selectLabel.place(relx=0.22,rely=0.2)
        self.genelLevelButton=Button(mGui,text='Gene Level  ',command=self.checkbutgene)
        self.genelLevelButton.place(relx=0.22,rely=0.27)
        self.sampleLevelButton=Button(mGui,text='Sample Level',command=self.checkbutsample)
        self.sampleLevelButton.place(relx=0.52,rely=0.27)        

    #help menu for grouping outliers#
    def helpmenu3(self):
        tkMessageBox.showwarning(title='Help',message=helpMG.helpscripts())

    #Checking which grouping outlier methods are selected #
    def checkbutsimMethods(self):
        if self.coexpress.get() and self.funsim.get() and self.keggPath.get():
            return
        
        if self.coexpress.get() and self.funsim.get():
            return

        if self.funsim.get() and self.keggPath.get():
            return

        if self.coexpress.get() and self.keggPath.get():
            return
        
        if self.funsim.get():
            return
        
        if self.keggPath.get():
            return

        if self.coexpress.get():
            return

    #FunSimMat method threshold value creating#
    def FunSimMatThreshold(self):
        self.funsimthresholdlabel=Label(mGui,text='Please enter threshold value')
        self.funsimthresholdlabel.place(relx=0.2,rely=0.4)
        self.funsimthresholdEntry=Entry(mGui,width=6)
        self.funsimthresholdEntry.place(relx=0.61,rely=0.4)

    #Coexpression method threshold value creating#
    def CoExpressionThreshold(self):
        self.coexpressthresholdlabel=Label(mGui,text='Please enter threshold value')
        self.coexpressthresholdlabel.place(relx=0.2,rely=0.75)
        self.coexpressthresholdEntry=Entry(mGui,width=6)
        self.coexpressthresholdEntry.place(relx=0.61,rely=0.75)

    #KEGG Pathway Participation method reference file selecting button #
    def KEGGPathwayButton(self):
        self.keggPathButton=Button(mGui,text='Select File',command=self.selectFileKegg)
        self.keggPathButton.place(relx=0.77,rely=0.53)
        self.keggPathEntry=Entry(mGui,width=30)
        self.keggPathEntry.place(relx=0.2,rely=0.53)
        self.KEGGthresholdlabel=Label(mGui,text='Please enter threshold value')
        self.KEGGthresholdlabel.place(relx=0.2,rely=0.63)
        self.KEGGthresholdEntry=Entry(mGui,width=6)
        self.KEGGthresholdEntry.place(relx=0.61,rely=0.63)
        thresholdValue='20'
        self.KEGGthresholdEntry.insert(0,thresholdValue)
        
    #Settings button after grouping methods have been selected #
    def similaritysetting(self):
        if self.coexpress.get() and self.funsim.get() and self.keggPath.get():
            self.FunSimMatThreshold() #Funsimmat method settings
            self.KEGGPathwayButton() #KEGG pathway participation method sttings
            self.CoExpressionThreshold() #Co-expression method setting
            self.similaritymethodsApplyButton=Button(mGui,text='Apply',command=self.similarityApply)
            self.similaritymethodsApplyButton.place(relx=0.44,rely=0.9)
            return
        
        if self.coexpress.get() and self.funsim.get():
            self.FunSimMatThreshold() #Funsimmat method settings
            self.CoExpressionThreshold() #Co-expression method setting
            self.similaritymethodsApplyButton=Button(mGui,text='Apply',command=self.similarityApply)
            self.similaritymethodsApplyButton.place(relx=0.44,rely=0.9)
            return

        if self.funsim.get() and self.keggPath.get():
            self.FunSimMatThreshold() #Funsimmat method settings
            self.KEGGPathwayButton() #KEGG pathway participation method sttings
            self.similaritymethodsApplyButton=Button(mGui,text='Apply',command=self.similarityApply)
            self.similaritymethodsApplyButton.place(relx=0.44,rely=0.9)
            return

        if self.coexpress.get() and self.keggPath.get():
            self.CoExpressionThreshold() #Co-expression method setting
            self.KEGGPathwayButton() #KEGG pathway participation method sttings
            self.similaritymethodsApplyButton=Button(mGui,text='Apply',command=self.similarityApply)
            self.similaritymethodsApplyButton.place(relx=0.44,rely=0.9)
            return
        
        if self.funsim.get():
            self.FunSimMatThreshold() #Funsimmat method settings
            self.similaritymethodsApplyButton=Button(mGui,text='Apply',command=self.similarityApply)
            self.similaritymethodsApplyButton.place(relx=0.44,rely=0.9)
            return
        
        if self.keggPath.get():
            self.KEGGPathwayButton() #KEGG pathway participation method sttings
            self.similaritymethodsApplyButton=Button(mGui,text='Apply',command=self.similarityApply)
            self.similaritymethodsApplyButton.place(relx=0.44,rely=0.9)
            return

        if self.coexpress.get():
            self.CoExpressionThreshold() #Co-expression method setting
            self.similaritymethodsApplyButton=Button(mGui,text='Apply',command=self.similarityApply)
            self.similaritymethodsApplyButton.place(relx=0.44,rely=0.9)
            return

    #Selecting KEGG file as a reference#
    def selectFileKegg(self):
        self.keggfilename=askopenfilename(parent=mGui)
        self.keggPathEntry.delete(0,END)
        self.keggPathEntry.insert(0,self.keggfilename)

    #Removing buttons and labels from grouping outlier section#
    def removeAndAddButtonAndLabel(self):
        #remove buttons and labels #
        self.selectmethodgeneLabel.destroy()
        self.checkbutton_funsim.destroy()
        self.checkbutton_keggPath.destroy()
        self.checkbutton_coexpress.destroy()
        self.settingsimButton.destroy()
        self.CancelsimButton.destroy()
        self.similaritymethodsApplyButton.destroy()
        try:
            self.coexpressthresholdlabel.destroy()
            self.coexpressthresholdEntry.destroy()
        except:
            pass
        try:
            self.keggPathButton.destroy()
            self.keggPathEntry.destroy()
            self.keggPathEntry.destroy()
            self.KEGGthresholdlabel.destroy()
            self.KEGGthresholdEntry.destroy()
        except:
            pass
        try:
            self.funsimthresholdlabel.destroy()
            self.funsimthresholdEntry.destroy()
        except:
            pass
        #adding labels that show the title of result page#
        self.similarityAppliedMethodLabel=Label(mGui,text='Applied Method',font='bold')
        self.similarityAppliedMethodLabel.place(relx=0.1,rely=0.3)
        self.similarityAppliedMethodNumberLabel=Label(mGui,text='Number of Similar Genes',font='bold')
        self.similarityAppliedMethodNumberLabel.place(relx=0.5,rely=0.3)

    #Adding buttons for next steps in grouping outlier result secction #
    def addbuttonsGroupingOutliers(self):
        self.GoSampleLevelButton=Button(mGui,text='Go Sample Level',command=self.goOutlierMenuBck)
        self.GoSampleLevelButton.place(relx=0.1,rely=0.7)

        self.removeOutlierButton=Button(mGui,text='Remove Outlier',command=self.removeOutlierMethod)
        self.removeOutlierButton.place(relx=0.45,rely=0.7)

        self.cancelsimilarityButton=Button(mGui,text='Quit',command=self.mQuit)
        self.cancelsimilarityButton.place(relx=0.8,rely=0.7)
        
    ### Applying Grouping methods ###
    def similarityApply(self):
        #when both normal and tumor samples applied and detected same gene as an outlier, removing same gene#
        if len(self.normalsets[0])>50:
            self.total_outlier_list=removeOutlierGene.removeGene(self.total_outlier_list)
        else:
            pass
        ### FunSimMat, Co-expression and KEGG Pathway Participation ###
        if self.coexpress.get() and self.funsim.get() and self.keggPath.get():
            FunSimMat.funsimmat(self.total_outlier_list,thresholdFSM=self.funsimthresholdEntry.get())
            CoExpression.co_express_calculation(self.cancerandnormaldatafull,self.total_outlier_list,threshold=self.coexpressthresholdEntry.get())
            KEGGpp.keggPathwaySim(self.total_outlier_list,self.keggPathEntry.get(),threshold=self.KEGGthresholdEntry.get())
            #merging outliers from both method#
            self.noneOutlier=[]
            self.temp=FunSimMat.genelist_fromuniprot + CoExpression.co_expressed + KEGGpp.kegg_outlier #outlier genes are kept for further studies
            for gene in self.temp:
                if gene in self.noneOutlier:
                    pass
                else:
                    self.noneOutlier.append(gene)
            #labeling outlier genes for removal #
            tempgene_list=[]
            for i in range(0,len(self.total_outlier_list)):
                if len(self.total_outlier_list[i][0]) > 1:
                    tempgene_list.append(self.total_outlier_list[i][0])
                else:
                    tempgene_list.append(self.total_outlier_list[i])
            self.RealOutlier=list(set(tempgene_list)-set(self.noneOutlier)) # outlier genes for removal
            #Removing and adding labels and buttons#
            self.removeAndAddButtonAndLabel()
            #REsults#
            self.similarityAppliedMethodText1=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText1.place(relx=0.1,rely=0.4)
            if len(FunSimMat.genelist_fromuniprot)==0:
                numberoffunctionalsim='Functional Similarity             None'
            else:
                numberoffunctionalsim='Functional Similarity            '+str(len(FunSimMat.genelist_fromuniprot))
            self.similarityAppliedMethodText1.insert(0.0,numberoffunctionalsim)

            self.similarityAppliedMethodText2=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText2.place(relx=0.1,rely=0.5)
            if len(KEGGpp.kegg_outlier)==0:
                numberofkeggpathway='Kegg Pathyway                     None'
            else:
                numberofkeggpathway='Kegg Pathyway                     '+str(len(KEGGpp.kegg_outlier))
            self.similarityAppliedMethodText2.insert(0.0,numberofkeggpathway)
            
            self.similarityAppliedMethodText3=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText3.place(relx=0.1,rely=0.6)
            numberofcoexpressin='Co-expression                     '+str(len(CoExpression.co_expressed))
            self.similarityAppliedMethodText3.insert(0.0,numberofcoexpressin)
            #Adding buttons #
            self.addbuttonsGroupingOutliers()
            self.exportkeggpathsButton=Button(mGui,text='Export Kegg Pathways',command=self.ExportKEGGfile)
            self.exportkeggpathsButton.place(relx=0.1,rely=0.8)
            return
        
        ### FunSimMat and Co-expression ###
        if self.coexpress.get() and self.funsim.get():
            FunSimMat.funsimmat(self.total_outlier_list,thresholdFSM=self.funsimthresholdEntry.get())
            CoExpression.co_express_calculation(self.cancerandnormaldatafull,self.total_outlier_list,threshold=self.coexpressthresholdEntry.get())
            #merging outliers from both method#
            self.noneOutlier=[]
            self.temp=FunSimMat.genelist_fromuniprot + CoExpression.co_expressed #outlier genes are kept for further studies
            for gene in self.temp:
                if gene in self.noneOutlier:
                    pass
                else:
                    self.noneOutlier.append(gene)
            #labeling outlier genes for removal #
            tempgene_list=[]
            for i in range(0,len(self.total_outlier_list)):
                if len(self.total_outlier_list[i][0]) > 1:
                    tempgene_list.append(self.total_outlier_list[i][0])
                else:
                    tempgene_list.append(self.total_outlier_list[i])
            self.RealOutlier=list(set(tempgene_list)-set(self.noneOutlier)) # outlier genes for removal
            #Removing and adding labels and buttons#
            self.removeAndAddButtonAndLabel()
            #REsults#
            self.similarityAppliedMethodText1=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText1.place(relx=0.1,rely=0.4)
            if len(FunSimMat.genelist_fromuniprot)==0:
                numberoffunctionalsim='Functional Similarity                   None'
            else:
                numberoffunctionalsim='Functional Similarity              '+str(len(FunSimMat.genelist_fromuniprot))
            self.similarityAppliedMethodText1.insert(0.0,numberoffunctionalsim)

            self.similarityAppliedMethodText2=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText2.place(relx=0.1,rely=0.5)
            numberofcoexpressin='Co-expression                   '+str(len(CoExpression.co_expressed))
            self.similarityAppliedMethodText2.insert(0.0,numberofcoexpressin)
            #Adding buttons #
            self.addbuttonsGroupingOutliers()
            return
        
        ### FunSimMat and KEGG Pathway Participation ####
        if self.funsim.get() and self.keggPath.get():
            FunSimMat.funsimmat(self.total_outlier_list,thresholdFSM=self.funsimthresholdEntry.get())
            KEGGpp.keggPathwaySim(self.total_outlier_list,self.keggPathEntry.get(),threshold=self.KEGGthresholdEntry.get())
            #merging outliers from both method#
            self.noneOutlier=[]
            self.temp=FunSimMat.genelist_fromuniprot + KEGGpp.kegg_outlier #outlier genes are kept for further studies
            for gene in self.temp:
                if gene in self.noneOutlier:
                    pass
                else:
                    self.noneOutlier.append(gene)
            #labeling outlier genes for removal #
            tempgene_list=[]
            for i in range(0,len(self.total_outlier_list)):
                if len(self.total_outlier_list[i][0]) > 1:
                    tempgene_list.append(self.total_outlier_list[i][0])
                else:
                    tempgene_list.append(self.total_outlier_list[i])
            self.RealOutlier=list(set(tempgene_list)-set(self.noneOutlier)) # outlier genes for removal
            #Removing and adding labels and buttons#
            self.removeAndAddButtonAndLabel()
            #REsults#
            self.similarityAppliedMethodText1=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText1.place(relx=0.1,rely=0.5)
            if len(KEGGpp.kegg_outlier)==0:
                numberofkeggpathway='Kegg Pathyway                   None'
            else:
                numberofkeggpathway='Kegg Pathyway              '+str(len(KEGGpp.kegg_outlier))
            self.similarityAppliedMethodText1.insert(0.0,numberofkeggpathway)

            self.similarityAppliedMethodText2=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText2.place(relx=0.1,rely=0.4)
            if len(FunSimMat.genelist_fromuniprot)==0:
                numberoffunctionalsim='Functional Similarity                   None'
            else:
                numberoffunctionalsim='Functional Similarity              '+str(len(FunSimMat.genelist_fromuniprot))
            self.similarityAppliedMethodText2.insert(0.0,numberoffunctionalsim)
            #Adding buttons #
            self.addbuttonsGroupingOutliers()
            self.exportkeggpathsButton=Button(mGui,text='Export Kegg Pathways',command=self.ExportKEGGfile)
            self.exportkeggpathsButton.place(relx=0.1,rely=0.8)
            return
        
        ### Co-expression and KEGG Pathway Participation ###
        if self.coexpress.get() and self.keggPath.get():
            CoExpression.co_express_calculation(self.cancerandnormaldatafull,self.total_outlier_list,threshold=self.coexpressthresholdEntry.get())
            KEGGpp.keggPathwaySim(self.total_outlier_list,self.keggPathEntry.get(),threshold=self.KEGGthresholdEntry.get())
            #merging outliers from both method#
            self.noneOutlier=[]
            self.temp=CoExpression.co_expressed + KEGGpp.kegg_outlier #outlier genes are kept for further studies
            for gene in self.temp:
                if gene in self.noneOutlier:
                    pass
                else:
                    self.noneOutlier.append(gene)
            #labeling outlier genes for removal #
            tempgene_list=[]
            for i in range(0,len(self.total_outlier_list)):
                if len(self.total_outlier_list[i][0]) > 1:
                    tempgene_list.append(self.total_outlier_list[i][0])
                else:
                    tempgene_list.append(self.total_outlier_list[i])
            self.RealOutlier=list(set(tempgene_list)-set(self.noneOutlier)) # outlier genes for removal
            #Removing and adding labels and buttons#
            self.removeAndAddButtonAndLabel()
            #REsults#
            self.similarityAppliedMethodText1=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText1.place(relx=0.1,rely=0.4)
            numberofcoexpressin='Co-expression                   '+str(len(CoExpression.co_expressed))
            self.similarityAppliedMethodText1.insert(0.0,numberofcoexpressin)

            self.similarityAppliedMethodText2=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText2.place(relx=0.1,rely=0.5)
            if len(KEGGpp.kegg_outlier)==0:
                numberofkeggpathway='Kegg Pathyway                   None'
            else:
                numberofkeggpathway='Kegg Pathyway                   '+str(len(KEGGpp.kegg_outlier))
            self.similarityAppliedMethodText2.insert(0.0,numberofkeggpathway)
            #Adding buttons #
            self.addbuttonsGroupingOutliers()
            self.exportkeggpathsButton=Button(mGui,text='Export Kegg Pathways',command=self.ExportKEGGfile)
            self.exportkeggpathsButton.place(relx=0.1,rely=0.8)
            return
            
        #### FunSimMat Method ####
        if self.funsim.get():
            FunSimMat.funsimmat(self.total_outlier_list,thresholdFSM=self.funsimthresholdEntry.get())
            self.noneOutlier=FunSimMat.genelist_fromuniprot #outlier genes are kept for further studies
            #labeling outlier genes for removal #
            tempgene_list=[]
            for i in range(0,len(self.total_outlier_list)):
                if len(self.total_outlier_list[i][0]) > 1:
                    tempgene_list.append(self.total_outlier_list[i][0])
                else:
                    tempgene_list.append(self.total_outlier_list[i])
            self.RealOutlier=list(set(tempgene_list)-set(self.noneOutlier)) # outlier genes for removal
            #Removing and adding labels and buttons#
            self.removeAndAddButtonAndLabel()
            #REsults#
            self.similarityAppliedMethodText=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText.place(relx=0.1,rely=0.4)
            if len(FunSimMat.genelist_fromuniprot)==0:
                numberoffunctionalsim='Functional Similarity                   None'
            else:
                numberoffunctionalsim='Functional Similarity              '+str(len(FunSimMat.genelist_fromuniprot))
            self.similarityAppliedMethodText.insert(0.0,numberoffunctionalsim)
            #Adding buttons #
            self.addbuttonsGroupingOutliers()
            return
        
        #### KEGG Pathway Participation ####
        if self.keggPath.get():
            KEGGpp.keggPathwaySim(self.total_outlier_list,self.keggPathEntry.get(),threshold=self.KEGGthresholdEntry.get())
            self.noneOutlier=KEGGpp.kegg_outlier # outlier genes are kept for further studies
            #labeling outlier genes for removal #
            tempgene_list=[]
            for i in range(0,len(self.total_outlier_list)):
                if len(self.total_outlier_list[i][0]) > 1:
                    tempgene_list.append(self.total_outlier_list[i][0])
                else:
                    tempgene_list.append(self.total_outlier_list[i])
            self.RealOutlier=list(set(tempgene_list)-set(self.noneOutlier)) # outlier genes for removal
            #Removing and adding labels and buttons#
            self.removeAndAddButtonAndLabel()
            #REsults#
            self.similarityAppliedMethodText=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText.place(relx=0.1,rely=0.4)
            if len(KEGGpp.kegg_outlier)==0:
                numberofkeggpathway='Kegg Pathyway                   None'
            else:
                numberofkeggpathway='Kegg Pathyway              '+str(len(KEGGpp.kegg_outlier))
            self.similarityAppliedMethodText.insert(0.0,numberofkeggpathway)
            #Adding buttons #
            self.addbuttonsGroupingOutliers()
            self.exportkeggpathsButton=Button(mGui,text='Export Kegg Pathways',command=self.ExportKEGGfile)
            self.exportkeggpathsButton.place(relx=0.1,rely=0.8)
            return
        
        #### Co-expression #####
        if self.coexpress.get():
            CoExpression.co_express_calculation(self.cancerandnormaldatafull,self.total_outlier_list,threshold=self.coexpressthresholdEntry.get())
            self.noneOutlier=CoExpression.co_expressed # outlier genes are kept for further studies
            #labeling outlier genes for removal #
            tempgene_list=[]
            for i in range(0,len(self.total_outlier_list)):
                if len(self.total_outlier_list[i][0]) > 1:
                    tempgene_list.append(self.total_outlier_list[i][0])
                else:
                    tempgene_list.append(self.total_outlier_list[i])
            self.RealOutlier=list(set(tempgene_list)-set(self.noneOutlier)) # outlier genes for removal
            #Removing and adding labels and buttons#
            self.removeAndAddButtonAndLabel()
            #Results #
            self.similarityAppliedMethodText=Text(mGui,height=1,width=60,wrap=None)
            self.similarityAppliedMethodText.place(relx=0.1,rely=0.4)
            numberofcoexpressin='Co-expression                   '+str(len(CoExpression.co_expressed))
            self.similarityAppliedMethodText.insert(0.0,numberofcoexpressin)
            #Adding buttons #
            self.addbuttonsGroupingOutliers()
            return

    #Exporting KEGG genes and pathways #
    def ExportKEGGfile(self):
        exportingKeggGenes.ExportKeggPathways(KEGGpp.finallist)
        tkMessageBox.showinfo(title='Info',message='File were successfully exported')

    #Returning main outlier menu after grouping outliers#
    def goOutlierMenuBck(self):
        #Removing labels and buttons#
        self.groupOutliergeneLabel.destroy()
        self.similarityAppliedMethodLabel.destroy()
        self.similarityAppliedMethodNumberLabel.destroy()
        self.helpmenumethodbutton.destroy()
        self.GoSampleLevelButton.destroy()
        self.removeOutlierButton.destroy()
        self.cancelsimilarityButton.destroy()
        try:
            self.similarityAppliedMethodText1.destroy()
        except:
            pass
        try:
            self.similarityAppliedMethodText2.destroy()
        except:
            pass
        try:
            self.similarityAppliedMethodText.destroy()
        except:
            pass
        try:
            self.exportkeggpathsButton.destroy()
        except:
            pass
        try:
            self.similarityAppliedMethodText3.destroy()
        except:
            pass
        #Adding the buttons and labels for the main outlier menu #
        self.selectLabel=Label(mGui,text='Please select outlier detection level type',font='Times 14 bold italic')
        self.selectLabel.place(relx=0.22,rely=0.2)
        self.genelLevelButton=Button(mGui,text='Gene Level  ',command=self.checkbutgene)
        self.genelLevelButton.place(relx=0.22,rely=0.27)
        self.sampleLevelButton=Button(mGui,text='Sample Level',command=self.checkbutsample)
        self.sampleLevelButton.place(relx=0.52,rely=0.27)

    #Removing out outliers after grouping outliers#
    def removeOutlierMethod(self):
        #removing sample level and gene level outliers#
        try:
            seperateSampleOutlier=self.sampleoutlierstring.split(',')
            self.sampleindex=[]
            # reading sample outliers from text widget and preparing for removal #
            for sample in seperateSampleOutlier:
                if self.selectedRbutton==1:
                    if sample in SampleLevelTCGA.newlist:
                        index=SampleLevelTCGA.newlist.index(sample)
                        self.sampleindex.append(index)
                    else:
                        pass
                else:
                    if sample in SampleLevelGEO.newlist:
                        index=SampleLevelGEO.newlist.index(sample)
                        self.sampleindex.append(index)
                    else:
                        pass
            reverseData=zip(*self.cancerandnormaldatafull)
            #deleting outlier samples #
            for sampleindex in sorted(self.sampleindex,reverse=True):
                del reverseData[int(sampleindex+1)]
            self.backcancerandnormaldatafull=zip(*reverseData)
            #deleting gene level outliers #
            cleaneddatagenes=GeneLevel.deleting_outliers(self.backcancerandnormaldatafull,self.RealOutlier)
            self.cleaneddatagenes=cleaneddatagenes
            #Removing labels and buttons#
            self.GoSampleLevelButton.destroy()
            self.removeOutlierButton.destroy()
            self.cancelsimilarityButton.destroy()
            self.helpmenumethodbutton.destroy()
            self.groupOutliergeneLabel.destroy()
            self.similarityAppliedMethodLabel.destroy()
            self.similarityAppliedMethodNumberLabel.destroy()
            self.similarityAppliedMethodText.destroy()
            try:
                self.exportkeggpathsButton.destroy()
            except:
                pass
            try:
                self.similarityAppliedMethodText1.destroy()
            except:
                pass
            try:
                self.similarityAppliedMethodText2.destroy()
            except:
                pass
            try:
                self.similarityAppliedMethodText3.destroy()
            except:
                pass
            #Results #
            self.dataCharacterwoutOutlierLabel=Label(mGui,text='Characteristics of Pure Data Set ',font='Times 18 bold italic')
            self.dataCharacterwoutOutlierLabel.place(relx=0.1,rely=0.2)

            self.numberofdeletedgenesText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofdeletedgenesText.place(relx=0.1,rely=0.3)
            numberofdeletedgenes='Number of Deleted Outlier Genes:     '+str(len(self.RealOutlier))
            self.numberofdeletedgenesText.insert(0.0,numberofdeletedgenes)

            self.numberofcleanedgenesText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcleanedgenesText.place(relx=0.1,rely=0.35)
            numberofcleangenes='Number of Genes of Data Set    :     '+str(len(self.cleaneddatagenes))
            self.numberofcleanedgenesText.insert(0.0,numberofcleangenes)

            self.numberofdeletedSampleText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofdeletedSampleText.place(relx=0.1,rely=0.4)
            numberofdeletedSample='Number of Deleted Sample       :     '+str(len(self.sampleindex))
            self.numberofdeletedSampleText.insert(0.0,numberofdeletedSample)
    
            self.numberofremainSampleText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofremainSampleText.place(relx=0.1,rely=0.45)
            numberofremainSample='Number of Sample of Data Set   :     '+str(len(reverseData)-1)
            self.numberofremainSampleText.insert(0.0,numberofremainSample)
            
            self.extractfileButton=Button(mGui,text='Extract Files',command=self.browseforextractgeneoutlier)
            self.extractfileButton.place(relx=0.2,rely=0.5)

            self.CancelextractfileButton=Button(mGui,text='Quit',command=self.mQuit)
            self.CancelextractfileButton.place(relx=0.6,rely=0.5)
        #when just gene level outlier is applied #
        except:
            #deleting gene level outliers #
            cleaneddatagenes=GeneLevel.deleting_outliers(self.cancerandnormaldatafull,self.RealOutlier)
            self.cleaneddatagenes=cleaneddatagenes
            #Removing labels and buttons#
            self.GoSampleLevelButton.destroy()
            self.removeOutlierButton.destroy()
            self.cancelsimilarityButton.destroy()
            self.helpmenumethodbutton.destroy()
            self.groupOutliergeneLabel.destroy()
            self.similarityAppliedMethodLabel.destroy()
            self.similarityAppliedMethodNumberLabel.destroy()
            try:
                self.similarityAppliedMethodText.destroy()
            except:
                pass
            try:
                self.exportkeggpathsButton.destroy()
            except:
                pass
            try:
                self.similarityAppliedMethodText1.destroy()
            except:
                pass
            try:
                self.similarityAppliedMethodText2.destroy()
            except:
                pass
            try:
                self.similarityAppliedMethodText3.destroy()
            except:
                pass
            #Results #
            self.dataCharacterwoutOutlierLabel=Label(mGui,text='Characteristics of Pure Data Set ',font='Times 18 bold italic')
            self.dataCharacterwoutOutlierLabel.place(relx=0.1,rely=0.2)

            self.numberofdeletedgenesText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofdeletedgenesText.place(relx=0.1,rely=0.3)
            numberofdeletedgenes='Number of Deleted Outlier Genes:          '+str(len(self.RealOutlier))
            self.numberofdeletedgenesText.insert(0.0,numberofdeletedgenes)

            self.numberofcleanedgenesText=Text(mGui,height=1,width=60,wrap=None)
            self.numberofcleanedgenesText.place(relx=0.1,rely=0.35)
            numberofcleangenes='Number of Genes of Data Set               '+str(len(self.cleaneddatagenes))
            self.numberofcleanedgenesText.insert(0.0,numberofcleangenes)

            self.extractfileButton=Button(mGui,text='Extract Files',command=self.browseforextractgeneoutlier)
            self.extractfileButton.place(relx=0.2,rely=0.45)

            self.CancelextractfileButton=Button(mGui,text='Quit',command=self.mQuit)
            self.CancelextractfileButton.place(relx=0.6,rely=0.45)

    #browsing path to extract files #
    def browseforextractgeneoutlier(self):
        self.browseEntry=Entry(mGui,width=26)
        self.browseEntry.place(relx=0.2,rely=0.65)
        self.browseButton=Button(mGui)
        self.browseButton.config(text="Browse",font=10,command=self.getdirectoryforextract)
        self.browseButton.place(relx=0.7,rely=0.65)

        self.extractfileButton.destroy()
        self.CancelextractfileButton.destroy()

    #directory settings for extracting #
    def getdirectoryforextract(self):
        mGui.fileNameExtract =  tkFileDialog.askdirectory(parent=mGui, title='Select Directory')
        pathname=mGui.fileNameExtract
        self.browseEntry.delete(0,END)
        self.browseEntry.insert(0,pathname)
        self.ExtractButton=Button(mGui,text="Extract",font=10,command=self.FinalExtract)
        self.ExtractButton.place(relx=0.3,rely=0.75)
        self.CancelExtractButton=Button(mGui,text="Quit",font=10,command=self.mQuit)
        self.CancelExtractButton.place(relx=0.6,rely=0.75)

    #extracting expression files #
    def FinalExtract(self):
        new_directory=self.browseEntry.get()
        os.chdir(new_directory)
        try:
            for sampleindex in sorted(self.sampleindex,reverse=True):
                del self.data_set[0][int(sampleindex+1)]

            for sampleindex in sorted(self.sampleindex,reverse=True):
                del self.data_set[1][int(sampleindex+1)]        
            #TCGA dataset#
            if self.selectedRbutton==1:
                final_data=[self.data_set[0]] + [self.data_set[1]] + self.cleaneddatagenes
                extractFile.Extracting_files_TCGA(final_data)
            #GEO dataset#
            else:
                final_data=[self.data_set[0]] + self.cleaneddatagenes
                extractFile.Extracting_files_GEO(final_data)
            tkMessageBox.showinfo(title='Info',message='Files were successfully extracted')
        except:
            #TCGA dataset#
            if self.selectedRbutton==1:
                final_data=[self.data_set[0]] + [self.data_set[1]] + self.cleaneddatagenes
                extractFile.Extracting_files_TCGA(final_data)
            #GEO dataset#
            else:
                final_data=[self.data_set[0]] + self.cleaneddatagenes
                extractFile.Extracting_files_GEO(final_data)
            tkMessageBox.showinfo(title='Info',message='Files were successfully extracted')

    # Ouit the program #
    def mQuit(self):
        mExit=tkMessageBox.askokcancel(title='quit',message='Do you want to quit?')
        if mExit>0:
            mGui.destroy()
            return           
mGui=Tk()
mGui.geometry("500x500+400+400")
mGui.title('outlier detection')
mlabel=Label(mGui,text="Outlier Detection Tool ",fg='white',bg='blue',height=3,font='Verdana 20 bold italic').pack(fill=BOTH)


app=Tool(mGui)
mGui.mainloop()
