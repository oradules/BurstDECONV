#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 17:38:18 2022

@author: mdouaihy
"""

import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import shutil
import timeit

#from joblib import Parallel, delayed

### imput: index of the file
class readDataInFile:
    def __init__(self,DataFilePath,outputFolder,data_type):
        dataFile='npzFile'+data_type
        PyFilePath = os.path.join(outputFolder,dataFile)
        ### where the images will be stored
        DataFilePath0 = os.path.join(outputFolder,'images')
        ### creating the folder
        if os.path.exists(PyFilePath): 
            shutil.rmtree(PyFilePath, ignore_errors = True)  
        os.mkdir(PyFilePath) 
        if os.path.exists(DataFilePath0): 
            shutil.rmtree(DataFilePath0, ignore_errors = True)  
        os.mkdir(DataFilePath0) 
        self.PyFilePath = PyFilePath
        self.DataFilePath0=DataFilePath0
        self.DataFilePath=DataFilePath
        self.lists, self.nexp = self.ListingFiles(self.DataFilePath)    
        self.file_name_list=self.lists
        for data_i in range(self.nexp):
            DataFileName=DataFilePath + self.file_name_list[data_i] 
            if '.~' in DataFileName:
                continue 
            self.a=self.read_data(data_i)
    def read_data(self,data_i):

        nexp=self.nexp
        file_name_list=self.file_name_list
        
        DataFileName = self.DataFilePath + file_name_list[data_i] #path to the specified file

        
        DataExp =pd.read_excel(DataFileName,usecols=np.arange(3,600),skiprows=0).to_numpy() #extracting the data from the excel sheet
        DataExp=DataExp[:,0:-1:3]
        
        Spot_names=pd.read_excel(DataFileName,usecols=np.arange(3,len(DataExp[0])*3+1),skipfooter =len(DataExp)+1).columns.tolist()
        lst1 = list(range(1, len(Spot_names), 3)) 
        lst2 = list(range(2, len(Spot_names), 3))
        Spot_names = np.delete(Spot_names, lst1+lst2)
        
        
        indx_to_remove = ~np.all(DataExp.transpose() == 0, axis=1) #avoid zero cells
        DataExp=DataExp.transpose()[indx_to_remove].transpose()
        Spot_names = np.array(Spot_names)[indx_to_remove]
        
        ## creating a folder to store the images
        dirname=file_name_list[data_i].replace('.xlsx','')

        WriteTo=self.DataFilePath0+'/'+dirname+'_image'
        if os.path.exists(WriteTo): 
            shutil.rmtree(WriteTo, ignore_errors = True)  
        os.mkdir(WriteTo)    

        ## file name and size of dataexp
        n=DataExp.shape
        fn='data_'+self.file_name_list[data_i].replace('_','').replace('.xlsx','')

        ## ploting the intensity spot for each nuclei in a 6*6 plot
        for iii in range(n[1]):

            ifig= math.floor((iii)/36)+1 
            h=plt.figure(ifig+data_i*10)

            plt.subplot(6,6, (iii%36+1))
            plt.plot(np.arange(1,len(DataExp)+1), DataExp[:,iii], color='black', linewidth=0.1)
            h.subplots_adjust(hspace = 0.5, wspace=0.5)

            plt.title(Spot_names[iii].replace('_','-'),fontsize=6, pad=2)
            plt.yticks(fontsize=5)

            if (iii+1)%36==0 or (iii+1)== n[1]:
                plt.xticks(fontsize=5)
                figfile=WriteTo+'/fig'+str(ifig)+'.pdf'
                h.savefig(figfile)
                
            else:                
                plt.xticks([], [])

        #saving the data
        sd=DataExp.shape      
        Frames=np.arange(0,sd[0])
        Samples=np.arange(1,sd[1])
        fname=self.PyFilePath+fn+'.npz'
        self.DataExp=DataExp
        self.Frames=Frames
        self.Samples=Samples

        np.savez(fname,DataExp=DataExp,Samples=Samples,Frames=Frames)
        self.fname=fname
        return self.DataFilePath0

    def ListingFiles(self,DataFilePath,extension='.xlsx'):
        self.DataFilePath=DataFilePath
        file_name_list = np.array(os.listdir(self.DataFilePath)) # list of subdirectories containing data from different genotypes
        nexp = len(file_name_list) # number of files inside he folder
        #  DataFilePath0='../images/' # place to store images
        #just to avoid erros in the opening of the files
        file_name_list_modified=np.array([])
        for i in range(nexp):
                if '.~' in file_name_list[i] :
                    file_name_list[i]=file_name_list[i].replace('._','',1)
                    
        for i in range(nexp):
            if not '._' in file_name_list[i] and extension in file_name_list[i]:
                file_name_list_modified=np.append(file_name_list_modified,file_name_list[i])
                    

        file_name_list=file_name_list_modified            
        file_name_list=np.unique(file_name_list)
        nexp=len(file_name_list)
        return file_name_list, nexp
            
        
            
