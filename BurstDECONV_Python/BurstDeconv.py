# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 15:28:13 2021

@author: mdouaihy
"""

import sys
sys.path.append('utilities/')

import multiprocessing as mp
import numpy as np
import os


print("Maximum number of threads available  ", mp.cpu_count())

""" ---------- USER PROVIDED INFORMATION ------------- """

DataFilePath = '../../artificial_data/sample data/'  ## path to the folder containing the data
outputFolder= '../../artificial_data/sample data/output/'
data_type = ''
extension ='.xlsx'

combined = 1 #if this is 1, combine all the result files in the same folder
        ############### if not, use each file separately########################
numberOfWorkers = 5 # it's recomended to use less than half of the total threads available 
colision = 1
retention = 0 # in seconds


npzFilePath=outputFolder
parameterFileName='drosoParameters'+ data_type #should always end with Parameters


##### parameters  

Intensity_for_1_Polym = 1 # calibration factor
Polym_speed = 25 # Polymerase speed'
TaillePreMarq = 1366 # length of mRNA until we reach the beginning of the MS2 sequence in bp
TailleSeqMarq = 1292 # length of MS2 sequence in bp
TaillePostMarq = 401 + 67*0 # length of mRNA from last loop of MS2 until the polymerase leaves the TS
EspaceInterPolyMin = 30 # minimal distance between polymerase
FrameLen = 4.64 # frame length in seconds or time resolution

FreqEchImg = 1/FrameLen 
DureeSignal = (TaillePreMarq + TailleSeqMarq + TaillePostMarq) / Polym_speed
FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed)

np.savez(npzFilePath+parameterFileName+'.npz', 
          Polym_speed = Polym_speed,  
          TaillePreMarq = TaillePreMarq,
          TailleSeqMarq = TailleSeqMarq,
          TaillePostMarq = TaillePostMarq,
          EspaceInterPolyMin = EspaceInterPolyMin,
          FrameLen = FrameLen,
          Intensity_for_1_Polym = Intensity_for_1_Polym,
          FreqEchImg = FreqEchImg,
          DureeSignal = DureeSignal,
          FreqEchSimu = FreqEchSimu,
          retention = retention
        )


##########################################################


""" Import the preprocessing package for Drosophila

#### please chose between these three type:
##   1- readDrosoData_2spot_1st_activation
##   2- readDrosoData_2spot_2nd_activation
##   3- readDrosoData_2spot_all_nuclei
##   4- readDrosoData_1spot

"""

from readDrosoData_1spot import readDataInFile

drosoSNAdata = readDataInFile(DataFilePath,outputFolder,data_type + '/', extension)


# import function to perform deconvolution
from deconvolveMyData import deconvolveMyData


#----------------- Run deconvolution----------------------#

DataFilePathToPreProcessedFile = os.path.join(outputFolder,'npzFile' +data_type + '/')



deconvolveMyData(DataFilePath = DataFilePathToPreProcessedFile, outputFolder=outputFolder, parameterFile=npzFilePath+parameterFileName+ '.npz',
                   number_of_workers=numberOfWorkers,data_type=data_type + '/', colision = colision)




from common_part_fitting import fit

# #----------------- Fit 2 state model-----------------------#


pathToDeconvolutionResultsFolder=os.path.join(outputFolder,'resultDec'+data_type + '/')
pardor=npzFilePath+parameterFileName+ '.npz'  #using the whole path
parameterPath=npzFilePath+parameterFileName+ '.npz'
FitResults=fit(pathToDeconvolutionResultsFolder,parameterPath,combined,outputFolder,data_type + '/')