#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 13:08:50 2023

@author: rachel
"""

from fit import fitLong

pathToDeconvolutionResultsFolder='/home/rachel/Downloads/3S_L+S/TR_3/resultDec/'
parameterFilePath='/home/rachel/Downloads/3S_L+S/hivParameters.npz'
fitLong(pathToDeconvolutionResultsFolder, parameterFilePath, combined=0, visualize=1, outputpath='/home/rachel/Downloads/3S_L+S/TR_3/', fit3exp=True)
