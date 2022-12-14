#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 00:34:49 2021

@author: mdouaihy
"""

import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
import math
import pandas as pd
from extentForPlot import extentForPlot
from scipy.io import loadmat
import seaborn as sns
from common_fit2_part import fit2
from common_fit3_part import fit3



class fit:
    def __init__(self,path,parameter,combined, outputpath,data_type, heterogenious_signal = 0):
        self.path=path
        self.parameterpath=parameter
        
        ### parameters used
        filecontent=np.load(self.parameterpath)
        Polym_speed=filecontent['Polym_speed']
        TaillePreMarq=filecontent['TaillePreMarq']
        TailleSeqMarq=filecontent['TailleSeqMarq']
        TaillePostMarq=filecontent['TaillePostMarq']
        EspaceInterPolyMin=filecontent['EspaceInterPolyMin']
        FreqEchImg=filecontent['FreqEchImg']
        FrameLen = filecontent['FrameLen']
        FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed) # how many interval(possible poly start position) in 1s
        self.FreqEchSimu =FreqEchSimu 
        
        ####### parameters for the plots
        fsz=16 #figure size
        
         ## function needed to set the the parameters for the color map
        cm_jet= plt.cm.get_cmap('jet') # set the colormap to jet array
        
        DataFilePath0 = outputpath+'Results'+data_type
        if os.path.exists(DataFilePath0):
            shutil.rmtree(DataFilePath0, ignore_errors = True)

        os.mkdir(DataFilePath0)


        ### creating of xls sheet for all kind of models
        
        ## 2 states
        
        xlsfilename2states = DataFilePath0 + '/fit2_results.xlsx'

        # Setting the Names of the Data output in the excel sheet
        xls_cont = pd.DataFrame({'Data': [],'k1p': [], 'k1m': [], 'k2': [],'p1': [], 'p2': [],
                                 'l1':[],'l2':[],'A1':[],'A2':[], 'Obj': [],'KS test':[],'Nuclei': [],'Frames': []})    

        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer2states = pd.ExcelWriter(xlsfilename2states, engine='xlsxwriter')

        # Convert the dataframe to an XlsxWriter Excel object.
        xls_cont.to_excel(writer2states, sheet_name='Sheet1', index=False)
        
        #############################
        
        ### 3 states M1
        xlsfilename3statesM1 = DataFilePath0 + '/fit3M1_results.xlsx'
         
        # Setting the Names of the Data output in the excel sheet
        xls_cont = pd.DataFrame({'Data': [],'k1p' : [],'k1m' : [],'k2p' : [],
                                 'k2m' : [],'k3' : [],'p1' : [],'p2': [],
                                 'p3' : [],'lambda1': [],'lambda2': [],
                                 'lambda3': [],'A1': [],'A2': [],'A3': [],'S1': [],
                                 'Obj': [],'KS test': [],'Nuclei': [],'Frames': []})  
        
        
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer3statesM1 = pd.ExcelWriter(xlsfilename3statesM1, engine='xlsxwriter')

        # Convert the dataframe to an XlsxWriter Excel object.
        xls_cont.to_excel(writer3statesM1, sheet_name='Sheet1', index=False)
        
        #########################
        
        ### 3 states M2
        xlsfilename3statesM2 = DataFilePath0 + '/fit3M2_results.xlsx'
         
        # Setting the Names of the Data output in the excel sheet
        xls_cont = pd.DataFrame({'Data': [],'k1p' : [],'k1m' : [],'k2p' : [],
                                 'k2m' : [],'k3' : [],'p1' : [],'p2': [],
                                 'p3' : [],'lambda1': [],'lambda2': [],
                                 'lambda3': [],'A1': [],'A2': [],'A3': [],'S1': [],
                                 'Obj': [],'KS test': [],'Nuclei': [],'Frames': []}) 
        
        
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer3statesM2 = pd.ExcelWriter(xlsfilename3statesM2, engine='xlsxwriter')

        # Convert the dataframe to an XlsxWriter Excel object.
        xls_cont.to_excel(writer3statesM2, sheet_name='Sheet1', index=False)
        
        #########################
        
        ####################################################################
        
        ## loading the result of the deconvolution
        NPZFilePath = path;  
        file_name_list = np.array(os.listdir(NPZFilePath)) # list of the data
        self.file_name_list=file_name_list
        nexp = len(file_name_list) # length of the list
        nfiles=nexp
        print(NPZFilePath)

        #######################################################################
        pooled =combined; #### if this is 1, pool all the result files from NPZfilePath
        ############### if not, use each file separately########################
    
        if not pooled:
            nfiles=nexp
        else:
            nfiles=1
            
        ### starting the fit for each file
        for ifile in range(nfiles):

            if pooled: 
                
            ### first compute max dimension
                nmax=1
                nmaxpos=1

                for iii in range(nexp):
                    fname=file_name_list[iii]
                    ffname=NPZFilePath+fname
                    if '.npz' in ffname:
                        fnameContent=np.load(ffname)
                    else:
                        fnameContent=loadmat(ffname)
                    DataPred = fnameContent['DataPred']
                    DataExp=fnameContent['DataExp']
                    PosPred=fnameContent['PosPred']
                    n2 =DataExp.shape
                    n3= PosPred.shape
            
                    if n2[0] >nmax:
                        nmax = n2[0]
            
                    if n3[0]> nmaxpos:
                        nmaxpos = n3[0]
               
                #### lump files, Procustes method 
                dataExp=np.empty((nmax, 0), int)
                dataPred=np.empty((nmax,0), int)
                posPred= np.empty((nmaxpos,0), int)
                tmax = np.empty((0, n2[1]), int)

        
            
                for iii in range(nexp): 
                    fname=file_name_list[iii]
                    ffname=NPZFilePath+fname
                    if '.npz' in ffname:
                        fnameContent=np.load(ffname)
                    else:
                        fnameContent=loadmat(ffname)
                    DataPred = fnameContent['DataPred']
                    DataExp=fnameContent['DataExp']
                    PosPred=fnameContent['PosPred']
                    n2 =DataExp.shape
                    n3= PosPred.shape
      
                    DataExp=np.append(DataExp,np.zeros((nmax-n2[0],n2[1])),axis=0)
                    DataPred=np.append(DataPred,np.zeros((nmax-n2[0],n2[1])),axis=0)
                    PosPred=np.append(PosPred,np.zeros((nmaxpos-n3[0],n3[1])),axis=0)
            
                    # we are adding all the data from different files together
                    dataExp = np.append(dataExp, DataExp, axis=1) 
                    dataPred = np.append(dataPred, DataPred, axis=1)
                    posPred = np.append(posPred, PosPred, axis=1)

                    tmax=np.append(tmax, n2[0]/FreqEchImg*np.ones(n2[1]))
            
        
                DataExp = dataExp.copy()
                DataPred=dataPred.copy()
                PosPred = posPred.copy()
            else:
                fname = file_name_list[ifile]
                #### full path file name
                ffname = NPZFilePath+ fname
                if '.npz' in ffname:
                    fnameContent=np.load(ffname)
                else:
                    
                    fnameContent=loadmat(ffname)
                DataPred = fnameContent['DataPred']
                DataExp=fnameContent['DataExp']
                PosPred=fnameContent['PosPred']
                n2=DataExp.shape
                tmax=n2[0]/FreqEchImg*np.ones(n2[1]) #### movie length, the same for all nuclei in a data sets 
  
            
            
            ### extract short name from result file name
            # iend=fname.index('_artificial')
            name=fname[7:-4]   
            self.name=name

            name = name.replace('CalibratedTraces','')
            

            ### where to write figure files 
            dirwrite = DataFilePath0+'/'+name+'_result'
            if os.path.exists(dirwrite):
                shutil.rmtree(dirwrite, ignore_errors = True)

            os.mkdir(dirwrite)

            n = DataExp.shape
            nexp = n[1]
            ## parameters
            frame_num = n[0]
   
            MIntensity = np.array([])
            T0 = np.array([])
            
            ## find first hit which is 1/5th of the max intensity
            for data_i in range(nexp):
                
                ifig= math.floor((data_i)/36)+1
                
                max_intensity=max(DataPred[:,data_i])
                MIntensity=np.append(MIntensity,max_intensity)
            
                ihit=np.where(DataPred[:,data_i]> max_intensity/5)[0]
            
                if len(ihit)== 0:
                    ihit=n[0]
                    t0o1 = (ihit)/FreqEchImg 
                else:
                    ihit=min(np.where(DataPred[:,data_i]> max_intensity/5)[0])
                    t0o1 = (ihit+1)/FreqEchImg 
            
                t0=t0o1   
                T0 = np.append(T0, t0) #stores the  first hit for each nuclei
                
                
                h = plt.figure(ifig+ ifile)#, figsize=[10,12]   
                plt.subplots_adjust(hspace=1,wspace=1)
                plt.subplot(6,6,(data_i%36+1))
                plt.fill_between(np.array([t0/60, tmax[data_i]/60]), np.array([150,150]), facecolors =np.array([0.9, 0.9, 0.9])) #this what is really analyzed
                plt.plot(np.arange(0, frame_num)/FreqEchImg/60, DataExp[:,data_i].T, color = 'k', linewidth = 0.1)
                plt.plot(np.arange(0, frame_num)/FreqEchImg/60, DataPred[:,data_i].T, color = 'r', linewidth = 0.1)
                plt.xlim(0,40)
                plt.ylim(0, 100)
                plt.xticks(fontsize=5)
                plt.yticks(fontsize=5)
                sns.despine()
                if data_i%36==35 or (data_i+1)==nexp:
                    figfile=dirwrite+'/figure'+str(ifig)+'.pdf'
                    h.savefig(figfile)
                    plt.close()
        
            ### Figure showing Data Signal Prediction
            h=plt.figure(40)
            sz= DataPred.shape
            Y_normal = np.arange(1,sz[1]+1)
            Y=Y_normal[::-1]
            X = np.arange(0, sz[0])/FreqEchImg/60
            plt.imshow(DataPred.T, cmap=cm_jet, extent=extentForPlot(X).result + extentForPlot(Y).result, aspect='auto', origin='upper')
            plt.xlabel('Time [min]', fontsize=12)
            plt.ylabel('Transcription site', fontsize=12)
            cb= plt.colorbar()
            cb.ax.tick_params(labelsize=fsz)
            figfile=dirwrite+'/DataPred_'+name+'.pdf'
            h.savefig(figfile, dpi=800) 
            plt.close()

            ### Figure showing Data Signal Experimental
            h = plt.figure(50)
            plt.imshow(DataExp.T, cmap=cm_jet, extent=extentForPlot(X).result + extentForPlot(Y).result, aspect='auto', origin='upper')
            plt.xlabel('Time [min]', fontsize=12)
            plt.ylabel('Transcription site', fontsize=12)
            plt.colorbar()
            figfile=dirwrite+'/DataExp_'+name+'.pdf'
            h.savefig(figfile, dpi=800) 
            plt.close()

            ### Figure showing Data Position Prediction

            time_window = 40
            nbr_frame_in_window = int(time_window/FrameLen)
            
            sd_pospred = np.array([PosPred.shape[0]//nbr_frame_in_window+1, nbr_frame_in_window])
            density_pospred = np.zeros((PosPred.shape[1], sd_pospred[0]))
            
            
            for ii in range(PosPred.shape[1]):
                PosPred_i = PosPred[:,ii]
                PosPred_i_reshape = np.resize(PosPred_i, sd_pospred[0]*sd_pospred[1]).reshape(sd_pospred[0],sd_pospred[1])
                PosPred_i_reshape[-1, -(sd_pospred[0]*sd_pospred[1]-len(PosPred_i)):] = np.zeros((sd_pospred[0]*sd_pospred[1]-len(PosPred_i)))
                
                density_pospred_i = np.sum(PosPred_i_reshape, axis =1)
                density_pospred[ii, :] = density_pospred_i
            
            
            h=plt.figure()
            
            X=np.arange(0,sd_pospred[0]*sd_pospred[1])*EspaceInterPolyMin/Polym_speed/60  -(TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed/60 ### time
            X = X[::nbr_frame_in_window]
            df = pd.DataFrame(density_pospred, columns=np.round(X,1))
            sns.heatmap(df, cmap = 'YlOrBr', xticklabels = round(len(X)/6), yticklabels = round(len(PosPred[0])/4)) #, vmin = 0, vmax = np.max(density_pospred)*2)
            plt.xlabel('Time [min]', fontsize=12)   
            plt.ylabel('Transcription site', fontsize=12)
            figfile=dirwrite+'/PosPred'+name+'.pdf'
            h.savefig(figfile, dpi=800) 
            plt.close()
            
            ### compute distribution of spacings
            nn=PosPred.shape

            dt=np.array([])
            dtc=np.array([])
            
            figfile=dirwrite +'/PosPred'+name+'.csv' ###  text file for pol positions 
            # fid = open(figfile,'w+')
            
            name_xls = np.arange(nn[1])
            name_xls = ['n '+ str(ls+1) for ls in name_xls] 
            
            
            times_to_xls=[]
            for i in range(nn[1]): #for all cells
                times = (np.where(PosPred[:,i]==1)[0]+1) / FreqEchSimu -(TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed 
                times_to_xls.append(( (times+(TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed  )/60).tolist())
                
                # fid.writelines([' \n'+ str(times/60)])

                if len(times) !=0:
                    dtimes = np.diff(times)

                    ### find first index larger than T0
                    to_determin_istart=np.where(times > T0[i] - (TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed )[0]
                    if len(to_determin_istart)==0:
                        istart=0
                    else:
                        istart= min(to_determin_istart)
                    dt=np.append(dt, dtimes[istart:])
                    if tmax[i]-times[-1]>0:
                        dtc = np.append(dtc, tmax[i]-times[-1]) ### the very last time

            # fid.close()
            
            xls_data=pd.DataFrame(times_to_xls, index = name_xls).T
            xls_data = xls_data.fillna(0)
            sort_xls_data = xls_data.sort_values(by = 0, axis = 1)
            sort_xls_data.to_csv(figfile, index=False)
            
            
            sd=DataExp.shape
            
            ########## save parameters results for 2 state model
            [res, resl, resh]=fit2(dirwrite,name,sd,dt,dtc)
            
            df1 = pd.DataFrame([res.tolist(), #best result
                                resl.tolist(), # low 
                                resh.tolist()]) # high


            df1.to_excel(writer2states,sheet_name='Sheet1', startrow=4*(ifile+1)-3, startcol=1, header=False, index=False)

            df2 = pd.DataFrame([name.replace('result_','')]) #filename
            
            df2.to_excel(writer2states,sheet_name='Sheet1', startrow=4*(1+ifile)-3, startcol=0, header=False, index=False)            
            
            ############################################################################
            
            
            ########## save parameters results for 3 state model
#            np.savez(dirwrite+'/fit3_comparison_'+name+'.npz',sd=sd,dt=dt,dtc=dtc)
            [resM1, reslM1, reshM1,resM2, reslM2, reshM2]=fit3(dirwrite,name,sd,dt,dtc)
            
            #### Model M1
            
            df1M1 = pd.DataFrame([resM1.tolist(), #best result
                    reslM1.tolist(), # low 
                    reshM1.tolist()]) # high
            
            df1M1.to_excel(writer3statesM1,sheet_name='Sheet1', startrow=4*(ifile+1)-3, startcol=1, header=False, index=False)
            df2M1 = pd.DataFrame([name.replace('result_','')]) #filename
            df2M1.to_excel(writer3statesM1,sheet_name='Sheet1', startrow=4*(ifile+1)-3, startcol=0, header=False, index=False)
            ########################################

            #### Model M1
            
            df1M2 = pd.DataFrame([resM2.tolist(), #best result
                    reslM2.tolist(), # low 
                    reshM2.tolist()]) # high

            df1M2.to_excel(writer3statesM2,sheet_name='Sheet1', startrow=4*(ifile+1)-3, startcol=1, header=False, index=False)
            df2M2 = pd.DataFrame([name.replace('result_','')]) #filename
            df2M2.to_excel(writer3statesM2,sheet_name='Sheet1', startrow=4*(ifile+1)-3, startcol=0, header=False, index=False)
            ########################################        

        writer2states.save()
        writer3statesM1.save()
        writer3statesM2.save()
