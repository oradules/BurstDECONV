#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 10:10:26 2021

@author: mdouaihy
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from ecdfEstimate import ecdf_bounds
from fit3M1Parameters import fit3M1Parameters
from fit3M2Parameters import fit3M2Parameters


def fit3(dirwrite,name,sd,dt,dtc):
    

    #Store: define matrices that store parameters and objective functions 
    # for each iteration of the least_square_function
    # to have better results of least square we ran the function for 100 iterations with random initial values for each iteration
    store = np.empty((0,6))

    if len(dt)!=0: 
        for cens in range(2):
            if cens:
                xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc),np.append( np.zeros(len(dt)), np.ones(len(dtc)) ) )
            else:
                xs,fs,flo,fup =ecdf_bounds(dt )

            ## fit distribution of spacings using combination of two exponentials

            xs = xs[:-1]
            fs = fs[:-1]
            flo = flo[:-1]
            fup = fup[:-1]

            ##############

            sN = np.sqrt(len(xs))

            def exp_fitness(k): #objective function
                return (np.log(k[3]*np.exp(k[0]*xs)+k[4]*np.exp(k[1]*xs)+(1-k[3]-k[4])*np.exp(k[2]*xs) ) -np.log(1-fs)) /sN # k: parameters

            k00 = np.array([-0.01, -0.01,-0.001, 0.25, 0.25])
            amp = np.array([np.log(100), np.log(100), np.log(100), np.log(1), np.log(1)])

            NbIterationinFit = 100

            test=0
            for mc in range(NbIterationinFit):
                print('iteration nbr ', mc)
                while test==0: #test is just to re-do the iteration until we encounter no error

                    ## change k00 which is the initial value
                    factor = np.exp(amp* (2* np.random.uniform(size=5)-1))
                    k0 = k00*factor
                    k0[3:5] = 2*np.random.uniform(size=2)-1

                    ## sort k0(1:2)
                    k0[0:3] = np.sort(k0[0:3])

                    ## impose constraints
                    if not (sum(k0[3:5])<1 and k0[0]*k0[3]+ k0[1] * k0[4] +k0[2]*(1-sum(k0[3:5])) <0):
                        while not (sum(k0[3:5])<1 and k0[0]*k0[3]+ k0[1] * k0[4] +k0[2]*(1-sum(k0[3:5])) <0):
                            k0[3:5] = 2*np.random.uniform(size=2)-1 # A1, A2 value

                    # Use the fcn lsqnonlin

                    try:
                        k = least_squares(exp_fitness, k0,bounds=(-np.inf,[0,0,0,np.inf,np.inf]), ftol = (1e-8), xtol= (1e-10), method='trf').x
                        obj = sum(exp_fitness(k)**2)
                        test = 1
                    except:
                        pass
                test=0

                # write down results


                ## sort k
                A = np.array([k[3],k[4], 1-k[3]-k[4]])  # A values before sorting 
                kk = np.sort(k[0:3])
                IX = np.argsort(k[0:3])
                k[:3]= kk
                A=A[IX]
                k[3:5]=A[0:2]

                store = np.append(store,np.append(k,np.array([obj])).reshape(1,len(np.append(k,np.array([obj])))),axis=0)
        
        [resM1,reslM1,reshM1]=fit3M1Parameters(store,dt,dtc,sd,dirwrite,name)
        [resM2,reslM2,reshM2]=fit3M2Parameters(store,dt,dtc,sd,dirwrite,name)




    else:
        
        resM1 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, sd[1], sd[0]]
        reslM1=np.zeros(len(8))
        reshM1=np.zeros(len(8))
        resM2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, sd[1], sd[0]]
        reslM2=np.zeros(len(8))
        reshM2=np.zeros(len(8))
    return [resM1, reslM1, reshM1,resM2, reslM2, reshM2]