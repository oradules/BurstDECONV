#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 00:39:16 2021

@author: mdouaihy
"""

import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from ecdfEstimate import ecdf_bounds
import math

def fit2(dirwrite,name,sd,dt,dtc):
        
    ####### parameters for the plots
    fsz=16 #figure size
    
    #Store: define matrices that store parameters and objective functions 
    # for each iteration of the least_square_function
    # to have better results of least square we ran the function for 100 iterations with random initial values for each iteration
    store = np.empty((0,5))
    
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
                return np.abs(np.log(k[2]*np.exp(k[0]*xs)+(1-k[2])*np.exp(k[1]*xs) ) -np.log(1-fs)) /sN # k: parameters

            k00 = np.array([-0.01, -0.001, 1])
            amp = np.array([np.log(100), np.log(100)])
            NbIterationinFit = 100

            test=0
            for mc in range(NbIterationinFit):
                print('iteration nbr ',mc)
                while test==0: #test is just to re-do the iteration until we encounter no error
                    ## change k00 which is the initial value
                    factor = np.exp(amp* (2* np.random.uniform(size=2)-1))
                    k0 = k00.copy()
                    k0[:2] = k0[0:2]*factor
                    k0[2]=2*np.random.uniform(size=1)-1

                    ## sort k0(1:2)
                    k0[0:2] = np.sort(k0[0:2])

                    ## impose constraints
                    if not (k0[0]*k0[2]+k0[1]*(1-k0[2])<0):
                        while not (k0[0]*k0[2]+k0[1]*(1-k0[2])<0):
                            k0[2] = 2*np.random.uniform(size=1)-1 # A1, A2 value

                    # Use the fcn lsqnonlin which is the least square function

                    try:
                        k = least_squares(exp_fitness, k0,bounds=(-np.inf,[0,0,np.inf]) ,ftol = (1e-8),max_nfev= 1e6, xtol= (1e-10)).x
                        obj = sum(exp_fitness(k)**2)
                        test = 1
                    except:
                        pass
                test=0

                # write down results


                ## sort k
                A = np.array([k[2],1-k[2]])  # A values before sorting 
                kk = np.sort(k[0:2])
                IX = np.argsort(k[0:2])
                k[:2]= kk
                A=A[IX]
                k[2]=A[0]
                k_obj=np.append(k,np.array([obj]))
                to_store=np.append(k_obj,cens)
                to_store=to_store.reshape(1,len(to_store))
                store = np.append(store,to_store,axis=0)

                
        # select optimal 

        # ind is index of the least square with real numbers
        ind = np.where(np.max(np.abs(store[:,0:3].imag )  ,axis=1) <1e-10)[0]
        # objmin is minimun of the above results (least square value of ind)
        objmin = np.min(store[ind,3])

        #ind is the index of the minimum value of the least square with real value
        indmin=np.argmin(store[ind,3])
        imin = ind[indmin]

        #overflow help us set the Uncertainty interval
        overflow = 1

        # ind index where the least square are real  and less than < (1+overflow)*objmin)
        ind = np.where( (store[:,3] < (1+overflow)*objmin) & (np.max(np.abs(store[:,0:3].imag   ),axis=1) <1e-10)) [0]
        ksel = store[ind,0:3].real
        kmin = store[imin,0:3].real
        censmin=store[imin,4].real

        if censmin:
            xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc),np.append( np.zeros(len(dt)), np.ones(len(dtc)) ) )
        else:
            xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc))


        ## Survival function
        h = plt.figure(70)
        plt.semilogy(xs, 1-fs, 'o', color='r', mfc = 'none', markersize=9, linestyle='', label = 'S funct from data') #empirical function
        plt.semilogy(xs, 1-flo, '--r') # lower confidence
        plt.semilogy(xs, 1-fup, '--r') # upper confidence
        pred = kmin[2]*np.exp(kmin[0]*xs)+(1-kmin[2])*np.exp(kmin[1]*xs)
        plt.semilogy(xs, pred, 'k', linewidth = 3, label = '2 States fit S funct') # predicted 2 exp
        plt.xlim(0,250)
        plt.ylim(1e-6, 1)
        plt.xlabel('Time [s]',fontsize=fsz)
        plt.ylabel('Survival function',fontsize = fsz)
        plt.legend()
        
        # compute 3 rates k1p,m k2 from the 3 parameters
        l1 = kmin[0]
        l2 = kmin[1]
        A1 = kmin[2]
        A2 = 1-A1
                
        S1 = A1*l1+(A2)*l2
        k2 = -S1
        S2 = A1*(l1)**2+(A2)*(l2)**2
        S3 = A1*(l1)**3+(A2)*(l2)**3
        k1m = S1-S2/S1
        k1p = (S3*S1-S2**2)/S1/(S1**2-S2)
        plt.title(r'$k_1^-$='+str("%.1g" % k1m)+' $k_1^+$='+str("%.1g" % k1p)+' $k_2$='+str("%.1g" % k2))

        figfile = dirwrite+'/Fit_2States_'+name+'.pdf'
        h.savefig(figfile)
        plt.close()

        ####### KS test

        thr= 20
        ind = np.where(xs>thr)[0]
        dist = np.max(np.abs( fs[ind] -1+pred[ind]))
        nsample=len(dt[dt>thr])
        N=10
        c=math.sqrt(nsample)*dist
        r=np.arange(1,N+1)
        aa=2*sum(np.exp(-2*c**2*r**2)*((-1)**(r-1)))
        
        ##########################

        ## optimal
        res = np.array([k1p,k1m,k2,k1m/(k1m+k1p),k1p/(k1m+k1p),l1,l2,A1,A2,objmin,aa,sd[1],sd[0]])
        
        # compute interval
        l1=ksel[:,0] 
        l1sel=l1.copy()
        l2=ksel[:,1] 
        l2sel=l2.copy()
        A1=ksel[:,2] 
        A1sel=A1.copy()
        A2=1-A1
        A2sel=A2.copy()
                
        S1 = ksel[:,2]*ksel[:,0]+(1-ksel[:,2])*ksel[:,1]
        K2 = -S1
        S2 = ksel[:,2]*ksel[:,0]**2+(1-ksel[:,2])*ksel[:,1]**2
        S3 = ksel[:,2]*ksel[:,0]**3+(1-ksel[:,2])*ksel[:,1]**3
        K1m = S1-S2/S1
        K1p = (S3*S1-S2**2)/S1/(S1**2-S2)
        #################

        
        P1=K1m/(K1m+K1p)
        P2=K1p/(K1m+K1p)

        resl= np.array([np.min(K1p),np.min(K1m),np.min(K2),np.min(P1),np.min(P2),np.min(l1sel),np.min(l2sel),np.min(A1sel),np.min(A2sel)])
        resl[0:5] = np.max(np.vstack([resl[0:5], np.zeros(5)]), axis=0)
        resh= np.array([np.max(K1p),np.max(K1m),np.max(K2),np.max(P1),np.max(P2),np.max(l1sel),np.max(l2sel),np.max(A1sel),np.max(A2sel)])
        
    else:

        res = [0, 0, 0, 0, 0, 0, sd[1],sd[0]]
        resh=np.zeros(5)
        resl=np.zeros(5)

    return [res, resl, resh]
