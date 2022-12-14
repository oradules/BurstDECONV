#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 18:14:04 2021

@author: mdouaihy
"""
import numpy as np
from ecdfEstimate import ecdf_bounds
import matplotlib.pyplot as plt


def fit3M1Parameters(store,dt,dtc,sd,dirwrite,name):
    
    fsz=16 #figure size
    
    
    ### finding the optimal finite parameters
    ind = np.where(np.max(np.abs(store[:,0:3].imag )  ,axis=1) <1e-10)[0] #takes the part where we have the imaginary part of lambda i almost 0
    
    ### this will reshape the matrix store in a way that we have an assending order of the last column objective 
    StoreAssirted=store[np.argsort(store[ind, -1])]
    
    
    objmin = np.min(store[ind,-1]).real #minimum of the objectives for all the previous lambdas
    
    indmin=np.argmin(store[ind,-1]) #finding the positions of the lambda where we have the lowest objective
    imin = ind[indmin] #finding the index of the positions of the lambda where we have the lowest objective
    
    kmin = store[imin,0:5].real# taking the lambda i's and A i's that are the fittest
    censmin=store[imin,-1]
    # compute 5 rates k1p,m k2p,m k3 from the 5 parameters
    l1 = kmin[0]
    l2 = kmin[1]
    l3 = kmin[2]
    A1 = kmin[3]
    A2 = kmin[4]
    A3 = 1-A1-A2

    
    L1=l1+l2+l3
    L2=l1*l2+l1*l3+l2*l3
    L3=l1*l2*l3
    S1=A1*l1+A2*l2+A3*l3
    S2=A1*l1**2+A2*l2**2+A3*l3**2
    S3=A1*l1**3+A2*l2**3+A3*l3**3
    
    # model M1
    k1p=-L3 *(S1  **2-S2) /(S2  **2-S1 *S3)   # k1p
    k2p=-(S2  **2-S1 *S3) /S1 /(S1  **2-S2)   # k2p
    k3=-S1   # k3
    k2m=(S1  **2-S2) /S1   # k2m
    k1m=-A1 *A2 *A3 *(l1-l2)  **2 *(l1-l3)  **2 *(l2-l3)  **2 *S1 /(S1  **2-S2) /(S2  **2-S1 *S3)   # k1m
    p1=k1m*k2m/(k1p*k2m+k1m*k2m+k1p*k2p) 
    p2=k1p*k2m/(k1p*k2m+k1m*k2m+k1p*k2p) 
    p3=k1p*k2p/(k1p*k2m+k1m*k2m+k1p*k2p) 
    
    
    
   
    #optimal 
    parameters=np.array([k1p,k1m,k2p,k2m,k3])
    indObjAssorted=1
    
    while np.isnan(parameters).any() or (parameters <= 0).any():
        try:
            print('ind of obj', indObjAssorted-1)
            print('parameters',parameters)
            
            objmin = StoreAssirted[indObjAssorted,-1] # np.min(store[ind,-1]).real #minimum of the objectives for all the previous lambdas
            imin =indObjAssorted# ind[indmin] #finding the index of the positions of the lambda where we have the lowest objective
            kmin = StoreAssirted[imin,0:5].real# taking the lambda i's and A i's that are the fittest
            censmin=StoreAssirted[imin,-1]
            
            # compute 5 rates k1p,m k2p,m k3 from the 5 parameters
            l1 = kmin[0]
            l2 = kmin[1]
            l3 = kmin[2]
            A1 = kmin[3]
            A2 = kmin[4]
            A3 = 1-A1-A2
            L1=l1+l2+l3
            L2=l1*l2+l1*l3+l2*l3
            L3=l1*l2*l3
            S1=A1*l1+A2*l2+A3*l3
            S2=A1*l1**2+A2*l2**2+A3*l3**2
            S3=A1*l1**3+A2*l2**3+A3*l3**3
        
            # model M1
            k1p=-L3 *(S1  **2-S2) /(S2  **2-S1 *S3)   # k1p
            k2p=-(S2  **2-S1 *S3) /S1 /(S1  **2-S2)   # k2p
            k3=-S1   # k3
            k2m=(S1  **2-S2) /S1   # k2m
            k1m=-A1 *A2 *A3 *(l1-l2)  **2 *(l1-l3)  **2 *(l2-l3)  **2 *S1 /(S1  **2-S2) /(S2  **2-S1 *S3)   # k1m
            p1=k1m*k2m/(k1p*k2m+k1m*k2m+k1p*k2p) 
            p2=k1p*k2m/(k1p*k2m+k1m*k2m+k1p*k2p) 
            p3=k1p*k2p/(k1p*k2m+k1m*k2m+k1p*k2p) 
            indObjAssorted=indObjAssorted+1
            parameters=np.array([k1p,k1m,k2p,k2m,k3])
            
        except:
            print('no parameters found')
            parameters=np.array([1000,1000,1000,1000,1000])
        
    #######################################################
    
    #### computing Uncertainty interval
    overflow = 1 #help us set the Uncertainty interval

    ind = np.where( (store[:,-1] < (1+overflow)*objmin) & (np.max(np.abs(store[:,0:3].imag   ),axis=1) <1e-10)) [0] #taking suboptimal objectives such that they are <2*optimal objective and they satisfy the fact that the imaginary part of lambda is almost 0
    ksel = store[ind,0:5].real #taking an array of lambda i's and Ai's where we have the the objective function less than <2*times the minimum objective function
    

    #### plot survival function
        
    if censmin:
        xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc),np.append( np.zeros(len(dt)), np.ones(len(dtc)) ) )
    else:
        xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc))


    h = plt.figure(70)
    plt.semilogy(xs, 1-fs, 'o', color='r', mfc = 'none', markersize=9, linestyle='', label = 'S funct from data') #empirical function
    plt.semilogy(xs, 1-flo, '--r') # lower confidence
    plt.semilogy(xs, 1-fup, '--r') # upper confidence
    pred = kmin[3]*np.exp(kmin[0]*xs)+kmin[4]*np.exp(kmin[1]*xs)+(1-kmin[3]-kmin[4])*np.exp(kmin[2]*xs) 
    plt.semilogy(xs, pred, 'k', linewidth = 2, label = '3 States fit S funct') # predicted 2 exp
    plt.xlim(0,250)
    plt.ylim(1e-6, 1)
    plt.xlabel('Time [s]',fontsize=fsz)
    plt.ylabel('Survival function',fontsize = fsz)
    plt.legend()
    plt.title(r'$k_1^-$='+str("%.1g" % k1m)+' $k_1^+$='+str("%.1g" % k1p) +' $k_2^-$='+str("%.1g" % k2m)+' $k_2^+$='+str("%.1g" % k2p)+' $k_3=$'+str("%.1g" % k3))
    figfile = dirwrite+'/Fit_3statesM1_'+name+'.pdf' 
    h.savefig(figfile)
    plt.close()
    
    ##############################################
    
    ## KS test
    thr =20
    ind = np.where(xs>thr)
    dist = np.max(np.abs(fs[ind]-1+pred[ind]))
    nsample=len(dt[dt>thr])
    N=10
    
    c=np.sqrt(nsample)*dist
    r=np.arange(1,N+1)
    aa = 2*sum(np.exp(-2*c**2*r**2)*((-1)**(r-1)))
    h=plt.figure(80)
    plt.plot(xs,fs,'kx')
    plt.plot(xs,1-pred,'ro', mfc = 'none')
    plt.xlabel('Time [s]',fontsize=fsz)
    plt.ylabel('CDF',fontsize=fsz)
    figfile=dirwrite+'/Fit_CDF_3statesM1_'+name+'.pdf'
    h.savefig(figfile)
    plt.close(h)
    
    #################



    #optimal 

    res= np.array([k1p,k1m,k2p,k2m,k3,p1,p2,p3,l1,l2,l3,A1,A2,A3,S1,objmin,aa,sd[1],sd[0]]) # optimum

    # compute interval
    ##################################################
    l1=ksel[:,0]
    l2=ksel[:,1] 
    l3=ksel[:,2] 
    A1=ksel[:,3] 
    A2=ksel[:,4] 
    A3=1-A1-A2 
    L1=l1+l2+l3 
    L2=l1 *l2+l1 *l3+l2 *l3 
    L3=l1 *l2 *l3 
    S1=A1 *l1+A2 *l2+A3 *l3 
    S2=A1 *l1  **2+A2 *l2  **2+A3 *l3  **2 
    S3=A1 *l1 **3+A2 *l2 **3+A3 *l3 **3 


    #### model M1
    K1p=-L3 *(S1 **2-S2) /(S2 **2-S1 *S3)  ### k1p
    K2p=-(S2 **2-S1 *S3) /S1 /(S1 **2-S2)  ### k2p
    K3=-S1.copy()  #### k3
    K2m=(S1 **2-S2) /S1  ### k2m
    K1m=-A1 *A2 *A3 *(l1-l2) **2 *(l1-l3) **2 *(l2-l3) **2 *S1 /(S1 **2-S2) /(S2 **2-S1 *S3)  ### k1m
    P1=K1m *K2m /(K1p *K2m+K1m *K2m+K1p *K2p) 
    P2=K1p *K2m /(K1p *K2m+K1m *K2m+K1p *K2p) 
    P3=K1p *K2p /(K1p *K2m+K1m *K2m+K1p *K2p) 

    resl= np.array([np.nanmin(K1p),np.nanmin(K1m),np.nanmin(K2p),np.nanmin(K2m),np.nanmin(K3),np.nanmin(P1),np.nanmin(P2),np.nanmin(P3),
                        np.nanmin(l1),np.nanmin(l2),np.nanmin(l3),np.nanmin(A1),np.nanmin(A2),np.nanmin(A3),np.nanmin(S1)])

    resl[0:8] = np.max(np.vstack([resl[0:8], np.zeros(8)]), axis=0)


    resh= np.array([np.nanmax(K1p),np.nanmax(K1m),np.nanmax(K2p),np.nanmax(K2m),np.nanmax(K3),np.nanmax(P1),np.nanmax(P2),np.nanmax(P3),
                    np.nanmax(l1),np.nanmax(l2),np.nanmax(l3),np.nanmax(A1),np.nanmax(A2),np.nanmax(A3),np.nanmax(S1)])


    return [res,resl,resh]
