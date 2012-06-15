# -*- coding: utf-8 -*-
"""
Created on Fri May 25 10:22:32 2012

@author: Maxim
"""

import numpy as ny

def BPOcoef(order):
    K = ny.ones([order+1,1],int)
    for ii in range(1,order):
        K[ii] = factorial(order)/(factorial(ii)*factorial(order-ii))
    return K

def factorial(n):
    n = ny.abs(int(n)) 
    if n < 1: n = 1 
    if n == 1: 
        return 1 
    else: 
        return n * factorial(n - 1)

def stretching(Npts=50,betaLE=4,betaTE=2):
    step = 1/(float(Npts-1))
    u = ny.arange(0,1+step,step)
    x = ny.zeros([Npts,1])
    ii = 0
    for xx in x:
        u[ii] = (ny.exp(betaLE*u[ii])-1)/(ny.exp(betaLE)-1)
        x[ii] = 1 - (ny.exp(betaTE*(1-u[ii]))-1)/(ny.exp(betaTE)-1)
        ii+=1
    return x

def CSTairfoil(Au,Al,N,Npts):
    Npts = int(Npts)
    Au = [float(Atmp) for Atmp in Au]
    Al = [float(Atmp) for Atmp in Al]
    N =  [float(Ntmp) for Ntmp in N]
    
    order = ny.shape(Au)[0]-1
    K = BPOcoef(order)
    Up = ny.zeros([Npts,2])
    Lo = ny.zeros([Npts,2])
    x = stretching(Npts)
    for ii in range(Npts):
        ysUp,ysLo = 0,0
        yc = x[ii]**N[0]*(1-x[ii])**N[1]
        for jj in range(order+1):
            ysUp = ysUp+x[ii]**(jj)*(1-x[ii])**(order-jj)*K[jj]*Au[jj]
            ysLo = ysLo+x[ii]**(jj)*(1-x[ii])**(order-jj)*K[jj]*Al[jj]
        Up[ii,:] = [x[ii],ysUp*yc]
        Lo[ii,:] = [x[ii],ysLo*yc]
    return Up,Lo