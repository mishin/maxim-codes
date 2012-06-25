# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 10:53:29 2012

@author: Maxim
"""

import numpy as ny
import matplotlib.pyplot as plt
import flapGeom2
import afAnalysis

def sensAnalysis(func, lb, ub, Npt,varNames=[]):
    step = (ub-lb)/(Npt-1)
    Xinit = (ub+lb)/2
    fval = ny.zeros([lb.shape[0],Npt])
    Xinc = ny.zeros([lb.shape[0],Npt])
    labels = ny.array(['pt0-X','pt1-Y','pt2-X','pt4-X','vect0','vect1','vect2','theta','overlap','gap'])
    
    for ii in range(lb.shape[0]):
        Xnew = ny.copy(Xinit)
        Xnew[ii] = lb[ii] - step[ii]
        for jj in range(Npt):
            Xnew[ii] = Xnew[ii] + step[ii]
            Xinc[ii,jj] = Xnew[ii]
            f,_,_ = func(Xnew)
            fval[ii,jj] = -f
        plt.figure(ii+1)
        plt.plot(Xinc[ii,:],fval[ii,:],'o-')
        plt.grid(True)
        plt.axis([lb[ii],ub[ii],2.,3.])
        plt.ylabel('Clmax')
        plt.xlabel(('X%d: %s'%(ii,labels[ii])))
    plt.show()


def costFcn(X):
    airfoilPath = 'GA37A315mod.dat'
    node1 = ny.zeros([3,2])
    node2 = ny.zeros([3,2])
    vectLenRatio1 = ny.zeros([3])
    
    node1[0,0] = X[0]
    node1[1,1] = X[1]
    node1[2,0] = X[2]
    node2[0,0] = X[3]
    vectLenRatio1[0] = X[4]
    vectLenRatio1[1] = X[5]
    vectLenRatio1[2] = X[6]
    theta       = X[7]
    overlap     = X[8]
    gap         = X[9]
    
    zTE = 0.002
    deflection = 35
    node1[1,0] = 0.7 #flap chord ratio
    Mach = 0.16
    Re = 4e6
    alphaStart = -20
    alphaEnd = 20
    alphaStep = 1
    try:
        airfoil = flapGeom2.getFlap(airfoilPath, node1,node2, vectLenRatio1, theta, gap, overlap, deflection,zTE)
        polar = afAnalysis.polar(airfoil)
        polar.calcJpolar(Mach,Re,alphaStart,alphaEnd,alphaStep,True)
        f = -polar.CLmax
        #g = 5. - polar.ClmaxAlpha
    except: f = 100
    h = []
    g = []
    return f,g,h
    


def testFcn(func):
    #lb = ny.array([0.75, 0.1, 0.71, 0.65, 0.1, 0.4, 0.1, 80 , -0.04, 0.005 ])
    #ub = ny.array([0.85 , 0.5, 0.85 , 0.7 , 0.5, 0.8, 0.5, 140,  0.04, 0.02])
    #0.827618	0.369286	0.805817	0.669824	0.333747	0.602765	0.273544	108.994754	0.030467	0.007191
    lb = ny.array([0.727618,0.269286,0.705817,0.569824,0.233747,0.502765,0.173544,88.994754,0.01,0.001])
    ub = ny.array([0.927618,0.469286,0.905817,0.769824,0.433747,0.702765,0.373544,128.994754,0.05,0.014191])
    sensAnalysis(costFcn,lb,ub,7)
    
testFcn(costFcn)