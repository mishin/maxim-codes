# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 10:53:29 2012

@author: Maxim
"""

import numpy as ny
import matplotlib.pyplot as plt
import flapGeom2
import afAnalysis

def sensAnalysis(func, lb, ub, Npt, dx, varNames=[]):
    step = (ub-lb)/(Npt-1)
    Xinit = (ub+lb)/2
    fval = ny.zeros([lb.shape[0],Npt])
    Xinc = ny.zeros([lb.shape[0],Npt])
    labels = ny.array(['pt0-X','pt1-Y','pt2-X','pt4-X','vect0','vect1','vect2','theta','overlap','gap'])
    grad = ny.zeros(lb.shape)
    for ii in range(lb.shape[0]):
        Xnew = ny.copy(Xinit)
        Xnew[ii] = lb[ii] - step[ii]
        for jj in range(Npt):
            Xnew[ii] = Xnew[ii] + step[ii]
            Xinc[ii,jj] = Xnew[ii]
            f,_,_ = func(Xnew)
            fval[ii,jj] = -f
        grad[ii] = (fval[ii,-1]-fval[ii,0]) / dx[ii]
        print grad
        plt.figure(ii+1)
        plt.plot(Xinc[ii,:],fval[ii,:],'o-')
        plt.grid(True)
        plt.axis([lb[ii],ub[ii],2.,3.])
        plt.ylabel('Clmax')
        plt.xlabel(('X%d: %s'%(ii,labels[ii])))
    plt.show()


def costFcn(X):
    airfoilPath = 'FinalOptimum20120720.dat'
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
    Mach = 0.0955
    Re =  2285283
    alphaStart = -10
    alphaEnd = 15
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
    xOpt = ny.array([0.843078,0.365108,0.804955,0.677721,0.527504,0.643172,0.788588,97.123078,0.036051,0.012782])
    lb = xOpt - ny.array([0.025,0.025,0.025,0.01,0.025,0.025,0.025,2,0.01, 0.01])
    ub = xOpt + ny.array([0.025,0.025,0.025,0.01,0.025,0.025,0.025,2,0.01, 0.01])
    #lb = xOpt - ny.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,1,0.005, 0.005])
    #ub = xOpt + ny.array([0.01,0.01,0.01,0.01,0.01,0.01,0.01,1,0.005, 0.005])
    #lb1 = ny.array([0.71, 0.1, 0.71, 0.65, 0.4, 0.4, 0.4, 60 , -0.04, 0.005])
    #ub1 = ny.array([0.9 , 0.9, 0.9 , 0.7 , 0.9, 0.9, 0.9, 150,  0.04, 0.05])
    dx = (ub-lb)
    sensAnalysis(costFcn,lb,ub,5,dx)
    
testFcn(costFcn)