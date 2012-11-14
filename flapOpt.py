# -*- coding: utf-8 -*-
"""
Created on Sat Jun 16 13:57:07 2012

@author: Maxim
"""

import numpy as ny
import matplotlib.pyplot as plt
import flapGeom2
import afAnalysis
import ga2

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
    
def final_result(X):
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
    airfoil = flapGeom2.getFlap(airfoilPath, node1,node2, vectLenRatio1, theta, gap, overlap, deflection,zTE)
    polar = afAnalysis.polar(airfoil)
    polar.calcJpolar(Mach,Re,alphaStart,alphaEnd,alphaStep,True)
    print polar.CLmax
    afAnalysis.afLib.writeFlap('optimized flap.dat', airfoil)
    plt.figure(2)
    plt.plot(airfoil.mainSec[:,0],airfoil.mainSec[:,1])
    plt.hold(True)
    plt.grid(True)
    plt.plot(airfoil.flap[:,0],airfoil.flap[:,1])
    plt.axis([0,1.2,-.5,.5])
    plt.show()


lb = ny.array([0.71, 0.1, 0.71, 0.65, 0.4, 0.4, 0.4, 60 , -0.04, 0.005])
ub = ny.array([0.9 , 0.9, 0.9 , 0.7 , 0.9, 0.9, 0.9, 150,  0.04, 0.05])
opt = ga2.gaOptions(lb,ub)
opt.MaxIterations(25,20)
opt.sigma = 1.
opt.HistFile = 'flapHistory2.dat'

fBest,xBest,fHist,xHist, Iter = ga2.gaMain(costFcn,opt)

final_result(xBest)