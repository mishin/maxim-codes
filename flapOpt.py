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
        g = 5. - polar.ClmaxAlpha
    except: f,g = 100,100
    h = []
    return f,g,h
    
lb = ny.array([0.71, 0.1, 0.71, 0.65, 0.1, 0.1, 0.1, 60 , -0.04, 0.005 ])
ub = ny.array([0.9 , 0.9, 0.9 , 0.7 , 0.9, 0.9, 0.9, 150,  0.04, 0.05])
opt = ga2.gaOptions(lb,ub)
opt.MaxIterations(400,200)
opt.sigma = 0.75
opt.HistFile = 'flapOptHist.dat'

fBest,xBest,fHist,xHist, Iter = ga2.gaMain(costFcn,opt)
