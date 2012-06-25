# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 21:19:41 2012

@author: Maxim
"""

import numpy as ny
import matplotlib.pyplot as plt
import flapGeom2
import afAnalysis
import ga2
import afLib

def costFcn(X,deflection):
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
    #deflection = 35
    node1[1,0] = 0.7 #flap chord ratio
    Mach = 0.16
    Re = 4e6
    alphaStart = -20
    alphaEnd = 10
    alphaStep = 1
    try:
        airfoil = flapGeom2.getFlap(airfoilPath, node1,node2, vectLenRatio1, theta, gap, overlap, deflection,zTE)
        polar = afAnalysis.polar(airfoil)
        #polar.calcJpolar(Mach,Re,alphaStart,alphaEnd,alphaStep,True)
        #f = -polar.CLmax
        #g = 5. - polar.ClmaxAlpha
    except: f = 100
    f = []
    h = []
    g = []
    
    wdir = r'D:\Documents\My Documents\1. Classes\5 Sem - Human Computer Interaction for MDO\term project\DoE'
    filename1 = ('flap overlap%.2f gap%.2f defl%.2f.txt'%(100*X[8],100*X[9],deflection))
    filename2 = ('flap overlap%.2f gap%.2f defl%.2f.iges'%(100*X[8],100*X[9],deflection))
    filepath1 = wdir + '\\' + filename1
    filepath2 = wdir + '\\' + filename2
    afLib.writeFlap(filepath1,airfoil)
    polar.saveAsIgs(filepath2,airfoil,True)
    
    plt.figure(1)
    plt.axis([0,1.2,-.5,.5])
    plt.grid(True)
    plt.hold(True)
    plt.plot(airfoil.mainSec[:,0],airfoil.mainSec[:,1])
    plt.plot(airfoil.flap[:,0],airfoil.flap[:,1])
    
    plt.show()
    
    return f,g,h
    
overlaps = ny.array([0.005,0.01,0.02])
gaps = ny.array([0.005, 0.01, 0.015])
deflections = ny.array([15.,25.,35])

for overlap in overlaps:
    for gap in gaps:
        for deflection in deflections:
            X = ny.array([0.839178,0.574781,0.823192,0.674157,0.291816,0.650752,0.297726,98.787865,overlap,gap])
            costFcn(X,deflection)