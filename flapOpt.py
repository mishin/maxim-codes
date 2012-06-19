# -*- coding: utf-8 -*-
"""
Created on Sat Jun 16 13:57:07 2012

@author: Maxim
"""

import numpy as ny
import matplotlib.pyplot as plt
import flapGeom
import afAnalysis
import ga2

def costFcn(X):
    airfoilPath = 'GA37A315mod.dat'
    nodes1 = ny.zeros([3,2])
    nodes2 = ny.zeros([3,2])
    vectLen1 = ny.zeros([3])
    
    nodes1[0,0] = X[0]
    
    nodes1[1,1] = X[1]
    nodes1[2,0] = X[2]
    nodes2[0,0] = X[3]
    vectLen1[0] = X[4]
    vectLen1[1] = X[5]
    vectLen1[2] = X[6]
    theta       = X[7]
    overlap     = X[8]
    gap         = X[9]
    
    teThickness = 0.002
    deflection = 35
    nodes1[1,0] = 0.7 #flap chord ratio
    Mach = 0.16
    Re = 4e6
    alphaStart = -20
    alphaEnd = 20
    alphaStep = 1
#    try:
    airfoil = flapGeom.getFlap(airfoilPath, nodes1,nodes2, vectLen1, theta, gap, overlap, deflection,teThickness)
    polar = afAnalysis.polar(airfoil)
    polar.calcJpolar(Mach,Re,alphaStart,alphaEnd,alphaStep,True)
    f = -polar.CLmax
    g = 5. - polar.ClmaxAlpha
#    except: f,g = 100,100
    h = []
    return f,g,h
    
    
lb = ny.array([0.71, 0.1, 0.71, 0.65, 0.01, 0.1, 0.01, 60, -0.03, 0.005 ])
ub = ny.array([0.9, 0.9, 0.85, 0.71,0.15,0.8,0.15, 150, 0.03,0.03])
opt = ga2.gaOptions(lb,ub)
opt.MaxIterations(400,200)
opt.sigma = 0.75
opt.HistFile = 'flapOptHist.dat'


fBest,xBest,fHist,xHist, Iter = ga2.gaMain(costFcn,opt)
