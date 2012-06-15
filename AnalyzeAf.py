# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 13:20:16 2012

@author: Maxim
"""

import CST
import afAnalysis
import numpy as ny
import matplotlib.pyplot as plt

xBest = ny.array([0.184069,0.173145,0.374482,0.311078,-0.178542,-0.089179,-0.224291])
    
    
Au = xBest[0:3]
Al = xBest[4:6]
Al = ny.hstack([-Au[0],Al])
Up,Lo = CST.CSTairfoil(Au,Al,ny.array([0.5,1.0]),80)

polar = AfAnalysis.polar(Up,Lo)
polar.calcXpolar(0.16,4.17e6,-10,22,0.5,50, True,True)

plt.figure(1)
plt.plot(polar.alpha, polar.CL)
plt.grid(True)
plt.xlabel('alpha')
plt.ylabel('lift coefficient')
plt.title('Optimum airfoil polar')
plt.figure(2)
plt.plot(polar.CD, polar.CL)
plt.grid(True)
plt.xlabel('drag coefficient')
plt.ylabel('lift coefficient')
plt.title('Optimum airfoil polar')
plt.show()