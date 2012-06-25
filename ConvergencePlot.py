# -*- coding: utf-8 -*-
"""
Created on Tue Jun 05 11:05:09 2012
convergence history plot
@author: Maxim
"""
import numpy as ny
import matplotlib.pyplot as plt
import time


HistFile = r'D:\Documents\My Documents\1. Classes\5 Sem - Human Computer Interaction for MDO\term project\flapOptHist.dat'
Xaxis = 0
Yaxis = 1

Hfile = open(HistFile,'rt')
Title = Hfile.readline()
Labels = Hfile.readline()
lines = Hfile.readlines()
Hfile.close()

Xlabel = Labels.split()[Xaxis]
Ylabel = Labels.split()[Yaxis]
Xval, Yval = [],[]

for line in lines:
    SegLine = line.split()
    Xval = ny.append(Xval,float(SegLine[Xaxis]))
    Yval = ny.append(Yval,-float(SegLine[Yaxis]))

plt.interactive(False)
plt.figure(1)
plt.plot(Xval,Yval)
plt.hold(True)
plt.title('Convergence history')
plt.xlabel('iterations')
plt.ylabel('c_l_max')
plt.axis([0,210,2.5,3.0])
plt.grid(True)
plt.show()

print 'done'