# -*- coding: utf-8 -*-
"""
Created on Tue Jun 05 11:05:09 2012
convergence history plot
@author: Maxim
"""
import numpy as ny
import matplotlib.pyplot as plt
import time


HistFile = r'D:\Documents\My Documents\Dropbox\cfd analysis\Flap Design Optimization\cl-history(16deg)(flap25deg)(X0Y-2.5)'
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
    Yval = ny.append(Yval,float(SegLine[Yaxis]))

plt.interactive(False)
plt.figure(1)
plt.plot(Xval,Yval)
plt.hold(True)
#plt.title(Title)
#plt.xlabel(Xlabel)
#plt.ylabel(Ylabel)
plt.draw()

#time.sleep(5)

print 'done'