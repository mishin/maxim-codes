# -*- coding: utf-8 -*-
"""
Created on Wed May 30 13:53:17 2012
af optimization test
@author: Maxim
"""

import numpy as ny
import ga2
import CST
import afAnalysis
import matplotlib.pyplot as plt
from datetime import datetime
    
def CostFcn(X):
    Au = X[0:3]
    Al = X[4:6]
    Al = ny.hstack([-Au[0],Al])
    Mach = 0.16
    Re = 4.17e6
    alphaStart = -10
    alphaEnd = 20
    alphaStep = 1.0
    clreq = 0.2
    
    Up,Lo = CST.CSTairfoil(Au,Al,ny.array([0.5,1.0]),80)
    tc = ny.max(Up[:,1]-Lo[:,1])
    
    g = ny.zeros([4])
    
    polar = AfAnalysis.polar(Up,Lo,'XfoilAnalysis')
    polar.calcXpolar(Mach,Re,alphaStart,alphaEnd,alphaStep,50,False)
    
    g[0] = 0.13 - tc
    g[1] = tc - 0.18
    g[2] = 1.6 - polar.CLmax
    f,alpha = polar.cdAtcl(clreq)
    g[3] = -alpha

    h = []
    return f, g,h

StartTime = datetime.now()


bestAfpath = r'D:\Documents\My Documents\1. Classes\5 Sem - Optimization\term project\result 20120610'
resultFile = bestAfpath + '\\result1.dat'
airfoilFile = bestAfpath + '\\airfoil1.dat'
historyFile = bestAfpath + '\\history.dat'

lb = ny.array([0.08,0.001,0.001,0.001,-0.4   ,-0.4   ,-0.4])
ub = ny.array([0.4 ,0.4  ,0.4  ,0.4  , -0.001, -0.001,-0.001])

opt = ga2.gaOptions(lb,ub)
#opt.PopSize = 5
opt.CorrNan = True
opt.MaxIterations(200,40)
opt.HistFile = historyFile
fBest,xBest,fHist,xHist,Iter = ga2.gaMain(CostFcn, opt)
print fBest, xBest

Au = xBest[0:3]
Al = xBest[4:6]
Al = ny.hstack([-Au[0],Al])
Up,Lo = CST.CSTairfoil(Au,Al,ny.array([0.5,1.0]),80)
f,g,h = CostFcn(xBest)

TimeConsumed = datetime.now()-StartTime

print fBest, xBest
print 'f:   ',f,'\ng:    ',g
print 'Time consumed:', TimeConsumed


AfAnalysis.AfLib.writeAirfoil(airfoilFile, 'optimized airfoil', Up,Lo)

result = open(resultFile,'wt')
result.write('xBest:\n')
for x in xBest:
    result.write('%.6f\t'%x)

result.write('\nf: %.6f\ng:'%f)
for gg in g:
    result.write('%.6f\t'%gg)
result.write('\n time:')
result.write(str(datetime.now()))
result.close()

plt.figure(1)
plt.plot(range(Iter),fHist,'o-')
plt.grid(True)
plt.xlabel('generation')
plt.ylabel('function value')
plt.title('convergence history')

plt.figure(2)
plt.plot(Up[:,0], Up[:,1])
plt.hold(True)
plt.plot(Lo[:,0], Lo[:,1])
plt.grid(True)
plt.axis([0,1,-0.4,0.4])
plt.show()

