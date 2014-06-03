# -*- coding: utf-8 -*-
"""
Created on Tue Jun 03 21:15:26 2014

@author: Maxim
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


def f(x):
    f=x[0]+x[1]
    return f
def g1(x):
    C1=x[0]**2*x[1]/20.0-1.
    return C1

def g2(x):
    C2=(x[0]+x[1]-5.0)**2/30.+(x[0]-x[1]-12.0)**2/120.0-1.
    return C2
    
def g3(x):
    C3=80.0/(x[0]**2+8.0*x[1]+5.0)-1
    return C3

def g2l(x):
    return 0.05*x[0]+x[1]-3.0

delta = 0.1
x = np.arange(0, 10, delta)
y = np.arange(0, 10, delta)
X, Y = np.meshgrid(x, y)

Z1 = np.array([g1([xx,yy]) for xx,yy in zip(X,Y)])
Z2 = np.array([g2([xx,yy]) for xx,yy in zip(X,Y)])
Z3 = np.array([g3([xx,yy]) for xx,yy in zip(X,Y)])
Z4 = np.array([g2l([xx,yy]) for xx,yy in zip(X,Y)])
Z1 = np.reshape(Z1,X.shape)
Z2 = np.reshape(Z2,X.shape)
Z3 = np.reshape(Z3,X.shape)
Z4 = np.reshape(Z4,X.shape)

cnstr = ({'type':'ineq','fun':g1},
         {'type':'ineq','fun':g2l},
         {'type':'ineq','fun':g3},)
rslt = minimize(f,[5,5],method='SLSQP',constraints=cnstr,bounds=((0,10),(0,10)))

plt.figure(1)
plt.hold(True)
plt.contour(X,Y,Z1,levels=[0])
plt.contour(X,Y,Z2,levels=[0],colors=['r'])
plt.contour(X,Y,Z3,levels=[0])
plt.contour(X,Y,Z4,levels=[0])
plt.plot(rslt.x[0],rslt.x[1],'ro')
plt.show()