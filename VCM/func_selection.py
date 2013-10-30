# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 20:15:08 2013

@author: Maxim
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
from scipy.optimize import minimize

def forrester_hi(x):
    return (5.0*x-2.0)**2.0*np.sin(12.*x-4.)
def forrester_low(x):
    A = 0.5
    B = 10
    C = -5.
    return A*forrester_hi(x)+B*(x-0.5)-C

def test_function_selection1():
    def lofi(x):
        return (x-2.2)**2.0/2+2.5
    def hifi(x):
        return np.sin(x*np.pi)+(x-2.0)**2.0/2 +2.0
    
    lofi = forrester_low
    hifi = forrester_hi
    rslt1 = minimize(lofi,0.1,method='SLSQP',bounds=[(0,1)],tol=1e-20)
    rslt2 = minimize(hifi,0.1,method='SLSQP',bounds=[(0,1)],tol=1e-20)
    print rslt1
    print rslt2
    
    x = np.linspace(0,1,50)
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax1.hold(True)
    ax1.plot(x,lofi(x),'b-')
    ax1.plot(x,hifi(x),'r-')
    ax1.plot(rslt1.x,lofi(rslt1.x),'bo')
    ax1.plot(rslt2.x,hifi(rslt2.x),'ro')
    plt.show()


def test_function_selection2():
    lb = np.array([0.1,0.1])
    ub = np.array([10.,10.])
    dx = 0.1
    #f = lambda x: (np.sin(2.*np.pi*x[0])**3 * np.sin(2.*np.pi*x[1]))/(x[0]**3*(x[0]+x[1]))
    f = lambda x: 4.*x[0]**2.+x[1]**3+x[0]*x[1]
    g1 = lambda x: x[0]**2-x[1]+1
    g2 = lambda x: 1-x[0]+(x[1]-4.)**2.
    
    x = y = np.arange(lb[0],ub[0],dx)
    X, Y = np.meshgrid(x,y)
    zs1 = np.array([f([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    zs2 = np.array([g1([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    zs3 = np.array([g2([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    Z1 = zs1.reshape(X.shape)
    Z2 = zs2.reshape(X.shape)
    Z3 = zs3.reshape(X.shape)
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.hold(True)
    cs1 = ax1.contour(X,Y,Z1,colors='r')
    cs2 = ax1.contour(X,Y,Z2,colors='y',levels=[0])
    cs3 = ax1.contour(X,Y,Z3,colors='g',levels=[0])
    plt.clabel(cs1, fontsize=9, inline=1)
    plt.clabel(cs2, fontsize=9, inline=1)
    plt.clabel(cs3, fontsize=9, inline=1)
    plt.show()

def test_function_selection3():
    def currin(x):
        val1 = 1.0-np.exp(-0.5/x[1])
        val2 = 2300.*x[0]**3 + 1900*x[0]**2 + 2092*x[0]+60
        val3 = 100.*x[0]**3+500*x[0]**2+4*x[0]+20.
        return val1*val2/val3
    
    def currin_low(x):
        x1 = x[0]
        x2 = x[1]
        maxarg = max([0, x2+0.05])
        val1 = currin([x1+1/20, x2+1/20])
        val2 = currin([x1+1/20, maxarg])
        val3 = currin([x1-1/20, x2+1/20])
        val4 = currin([x1-1/20, maxarg])
        return 0.25*(val1+val2+val3+val4)

    lb = [0.05,0.05]
    ub = [1,1]
    dx=0.02
    x = y = np.arange(lb[0],ub[0],dx)
    X, Y = np.meshgrid(x,y)
    zs1 = np.array([currin([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    zs2 = np.array([currin_low([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    Z1 = zs1.reshape(X.shape)
    Z2 = zs2.reshape(X.shape)
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.hold(True)
    cs1 = ax1.contour(X,Y,Z1,colors='r')
    cs2 = ax1.contour(X,Y,Z2,colors='k')
    plt.clabel(cs1, fontsize=9, inline=1)
    plt.clabel(cs2, fontsize=9, inline=1)
    plt.show()

def test_function_selection4():
    f = lambda x: np.exp(x[0]/3)+np.exp(x[1]/5)-x[0]
    #g2 = lambda x: g1(x)-0.2*(x[0]*x[1])
    g1 = lambda x: x[0]**2.*x[1]/20.-1
    #g1 = lambda x: x[0]+x[1]
    g2 = lambda x: (x[0]+x[1]-5)**2./30 + (x[0]-x[1]-12)**2./120 - 1
    g3 = lambda x: (x[0])**2.*(x[1])/20+(x[1]-x[0])/5.-2
    #g3 = lambda x: 80./(x[0]**2.+8*x[1]+5)-1
    lb = [1.0, 10.0]
    ub = [10., 10.]
    dx = 0.5
    rslt = minimize(f,[5,5],bounds=((0,10),(0,10)),constraints=({'type':'ineq','fun':g2},{'type':'ineq','fun':g1}),method='SLSQP')
    print rslt
    x = y = np.arange(lb[0],ub[0],dx)
    X, Y = np.meshgrid(x,y)
    zs1 = np.array([f([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    zs2 = np.array([g1([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    zs3 = np.array([g2([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    zs4 = np.array([g3([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    Z1 = zs1.reshape(X.shape)
    Z2 = zs2.reshape(X.shape)
    Z3 = zs3.reshape(X.shape)
    Z4 = zs4.reshape(X.shape)
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
#    ax1 = Axes3D(fig)
    ax1.hold(True)
#    ax1.plot_wireframe(X,Y,Z1,color='b')
#    ax1.plot_wireframe(X,Y,Z2,color='r')
#    ax1.plot_wireframe(X,Y,Z3,color='y')
    cs1 = ax1.contour(X,Y,Z1,colors='b')
    cs2 = ax1.contour(X,Y,Z2,colors='r',levels=[0])
    cs3 = ax1.contour(X,Y,Z3,colors='y',levels=[0])
    cs4 = ax1.contour(X,Y,Z4,colors='g')
    plt.clabel(cs1, fontsize=9, inline=1)
    plt.clabel(cs2, fontsize=9, inline=1)
    plt.clabel(cs3, fontsize=9, inline=1)
    ax1.plot(rslt.x[0],rslt.x[1],'ro')
    plt.clabel(cs4, fontsize=9, inline=1)
    ax1.axis([0,10,0,10])
    plt.show()

def test_function_selection5():
    forrester = lambda x: (5.0*x-2.0)**2.0*np.sin(12.*x-4.)
    def f(x):
        x = np.array(x)*0.1
        return forrester(np.linalg.norm(x)/1.2)+5.*(x[0]+x[1])
    f2 = lambda x: np.exp(x[0]/3)+np.exp(x[1]/5)-x[0]
    lb = [0, 0.0]
    ub = [10., 10.]
    dx = 0.5
    x = y = np.arange(lb[0],ub[0],dx)
    X, Y = np.meshgrid(x,y)
    zs1 = np.array([f([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    Z1 = zs1.reshape(X.shape)
    zs2 = np.array([f2([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    Z2 = zs2.reshape(X.shape)
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    #ax1 = Axes3D(fig)
    #ax1.plot_wireframe(X,Y,Z1)
    cs1 = ax1.contour(X,Y,Z1,levels=np.arange(0,20,1))
    plt.clabel(cs1)
    
    fig2 = plt.figure(2)
    ax2 = Axes3D(fig2)
    ax2.hold(True)
    ax2.plot_wireframe(X,Y,Z1,color='y')
    ax2.plot_wireframe(X,Y,Z2,color='r')
    plt.show()

if __name__=="__main__":
    test_function_selection5()
    #test_function_selection2()