# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 17:24:53 2013

@author: Maxim
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class Kriging1D:
    def __init__(self,x,y):
        self.x = x
        self.y = y
        self._n = len(self.x)
        self._create_model()
    
    def _create_model(self):
        self.k = self._n*(self._n-1)/2
        rho = np.zeros(self.k)
        nu = np.zeros(self.k)
        i,k = 1,0
        for xi in self.x[1:]:
            j = 0
            for xj in self.x[:i]:
                rho[k] = np.linalg.norm([xi-xj])
                nu[k] = (self.y[i]-self.y[j])**2.0
                j += 1                
                k += 1
            i += 1
        popt,pcov = curve_fit(self._dispersion_func,rho,nu)
        beta, gamma = popt[0], popt[1]
        self._beta = popt[0]
        self._gamma = popt[1]
        self._calc_model_matrix()
        plt.figure(2)
        plt.hold(True)
        plt.plot(rho,nu,'go')
        _r = np.linspace(0,1,50)
        plt.plot(_r,self._dispersion_func(_r,beta,gamma),'r-')
        plt.plot()

    def _d(self,a,b):
        rho = np.linalg.norm([a-b])
        return self._dispersion_func(rho,self._beta,self._gamma)

    def _dispersion_func(self,rho,beta,gamma):
        return beta*(1.0-np.exp(-gamma*rho**2.0))
    
    def _calc_model_matrix(self):
        A = np.ones([self._n+1,self._n+1])
        A[-1,-1] = 0.0
        for i in range(self._n):
            A[i,i] = self._dispersion_func(0,self._beta,self._gamma)
        i,k = 1,0
        for xi in self.x[1:]:
            j = 0
            for xj in self.x[:i]:
                A[i,j] = self._d(xi,xj)
                A[j,i] = self._d(xi,xj)
                j += 1                
                k += 1
            i += 1
        self.modelMatrix = A
    
    def __call__(self,x):
        d = np.ones(self._n+1)
        for i,_x in enumerate(self.x):
            d[i] = self._d(_x,x)
        weight = np.linalg.solve(self.modelMatrix,d)
        y = 0
        for w,_y in zip(weight[:-1],self.y):
            y += w*_y
        D = np.dot(weight,d)
        return y,D





def func(x):
    return (6.*x-2.)**2*np.sin(12.*x-4)

def run_test():
    #x = np.array([0.0,.1,.2,.4,.75,.9,1.])
    x = np.array([0.0,0.2,0.8,1.0])
    y = np.array([func(_x) for _x in x])
    xp = np.linspace(0,1,100)
    yp = func(xp)

    krig = Kriging1D(x,y)
    xk = np.linspace(0,1,100)
    yk = np.zeros(len(xk))
    dk = np.zeros(len(xk))
    for i,xx in enumerate(xk):
        yk[i],dk[i] = krig(xx)
    
    print sum(dk)
    plt.figure(1)
    plt.title('Kriging test')
    plt.hold(True)
    plt.plot(x,y,'ro')
    plt.plot(xp,yp,'b-')
    plt.plot(xk,yk,'r-')
    plt.fill_between(xk,yk+0.5*dk,yk-0.5*dk,color='#cccccc')
    plt.grid(True)
    plt.axis([0,1,-10,20])
    plt.legend(['sample point','exact function','kriging','confidence interval'],'upper left')
    plt.show()


if __name__=="__main__":
    run_test()