# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 20:15:08 2013

@author: Maxim
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

def test_function_selection():
    def lofi(x):
        return (x-2)**2.0/2+2.0
    def hifi(x):
        return np.sin(x*np.pi/2.0)+(x-2.0)**2.0/2 +2.0
    
    
    rslt1 = minimize(lofi,3.0,method='SLSQP',tol=1e-20)
    rslt2 = minimize(hifi,3.0,method='SLSQP',tol=1e-20)
    print rslt1
    print rslt2
    
    x = np.linspace(0,5,50)
    
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax1.hold(True)
    ax1.plot(x,lofi(x),'b-')
    ax1.plot(x,hifi(x),'r-')
    ax1.plot(rslt1.x,lofi(rslt1.x),'bo')
    ax1.plot(rslt2.x,hifi(rslt2.x),'ro')
    plt.show()

if __name__=="__main__":
    test_function_selection()