# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 12:16:22 2013

@author: Maxim
"""

import matplotlib.pyplot as plt
import numpy as np
from pyDOE import lhs, fullfact

def run_test_1():
    x = lhs(2,2)
    print x

    plt.figure(1)
    plt.plot(x[:,0],x[:,1],'ro')
    #plt.axis([-.5,1.5,-.5,1.5])
    plt.grid(True)
    plt.show()

if __name__=="__main__":
    run_test_1()