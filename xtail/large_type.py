# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 23:13:19 2013

@author: Maxim
"""
from numpy import array, linspace, cos, sin, pi, hstack, vstack, flipud, transpose, concatenate
import matplotlib.pyplot as plt

class Fuselage():
    def __init__(self):
        self.xprofile = list()
        self.yprofile = list()
        self.length = 0.0
        self.diameter = 0.0
    
    def _naca_profile(self,length,diameter,trailingEdgeGap=0.0,npts=30):
        tc = diameter / 2./length
        x = flipud(self._cos_distribution(linspace(0,1,npts)))
        y = tc/0.2*(0.2969*(x)**0.5-0.1281*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)
        self.xprofile = x*length
        self.yprofile = y*length + x*trailingEdgeGap
        xcrd = hstack([flipud(self.xprofile), self.xprofile[1:]])
        ycrd = hstack([flipud(self.yprofile),-self.yprofile[1:]])
        self.coord = transpose(vstack([xcrd,ycrd]))
    def _cos_distribution(self,x):
        def cos_curve(x):
            return (cos(x*pi)+1.0)/2.0
        return array([cos_curve(xx) for xx in x])

def test():
    fus = Fuselage()
    fus._naca_profile(2.,0.5,0.1,20)
    plt.figure(1)
    plt.plot(fus.coord[:,0],fus.coord[:,1],'bo-')
    plt.axis('equal')
    plt.grid(True)
    plt.show()

if __name__=="__main__":
    test()