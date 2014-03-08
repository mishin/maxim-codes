# -*- coding: utf-8 -*-
"""
Created on Sat Mar 08 13:58:36 2014

@author: Maxim
"""

from airfoil import *

def xfoil_convergence():
    numberOfPoints = np.linspace(40,200,10)
    #af = read_txt('GA37A315mod_reworked.dat')
    af = read_txt('pb092_closed.txt')
    
    n = len(numberOfPoints)
    clmax = np.zeros(n)
    cd0 = np.zeros(n)
    alphaClmax = np.zeros(n)
    alphaCdmin = np.zeros(n)
    stall = 'calcfoil'
    trans = 'drelaAfter1991'
    surface = 'smooth'
    distribution = 'cos'
    M = 0.16
    Re = 3.3e6
    
    
    for i,pts in enumerate(numberOfPoints):
        af.redim(pts,True,distribution)
        af.display()
        try:
            #pol = af.get_jfoil_polar(M,Re,[-20,20,1.0],stall,trans,surface)
            pol = af.get_xfoil_polar(M,Re,[-20,20,1.0])
            clmax[i] = pol.clmax
            cd0[i] = pol.cdmin
            alphaClmax[i] = pol.alphaClmax
            alphaCdmin[i] = pol.alphaCdmin
        except:
            clmax[i] = None
            cd0[i] = None
            alphaClmax[i] = None
            alphaCdmin[i] = None
    
    xmin = numberOfPoints[0]-10
    xmax = numberOfPoints[-1]
    fig = plt.figure(1)
    fig.suptitle('%s distribution'%distribution)
    ax1 = fig.add_subplot(411)
    ax1.grid(True)
    ax1.plot(numberOfPoints,clmax,'ro-')
    ax1.axis([xmin, xmax, 1.2,1.8])
    ax1.set_ylabel('clmax')
    
    ax2 = fig.add_subplot(412)
    ax2.grid(True)
    ax2.plot(numberOfPoints,alphaClmax,'ro-')
    ax2.axis([xmin,xmax,12,22])
    ax2.set_ylabel('alpha clmax')
    
    ax3 = fig.add_subplot(413)
    ax3.grid(True)
    ax3.plot(numberOfPoints,cd0*1e3,'ko-')
    ax3.axis([xmin,xmax,3,10])
    ax3.set_ylabel('cdmin')
    
    ax4 = fig.add_subplot(414)
    ax4.grid(True)
    ax4.plot(numberOfPoints,alphaCdmin,'ko-')
    ax4.axis([xmin,xmax,0,2])
    ax4.set_ylabel('alpha cdmin')
    ax4.set_xlabel('number of points')
    plt.show()


if __name__=="__main__":
    xfoil_convergence()