# -*- coding: utf-8 -*-
"""
Created on Sat Jun 28 18:43:12 2014

@author: Maxim
"""

from geometry import BezierCurve, get_sine_distribution, PtsArray
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize



def get_bezier_airfoil(x=None,nPts=30,show=False):
    if x==None:
        x = np.array([0.03,0.04,0.06,0.04,  0.01,0.01,-0.03])
    nodesP = np.zeros([6,2])
    nodesC = np.zeros([5,2])
    
    nodesP[1:,0] = np.array([0.0,0.1,0.4,0.75,1.0])
    nodesP[1:-1,1] = x[:4]
    
    nodesC[1:,0] = np.array([0.3,0.6,0.9,1.0])
    nodesC[1:-1,1] = x[4:]

    curveUp = BezierCurve(nodesP)
    curveCamber = BezierCurve(nodesC)
    t = get_sine_distribution(nPts)
    
    
    ptsCamber = curveCamber(np.linspace(0,1,nPts))
    curveCamber2 = interp1d(ptsCamber[:,0],ptsCamber[:,1],'cubic')
    
    _ptsUp = curveUp(t)
    ptsUpX = _ptsUp[:,0]
    ptsUpY = _ptsUp[:,1] + curveCamber2(_ptsUp[:,0])
    ptsLoY = -_ptsUp[:,1] + curveCamber2(_ptsUp[:,0])
    ptsUp = np.transpose(np.vstack([ptsUpX,ptsUpY]))
    ptsLo = np.transpose(np.vstack([ptsUpX,ptsLoY]))

    if show:
        plt.figure()
        plt.hold(True)
        plt.plot(nodesP[:,0],nodesP[:,1],'rs--')
        plt.plot(nodesC[:,0],nodesC[:,1],'bs--')
        plt.plot(ptsUpX,ptsUpY,'k-')
        plt.plot(ptsUpX,ptsLoY,'k-')
        plt.plot(ptsCamber[:,0],ptsCamber[:,1],'g-')
        plt.axis('equal')
        plt.show()
    
    return ptsUp, ptsLo


def fit_airfoil():
    import airfoil
    af = airfoil.load('GA37A315mod')
    curveUp = interp1d(af.ptsUp[:,0],af.ptsUp[:,1],'cubic')
    curveLo = interp1d(af.ptsLo[:,0],af.ptsLo[:,1],'cubic')
    args = (curveUp,curveLo,)
    def func(x,curveUp,curveLo):
        ptsUp,ptsLo = get_bezier_airfoil(x)
        err = 0.0
        for pt in ptsUp:
            err += (curveUp(pt[0])-pt[1])**2.0
        for pt in ptsLo:
            err += (curveLo(pt[0])-pt[1])**2.0
        return err
    
    x0 = np.array([0.03,0.04,0.06,0.04,  0.01,0.01,-0.03])
    rslt = minimize(func,x0,method='SLSQP',args=args)
    ptsUp,ptsLo = get_bezier_airfoil(rslt.x,30)
    
    plt.figure()
    plt.hold(True)
    plt.plot(ptsUp[:,0],ptsUp[:,1],'r-')
    plt.plot(ptsLo[:,0],ptsLo[:,1],'r-')
    plt.plot(af.pts[:,0],af.pts[:,1],'k--')
    plt.axis('equal')
    plt.show()
        

if __name__=="__main__":
    fit_airfoil()