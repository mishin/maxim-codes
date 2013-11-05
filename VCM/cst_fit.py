# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 19:48:36 2013

@author: Maxim
"""

def cst_fit():
    from numpy import array, zeros, ones
    from scipy.optimize import minimize
    import airfoil
    import matplotlib.pyplot as plt
    
    af = airfoil.Airfoil()
    af.read_txt(r'C:\Users\Maxim\Dropbox\2. projects\VCM\transonic airfoil\rae2822.txt',afType=2)
    print af.coord
    #af.naca4(nPts=30)
    #af.naca4()
    up, lo = af.upPts, af.loPts
    args1 = up,lo

    def _obj(x,*args):
        up,lo = args
        Au = array([x[0],x[1],x[2],x[3]])
        Al = array([-x[0],x[4],x[5],x[6]])
        afcst = airfoil.cst(Au,Al)
        sqErr = 0.0
        for i,pt in enumerate(up[:,0]):
            sqErr += (afcst.upCurve(pt) - up[i,1])**2.0
        for i,pt in enumerate(lo[:,0]):
            sqErr += (afcst.loCurve(pt) - lo[i,1])**2.0
        return sqErr
    
    bnds = ((0.0,1.0),(0.0,1.0),(0.0,1.0),(0.0,1.0),(-1.0,0.0),(-1.0,0.0),(-1.0,0.0))
    rslt = minimize(_obj,x0=0.1*ones(7),bounds=bnds,method='SLSQP',args=args1)
    xnew = rslt.x
    print rslt
    
    #xnew[:4] += +0.075
    #xnew[4:] += -0.075
    print xnew
    afnew = airfoil.cst(xnew[0:4],[-xnew[0],xnew[4],xnew[5],xnew[6]])
    
    plt.figure(1)
    plt.plot(af.coord[:,0],af.coord[:,1],'bo-')
    plt.grid(True)
    plt.hold(True)
    plt.axis('equal')
    plt.plot(afnew.coord[:,0],afnew.coord[:,1],'ro-')
    plt.show()


if __name__=="__main__":
    cst_fit()