# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 19:48:36 2013

@author: Maxim
"""
import airfoil
from numpy import array, zeros, ones
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def cst_fit():

    
    af = airfoil.Airfoil()
    #af.read_txt(r'C:\Users\Maxim\Dropbox\2. projects\VCM\transonic airfoil\rae2822.txt',afType=2)
    af.read_txt(r'D:\lab_sync\1. Projects\RENN\airfoil design\AG24new.txt')
    #af.naca4()
    print af.coord
    #af.naca4(nPts=30)
    #af.naca4()
    up, lo = af.upPts, af.loPts
    args1 = up,lo

    def _obj(x,*args):
        up,lo = args
        Au = array([x[0],x[1],x[2],x[3]])
        #Al = array([-x[0],x[4],x[5],x[6]])
        Al = array([x[4],x[5],x[6],x[7]])
        afcst = airfoil.cst(Au,Al)
        sqErr = 0.0
        for i,pt in enumerate(up[:,0]):
            sqErr += (afcst.upCurve(pt) - up[i,1])**2.0
        for i,pt in enumerate(lo[:,0]):
            sqErr += (afcst.loCurve(pt) - lo[i,1])**2.0
        return sqErr
    
    def _err(x,*args):
        up,lo = args
        Au = array([x[0],x[1],x[2],x[3]])
        #Al = array([-x[0],x[4],x[5],x[6]])
        Al = array([x[4],x[5],x[6],x[7]])
        afcst = airfoil.cst(Au,Al)
        upErr = zeros(len(up))
        loErr = zeros(len(lo))
        for i,pt in enumerate(up[:,0]):
            upErr[i] = afcst.upCurve(pt) - up[i,1]
        for i,pt in enumerate(lo[:,0]):
            loErr[i] = (afcst.loCurve(pt) - lo[i,1])
        return upErr, loErr
    
    bnds = ((0.0,1.0),(-0.3,1.0),(-0.3,1.0),(-0.3,1.0),(-1.0,0.0),(-1.0,0.3),(-1.0,0.3),(-1.0,0.3))
    rslt = minimize(_obj,x0=0.1*ones(8),bounds=bnds,method='SLSQP',args=args1)
    xnew = rslt.x
    upErr, lowErr = _err(xnew,up,lo)
    
    #xnew[:4] += +0.075
    #xnew[4:] += -0.075
    print xnew
    #xnew += 0.05
    afnew = airfoil.cst(xnew[0:4],xnew[4:8])
    print afnew.thickness
    fig = plt.figure(1)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.hold(True)
    ax1.set_title('AG24 a')
    ax2.set_title('Fit error')
    ax1.plot(af.coord[:,0],af.coord[:,1],'k-')
    ax1.plot(afnew.coord[:,0],afnew.coord[:,1],'r.-')
    ax1.axis('equal')
    ax1.legend(['Original Airfoil','CST fit'])
    ax2.hold(True)
    ax2.plot(up[:,0],upErr,'ks-',markevery=3)
    ax2.plot(lo[1:,0],lowErr[1:],'k^-',markevery=3)
    ax2.set_ylabel('error')
    ax2.legend(['Upper curve','Lower curve'],'upper right')
    plt.show()

def tmp_cst_plot():
    	
    Au = [0.120718,0.135291,0.174915,0.142468]
    Al = [-0.120718,-0.150704,-0.176010,-0.040013]
    af = airfoil.cst(Au,Al)
    print af.thickness
    af.plot()

if __name__=="__main__":
    cst_fit()