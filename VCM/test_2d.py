# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 18:03:14 2013

@author: Maxim
"""
import numpy as np
from vcm import *
from pyDOE import lhs, ff2n
from tmp_doe import read_samples
import matplotlib.pyplot as plt

"""
exact solution: x1=0.8842, x2=1.1507, f=5.6683
"""

def flow(x):
    return 4.*(x[0]+0.1)**2.+(x[1]-0.1)**3.+x[0]*x[1]+0.1
def fhigh(x):
    return 4.*x[0]**2.+x[1]**3.+x[0]*x[1]
def glow(x):
    return -(1./x[0]+1./(x[1]+0.1)-2.001)
def ghigh(x):
    return -(1./x[0]+1./x[1]-2.)


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
        
def run_2d_example():
    tol  = 1.0e-4
    lb   = np.array([1.0,1.0])
    ub   = np.array([10.,10.])
    x0   = np.array([5.0,5.0])
    iterMax = 50
    
    xDoe = list()
    xDoe.append([lb[0],lb[1]])
    xDoe.append([ub[0],ub[1]])
    xDoe.append([lb[0],ub[1]])
    xDoe.append([ub[0],lb[1]])
    xDoe.append([(lb[0]+ub[0])/2.,(lb[1]+ub[1])/2.])

    delta = min([min(xu-x,x-xl) for x,xu,xl in zip(x0,ub,lb)])
    trustRegion = TrustRegionManagement(delta, 0.25, 0.75, 1.25, 0.3, 2.0)
    err = tol + 1
    niter = 0
    gNotConverged=True
    
    fscaled = ScaledFunction(flow,fhigh,1,'mult')
    gscaled = ScaledFunction(glow,ghigh,1,'mult')
    
    fscaled._initialize_by_doe_points(xDoe)
    gscaled._initialize_by_doe_points(xDoe)
    
    while err>tol and iterMax>niter and gNotConverged:
        fscaled.construct_scaling_model(x0)
        gscaled.construct_scaling_model(x0)
        bnds = [(x0[0]-delta,x0[0]+delta),(x0[1]-delta,x0[1]+delta)]
        cnstr = {'type':'ineq','fun':gscaled}
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,tol=1e-6,
                        constraints=cnstr)
        xnew = rslt.x
        fnew = rslt.fun
        rho1,fHiNew = fscaled.get_trust_region_ratio(xnew)
        rho2,gHiNew = gscaled.get_trust_region_ratio(xnew)
        rho = min([rho1,rho2])
        err = np.linalg.norm([x0-xnew])
        delta = trustRegion.adjust(rho,err)
        x0 = xnew
        niter +=1
        print '%.4f\t%.4f\t'%(x0[0],x0[1]),'%.2f\t%.2e\t'%(rho,delta),'%.2e'%err,'%.4f\t%.4f\t'%(glow(xnew),ghigh(xnew))
    gscaled.funcHi.display()
    fscaled.funcHi.display()
    gscaled.funcLo.display()
    fscaled.funcLo.display()
    print niter

def currin(x):
    val1 = 1.0-np.exp(-0.5/x[1])
    val2 = 2300.*x[0]**3 + 1900*x[0]**2 + 2092*x[0]+60
    val3 = 100.*x[0]**3+500*x[0]**2+4*x[0]+20.
    return -val1*val2/val3

def currin_low(x):
    x1 = x[0]
    x2 = x[1]
    maxarg = max([0, x2+0.05])
    val1 = currin([x1+1/20, x2+1/20])
    val2 = currin([x1+1/20, maxarg])
    val3 = currin([x1-1/20, x2+1/20])
    val4 = currin([x1-1/20, maxarg])
    return -0.25*(val1+val2+val3+val4)
    
def run_2d_unc():
    tol  = 1.0e-4
    lb   = np.array([0.1,0.1])
    ub   = np.array([1.,1.])
    x0   = np.array([.5,.5])
    iterMax = 50
    
    xDoe = list()
    xDoe.append([lb[0],lb[1]])
    xDoe.append([ub[0],ub[1]])
    xDoe.append([lb[0],ub[1]])
    xDoe.append([ub[0],lb[1]])
    xDoe.append([(lb[0]+ub[0])/2.,(lb[1]+ub[1])/2.])

    delta = min([min(xu-x,x-xl) for x,xu,xl in zip(x0,ub,lb)])
    trustRegion = TrustRegionManagement(delta, 0.25, 0.75, 1.25, 0.3, 2.0)
    err = tol + 1
    niter = 0
    
    fscaled = ScaledFunction(currin_low,currin,1,'mult')
    fscaled._initialize_by_doe_points(xDoe)
    
    while err>tol and iterMax>niter:
        fscaled.construct_scaling_model(x0)
        xlb = [max(lb[i],x0[i]-delta) for i in range(len(x0))]
        xub = [min(ub[i],x0[i]+delta) for i in range(len(x0))]
        bnds = [(xlb[0],xub[0]),(xlb[1],xub[1])]
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,tol=1e-6)
        xnew = rslt.x
        fnew = rslt.fun
        rho,fhiNew = fscaled.get_trust_region_ratio(xnew)
        err = np.linalg.norm([x0-xnew])
        delta = trustRegion.adjust(rho,err)
        x0 = xnew
        niter +=1
        if niter%5: print_header()
        print '%.4f\t%.4f\t'%(x0[0],x0[1]),'%.2f\t%.2e\t'%(rho,delta),'%.2e'%err
    fscaled.funcHi.display()
    fscaled.funcLo.display()
    print niter

def exact_solution():
    rslt = minimize(currin,[.5,.5],method='SLSQP',bounds=((0.0,1.),(.0,1.)))
    print rslt

def get_bounds(x0,delta,lb,ub):
    bnds = list()
    for i,xx in enumerate(x0):
        bnd = (max([lb[i],xx-delta]),min([ub[i],xx+delta]))
        bnds.append(bnd)
    return bnds

def ForresterFunc(x):
    return (5.0*x-2.)**2*sin(12.*x-4)

def numerical_example_2():
#    fHigh  = lambda x: np.exp(x[0]/3) + np.exp(x[1]/5) - x[0]
#    fLow   = lambda x: x[0] + x[1]

    
    def fHigh(x):
        x = 0.1*x
        return ForresterFunc(np.linalg.norm(x)/1.2)+5.*(x[0]+x[1])
    fLow  = lambda x: np.exp(x[0]/3) + np.exp(x[1]/5) - x[0]
    g1High = lambda x: x[0]**2*x[1]/20. - 1
    g1Low  = lambda x: x[0]**2*x[1]/20. + (x[0]+x[1])/5. - 2
    g2     = lambda x: (x[0]+x[1]-5)**2/30 + (x[0]-x[1]-12)**2/120 - 1
    
    tol  = 1.0e-4
    gtol = 1.0e-4
    maxIter = 50
    lb   = np.array([0.0,0.0])
    ub   = np.array([10.,10.])
    x0   = np.array([5.0,5.0])
    
    dx = 0.1
    xmsh = ymsh = np.arange(-0.5,10.5,dx)
    Xmsh, Ymsh = np.meshgrid(xmsh,ymsh)
    zs1 = np.array([fHigh(np.array([x,y])) for x,y in zip(np.ravel(Xmsh),np.ravel(Ymsh))])
    zs2 = np.array([g1High(np.array([x,y])) for x,y in zip(np.ravel(Xmsh),np.ravel(Ymsh))])
    zs3 = np.array([g2(np.array([x,y])) for x,y in zip(np.ravel(Xmsh),np.ravel(Ymsh))])
    Z1 = zs1.reshape(Xmsh.shape)
    Z2 = zs2.reshape(Xmsh.shape)
    Z3 = zs3.reshape(Xmsh.shape)
    

    #plt.show()
    nDoe = 4 # nDoe>20 creates numerical noise and fails to converge
    
    err = tol+1.
    niter = 0
    gConverged = False
    xConverged = False
    delta = min([min(xu-x,x-xl) for x,xu,xl in zip(x0,ub,lb)])
    trustRegion = TrustRegionManagement(delta, 0.25, 0.75, 1.25, 0.3, 2.0)

    #xDoe = lhs(len(x0),nDoe)
    #xDoe = ff2n(2)
    xDoe = read_samples()
    xDoe = denormalize(xDoe[:nDoe+1],lb,ub,1)
    print xDoe
    fscaled = ScaledFunction(fLow, fHigh,0,'add')
    gscaled = ScaledFunction(g1Low, g1High,0,'add')

    fscaled._initialize_by_doe_points(xDoe)
    gscaled._initialize_by_doe_points(xDoe)
    def print_header():
        print 'x1\tx2\tf\trho\tdelta\terr\tgScaled\tgHigh'
    while xConverged==False or gConverged==False:
        fscaled.construct_scaling_model(x0)
        gscaled.construct_scaling_model(x0)
        bnds = get_bounds(x0,delta,lb,ub)
        cnstr = ({'type':'ineq','fun':gscaled},{'type':'ineq','fun':g2})
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,constraints=cnstr,tol=1e-10)
        xnew = rslt.x
        fnew = rslt.fun
        gScaledNew = gscaled(xnew)
        rho1, fHighNew = fscaled.get_trust_region_ratio(xnew)
        rho2, gHighNew = gscaled.get_trust_region_ratio(xnew)
        rho = min([rho1,rho2])
        err = np.linalg.norm([x0-xnew])

        delta = trustRegion.adjust(rho,err)
        x0 = xnew
        if niter%5==0:
            print_header()
        print '%.4f\t%.4f\t%.4f\t'%(x0[0],x0[1],fnew),'%.2f\t%.2e\t'%(rho,delta),'%.2e'%err,'%.4f\t%.4f\t'%(gScaledNew,gHighNew)
        niter += 1
        
        if err<=tol or niter>maxIter:
            xConverged=True
        else:
            xConverged=False

        if (-gtol<=gHighNew<=gtol and gScaledNew<=gtol) or (gHighNew>=0.0 and gScaledNew>gtol):
            gConverged = True
        else:
            gConverged = False
        xbounds = np.array([bnds[0][0],bnds[0][0],bnds[0][1],bnds[0][1],bnds[0][0]])
        ybounds = np.array([bnds[1][0],bnds[1][1],bnds[1][1],bnds[1][0],bnds[1][0]])
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.hold(True)
        zs4 = np.array([fscaled(np.array([x,y])) for x,y in zip(np.ravel(Xmsh),np.ravel(Ymsh))])
        zs5 = np.array([gscaled(np.array([x,y])) for x,y in zip(np.ravel(Xmsh),np.ravel(Ymsh))])
        Z4 = zs4.reshape(Xmsh.shape)
        Z5 = zs5.reshape(Xmsh.shape)
        cs1 = ax.contour(Xmsh,Ymsh,Z1,colors='m',levels=[0,2,4,6,9,12,15])
        ax.contour(Xmsh,Ymsh,Z2,colors='r',levels=[0])
        ax.contour(Xmsh,Ymsh,Z3,colors='k',levels=[0])
        cs2 = ax.contour(Xmsh,Ymsh,Z4,colors='g',levels=[0,2,4,6,9,12,15])
        ax.contour(Xmsh,Ymsh,Z5,colors='b',levels=[0])
        tmpx, tmpy = [-5,-5],[1,2]
        ax.plot(tmpx,tmpy,'m')
        ax.plot(tmpx,tmpy,'g')
        ax.plot(tmpx,tmpy,'r')
        ax.plot(tmpx,tmpy,'b')
        ax.plot(tmpx,tmpy,'k')
        ax.plot(np.array(fscaled.xPrev)[:,0],np.array(fscaled.xPrev)[:,1],'ko')
        ax.plot(xbounds,ybounds,'k--')
        ax.axis([-0.5,10.5,-0.5,10.5])
        plt.clabel(cs1, fontsize=9, inline=1)
        plt.clabel(cs2, fontsize=9, inline=1)
        
        plt.legend(['Hifi objective','Lowfi objective','Hifi constraint 1','Lofi constraint 1','constraint 2','Hifi sample points','Bounds'])
        plt.show()
        plt.cla()
        #print xConverged, gConverged, xConverged and gConverged
    fscaled.funcHi.display()
    fscaled.funcLo.display()
    gscaled.funcHi.display()
    gscaled.funcLo.display()
    print niter
    
    xexact = np.array([ 3.11388573,  2.06264578])
    print np.linalg.norm(xexact-xnew)
    print fnew-2.5367139107528409
if __name__=="__main__":
    #run_2d_example()
    numerical_example_2()