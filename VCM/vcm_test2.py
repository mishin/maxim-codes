# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 16:46:59 2013

@author: Maxim
"""
from vcm import *
from numpy import linalg

def _fhigh(x):
    return 4.0*x[0]**2 + x[1]**3 + x[0]*x[1]
def _ghigh(x):
    return 1.0/x[0] + 1.0/x[1] - 2.0
def _flow(x):
    return 4.0*(x[0]+0.1)**2 + (x[1]-0.1)**3 + x[0]*x[1]+0.1
def _glow(x):
    return 1.0/x[0] + 1.0/(x[1]+0.1) - 2.0 - 0.001
def _phigh(x):
    return _fhigh(x) + 100000.0*(max([0,_ghigh(x)]))**2.0
def _plow(x):
    return _flow(x) + 100000.0*(max([0,_glow(x)]))**2.0
def _fhigh1(x):
    return x*x - 4.0*x + 2.0
def _flow1(x):
    return x+1.0

def vcm_test_2d():
    x0 = array([1.5,1.5])
    delta = 1.0
    eta1 = 0.25
    eta2 = 0.75
    eta3 = 1.25
    c1 = 0.25
    c2 = 2.0
    tol = 1e-3
    err = tol+1
    fscaled = ScaledFunction(_plow,_phigh)
    plt.figure(1)
    plt.grid(True)
    plt.hold(True)
    while err>tol:
        fscaled.construct_scaling_model(x0)
        bnds = [(x0[0]-delta, x0[0]+delta),(x0[1]-delta, x0[1]+delta)]
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds)
        xnew = rslt.x
        fnew = rslt.fun
        rho = fscaled.get_thrust_region_ratio(xnew)
        err = linalg.norm(xnew-x0)
        if rho<eta1 or rho>=eta3:
            delta *= c1
        elif eta2<rho<eta3:
            if err==delta:
                gamma = c2
            else:
                gamma = 1.0
            delta *= gamma
        x0 = xnew
        print '%.4f\t'%rho, xnew, '%.4f\t'%fnew, '%.4f\t'%delta, '%.4f\t'%err
        plt.plot(xnew[0],xnew[1],'ro')
    fscaled.funcHi.display()
    fscaled.funcLo.display()
    plt.show()

def vcm_test_2dwire():
    x0 = array([3.0,3.0])
    func = TestFunction(_fhigh,'test function')
    ts = func.get_taylor(x0)
    
    x = linspace(0.5,10,50)
    y = linspace(0.5,10,50)
    X,Y = meshgrid(x,y)
    zs1 = array([func([xx,yy]) for xx,yy in zip(ravel(X), ravel(Y))])
    zs2 = array([ts([xx,yy]) for xx,yy in zip(ravel(X), ravel(Y))])
    Z1 = zs1.reshape(X.shape)
    Z2 = zs2.reshape(X.shape)
    
    xrbf = array([0.5,0.5,10,10,5])
    yrbf = array([0.5,10,0.5,10,5])
    zrbf = func([xrbf,yrbf])
    
    args = tuple()
    args = args + (xrbf,)
    args = args + (yrbf,)
    args = args + (zrbf,)
    rbf = Rbf(*args)
    
    zs3 = array([rbf(xx,yy) for xx,yy in zip(ravel(X), ravel(Y))])
    Z3 = zs3.reshape(X.shape)
    
    fig = plt.figure(1)
    ax = Axes3D(fig)
    ax.hold(True)
    ax.plot_wireframe(X,Y,Z1,color='y')
    ax.plot_wireframe(X,Y,Z2,color='b')
    ax.plot_wireframe(X,Y,Z3,color='r')
    ax.plot([x0[0]],[x0[1]],[func(x0)],'ro')
    plt.show()

def vcm_test_2d_2():
    vcmSolver = VCMoptimization()
    vcmSolver.set_objective_vcm(_flow1,_fhigh1,3.0)
    #vcmSolver.add_constraint_vcm(_glow,_ghigh)
    vcmSolver.solve()

def exact_solution():
    x0 = array([1.5,1.5])
    bnds = [(0.1,10),(0.1,10)]
    g = lambda x: -_ghigh(x)
    cnstr = ({'type':'ineq','fun':g},)
    rslt = minimize(_fhigh,x0,bounds=bnds, constraints=cnstr, method='SLSQP')
    print rslt

if __name__=="__main__":
    vcm_test_2d_2()