# -*- coding: utf-8 -*-
"""
Created on Sat Dec 07 14:33:53 2013

@author: Maxim
"""

from ScaledFunctions import *
from scipy.optimize import minimize, rosen
from write_data import FileOutput
from mpl_toolkits.mplot3d import Axes3D

fhigh1d = lambda x: (6.0*x-2.0)**2.0*np.sin(12.*x-4.)

def fhigh(x):
    A = 0.07
    B = 0.9
    C = 0.2
    return fhigh1d(np.linalg.norm(A*x)) + B*np.sum(x) + C

def flow(x):
    A = 0.035
    B = 0.5
    C = 0.1
    return fhigh1d(np.linalg.norm(A*x)) + B*np.sum(x) + C

fhigh = lambda x: rosen(x)
flow = lambda x: 0.5*rosen(x) + 1.0*x[0] + .1*x[1] -1.0
g1high   = lambda x: x[0]**2*x[1]/20. -.1
g1low  = lambda x: x[0]**2*x[1]/20. + (x[0]+x[1])/5. -2
g2     = lambda x: (x[0]+x[1]-5)**2/30 + (x[0]-x[1]-12)**2/120 - 1


#fhigh = lambda x: 4.*x[0]**2. + x[1]**3. + x[0]*x[1]
#flow = lambda x: 4.*(x[0]+0.1)**2. + (x[1]-0.1)**3. + x[0]*x[1]+0.1
#
#g1high = lambda x: -(1./x[0]+ 1/x[1] -2.)
#g1low = lambda x: -(1./x[0] + 1./(x[1]+0.1)-2.001)


def problem_2d():
    cnstr = ({'type':'ineq','fun':g1high},{'type':'ineq','fun':g2})
    #cnstr = ({'type':'ineq','fun':g1high},)
    print minimize(fhigh,[3.,3.],method='SLSQP',bounds=((0,10),(0,10),),constraints=cnstr)
    xL = np.array([-2,-2])
    xU = np.array([2,2])
    x0 = np.array([1,1])
    delta = 1.0
    nw = 3
    tol = 1e-4
    err = tol+1
    iterMax = 50
    nIter = 0

    prefix = 'AVCM'
    fscaled = AdditiveScaling(fhigh,flow,nw)
    gscaled = AdditiveScaling(g1high,g1low,nw)
    #xDOE = np.array([[10.0,3.3333],[0.0,6.6667],[3.3333,0.0],[6.6667,10.0]])
    #fscaled.initialize_by_points(xDOE)
    #gscaled.initialize_by_points(xDOE)

    trReg = TrustRegionManagement(delta)

    report = FileOutput('twodim_results_20140717.py')
    _x10 = list()
    _x20 = list()
    _x1new = list()
    _x2new = list()
    _fh = list()
    _fl = list()
    _fs = list()
    _gh = list()
    _gl = list()
    _gs = list()
    _delta = list()
    _rho = list()
    _err = list()

    while err>tol and nIter<iterMax:
        fscaled.construct_scaling_model(x0)
        gscaled.construct_scaling_model(x0)
        bnds = get_bounds(x0,delta,xL,xU)
        cnstr = ({'type':'ineq','fun':gscaled},{'type':'ineq','fun':g2})
        #cnstr = ({'type':'ineq','fun':gscaled},)
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,constraints=cnstr)
        xnew = rslt.x
        fnew = rslt.fun
        rho1 = fscaled.get_trust_region_ratio(xnew)
        rho2 = gscaled.get_trust_region_ratio(xnew)
        err = np.linalg.norm([x0-xnew])
        delta1,x0 = trReg.adjust(rho1,x0,xnew)
        delta2,tmp = trReg.adjust(rho2,x0,xnew)
        delta = min(delta1,delta2)
        nIter += 1
        print '%.4f\t%.4f\t%.4f\t%.4f\t%.2e\t%.4f\t%.4f'%(xnew[0],xnew[1],x0[0],x0[1], err, fnew, rho1)

        _x10.append(x0[0])
        _x20.append(x0[1])
        _x1new.append(xnew[0])
        _x2new.append(xnew[1])
        _fh.append(fhigh(xnew))
        _fl.append(flow(xnew))
        _fs.append(fnew)
        _gh.append(g1high(xnew))
        _gl.append(g1low(xnew))
        _gs.append(gscaled(xnew))
        _delta.append(delta)
        _rho.append(rho1)
        _err.append(err)
    
    report.write_string('#===== %s ====='%prefix)
    report.write_string(prefix+'fhEval = %d'%fscaled.fHigh._nEval)
    report.write_string(prefix+'flEval = %d'%fscaled.fLow._nEval)
    report.write_string(prefix+'ghEval = %d'%gscaled.fHigh._nEval)
    report.write_string(prefix+'glEval = %d'%gscaled.fLow._nEval)
    report.write_string(prefix+'nIter = %d'%nIter)
    report.write_array(_x10,prefix+'x10')
    report.write_array(_x20,prefix+'x20')
    report.write_array(_x1new,prefix+'x1new')
    report.write_array(_x2new,prefix+'x2new')
    report.write_array(_fh,prefix+'fh')
    report.write_array(_fl,prefix+'fl')
    report.write_array(_fs,prefix+'fs')
    report.write_array(_gh,prefix+'gh')
    report.write_array(_gl,prefix+'gl')
    report.write_array(_gs,prefix+'gs')
    report.write_array(_delta,prefix+'delta')
    report.write_array(_rho,prefix+'rho')
    report.write_array(_err,prefix+'err')
    


def plot_functions():
    x = y = np.arange(0,10,0.25)
    X, Y = np.meshgrid(x,y)
    Z1 = np.array([fhigh(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z2 = np.array([flow(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z1 = np.reshape(Z1,X.shape)
    Z2 = np.reshape(Z2,X.shape)
    Z3 = np.array([g1high(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z4 = np.array([g2(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z5 = np.array([g1low(np.array([_x,_y])) for _x,_y in zip(np.ravel(X),np.ravel(Y))])
    Z3 = np.reshape(Z3,X.shape)
    Z4 = np.reshape(Z4,X.shape)
    Z5 = np.reshape(Z5,X.shape)
    
    fig = plt.figure(1)
    ax = Axes3D(fig)
    #ax = fig.add_subplot(111)
    ax.hold(True)
    ax.plot_wireframe(X,Y,Z1,color='r',label='High fidelity function')
    ax.plot_wireframe(X,Y,Z2,color='b',label='Low fidelity function')
#    ax.plot_wireframe(X,Y,Z3,color='y')
#    ax.plot_wireframe(X,Y,Z4,color='g')
#    ax.plot([r.x[0]],[r.x[1]],[r.fun],'ro')
    ax.legend()
    fig2 = plt.figure(2)
    ax = fig2.add_subplot(111)
    ax.hold(True)
    ax.contour(X,Y,Z1,levels=[7.5])
    #ax.contour(X,Y,Z2,levels=[0],color='r')
    ax.contour(X,Y,Z5,levels=[0],colors='r')
    ax.contour(X,Y,Z3,levels=[0],colors='k')
    ax.contour(X,Y,Z4,levels=[0],colors='b')
    plt.show()

def solve_exact():
    x0 = np.array([3.,3.])
    bnds = ((0.,10),(0,10))
    cnstr = ({'type':'ineq','fun':g1high},{'type':'ineq','fun':g2})
    rslt = minimize(fhigh,x0,method='SLSQP',bounds=bnds,constraints=cnstr)
    xVCM = np.array([6.550791,6.550791])
    xAVCM = np.array([6.550898,6.550898])
    xGVFM = np.array([3.113886,2.062646])
    ref = 200.**0.5
    print np.linalg.norm(xVCM-rslt.x)/ref*100.
    print np.linalg.norm(xAVCM-rslt.x)/ref*100.
    print np.linalg.norm(xGVFM-rslt.x)/ref*100.
if __name__=="__main__":
    problem_2d()
    #plot_functions()
