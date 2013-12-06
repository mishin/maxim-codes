# -*- coding: utf-8 -*-
"""
Created on Thu Dec 05 20:49:28 2013

@author: Maxim
"""

from ScaledFunctions import *
from scipy.optimize import minimize
from write_data import FileOutput

fhigh = lambda x: (6.0*x-2.0)**2.0*np.sin(12.*x-4.)
flow  = lambda x: 0.5*fhigh(x) + 8.*(x-.5)+5.

def run_test_additive():
    xL, x0, xU = 0.0, 0.5, 1.0
    delta = 0.5
    xDOE = np.array([0.2,0.55,0.8])
    prefix = 'GVFMhybr'
    fsc = HybridScaling(fhigh,flow,0,1e-6,0.)
    report = FileOutput('onedim_results_20131206_2.py')
    fsc.initialize_by_points(xDOE)
    
    trReg = TrustRegionManagement(delta)
    tol = 1e-6
    err = tol+1
    iterMax = 50
    nIter = 0
    x = np.linspace(0,1,100)
    yh = np.array([fhigh(_x) for _x in x])
    yl = np.array([flow(_x) for _x in x])
    
    _x0 = list()
    _fh = list()
    _fl = list()
    _fs = list()
    _rho = list()
    _err = list()
    _delta = list()
    while err>tol and nIter<iterMax:
        fsc.construct_scaling_model(x0)
        bnds = get_bounds(x0,delta,xL,xU)
        rslt = minimize(lambda _x: fsc(_x,False),x0,method='SLSQP',bounds=bnds)
        xnew = rslt.x[0]
        fnew = rslt.fun
        rho = fsc.get_trust_region_ratio(xnew)
        err = abs(x0-xnew)
        delta,x0 = trReg.adjust(rho,x0,xnew)
        nIter += 1
        print '%.4f\t%.4f\t%.2e\t%.4f\t%.4f'%(xnew,x0, err, fnew, rho)
        
        _x0.append(xnew)
        _fh.append(fhigh(xnew))
        _fl.append(flow(xnew))
        _fs.append(fnew)
        _rho.append(rho)
        _err.append(err)
        _delta.append(delta)
    
    print nIter
    print fsc._scalingAdd.fHigh._nEval
    print fsc._scalingAdd.fLow._nEval

    report.write_string('#===== %s ====='%prefix)
    report.write_string(prefix+'fHighEval = %d'%fsc._scalingAdd.fHigh._nEval)
    report.write_string(prefix+'fLowEval = %d'%fsc._scalingAdd.fLow._nEval)
    report.write_string(prefix+'nIter = %d'%nIter)
    report.write_array(_x0,prefix+'x0')
    report.write_array(_fh,prefix+'fh')
    report.write_array(_fl,prefix+'fl')
    report.write_array(_fs,prefix+'fs')
    report.write_array(_rho,prefix+'rho')
    report.write_array(_err,prefix+'err')
    report.write_array(_delta,prefix+'delta')
    report.write_array(yh,'fhigh')
    report.write_array(yl,'flow')
    fSc = np.array([fsc(_x) for _x in x])
    report.write_array(fSc,prefix+'fscaled')
    report.write_array(fsc._scalingAdd.fHigh._histXpart,prefix+'xHist')
    report.write_array(fsc._scalingAdd.fHigh._histFpart,prefix+'fHist')

    fsc._scalingAdd.fHigh.write_history('history.txt')
    plt.figure(1)
    plt.hold(True)
    plt.plot(x,yh,'b')
    plt.plot(x,yl,'r')
    plt.plot(x,fSc,'k')
    plt.plot(fsc._scalingAdd.fHigh._histXpart,fsc._scalingAdd.fHigh._histFpart,'go')
    plt.plot(x0,fsc.fHigh(x0,False),'ro')
    plt.axis([0,1,-10,20])
    plt.show()
    plt.cla()

if __name__=="__main__":
    run_test_additive()