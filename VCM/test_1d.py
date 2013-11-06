# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 16:51:31 2013

@author: Maxim
"""

from vcm import *

class OutputData():
    def __init__(self,path,mode='wt'):
        self.path = path
        self.fid = open(path,mode)
    def write_string(self,string):
        self.fid.write(string)
    def write_array(self,name,data):
        self.write_string('%s = array(['%name)
        for val in data[:-1]:
            self.write_string('%.4f, '%val)
        self.write_string('%.4f])\n\n'%data[-1])
    def close(self):
        self.fid.close()

def vcm_test_1d():
#    def _flow(x):
#        return x+1.0
#    def _fhi(x):
#        return x*x - 4.0*x + 2.0
#    def _flow(x):
#        return (x-2.15)**2.0/2+1.25
#    def _fhi(x):
#        return sin(x*pi/2.0)+(x-2.0)**2.0/2
#    def _flow(x):
#        return (x-2.2)**2.0/2+2.5
#    def _fhi(x):
#        return sin(x*pi)+(x-2.0)**2.0/2 +2.0
    def forrester(x):
        return (5.*x-2.)**2.*sin(12.*x-4.)
    def forrester_low(x):
        A, B, C = 0.5, 10., -5.
        return A*forrester(x)+B*(x-0.5)-C
    output = OutputData('1D_example_results.py','a')
    #output.write_string('from numpy import array')
    _fhi = forrester
    _flow = forrester_low
    eta1 = 0.25
    eta2 = 0.75
    eta3 = 1.25
    c1 = 0.3
    c2 = 2.0
    tol = 1.0e-3
    err = tol + 1
    niter = 0
    fscaled = HybridScaledFunction(_flow,_fhi,0,weight=0.5)
    #fscaled = ScaledFunction(_flow,_fhi,3,'mult')
    lb = 0.0
    ub = 1.0
    x0 = 0.5
    delta = min([abs(x0-lb),abs(x0-ub)])*0.5
    dSpace = ub-lb
    x = linspace(lb-0.1*dSpace,ub+0.1*dSpace,75)
    #output.write_array('X',x)
    #output.write_array('Yhigh',forrester(x))
    #output.write_array('Ylow',forrester_low(x))
    xdoe = [x0-delta,x0+delta]
    xdoe = [lb,ub]
    fscaled._initialize_by_doe_points(xdoe)
    fhiNew = fscaled.funcHi(x0)
    while err>tol:
        fscaled.construct_scaling_model(x0,fhiNew)
        x0plt = x0
#        plt.figure(1)
#        plt.hold(True)
#        plt.plot(x,fscaled.funcLo(x),'r.-',linewidth=1.5)
#        plt.plot(x,fscaled.funcHi(x),'b-',linewidth=1.5)
#        fsc = [fscaled(xsc) for xsc in x]
#        plt.plot(x,fsc,'k-',linewidth=1.5)
#        plt.plot(fscaled.xPrev,fscaled.fPrev,'ro',ms=7,mew=0)
#        plt.plot(array([x0plt,x0plt])+delta,[-50,50],'k--')
#        plt.plot(array([x0plt,x0plt])-delta,[-50,50],'k--')
#        plt.legend(['Low-fidelity','High-fidelity','Scaled','Sample point','Trust region'],'upper left')
#        #plt.axis([lb-0.5,ub+0.5,-2,6])
#        plt.axis([lb-0.1*dSpace,ub+0.1*dSpace,-10,20])
#        plt.show()
#        plt.cla()
        
        bnds = [(max([lb,x0-delta]), min([x0+delta,ub]))]
        rslt = minimize(fscaled,x0,method='SLSQP',bounds=bnds,tol=1e-6)
        xnew = rslt.x[0]
        fnew = rslt.fun
        rho,fhiNew = fscaled.get_trust_region_ratio(xnew)
        if rho<=eta1 or rho>=eta3:
            delta *= c1
        elif eta2<rho<eta3:
            if abs(x0-xnew)==delta:
                gamma = c2
            else:
                gamma = 1.0
            delta *= gamma
        err = abs(x0-xnew)
        
        x0 = xnew
        print 'rho:%.4f\tx:%.4f\tf:%.4f\tdelta:%.4e'%(rho, xnew, fnew, delta)
        niter += 1
    y = array([fscaled(_x) for _x in x])
    output.write_array('YgvfmHybr',y)
    output.close()


#    fscaled.funcHi.display()
#    fscaled.funcLo.display()
    fscaled.fAdd.funcHi.display()
    fscaled.fAdd.funcLo.display()
    print niter

if __name__=="__main__":
    vcm_test_1d()