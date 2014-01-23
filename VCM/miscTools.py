# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:55:31 2012
Contains an assortment of tools for basic math and geometry
@author: daniel
"""
from datetime import datetime
from numpy import zeros

def read_tabulated_data_txt(path,idxHeader=0):
    """
    Read tabulated data from text file.
    
    Parameters
    ----------
    
    path : string
        path to text file
    idxHeader : integer
        index of line where header (variable names) is.
    """
    fid = open(path,'rt')
    lines = fid.readlines()[idxHeader:]
    fid.close()
    nData = len(lines)-1
    nVars = len(lines[0].split())
    values = zeros([nVars,nData])
    for i,line in enumerate(lines[1:]):
        seg = line.split()
        for j,val in enumerate(seg):
            values[j,i] = float(val)
    data = {}
    for k,varName in enumerate(lines[0].split()):
        data[varName] = values[k]
    return data

def normalize(X,lb,ub):
    """
    Normalize variable or vector to [-1,1]
    
    Parameters
    ----------
    
    X : float, array
        array of variable values
    lb : float, array
        lower bound
    ub : floar, array
        upper bound
    """
    #return (X-lb) / (ub-lb)
    return (2.*X-(lb+ub)) / (ub-lb)

def denormalize(X,lb,ub):
    """
    Denormalize variable or vector from [-1,1]
    
    Parameters
    ----------
    
    X : float, array
        array of normalized variable values
    lb : float, array
        lower bound
    ub : floar, array
        upper bound
    """
    #return (X+1.0)*(ub-lb)/2.0+lb
    return 0.5* (X*(ub-lb) + (lb+ub))

def integrateRomberg(f,a,b,tolerance):
    """    
    Romberg integration integrates a funciton with increasing numbers of
    intervals until a desired tolerance is reached.
    """
    import numpy
    iMax=100
    errorVal=numpy.inf
    R=numpy.zeros([iMax,iMax])
    i=0
    R[i,0]=(b-a)/2.0*(f(a)+f(b))    
    while errorVal>tolerance and i<iMax-1:
        n=i+1
        hn=(b-a)/(2.0**n)
        R1=0.0
        for k in range(1,2**(n-1)+1,1):
            R1=R1+f(a+(2*k-1)*hn)
        R[i+1,0]=R[i,0]/2+hn*R1
        for j in range(0,i+1,1):
            m=j+1
            R[i+1,j+1]=1.0/(4**m-1)*(4**m*R[i+1,j]-R[i,j])
        errorVal=numpy.abs(R[i+1,j+1]-R[i+1,j])
        i=i+1
    return R[i,i]
    
def frange(start, stop, n):
    """    
    Matlab-like increment generator
    """
    L = [0.0] * n
    nm1 = n - 1
    nm1inv = 1.0 / nm1
    for i in range(n):
        L[i] = nm1inv * (start*(nm1 - i) + stop*i)
    return L

def polyCentroid(x,y):
    """    
    Calculates the centroid of any closed, non-self-intersecting polygon
    """
    assert len(x)==len(y)
    Cx=0.0
    Cy=0.0
    A =0.0
    for i in range(len(x)-1):
        A+=0.5*(x[i]*y[i+1]-x[i+1]*y[i])
    for i in range(len(x)-1):
        Cx+=1./6./A*(x[i]+x[i+1])*(x[i]*y[i+1]-x[i+1]*y[i])
        Cy+=1./6./A*(y[i]+y[i+1])*(x[i]*y[i+1]-x[i+1]*y[i])   
    return Cx,Cy,A
    
class Timer:
    def __init__(self):
        self.timeStart = 0.0
        self.timeStop   = 0.0
        self.timeLoop  = list()
        self.durationTotal = 0.0
        self.durationLoop = list()
        self.loopName = list()
    def start(self):
        self.timeStart = datetime.now()
    def finish(self,disp=True):
        self.timeStop = datetime.now()
        self.durationTotal = self.timeStop - self.timeStart
        self._display(self.durationTotal)
    def _display(self,timedelta):
        d = timedelta.days
        s = timedelta.seconds
        out = ''
        if d!=0:
            out += '%d days  '%d
        out += '%d seconds'%s
        print out
    def loop(self,disp=True):
        pass

def testFcn(x):
    import numpy
    return 2.0/numpy.pi**.5*numpy.exp(-x**2)
def testRomberg():
    print integrateRomberg(testFcn,0,1,10**-9)

def testPolyCentroid():
    import numpy
    x=numpy.array([0,0,2,2,0])
    y=numpy.array([0,1,1,0,0])
    print polyCentroid(x,y)

def test_read():
    path = r'D:\polar_40.txt'
    data = read_tabulated_data_txt(path,1)
    print data['alpha']
    print data['cl']

if __name__=="__main__":
    test_read()
    #testRomberg()
    #testPolyCentroid()