# -*- coding: utf-8 -*-
"""
Created on Thu Jan 09 11:50:25 2014

@author: Maxim
"""

from datetime import datetime
import numpy as np
from scipy.interpolate import Rbf


class RbfMod():
    def __init__(self,x,y):
        if x.ndim>1:
            x = np.transpose(x)
            x = self._get_tuple(x)
            args = x + (y,)
        else:
            args = (x,y)
        self.rbf = Rbf(*args)

    def __call__(self,x):
        if hasattr(x,'__iter__'):
            x = self._get_tuple(x)
            return self.rbf(*x)
        else:
            return self.rbf(x)

    def _get_tuple(self,xArray):
        xTuple = tuple()
        for x in xArray:
            xTuple += (x,)
        return xTuple


class SaveTextData():
    def __init__(self,path):
        self.fid = open(path,'wt')
        self.ndecimal = 8
        self.code = '%e'
        self.fid.write('import numpy as np\n\n')

    def write_array(self,inputArray,name):
        self.fid.write('%s = np.array(['%name)
        self._write_array(inputArray)
        self.fid.write('])\n')
    
    def _write_array(self,inputArray):
        for val in inputArray[:-1]:
            self.fid.write(self.code%val)
            self.fid.write(', ')
        self.fid.write(self.code%inputArray[-1])
        
    def close(self):
        self.fid.close()


def read_tabulated_data(path,header=False,firstRow=0):
    if header:
        return read_tabulated_data_with_header(path,firstRow)
    else:
        return read_tabulated_data_without_header(path,firstRow)
            
def read_tabulated_data_with_header(path,firstRow=0):
    fid = open(path,'rt')
    data = {}
    segHeader = fid.readline().split()
    ncol = len(segHeader)
    lines = fid.readlines()
    nrow = len(lines)
    values = np.zeros([nrow,ncol])
    for i,line in enumerate(lines):
        seg = line.split()
        for j,val in enumerate(seg):
            values[i,j] = float(val)
    for i,varname in enumerate(segHeader):
        data[str(varname)] = values[:,i]
    return data

def read_tabulated_data_without_header(path,firstRow=0):
    fid = open(path,'rt')
    if not firstRow==0:
        for i in range(firstRow):
            fid.readline()
    lines = fid.readlines()
    nrow = len(lines)
    ncol = len(lines[0].split())
    values = np.zeros([nrow,ncol])
    for i,line in enumerate(lines):
        seg = line.split()
        for j,val in enumerate(seg):
            values[i,j] = float(val)
    return values

def write_tabulated_data_without_header(path,data):
    fid = open(path,'wt')
    for line in data:
        for val in line:
            fid.write('%.6f\t'%val)
        fid.write('\n')
    fid.close()

class Normalization:
    def __init__(self,lb,ub,xMin=-1.,xMax=1.):
        self.lb = np.array(lb,dtype=float)
        self.ub = np.array(ub,dtype=float)
        self.xMin = float(xMin)
        self.xMax = float(xMax)
        self._delta1 = self.ub - self.lb
        self._delta2 = self.xMax - self.xMin
        self._ratio1 = self._delta2/self._delta1

    def normalize(self,x):
        x = np.array(x,dtype=float)
        return (x-self.lb)*self._ratio1+self.xMin

    def denormalize(self,xnorm):
        xnorm = np.array(xnorm,dtype=float)
        return (xnorm-self.xMin)/self._ratio1 + self.lb
        
        

class Timer:
    def __init__(self,header=''):
        self.timeStart     = 0.0
        self.timeFinish    = 0.0
        self.timeLapStart  = 0.0
        self.timeLapFinish = 0.0
        self.timeLastLap   = 0.0
        self.durationTotal = 0.0
        self.durationLap   = list()
        self.headersLap    = list()
        self._nLaps        = 0
        self._header       = header
        self.start(header)
    
    def start(self,header=''):
        self.timeStart     = datetime.now()
        self.timeLapStart  = self.timeStart
        self.timeLoopStart = self.timeStart

    def lap(self,header=None,display=True):
        self._nLaps += 1
        self.timeLapFinish = datetime.now()
        duration = self.timeLapFinish - self.timeLapStart
        self.timeLapStart = self.timeLapFinish
        self.durationLap.append(duration.total_seconds())
        self.durationTotal = sum(self.durationLap)
        if header==None:
            header = 'task %d'%self._nLaps
        self.headersLap.append(header)
        if display:
            self._disp_lap(-1)
    
    def stop(self,header='',detailed=False):
        self.lap(header, detailed==False)
        print 'Timer =>',self._header
        self.display(detailed)

    def display(self,detailed=False):
        if detailed:
            self._disp_detailed()
            self._disp_total()
        else:
            self._disp_total()

    def _disp_detailed(self):
        for i in range(self._nLaps):
            self._disp_lap(i)

    def _disp_total(self):
        print '{:<14} : '.format('TOTAL'), self.durationTotal, 'sec'

    def _disp_lap(self,idx):
        label = self.headersLap[idx]
        duration = self.durationLap[idx]
        print '{:<14} : '.format(label), duration, 'sec'


# --- debug ---
def run_test1():
    timer = Timer('test task')
    a = np.arange(0,1e5,1)
    for j in range(15):
        for aa in a:
            b = np.sin(aa)
    timer.lap('sin')
    for aa in a:
        b = np.sin(aa)
    timer.stop('cosine',True)

def run_test2():
    lb = np.array([40, 40, 6.00, 3.00, 0.5, 1.00, 3.000, -4, -4])
    ub = np.array([60, 60, 7.50, 5.25, 1.8, 1.80, 3.200,  0,  0])
    x0 = np.array([55, 55, 6.91, 4.15, 1.1, 1.44, 3.115,  0, -3])
    n = Normalization(lb,ub,-1,1)
    xnorm1 = n.normalize(x0)
    xorig = n.denormalize(xnorm1)
    print xnorm1, xorig

def run_test3():
    data = read_tabulated_data('test_file2.txt',header=False)
    print data['col1']

def run_test4():
    a = np.linspace(0,20,3)
    fid = SaveTextData('new.txt')
    fid.write_array(a,'a')
    fid.close()

if __name__=="__main__":
    run_test2()