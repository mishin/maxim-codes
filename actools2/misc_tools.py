# -*- coding: utf-8 -*-
"""
Created on Thu Jan 09 11:50:25 2014

@author: Maxim
"""

from datetime import datetime
import numpy as np

def normalize(X,lb,ub,normLb=-1.,normUb=1.):
    delta = normUb - normLb
    #TODO

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


if __name__=="__main__":
    run_test1()