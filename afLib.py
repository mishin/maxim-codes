# -*- coding: utf-8 -*-
"""
Created on Wed May 09 13:45:19 2012

@author: Maxim
"""
import numpy as ny
from scipy.interpolate import interp1d

class airfoil:
    def __init__(self,up=[],lo=[],mainSec=[],flap=[],name = 'airfoil'):
        self.up = up
        self.lo = lo
        self.name = name
        self.hasFlap = False
        if up!=[] and lo!=[]: self.__interp__()
        if mainSec!=[] and flap!=[]:
            self.hasFlap = False
            self.mainSec = mainSec
            self.flap = flap

    def calcGeom(self):
        self.thickness = ny.max(self.up-self.lo)
        
    def getThicknessAtX(self,x):
        return self.upCurve(x) - self.loCurve(x)
        
    def __interp__(self):
        self.upCurve = interp1d(self.up[:,0],self.up[:,1],'cubic')
        self.loCurve = interp1d(self.lo[:,0],self.lo[:,1],'cubic')

def readAirfoil(filePath):
    file = open(filePath,'r')
    AirfoilName = file.readline()
    lines = file.readlines()
    Npts = len(lines)
    coord = ny.zeros((Npts,2))
    ii = 0
    for line in lines:
        segLine = line.split()
        coord[ii][0],coord[ii][1] = float(segLine[0]),float(segLine[1])
        if coord[ii][0]==0 and coord[ii][1]==0:
            zeroPt = ii
        ii = ii+1
    Up = ny.zeros((zeroPt+1,2))
    Lo = ny.zeros((Npts-zeroPt,2))
    for ii in range(zeroPt+1):
        Up[ii][0],Up[ii][1]=coord[zeroPt-ii][0],coord[zeroPt-ii][1]
    for ii in range(zeroPt,Npts):
        Lo[ii-zeroPt][0],Lo[ii-zeroPt][1] = coord[ii][0],coord[ii][1]
    file.close()
    return airfoil(Up,Lo,AirfoilName)

def writeAirfoil(filePath,airfoil):
    file = open(filePath,'wt')
    file.write(airfoil.name)
    file.write('\n')
    Up = ny.flipud(airfoil.up)
    for ii in range(len(Up)):
        file.write(('%.6f    %.6f\n' %(Up[ii][0],Up[ii][1])))
    for ii in range(1,len(airfoil.lo)):
        file.write(('%.6f    %.6f\n' %(airfoil.lo[ii][0],airfoil.lo[ii][1])))
    file.close()

def writeFlap(filePath, airfoil):
    file = open(filePath,'wt')
    file.write(airfoil.name + '\n')
    for mainSecPt in airfoil.mainSec:
        file.write(('%.6f    %.6f\n' %(mainSecPt[0],mainSecPt[1])))
    file.write('9999.9    9999.9\n')
    for flapPt in airfoil.flap:
        file.write(('%.6f    %.6f\n' %(flapPt[0],flapPt[1])))
    file.close()