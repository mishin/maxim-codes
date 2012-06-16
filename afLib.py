# -*- coding: utf-8 -*-
"""
Created on Wed May 09 13:45:19 2012

@author: Maxim
"""
import numpy as ny
from scipy.interpolate import interp1d

class airfoil:
    def __init__(self,up=[],lo=[],name = 'airfoil'):
        self.up = up
        self.lo = lo
        self.name = name
        if up!=[] and lo!=[]:
            self.upCurve = interp1d(up[:,0],up[:,1],'cubic')
            self.loCurve = interp1d(lo[:,0],lo[:,1],'cubic')
    def calcGeom(self):
        self.thickness = ny.max(self.up-self.lo)
    def getThicknessAtX(self,x):
        return self.upCurve(x) - self.loCurve(x)

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
    file.write(airfoil.name + '\n')
    Up = airfoil.up[ ::-1,:]
    for ii in range(len(Up)):
        file.write(('%.6f    %.6f\n' %(Up[ii][0],Up[ii][1])))
    for ii in range(1,len(airfoil.lo)):
        file.write(('%.6f    %.6f\n' %(airfoil.lo[ii][0],airfoil.lo[ii][1])))
    file.close()

def writeFlap(filePath,body, flap,name='flap'):
    file = open(filePath,'wt')
    file.write(name + '\n')
    for bodyPt in body:
        file.write(('%.6f    %.6f\n' %(bodyPt[0],bodyPt[1])))
    file.write('9999.9    9999.9\n')
    for flapPt in flap:
        file.write(('%.6f    %.6f\n' %(flapPt[0],flapPt[1])))
    file.close()