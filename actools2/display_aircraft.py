# -*- coding: utf-8 -*-
"""
Created on Wed Apr 09 22:16:43 2014

@author: Maxim
"""
import numpy as np
from mayavi import mlab

from geometry import rotate_2d

def flying_wing_display(aircraft,showAxes=False):
    adisp = AircraftDisplay()
    adisp.add_wing(aircraft.wing)
    adisp.display(showAxes=showAxes)

class AircraftDisplay(object):
    def __init__(self):
        self.xMesh = list()
        self.yMesh = list()
        self.zMesh = list()
        self.meshOpacity = list()
        self.meshColor = list()
        self.points = list()
        self.text = list()
    
    def add_wing(self,wing,nSurfPts=25,opacity=0.7,color=(0.49,1,0.83)):
        nSeg = wing.nSeg
        nSec2 = 2*nSeg+1
        xMesh = np.zeros([nSec2,nSurfPts])
        yMesh = np.zeros([nSec2,nSurfPts])
        zMesh = np.zeros([nSec2,nSurfPts])
        for i in range(nSeg+1):
            afPts = wing.airfoils[i].redim(nSurfPts,overwrite=False)
            section = afPts*wing.chords[i]
            axis1 = np.array([wing.chords[i]*0.25, 0.0])
            section = rotate_2d(section,axis1,wing.secAngles[i])
            secX = section[:,0] + wing.secApex[i,0]
            secY = np.zeros(nSurfPts) + wing.secApex[i,1]
            secZ = section[:,1] + wing.secApex[i,2]
            if i==0:
                xMesh[nSeg] = np.copy(secX)
                yMesh[nSeg] = np.copy(secY)
                zMesh[nSeg] = np.copy(secZ)
            else:
                xMesh[nSeg-i] = np.copy(secX)
                xMesh[nSeg+i] = np.copy(secX)
                yMesh[nSeg-i] = np.copy(-secY)
                yMesh[nSeg+i] = np.copy(secY)
                zMesh[nSeg-i] = np.copy(secZ)
                zMesh[nSeg+i] = np.copy(secZ)
        self.xMesh.append(xMesh)
        self.yMesh.append(yMesh)
        self.zMesh.append(zMesh)
        self.meshOpacity.append(opacity)
        self.meshColor.append(color)
            
    def add_points(self,x,y,z,text=None):
        pass
    
    def display(self,bgcolor=None,showAxes=False):
        mlab.figure(bgcolor=bgcolor)
        for x,y,z,op,col in zip(self.xMesh,self.yMesh,self.zMesh,self.meshOpacity,self.meshColor):
            mlab.mesh(x,y,z,opacity=op,color=col)
        if showAxes:
            mlab.axes()
        mlab.show()