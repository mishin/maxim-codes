# -*- coding: utf-8 -*-
"""
Created on Wed Apr 09 22:16:43 2014

@author: Maxim
"""
import numpy as np
from mayavi import mlab

class AircraftDisplay(object):
    def __init__(self):
        self.meshes = list()
        self.points = list()
        self.text = list()
    
    def add_wing(self,wing,nSurfPts=30,opacity=0.7):
        nSeg = wing.nSeg
        xMesh = None
        yMesh = None
        zMesh = None
    
    def add_points(self,x,y,z,text=None):
        pass
    
    def display(self):
        pass