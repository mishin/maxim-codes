# -*- coding: utf-8 -*-
"""
Created on Wed May 09 13:45:19 2012

@author: Maxim
"""
import numpy as ny

class airfoil:
    def __init__(self,Up=[],Lo=[],name = 'airfoil'):
        self.Up = Up
        self.Lo = Lo
        self.Name = name
        if Lo!=[]:
            self.thickness = ny.max(Up[:,1]-Lo[:,1])
            thicknessPt = ny.argmax(Up[:,1]-Lo[:,1])
            self.camber = Up[thicknessPt,1] - self.thickness/2

    def readAirfoil(self,filePath):
        file = open(filePath,'r')
        self.Name = file.readline()
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
        self.Up, self.Lo = Up, Lo
    
    def writeAirfoil(self,filePath):
        file = open(filePath,'wt')
        file.write(self.Name + '\n')
        Up = self.Up[ ::-1,:]
        Lo = self.Lo
        for ii in range(len(Up)):
            file.write(('%.6f    %.6f\n' %(Up[ii][0],Up[ii][1])))
        for ii in range(1,len(Lo)):
            file.write(('%.6f    %.6f\n' %(Lo[ii][0],Lo[ii][1])))
        file.close()
