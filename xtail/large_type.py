# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 23:13:19 2013

@author: Maxim
"""
from numpy import array, linspace, cos, sin, pi, hstack, vstack, flipud, transpose, concatenate, zeros
import matplotlib.pyplot as plt
from scipy.integrate import simps
import dbTools

class LargeTypeASM:
    def __init__(self):
        """
        Large type air to surface missile configuration
        """
        self.body = Fuselage()
        self.wing = Wing()
        self.tail = Tail()
        self.CG = zeros(3)
        self.mass = 0.0
    
    def read_from_db(self,name,dbFile=None):
        if dbFile==None:
            dbFile = 'LargeTypeASM.xls'
        wb = dbTools.loadDB(dbFile,'r')
        sh = dbTools.readDB(wb.selectByName(name))
        
        i = sh.findHeader('SECTION: BODY')
        self.body.xprofile = sh.readRow(i+1,1,True)
        self.body.yprofile = sh.readRow(i+2,1,True)
        self.body.set_length(sh.readRow(i+3,1))
        self.body.set_diameter(sh.readRow(i+4,1))
        
        i = sh.findHeader('SECTION: WING')
        self.wing.chords         = sh.readRow(i+1,1,True)
        self.wing.incidenceAngle = sh.readRow(i+2,1,True)
        self.wing.sectionSpan    = sh.readRow(i+3,1)
        self.wing.leSweep        = sh.readRow(i+4,1)
        self.wing.dihedral       = sh.readRow(i+5,1)
        self.wing.airfoilPath    = sh.readRow(i+6,1)
        self.wing.location       = sh.readRow(i+7,1)
        
        i = sh.findHeader('SECTION: TAIL')
        self.tail.chords             = sh.readRow(i+1,1,True)
        self.tail.incidenceAngle     = sh.readRow(i+2,1,True)
        self.tail.sectionSpan        = sh.readRow(i+3,1)
        self.tail.leSweep            = sh.readRow(i+4,1)
        self.tail.dihedral           = sh.readRow(i+5,1)
        self.tail.airfoilPath        = sh.readRow(i+6,1)
        self.tail.location           = sh.readRow(i+7,1)
        self.tail.xangle             = sh.readRow(i+8,1)
        self.tail.centerOffset       = sh.readRow(i+9,1)
        self.tail.elevatorChordRatio = sh.readRow(i+10,1)
    
    def write_to_db(self,name,dbFile=None):
        if dbFile==None:
            dbFile = 'LargeTypeASM.xls'
        db = dbTools.loadDB(dbFile,mode='w')
        sh = dbTools.writeDB(db.add_sheet(name))
        sh.writeRow('SECTION: BODY')
        sh.writeRow('x profile',self.body.xprofile)
        sh.writeRow('y profile',self.body.yprofile)
        sh.writeRow('length',self.body.length)
        sh.writeRow('diameter',self.body.diameter)
        sh.writeRow('SECTION: WING')
        sh.writeRow('chords',self.wing.chords)
        sh.writeRow('incidence',self.wing.incidenceAngle)
        sh.writeRow('span',self.wing.sectionSpan)
        sh.writeRow('leading edge sweep',self.wing.leSweep)
        sh.writeRow('dihedral',self.wing.dihedral)
        sh.writeRow('airfoil',self.wing.airfoilPath)
        sh.writeRow('location X,Z',self.wing.location)
        sh.writeRow('SECTION: TAIL')
        sh.writeRow('chords',self.tail.chords)
        sh.writeRow('incidence',self.tail.incidenceAngle)
        sh.writeRow('span',self.tail.sectionSpan)
        sh.writeRow('leading edge sweep',self.tail.leSweep)
        sh.writeRow('dihedral',self.tail.dihedral)
        sh.writeRow('airfoil',self.tail.airfoilPath)
        sh.writeRow('location X,Z',self.tail.location)
        sh.writeRow('tail X-angle',self.tail.xangle)
        sh.writeRow('tail center offset', self.tail.centerOffset)
        sh.writeRow('elevator chord ratio',self.tail.elevatorChordRatio)
        db.save_db()

class Fuselage():
    def __init__(self):
        self.xprofile = list()
        self.yprofile = list()
        self.length = 0.0
        self.diameter = 0.0
        self.wettedArea = 0.0
        self.frontalArea = 0.0

    def _analyze_geometry(self):
        self._diameter = 2.0*max(self.yprofile)
        self._length = self.xprofile[-1] - self.xprofile[0]
        self.frontalArea = self.diameter**2*pi/4
        _c = array([2.0*pi*y for y in self.yprofile])
        self.wettedArea = simps(_c,self.xprofile)
    
    def set_length(self,newLength):
        if self.length==0.0:
            self._analyze_geometry()
            self.length = self._length
        ratio = newLength/self.length
        self.length = newLength
        self.xprofile *= ratio
    
    def set_diameter(self,newDiameter):
        if self.diameter==0.0:
            self._analyze_geometry()
            self.diameter = self._diameter
        ratio = newDiameter/self.diameter
        self.diameter = newDiameter
        self.yprofile *= ratio

    def set_naca_profile(self,length,diameter,trailingEdgeGap=0.0,npts=30):
        """
        Creates x, y profile points using NACA 4 series equation.
        
        Parameters
        ----------
        
        length : float, m
            length of the fuselage
        diameter : float, m
            diameter of the fuselage. Note that final diameter may be different 
            due to finite trailing edge gap.
        trailingEdgeGap : float, m
            trailing edge gap of the profile
        nPts : int
            number of points
        """
        tc = diameter / length
        x = flipud(self._cos_distribution(linspace(0,1,npts)))
        y = tc/0.2*(0.2969*(x)**0.5-0.1281*x-0.3516*x**2+0.2843*x**3-0.1015*x**4)
        self.xprofile = x*length
        self.yprofile = y*length + x*trailingEdgeGap
        xcrd = hstack([flipud(self.xprofile), self.xprofile[1:]])
        ycrd = hstack([flipud(self.yprofile),-self.yprofile[1:]])
        self.coord = transpose(vstack([xcrd,ycrd]))
    def _cos_distribution(self,x):
        def cos_curve(x):
            return (cos(x*pi)+1.0)/2.0
        return array([cos_curve(xx) for xx in x])


class Wing:
    def __init__(self, nSec=2):
        self.airfoil        = object()
        self.airfoilPath    = ''
        self.chords         = zeros(nSec)
        self.sectionSpan    = zeros(nSec-1)
        self.leSweep        = zeros(nSec)
        self.incidenceAngle = zeros(nSec)
        self.location       = zeros(2)
        self.dihedral       = zeros(nSec-1)
    
    def _analyze_geometry(self):
        self.mac = 0.0
        self.span = 0.0
        self.area = 0.0
        for i,b in enumerate(self.sectionSpan):
            self.area += (self.chords[i]+self.chords[i+1])*b / 2.0
            self.span += b
        self.area *= 2.0
        self.span *= 2.0
        

class Tail(Wing):
    def __init__(self):
        Wing.__init__(self)
        self.numberOfTails = 4
        self.xangle = 45.0
        self.centerOffset = 0.1
        self.elevatorChordRatio = 0.5

def test():
    gb = LargeTypeASM()
    gb.read_from_db('example1')
    gb.body.set_naca_profile(2,0.5)
    gb.wing.sectionSpan
    gb.write_to_db('new')
    

if __name__=="__main__":
    test()