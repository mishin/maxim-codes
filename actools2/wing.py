# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 22:29:53 2014

@author: Maxim
"""
import numpy as np
from scipy.interpolate import interp1d
import airfoil
from paths import MyPaths

path = MyPaths()


class ControlSurface:
    def __init__(self,location,wingChords,wingSpans,inverted=False,symmetric=True):
        self.location    = np.asarray(location)
        self.inverted    = bool(inverted)
        self.symmetric    = bool(symmetric)
        self._wingChords = np.asarray(wingChords)
        self._wingSpans  = np.asarray(wingSpans)
        self.type = 0 # option for flap only
        self._process_data()

    def _process_data(self):
        """ Note area - cs area per side """
        self.chords = self._wingChords * (1-self.location)
        secStart = np.argmax(self.location<1)
        secEnd = len(self.location) - np.argmax(self.location[::-1]<1)-1
        area = 0.0
        for i in range(secStart,secEnd):
            area += (self.chords[i]+self.chords[i+1])*self._wingSpans[i]
        self.areaPerSide = area
        if self.symmetric:
            self.area = 2.0*self.areaPerSide
        else:
            self.area = self.areaPerSide


class EquivalentWing(object):
    def __init__(self,wing):
        pass
    
    def _get_linear_equivalent(self):
        pass
    def _get_angular_equivalent(self):
        pass

class Wing(object):
    def __init__(self):
        self.locationLE        = np.zeros(3)
        self.segSpans          = None
        self.chords            = None
        self.airfoilNames      = None
        self.airfoils          = list()
        self.secOffset         = None
        self.segDihedral       = None
        self.segAreas          = None
        self.secTwist          = None
        self.incidence         = 0.0
        self.secAngles         = None
        self.leadingEdge       = np.zeros(3)
        self.material          = None
        self.MAC               = 0.0
        self.MAClocation       = np.zeros(2) # x,y
        self.nSeg              = 0
        self.secThickness      = None
        self.csArea            = 0
    
    def locate_on_wing(self,chordRatio,spanRatio):
        """
        returns coordinate of the point that is located on the wing with 
        given span and chord ratio
        """
        y = self.span/2.0 * spanRatio
        chords = interp1d(self.secApex[:,1],self.chords,'linear')
        xLE = interp1d(self.secApex[:,1],self.secApex[:,0],'linear')
        zLE = interp1d(self.secApex[:,1],self.secApex[:,2],'linear')
        x = xLE(y) + chords(y)*chordRatio
        z = zLE(y)
        return np.array([x,y,z])
        
    def _process_data(self,updateAirfoils=False):
        self.nSeg = len(self.segSpans)
        self.nSec = self.nSeg+1
        self.segDihedralRad = np.radians(self.segDihedral)
        if updateAirfoils:
            self._load_airfoils()
        self._calc_geometry_data()
        self._calc_2d_points_data()
    
    def _calc_2d_points_data(self):
        npts = 5*self.nSeg
        x = np.zeros([npts])
        y = np.zeros([npts])
        for i in range(self.nSeg):
            _n = 5*i
            x[_n] = None
            y[_n] = None
            x[_n+1] = self.secApex[i,0]
            y[_n+1] = self.secApex[i,1]
            x[_n+2] = self.secApex[i+1,0]
            y[_n+2] = self.secApex[i+1,1]
            x[_n+3] = self.secApex[i+1,0] + self.chords[i+1]
            y[_n+3] = self.secApex[i+1,1]
            x[_n+4] = self.secApex[i,0] + self.chords[i]
            y[_n+4] = self.secApex[i,1]
        ynew = np.hstack([np.flipud(x[1:]), x])
        self.y2d = -ynew
        self.x2d = np.hstack([-np.flipud(y[1:]), y])
    
    def set_elevon(self,ailLocation):
        self.elevon = ControlSurface(ailLocation,self.chords,self.segSpans,False)
    
    def set_flap(self,flapLocation):
        self.flap = ControlSurface(flapLocation,self.chords,self.segSpans,False)

    def _load_airfoils(self,xlsPath=None):
        if xlsPath==None:
            xlsPath = path.db.airfoil
        self.secThickness = np.zeros(self.nSec)
        for i,name in enumerate(self.airfoilNames):
            self.airfoils.append(airfoil.load(name))
            self.secThickness[i] = self.airfoils[i].thickness

    def _calc_geometry_data(self):
        self._calc_apex()
        self._calc_segment_data()
        self._calc_mac()
        self._calc_span()
        self._calc_wetted_area()
        self._calc_angles()
        self._calc_elastic_axis_sweep()
        self._calc_controls()
        self._calc_equiv_trapz_wing()
        self.taperRatio = self.chords[-1]/self.chords[0]
    
    def _calc_controls(self):
        self.csArea = self.elevon.area + self.flap.area
    
    def _get_cs_area(self,location):
        startSec = np.argmax(location<1)
        endSec = len(location) - np.argmax(location[::-1]<1)-1
        area = 0.0
        for i in range(startSec,endSec):
            c1 = self.chords[i]*(1-location[i])
            c2 = self.chords[i+1]*(1-location[i+1])
            area += (c1 + c2)*self.segSpans[i]
        return area

    def _calc_span(self):
        self.span = 2.0*np.sum(self.segSpans)
        self.aspectRatio = self.span**2.0/self.area

    def _calc_apex(self):
        # calculate leading edge point of each section - section incidence is 
        # not considered
        secApex = np.zeros([3,self.nSec])
        secApex[0,1:] = self.secOffset.cumsum()
        secApex[1,1:] = self.segSpans.cumsum()
        secApex[2,1:] = secApex[0,1:]*np.tan(np.radians(self.segDihedral))
        self.secApex = np.transpose(secApex)
    
    def _calc_angles(self):
        # calculates angle of each section
        self.secAngles = np.zeros(self.nSec)
        self.secAngles[0] = self.incidence
        for i,twist in enumerate(self.secTwist):
            self.secAngles[i+1] = self.secAngles[i]+twist
    
    def _calc_elastic_axis_sweep(self,FSfrontSpar=0.3,FSaftSpar=0.7):
        """
        Assumed values are location of front and aft spars. Typically the elastic 
        axis is located between two spars. 
        NOTE: This calculation is very approximate.
        """
        eaLoc = (FSaftSpar-FSfrontSpar)/2.0
        xr = self.secApex[0,0] + self.chords[0]*eaLoc
        xt = self.secApex[-1,0] + self.chords[-1]*eaLoc
        yr = self.secApex[0,1]
        yt = self.secApex[-1,1]
        self.sweepElasticRad = np.arctan((yt-yr)/(xt-xr))
        self.sweepElasticDeg = np.degrees(self.sweepElasticRad)
        dx, dy = xt-xr, yt-yr
        self.structuralSpan = (dx*dx+dy*dy)**0.5

    def _calc_mac(self):
        """ calculates mean aerodynamic chord and it's location"""
        self.segMAC = np.zeros(self.nSeg)
        self.segAreas = np.zeros(self.nSeg)
        def mac_of_one_segment(c1, c2, x1, x2, y1, y2):
            span = y2 - y1
            area = 0.5*(c1+c2)*span
            taper = c2/c1
            taper1 = taper + 1.0
            frac = (taper+taper1)/(3.0*taper1)
            mac = c1*(taper*taper + taper + 1.0) / (1.5*taper1)
            ymac = y1 + frac*span
            xlemac = x1 + frac*(x2-x1)
            return area, mac, ymac, xlemac
        area, mac, ymac, xlemac = 0.0, 0.0, 0.0, 0.0
        for i in range(self.nSeg):
            c1,c2 = self.chords[i], self.chords[i+1]
            x1,x2 = self.secApex[i,0], self.secApex[i+1,0]
            y1,y2 = self.secApex[i,1], self.secApex[i+1,1]
            aseg, macseg, yseg, xseg = mac_of_one_segment(c1,c2,x1,x2,y1,y2)
            self.segMAC[i] = macseg
            self.segAreas[i] = aseg
            area   += aseg
            mac    += macseg*aseg
            ymac   += yseg*aseg
            xlemac += xseg*aseg
        self.MAC = mac/area
        self.MAClocation = np.array([xlemac/area, ymac/area])
        self.area = 2.0*area

    def _calc_wetted_area(self):
        wettedArea = 0.0
        for i in range(len(self.segSpans)-1):
            root = self.airfoils[i].length*self.chords[i]
            tip = self.airfoils[i+1].length*self.chords[i+1]
            wettedArea += (root+tip)*self.segSpans[i+1]
        self.wettedArea = wettedArea
    
    def _calc_equiv_trapz_wing(self):
        #FIXME: calculate in same routine with segment sweep angles
        tc = 0.0
        tcLoc = 0.0
        camber = 0.0
        tanLambdaLE = 0.0
        cosLambdaC2 = 0.0
        cosLambdaC4 = 0.0
        tanLambdaTE = 0.0
        for i,segArea in enumerate(self.segAreas):
            tcAvg = 0.5 *(self.airfoils[i].thickness + self.airfoils[i+1].thickness)
            camberAvg = 0.5 * (self.airfoils[i].camber + self.airfoils[i+1].camber)
            tcLocAvg = 0.5 * (self.airfoils[i].thicknessLocation + self.airfoils[i+1].thicknessLocation)
            tc += tcAvg**2.*segArea
            camber += camberAvg**2. *segArea
            tcLoc += tcLocAvg**2. *segArea
            tanLambdaLE += np.tan(self.segSweepLErad[i])*segArea
            cosLambdaC2 += np.cos(self.segSweepC2rad[i])*segArea
            cosLambdaC4 += np.cos(self.segSweepC4rad[i])*segArea
            tanLambdaTE += np.tan(self.segSweepTErad[i])*segArea
        self.equivThickness = (2.*tc/self.area)**0.5
        self.equivCamber = (2.*camber/self.area)**0.5
        self.equivThicknessLoc = (2.*tcLoc/self.area)**0.5
        self.equivSweepLErad = np.arctan(2.*tanLambdaLE/self.area)
        self.equivSweepLEdeg = np.degrees(self.equivSweepLErad)
        self.equivSweepC2rad = np.arccos(2.*cosLambdaC2/self.area)
        self.equivSweepC2deg = np.degrees(self.equivSweepC2rad)
        self.equivSweepC4rad = np.arccos(2.*cosLambdaC4/self.area)
        self.equivSweepC4deg = np.degrees(self.equivSweepC4rad)
        self.equivSweepTErad = np.arctan(2.*tanLambdaTE/self.area)
        self.equivSweepTEdeg = np.degrees(self.equivSweepTErad)
        self.equivLEradius = 1.1019*self.equivThickness**2.0
    
    def _calc_segment_data(self):
        self.segSweepLErad = np.arctan(self.secOffset/self.segSpans)
        self.segSweepC2rad = np.zeros(self.nSeg)
        self.segSweepC4rad = np.zeros(self.nSeg)
        self.segSweepTErad = np.zeros(self.nSeg)
        for i in range(self.nSeg):
            dApex = self.secApex[i+1,0] - self.secApex[i,0]
            dc = self.chords[i+1] - self.chords[i]
            dx1 = dApex + dc*0.5
            dx2 = dApex + dc*0.25
            dx3 = dApex + dc
            self.segSweepC2rad[i] = np.arctan(dx1/self.segSpans[i])
            self.segSweepC4rad[i] = np.arctan(dx2/self.segSpans[i])
            self.segSweepTErad[i] = np.arctan(dx3/self.segSpans[i])
        self.segSweepLEdeg = np.degrees(self.segSweepLErad)
        self.segSweepC2deg = np.degrees(self.segSweepC2rad)
        self.segSweepC4deg = np.degrees(self.segSweepC4rad)
        self.segSweepTEdeg = np.degrees(self.segSweepTErad)
    
    def get_max_segment_length(self,width):
        """
        Calculates maximum length of payload (assumed to be rectangular) 
        that can be located inside wing central part (1st segment). 
        Method returns front fuselage station (x-coordinate) of payload and 
        maximum length.
        
        Parameters
        ----------
        
        width : float, m
            payload width
        """
        width *= 0.5

        if width>self.segSpans[0]:
            return 0.0, 0.0
        
        fsFront = width * np.tan(self.segSweepLErad[0])
        if self.segSweepTEdeg[0]>=0:
            fsAft = self.secApex[0,0] + self.chords[0]
        else:
            fsAft = self.chords[0] + width *np.tan(self.segSweepTErad[0])
        length = fsAft - fsFront
        return fsFront, length