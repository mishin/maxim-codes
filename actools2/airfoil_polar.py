# -*- coding: utf-8 -*-
"""
Created on Thu Jan 09 21:39:53 2014

@author: Maxim
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d
from scipy.optimize import fminbound


class Interp2D:
    """
    2D interpolation using rectangular bivariate spline
    scipy.interpolate.RectBivariateSpline
    
    Returns interpolated value if requested value is between bounds otherwise 
    returns value at the boundary.
    """
    def __init__(self,x,y,z):
        self._model = RectBivariateSpline(x,y,z)
        self.x = x
        self.xmax = max(x)
        self.xmin = min(x)
        self.y = y
        self.ymax = max(y)
        self.ymin = min(y)
    def __call__(self,x,y):
        if x>self.xmax: x = self.xmax
        if x<self.xmin: x = self.xmin
        if y>self.ymax: y = self.ymax
        if y<self.ymin: y = self.ymin
        #FIXME: returns 2D array
        return self._model(x,y)[0,0]


class AirfoilPolar:
    def __init__(self):
        self.source = None
        self.cl = None
        self.cd = None
        self.cm = None
        self.alpha = None
        self.Mach = None
        self.Re  = None
        self.cdp = None
        self.TU = None
        self.TL = None
        self.SU = None
        self.SL = None
        self.LD = None
        self.AC = None
        self.CP = None

    def _calc_clmax(self,simple=True):
        #FIXME: clmax trough optimization is not working
        """ calculate maximum lift coefficient for all mach numbers stored """
        if hasattr(self.Mach, '__iter__'):
            n = len(self.Mach)
            self.clmax = np.zeros(n)
            self.alphaClmax = np.zeros(n)
            for i in range(n):
                _cl = self.cl[i]
                self.clmax[i],self.alphaClmax = self._get_clmax(_cl,self.alpha,simple)
        else:
            self.clmax,self.alphaClmax = self._get_clmax(self.cl,self.alpha,simple)

    def _get_clmax(self,cl,alpha,simple=True):
        if simple:
            idx = np.argmax(cl)
            return cl[idx], alpha[idx]
        else:
            clCurve = interp1d(-cl,alpha,'cubic')
            alphaClmax = fminbound(clCurve,alpha[1],alpha[-2])
            clmax = -clCurve(alphaClmax)
            return clmax, alphaClmax

    def _create_splines(self):
        """
        creates 1D or 2D spline of lift, drag and moment coefficient using Interp2D
        
        Parameters
        ----------
        
        mach : bool
            if array of Mach numbers exists then generates 2D spline of aero        coefficients vs. angle of attack and Mach number, otherwise generates 1D spline 
            for angle of attack at Mach number stored in self.Mach[0]
        """
        if hasattr(self.Mach,'__iter__') and len(self.Mach)>1:
            self.clAlpha = Interp2D(self.Mach,self.alpha,self.cl)
            self.cdAlpha = Interp2D(self.Mach,self.alpha,self.cd)
            self.cmAlpha = Interp2D(self.Mach,self.alpha,self.cm)
        else:
            self.clAlpha = interp1d(self.alpha,self.cl,'cubic')
            self.cdAlpha = interp1d(self.alpha,self.cd,'cubic')
            self.cmAlpha = interp1d(self.alpha,self.cm,'cubic')