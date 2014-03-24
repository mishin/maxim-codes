# -*- coding: utf-8 -*-
"""
Created on Thu Jan 09 21:39:53 2014

@author: Maxim
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d
from scipy.optimize import fminbound, bisect
import matplotlib.pyplot as plt

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


class AirfoilPolar1D:
    def __init__(self):
        self.name = None
        self.source = None
        self.cl = list()
        self.cd = list()
        self.cm = list()
        self.LD = list()
        self.alpha = list()
        self.cdp = list()
        self.Mach = None
        self.Re = None
        self.clmax = None
        self.cdmin = None
        self.alphaClmax = None
        self.alphaCdmin = None
        self._alphaCl = None
        self._alphaCd = None
        self._alphaCm = None
    
    def _create_splines(self):
        self._alphaCl = interp1d(self.alpha,self.cl,'cubic')
        self._alphaCd = interp1d(self.alpha,self.cd,'cubic')
        self._alphaCm = interp1d(self.alpha,self.cm,'cubic')
    
    def _calculate_data(self):
        self._calc_clmax()
        self._calc_cdmin()
        
    def _calc_clmax(self):
        if self._alphaCl==None:
            self._create_splines()
        f = lambda alpha: -self._alphaCl(alpha)
        alphaClmax = fminbound(f, self.alpha[0], self.alpha[-1])
        self.alphaClmax = alphaClmax
        self.clmax = self._alphaCl(alphaClmax)
    
    def _calc_cdmin(self):
        if self._alphaCd==None:
            self._create_splines()
        alphaCdmin = fminbound(self._alphaCd,self.alpha[0],self.alpha[-1])
        self.alphaCdmin = alphaCdmin
        self.cdmin = self._alphaCd(alphaCdmin)
        
    def get_clmax(self):
        if self.clmax==None:
            self._calc_clmax()
        return self.clmax
        
    def get_alpha_clmax(self):
        if self.alphaClmax==None:
            self._calc_clmax()
        return self.alphaClmax
    
    def get_cdmin(self):
        if self.cdmin==None:
            self._calc_cdmin()
        return self.cdmin
    
    def get_alpha_cdmin(self):
        if self.alphaCdmin==None:
            self._calc_cdmin()
        return self.alphaCdmin
    
    def get_cd_at_cl(self,cl):
        f = lambda alpha: self._alphaCl(alpha)-cl
        alpha = bisect(f,self.alpha[0],self.alpha[-1],xtol=1e-4)
        return self._alphaCd(alpha)
    
    
    def display(self):
        plt.figure(1)
        plt.grid(True)
        plt.plot(self.alpha, self.cl,'ko-')
        plt.xlabel('Angle of Attack,deg')
        plt.ylabel('Lift coefficient')
        plt.figure(2)
        plt.grid(True)
        plt.plot(self.cd, self.cl, 'bo-')
        plt.xlabel('Drag coefficient')
        plt.ylabel('Lift coefficient')
        plt.show()
    
    def __repr__(self):
        out = 'airfoil: %s\n'%self.name
        out += 'analysis: %s\n'%self.source
        out += 'Mach = %.4f\n'%self.Mach
        out += 'Re = %.4e\n'%self.Re
        out += 'alpha\tcl\tcd\tcm\n'
        for a,l,d,m in zip(self.alpha,self.cl,self.cd,self.cm):
            out += '%.2f\t%.4f\t%.8f\t%.4f\n'%(a,l,d,m)
        return out
    
    def save(self,path):
        fid = open(path,'wt')
        fid.write(self.__repr__())
        fid.close()
    

class AirfoilPolar:
    def __init__(self):
        self.source = None
        self.cl = list()
        self.cd = list()
        self.cm = list()
        self.alpha = list()
        self.Mach = None
        self.Re  = None
        self.cdp = list()
        self.TU = list()
        self.TL = list()
        self.SU = list()
        self.SL = list()
        self.LD = list()
        self.AC = list()
        self.CP = list()

    def _calc_clmax(self,simple=True):
        #FIXME: clmax trough optimization is not working
        """ calculate maximum lift coefficient for all mach numbers stored """
        if hasattr(self.Mach, '__iter__'):
            n = len(self.Mach)
            self.clmax = np.zeros(n)
            self.alphaClmax = np.zeros(n)
            for i in range(n):
                _cl = self.cl[i]
                self.clmax[i],self.alphaClmax[i] = self._get_clmax(_cl,self.alpha,simple)
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