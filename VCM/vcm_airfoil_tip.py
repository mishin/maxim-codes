# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 22:15:13 2014

@author: Maxim
"""
from airfoil import *
import geometry as geom
from scipy.optimize import minimize
from flight_conditions import FlightConditions
import sys
from misc_tools import read_tabulated_data_without_header, Normalization
import os
from CFDpolar import CFDsolver,simple_cfd_analysis


def get_tangency_vector(curve,x,dx=1e-5):
    dy = curve(x+dx)-curve(x)
    l = (dy*dy  + dx*dx)**0.5
    return np.array([dx/l, dy/l])

def split_points(pts,xmax=None,xmin=None):
    ptsNew = np.zeros([1,pts.shape[1]])
    if xmax==None:
        xmax = pts[-1,0]
    if xmin==None:
        xmin = pts[0,0]
    for pt in pts:
        if xmin<=pt[0]<=xmax:
            ptsNew = np.vstack([ptsNew,pt])
    return ptsNew[1:]

def get_new_airfoil(x,*args):
    dx = 1e-5
    af0 = args[0]
    afNew = Airfoil()
    x1, x2, y2, x3, l1, l2, l3 = x[0],x[1],x[2],x[3],x[4],x[5],x[6]

    nodesUp = np.zeros([4,2])
    nodesLo = np.zeros([4,2])
    af0._create_splines()
    y1 = af0._curveUp2(x1)
    y3 = af0._curveLo2(x3)
    nodesUp[0] = np.array([x1,y1])
    nodesUp[2] = np.array([x2,y2+l2])
    nodesUp[3] = np.array([x2,y2])
    nodesLo[0] = np.array([x2,y2])
    nodesLo[1] = np.array([x2,y2-l2])
    nodesLo[3] = np.array([x3,y3])
    tan1 = get_tangency_vector(af0._curveUp2,x1,dx)
    tan2 = get_tangency_vector(af0._curveLo2,x3,dx)
    nodesUp[1] = nodesUp[0] - l1*tan1
    nodesLo[2] = nodesLo[3] - l3*tan2
    
    _nodesUp = geom.PtsArray()
    _nodesLo = geom.PtsArray()
    _nodesUp.xy(nodesUp[:,0],nodesUp[:,1])
    _nodesLo.xy(nodesLo[:,0],nodesLo[:,1])
    curve1 = geom.BezierCurve(_nodesUp)
    curve2 = geom.BezierCurve(_nodesLo)
    t1 = np.linspace(0,1,30)
    curve1pts = curve1(t1)
    curve2pts = curve2(t1)
    afUp = split_points(af0.ptsUp,xmin=x1)
    afLo = split_points(af0.ptsLo,xmin=x3)
    newPts = np.vstack([np.flipud(afUp),curve1pts[:,0:2],curve2pts[:,0:2],afLo])
    afNew.pts = newPts
    idx = np.argmin(newPts[:,0])
    afNew.ptsUp = np.flipud(newPts[0:idx])
    afNew.ptsLo = newPts[idx:]
    #nodesAll = np.vstack([nodesUp,nodesLo])
    return afNew


class AirfoilAnalysis:
    def __init__(self):
        self.gclmax = 1.8
        self.galphaStall = 18.0
        self.Mach = 0.0
        self.Re = 2.4e6
        self.cruise = FlightConditions(self.Mach,0)
        self.cruise.Re = self.Re
        self.alphaSeq = [-10,20,0.5]
        self.clCruise = [0.2,0.3]
        self.af0 = read_txt('pb092_geometry.txt')
        self.af = None
        self.dx = 1e-4

    def set_airfoil(self,x):
        self.af = get_new_airfoil(x,*(self.af0,))

    def get_J_polar(self):
        return self.af.get_J_polar(self.Mach,self.Re,self.alphaSeq)

    def get_cfd_polar(self):
        pol = simple_cfd_analysis(self.af,self.cruise)
        pol._calc_clmax()
        pol._calc_cdmin()
        return pol

    def f(self,x):
        self.set_airfoil(x)
        pol = self.get_J_polar()
        cd = array([pol.get_cd_at_cl(cl) for cl in self.clCruise])
        return cd.mean()
    
    def g1high(self,x):
        self.set_airfoil(x)
        pol = self.get_cfd_polar()
        return pol.clmax - self.gclmax
        
    def g1low(self,x):
        self.set_airfoil(x)
        pol = self.get_J_polar()
        return pol.clmax - self.gclmax

    def g2high(self,x):
        self.set_airfoil(x)
        pol = self.get_cfd_polar()
        return pol.alphaClmax - self.galphaStall
        
    def g2low(self,x):
        self.set_airfoil(x)
        pol = self.get_J_polar()
        return pol.alphaClmax - self.galphaStall
    
    def jac(self,x,fun,dx=None):
        if dx==None:
            dx = self.dx
        grad = np.zeros(len(x))
        fval = self.f(x)
        for i,xx in enumerate(x):
            X = np.copy(x)
            X[i] = X[i]+dx
            grad[i] = (func(X)-fval)/dx
        return grad

    def df(self,x):
        return self.jac(x,self.f)


def vfm_tip_airfoil_design():
    aa = AirfoilAnalysis()