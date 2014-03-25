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
        self.af0 = read_txt('pb092_geometry.txt')
        self.af = None
    
    def set_airfoil(self,x):
        self.af = get_new_airfoil(x,*(self.af0,))
        
    def f(self,x):
        self.set_airfoil(x)
    
    def g1high(self,x):
        self.set_airfoil(x)
    def g1low(self,x):
        self.set_airfoil(x)
    def g2high(self,x):
        self.set_airfoil(x)
    def g2low(self,x):
        self.set_airfoil(x)
    
    def jac(self,x,dx=1e-4):
        pass
    def df(self,x,dx=1e-4):
        pass
    
    def run_cfd(self,x):
        pass