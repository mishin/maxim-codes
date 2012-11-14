# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 20:42:05 2012

@author: Maxim
"""

import numpy as ny
import math
import afLib
import scipy.optimize as root
import curves
import matplotlib.pyplot as plt
from scipy.optimize import minimize, minimize_scalar

def getFlap(airfoilPath, node1,node2, vectLenRatio1, theta, gap, overlap, deflection,zTE):
    dx = 1e-4
    dt = 1e-4
    Npts = ny.array([25,10]) #flap, main section
    tmpAxis = ny.array([.8,-.5])
    theta = math.radians(theta)
    if overlap < -0.9*gap: gap = ny.abs(overlap)/0.9
    
    airfoil = afLib.readAirfoil(airfoilPath)
    
    node1[0,1] = airfoil.upCurve(node1[0,0])
    node1[1,1] = airfoil.getThicknessAtX(node1[1,0])*node1[1,1] + airfoil.loCurve(node1[1,0])
    node1[2,1] = airfoil.loCurve(node1[2,0])
    
    tan1 = ny.zeros([3,2]) #tangency of vectors forming flap curve
    tan1[0,0] = -dx
    tan1[0,1] = -(airfoil.upCurve(node1[0,0]+dx) - airfoil.upCurve(node1[0,0]))
    tan1[1,1] = -1.0
    tan1[2,0] = dx
    tan1[2,1] = airfoil.loCurve(node1[2,0]+dx) - airfoil.loCurve(node1[2,0])
    
    tan2 = ny.zeros([2,2]) #tangency of vectors forming main section lower curve
    tan2[0,0] = dx
    tan2[0,1] = airfoil.loCurve(node2[0,0]+dx) - airfoil.loCurve(node2[0,0])
    tan2[1,0] = math.cos(theta)
    tan2[1,1] = math.sin(theta)
    
    for ii in range(3): tan1[ii,:] = curves.normVect(tan1[ii,:])
    for ii in range(2): tan2[ii,:] = curves.normVect(tan2[ii,:])
    
    line1 = curves.linePtDir2D(node1[0,:],tan1[0,:])
    line2 = curves.linePtDir2D(node1[1,:],tan1[1,:])
    line3 = curves.linePtDir2D(node1[2,:],tan1[2,:])
    
    intersect1 = curves.lineIntersect(line1,line2)
    intersect2 = curves.lineIntersect(line3,line2)
    
    vectLenMax1 = ny.zeros(3)
    vectLenMax1[0] = curves.dist2pts(intersect1,node1[0,:])
    vectLenMax1[2] = curves.dist2pts(intersect2,node1[2,:])
    tmp = ny.array([airfoil.upCurve(node1[1,0])-node1[1,1] , node1[1,1] - airfoil.loCurve(node1[1,0])])
    vectLenMax1[1] = ny.min(tmp)
    del tmp
    
    vectLen1 = vectLenRatio1 * vectLenMax1
    vect1 = ny.zeros([3,2])
    for ii in range(3): vect1[ii,:] = tan1[ii,:] * vectLen1[ii]
    
    xytCurve1,_ = curves.pwBezier(node1,vect1,Npts[0])
    xytCurve2 = curves.rotate2D(xytCurve1,ny.array([0.7,0]),math.degrees(theta))
    splitPt1 = xytCurve2[ny.argmin(xytCurve2[:,0]),:]
    splitPt2 = xytCurve2[ny.argmax(xytCurve2[:,0]),:]
    xytCurveInt2 = curves.xytCurveInterp(xytCurve2)
    def get_dydt(t):
        if t == xytCurve2[-1,2]: t -=dt
        pt1 = xytCurveInt2.getCoord(t)
        pt2 = xytCurveInt2.getCoord(t+dt)
        return (pt2[1]-pt1[1])/dt
    t_theta = root.bisect(get_dydt,splitPt1[2],splitPt2[2])
    
    xytCurveInt1 = curves.xytCurveInterp(xytCurve1)
    node2[0,1] = airfoil.loCurve(node2[0,0])
    node2[2,:] = xytCurveInt1.getCoord(t_theta)[0:2]
    line4 = curves.linePtDir2D(node2[0,:],tan2[0,:])
    line5 = curves.linePtDir2D(node2[2,:],tan2[1,:])
    node2[1,:] = curves.lineIntersect(line4,line5)

    xytCurve3 = curves.BezierCurve(node2,Npts[1])

    def getTEthickness(t):
        pt1 = xytCurveInt1.getCoord(t)
        yy = airfoil.upCurve(pt1[0])
        thickness = yy-pt1[1]
        return thickness-zTE

    t_split = root.bisect(getTEthickness,0,1.)
    x_split = xytCurveInt1.getCoord(t_split)[0]
    
    split1 = curves.xyCurveSplit(airfoil.up,[],x_split)
    split2 = curves.xyCurveSplit(airfoil.up,node1[0,0])
    split3 = curves.xyCurveSplit(airfoil.lo,[],node2[0,0])
    split4 = curves.xyCurveSplit(airfoil.lo,node1[-1,0])
    split5 = curves.xytCurveSplit(xytCurve1,t_split,t_theta)
    
    mainSecCurve = ny.array([split1,split3,xytCurve3[:,0:2],split5[:,0:2]])
    mainSecCurve = curves.join(mainSecCurve,[True,False,False,True])
    
    flapCurve = ny.array([split2,xytCurve1[:,0:2],split4])
    flapCurve = curves.join(flapCurve,[True,False,False])
    flapCurve = curves.rotate2D(flapCurve,tmpAxis,deflection)
    
    tmpFlap = curves.rotate2D(xytCurve1,tmpAxis,deflection)
    xytCurveInt3 = curves.xytCurveInterp(tmpFlap)
    refPt = mainSecCurve[-1,:]
    
    def getFlapX(t):
        return xytCurveInt3.getCoord(t)[0]
    dist = minimize(getFlapX,1)
    dx = (refPt[0]-overlap) - dist.fun
    tmpFlap[:,0] += dx
    xytCurveInt2 = curves.xytCurveInterp(tmpFlap)

    def getGap(dy):
        tmpFlap2 = ny.copy(tmpFlap)
        tmpFlap2[:,1] = tmpFlap2[:,1] + dy
        xytCurveInt4 = curves.xytCurveInterp(tmpFlap2)
        def getDist(t):
            pt = xytCurveInt4.getCoord(t)
            return curves.dist2pts(pt,refPt)
        currentGap = minimize_scalar(getDist,bounds=(0.,2.),method='bounded')
        return (currentGap.fun - gap)
    
    dy = root.fsolve(getGap,0)
    flapCurve[:,0] = flapCurve[:,0] + dx
    flapCurve[:,1] = flapCurve[:,1] + dy
    
    airfoil.mainSec = mainSecCurve
    airfoil.flap = flapCurve
    airfoil.hasFlap = True
    
    return airfoil
    
    
    
node1 = ny.zeros([3,2])
node2 = ny.zeros([3,2])
vectLenRatio1 = ny.zeros(3)

node1[0,0] = 0.9
node1[1,0] = 0.7
node1[1,1] = 0.5 #height ratio
node1[2,0] = 0.9

node2[0,0] = 0.65

vectLenRatio1[0] = 0.8
vectLenRatio1[1] = 0.1
vectLenRatio1[2] = 0.9

theta = 90
overlap = -0.02
gap = 0.05

zTE = 0.002 #trailing edge thickness
deflection = 35
airfoilPath = 'GA37A315mod.dat'

#getFlap(airfoilPath, node1,node2, vectLenRatio1, theta, gap, overlap, deflection,zTE)
