# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 14:22:32 2012

@author: Maxim
"""

import numpy as ny
import math
import afLib
import scipy.optimize as root
import curves
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def getFlap(airfoilPath, nodes1,nodes2, vectLen1, theta, gap, overlap, deflection,teThickness):
    dx = 1e-4
    
    flapNpts = 15
    mainSecTEpts = 8
    tmpAxis = ny.array([.8,-.5])
    if overlap < -0.9*gap: gap = ny.abs(overlap)/0.9
    af = afLib.readAirfoil(airfoilPath)
    xyCurve1 = af.upCurve
    xyCurve2 = af.loCurve
    
    nodes1[0,1] = xyCurve1(nodes1[0,0])
    nodes1[2,1] = xyCurve2(nodes1[2,0])
    nodes2[0,1] = xyCurve2(nodes2[0,0])
    
    nodes1[1,1] = nodes1[1,1]*af.getThicknessAtX(nodes1[1,0])+xyCurve2(nodes1[1,0])
    tmpVecLen   = ny.array([xyCurve1(nodes1[1,0])-nodes1[1,1], nodes1[1,1]-xyCurve2(nodes1[1,0])])
    vectLen1[1] = vectLen1[1]*ny.min(tmpVecLen)
    del tmpVecLen
    
    tangency1 = ny.zeros([3,2])
    tangency2 = ny.zeros([2,2])
    tangency1[0,0] = dx
    tangency1[0,1] = xyCurve1(nodes1[0,0]+dx) - xyCurve1(nodes1[0,0])
    tangency1[1,1] = -1.
    tangency1[2,0] = dx
    tangency1[2,1] = xyCurve2(nodes1[2,0]+dx) - xyCurve2(nodes1[2,0])
    tangency1[0,:] = - tangency1[0,:]
    for ii in range(3): tangency1[ii,:] = curves.normVect(tangency1[ii,:])
    vect1 = ny.zeros([3,2])
    for ii in range(3): vect1[ii,:] = tangency1[ii,:]*vectLen1[ii]
    
    tangency2[0,0] = dx
    tangency2[0,1] = xyCurve2(nodes2[0,0]+dx) - xyCurve2(nodes2[0,0])
    tangency2[1,0] = math.cos(math.radians(theta))
    tangency2[1,1] = math.sin(math.radians(theta))
    
    xytCurve1,_ = curves.pwBezier(nodes1,vect1,flapNpts) #flap curve
    xytCurveInt1 = curves.xytCurveInterp(xytCurve1)

    def getTangency(t):
        dt = 1e-4
        pt1 = xytCurveInt1.getCoord(t)
        pt2 = xytCurveInt1.getCoord(t+dt)
        ctan = (pt1[0]-pt2[0])/(pt1[1]-pt2[1])
        ctanTheta = math.cos(math.radians(theta))/math.sin(math.radians(theta))
        return ctanTheta - ctan

    t_theta = root.bisect(getTangency,0.5,1.5)
    nodes2[2,:] = xytCurveInt1.getCoord(t_theta)[0:2]
    line1 = curves.linePtDir2D(nodes2[0,:],tangency2[0,:])
    line2 = curves.linePtDir2D(nodes2[2,:],tangency2[1,:])
    nodes2[1,:] = curves.lineIntersect(line1,line2)
    
    xytCurve2 = curves.BezierCurve(nodes2,mainSecTEpts)
    
    def getTEthickness(t):
        pt1 = xytCurveInt1.getCoord(t)
        yy = xyCurve1(pt1[0])
        thickness = yy-pt1[1]
        return thickness-teThickness
        
    t_split = root.bisect(getTEthickness,0,1.0)
    x_split = xytCurveInt1.getCoord(t_split)[0]
    
    split1 = curves.xyCurveSplit(af.up,[],x_split)
    split2 = curves.xyCurveSplit(af.up,nodes1[0,0])
    split3 = curves.xyCurveSplit(af.lo,[],nodes2[0,0])
    split4 = curves.xyCurveSplit(af.lo,nodes1[-1,0])
    split5 = curves.xytCurveSplit(xytCurve1,t_split,t_theta)
    
    mainSecCurve = ny.array([split1,split3,xytCurve2[:,0:2],split5[:,0:2]])
    mainSecCurve = curves.join(mainSecCurve,[True,False,False,True])
    
    flapCurve = ny.array([split2,xytCurve1[:,0:2],split4])
    flapCurve = curves.join(flapCurve,[True,False,False])
    flapCurve = curves.rotate2D(flapCurve,tmpAxis,deflection)
    
    tmpFlap = curves.rotate2D(xytCurve1,tmpAxis,deflection)
    xytCurveInt2 = curves.xytCurveInterp(tmpFlap)
    refPt = mainSecCurve[-1,:]
    
    def getFlapX(t):
        return xytCurveInt2.getCoord(t)[0]
    dist = minimize(getFlapX,1)
    dx = (refPt[0]-overlap) - dist.fun
    
    tmpFlap[:,0] += dx
    xytCurveInt2 = curves.xytCurveInterp(tmpFlap)
    
    def getGap(dy):
        tmpCurve = ny.copy(tmpFlap)
        tmpCurve[:,1] += dy
        xytCurveInt3 = curves.xytCurveInterp(tmpCurve)
        def getDist(t):
            pt = xytCurveInt3.getCoord(t)
            return curves.dist2pts(pt,refPt)
        CurGap = minimize(getDist,1.).fun
        return CurGap-gap

    dy = root.fsolve(getGap,0)
    
    flapCurve[:,0] = flapCurve[:,0] + dx
    flapCurve[:,1] = flapCurve[:,1] + dy
    
    af.mainSec = mainSecCurve
    af.flap = flapCurve
    af.hasFlap = True
    
#    plt.plot(mainSecCurve[:,0],mainSecCurve[:,1])
#    plt.hold(True)
#    plt.axis([0,1.2,-0.5,0.5])
#    plt.grid(True)
#    plt.plot(flapCurve[:,0],flapCurve[:,1])
#    plt.show()
    
    return af

airfoilPath = 'GA37A315mod.dat'
nodes1 = ny.zeros([3,2])
nodes2 = ny.zeros([3,2])
vectLen1 = ny.zeros([3])
nodes1[0,0] = 0.85
nodes1[1,0] = 0.7 #flap chord ratio
nodes1[1,1] = 0.3
nodes1[2,0] = 0.8
nodes2[0,0] = 0.65
vectLen1[0] = 0.1
vectLen1[1] = 0.5
vectLen1[2] = 0.08
theta = 120
teThickness = 0.002
deflection = 35
overlap = -0.03
gap = 0.05

#getFlap(airfoilPath, nodes1,nodes2, vectLen1, theta, gap, overlap, deflection,teThickness)
