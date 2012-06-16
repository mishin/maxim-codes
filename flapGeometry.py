# -*- coding: utf-8 -*-
"""
Created on Thu Jun 07 15:08:44 2012

@author: Maxim
"""
import numpy as ny
import math
import afLib
import scipy.optimize as root
import scipy.interpolate as interp
import curves
import matplotlib.pyplot as plt
from scipy.optimize import minimize

pt1X = 0.85
pt2Y = 0.3
pt2X = 0.7  #flap chord ratio
pt3X = 0.8
V1 = 0.1
V2 = 0.5
V3 = 0.08
theta = 90
pt4X = 0.65
BteThickness = 0.002
gap = 0.02
overlap = -0.02
deflection = 35

AfFile = 'GA37A315mod.dat'
dx = 1e-3

FNode = ny.zeros([3,2])
FVect = ny.zeros([3,2])
FTan  = ny.zeros([3,2])
BNode = ny.zeros([3,2])
BTan  = ny.zeros([2,2])

FNode[0,0] = pt1X
FNode[1,0] = pt2X
FNode[2,0] = pt3X
BNode[2,0] = pt4X

FTan[1,1] = 1.0
BTan[0,0] = math.cos(math.radians(theta))
BTan[0,1] = math.sin(math.radians(theta))

Up, Lo, _ = afLib.readAirfoil(AfFile)
AfUpCurve = interp.interp1d(Up[:,0],Up[:,1],'cubic')
AfLoCurve = interp.interp1d(Lo[:,0],Lo[:,1],'cubic')
#flap curve
pt2Ylb = AfLoCurve(pt2X)
pt2Yub = AfUpCurve(pt2X)
FNode[1,1] = (pt2Yub - pt2Ylb)*pt2Y + pt2Ylb
V2allowed = ny.min([FNode[1,1]-pt2Ylb, pt2Yub-FNode[1,1]])

FNode[0,1] = AfUpCurve(FNode[0,0])
dy = AfUpCurve(FNode[0,0]+dx) - FNode[0,1]
FTan[0,:] = -curves.normVect( ny.array([dx,dy]) ) * V1
FTan[1,:] = -FTan[1,:] * V2 * V2allowed
FNode[2,1] = AfLoCurve(FNode[2,0])
dy = AfLoCurve(FNode[2,0]+dx) - FNode[2,1]
FTan[2,:] = curves.normVect( ny.array([dx,dy]) ) * V3

FCurve,_ = curves.PwBezier(FNode,FTan,20)
#body curve
XFlap = interp.interp1d(FCurve[:,2],FCurve[:,0],'cubic')
YFlap = interp.interp1d(FCurve[:,2],FCurve[:,1],'cubic')

def Ctan5(t):
    dt = 1e-4
    dy = YFlap(t+dt) - YFlap(t)
    dx = XFlap(t+dt) - XFlap(t)
    CtanTheta = math.cos(math.radians(theta))/math.sin(math.radians(theta))
    return CtanTheta - dx/dy

pt5T = root.bisect(Ctan5,0.5,1.5)
XFlap = interp.interp1d(FCurve[:,2],FCurve[:,0],'cubic')
YFlap = interp.interp1d(FCurve[:,2],FCurve[:,1],'cubic')
BNode[0,0] = XFlap(pt5T)
BNode[0,1] = YFlap(pt5T)
BNode[2,1] = AfLoCurve(BNode[2,0])

dy = AfLoCurve(BNode[2,0]+dx) - BNode[2,1]
BTan[1,:] = curves.normVect( ny.array([dx,dy]) )

Line1 = curves.linePtDir2D(BNode[0,:],BTan[0,:])
Line2 = curves.linePtDir2D(BNode[2,:],BTan[1,:])
BNode[1,:] = curves.lineIntersect(Line1,Line2)

BCurve = curves.BezierCurve(BNode,11)

#trailing edge thickness of body part
Fcurve1 = interp.interp1d(FCurve[:,2],FCurve[:,0])
Fcurve2 = interp.interp1d(FCurve[:,2],FCurve[:,1])

def findTE(t):
    xx = Fcurve1(t)
    yFlap = Fcurve2(t)
    yAf = AfUpCurve(xx)
    return (yAf-yFlap)-BteThickness
    
t_teThickness = root.bisect(findTE,0,0.5)

FCurve2 = curves.curvesplit2(FCurve,t_teThickness,0)
BCurveLo = curves.curvesplit1(Lo,0,BNode[2,0])
BCurveUp = curves.curvesplit1(Up,0,Fcurve1(t_teThickness))
BCurveLo2 = BCurve[:,0:2]
BCurveUp2 = curves.curvesplit2(FCurve2,0,pt5T)[:,0:2]
BodyCurve = ny.vstack([ny.flipud(BCurveUp),BCurveLo,ny.flipud(BCurveLo2),ny.flipud(BCurveUp2)])

#flap curve
FCurveUp = curves.curvesplit1(Up,FNode[0,0],1.)
FCurveLo = curves.curvesplit1(Lo,FNode[2,0],1.)
FlapCurve = ny.vstack([ny.flipud(FCurveUp),FCurve[:,0:2],FCurveLo])

#flap position
tmpAxis = ny.array([.7,-0.5])
FlapCurve = curves.rotate2D(FlapCurve,tmpAxis,deflection)
refPt = BodyCurve[-1,:] 
currentOverlap = refPt[0] - ny.min(FlapCurve[:,0])
FlapCurve[:,0] = FlapCurve[:,0] + (currentOverlap-overlap)

def getGap(dy,curve,refPt,reqGap):
    curve[:,1] = curve[:,1] + dy
    curve = curve - refPt
    gaps = (curve[:,0]**2+curve[:,1]**2)**0.5
    return ny.min(gaps)-reqGap

dy = getGap(dy,FlapCurve,refPt,gap)
FlapCurve[:,1] = FlapCurve[:,1] + dy
print getGap(0,FlapCurve,refPt,gap)

plt.plot(FCurve[:,0],FCurve[:,1],'-')
plt.hold(True)
plt.plot(Up[:,0],Up[:,1],'g')
plt.plot(Lo[:,0],Lo[:,1],'g')
plt.plot(BCurve[:,0],BCurve[:,1],'r')
plt.axis([0,1,-0.4,0.4])

plt.figure(2)
plt.plot(BodyCurve[:,0],BodyCurve[:,1],'-')
plt.axis([0,1,-0.4,0.4])
plt.hold(True)
plt.plot(FlapCurve[:,0],FlapCurve[:,1],'-')
plt.show()

