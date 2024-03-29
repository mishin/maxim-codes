# -*- coding: utf-8 -*-
"""
Created on Thu Jun 07 10:25:47 2012
Bezier curves
@author: Maxim
"""
import numpy as ny
import matplotlib.pyplot as plt
import sys
import scipy.interpolate as interp
import math

def factorial(n):
    n = ny.abs(int(n)) 
    if n < 1: n = 1 
    if n == 1: 
        return 1 
    else: 
        return n * factorial(n - 1)

def BPOcoef(order):
    K = ny.ones([order+1,1],int)
    for ii in range(1,order):
        K[ii] = factorial(order)/(factorial(ii)*factorial(order-ii))
    return K

def BezierCurve(Pts, opt):
    order, dimm = Pts.shape
    order -=1
    t = 0
    K = BPOcoef(order)
    if opt>1:
        Npts = int(opt)
        Curve = ny.zeros([Npts,dimm+1])
        Step = 1./(Npts-1)
        
        for ii in range(Npts):
            Curve[ii,-1] = t
            for jj in range(dimm):
                for kk in range(order+1):
                    Curve[ii,jj] += K[kk]*Pts[kk,jj]*(1-t)**(order-kk)*t**kk
            t +=Step
    elif opt>=0 and opt<=1:
        Curve = ny.zeros([1,dimm+1])
        t = float(opt)
        Curve[0,-1] = t
        for jj in range(dimm):
            for kk in range(order+1):
                Curve[0,jj] += K[kk]*Pts[kk,jj]*(1-t)**(order-kk)*t**kk
    return Curve

def pwBezier(Node, Vect, Npts):
    if Node.shape != Vect.shape: sys.exit('Node and Vector dimmensions should match')
    Nseg,dimm = Node.shape
    Nseg -=1
    #if hasattr(Npts, "__iter__") == False:
    #    tmp = ny.zeros([Nseg])
    #    tmp[:] = Npts
    #    Npts = tmp
    SegCurve = ny.zeros([Nseg,Npts,dimm+1])
    ControlPts = ny.zeros([4,dimm])
    for ii in range(Nseg):
        ControlPts[0,:] = Node[ii,:]
        ControlPts[1,:] = Node[ii,:]   + Vect[ii,:]
        ControlPts[2,:] = Node[ii+1,:] - Vect[ii+1,:]
        ControlPts[3,:] = Node[ii+1,:]
        SegCurve[ii,:] = BezierCurve(ControlPts,Npts)
        if ii == 0:
            Curve = SegCurve[0,:]
        else:
            tmp = SegCurve[ii,:]
            tmp[:,-1] = tmp[:,-1]+ii
            tmp = ny.delete(tmp,0,0)
            Curve = ny.vstack([Curve, tmp])
    return Curve, SegCurve

class xytCurveInterp:
    def __init__(self,curve):
        self.points = curve
        self.xt = interp.interp1d(curve[:,2],curve[:,0],'cubic')
        self.yt = interp.interp1d(curve[:,2],curve[:,1],'cubic')
    def getCoord(self,t):
        point = ny.array([self.xt(t),self.yt(t),t])
        return point

def normVect(Vector):
    length = 0
    for VectComp in Vector:
        length += VectComp**2
    length = length**0.5
    return Vector/length

class linePtDir2D:
    def __init__(self, StartPt, Direction, Length=1):
        self.Direction = normVect(Direction)
        self.StartPt = StartPt
        self.EndPt = StartPt + Length*self.Direction
        self.Length = Length
        sp = self.StartPt
        ep = self.EndPt
        #Ax+By+C=0
        self.A = sp[1]-ep[1]
        self.B = ep[0]-sp[0]
        self.C = (sp[0]*ep[1]-ep[0]*sp[1])

def lineIntersect(line1,line2):
    Xint = -(line1.C*line2.B - line2.C*line1.B)/(line1.A*line2.B - line2.A*line1.B)
    Yint = -(line1.A*line2.C - line2.A*line1.C)/(line1.A*line2.B - line2.A*line1.B)
    return ny.array([Xint,Yint])

def xyCurveSplit(curve,startPt,endPt=[]):
    NewCurve = ny.zeros([1,2])
    curveInt = interp.interp1d(curve[:,0],curve[:,1],'cubic')
    if startPt == []:
        for ii in range(curve.shape[0]):
            if curve[ii,0] < endPt: NewCurve = ny.vstack([NewCurve,curve[ii,:]])
        NewCurve = ny.vstack([NewCurve,[endPt,curveInt(endPt)]])
    elif endPt == []:
        NewCurve = ny.vstack([NewCurve,[startPt,curveInt(startPt)]])
        for ii in range(curve.shape[0]):
            if curve[ii,0] > startPt: NewCurve = ny.vstack([NewCurve,curve[ii,:]])
    else:
        NewCurve = ny.vstack([NewCurve,[startPt,curveInt(startPt)]])
        for ii in range(curve.shape[0]):
            if endPt > curve[ii,0] > startPt: NewCurve = ny.vstack([NewCurve,curve[ii,:]])
        NewCurve = ny.vstack([NewCurve,[endPt,curveInt(endPt)]])
    NewCurve = ny.delete(NewCurve,0,0)
    return NewCurve
        
def xytCurveSplit(curve,startPt,endPt=[]):
    NewCurve = ny.zeros([1,3])
    curveInt = xytCurveInterp(curve).getCoord
    if startPt == []:
        for ii in range(curve.shape[0]):
            if curve[ii,2] < endPt: NewCurve = ny.vstack([NewCurve,curve[ii,:]])
        NewCurve = ny.vstack([NewCurve,curveInt(endPt)])
    elif endPt == []:
        NewCurve = ny.vstack([NewCurve,curveInt(startPt)])
        for ii in range(curve.shape[0]):
            if curve[ii,2] > startPt: NewCurve = ny.vstack([NewCurve,curve[ii,:]])
    else:
        NewCurve = ny.vstack([NewCurve,curveInt(startPt)])
        for ii in range(curve.shape[0]):
            if endPt > curve[ii,2] > startPt: NewCurve = ny.vstack([NewCurve,curve[ii,:]])
        NewCurve = ny.vstack([NewCurve,curveInt(endPt)])
    NewCurve = ny.delete(NewCurve,0,0)
    return NewCurve

def rotate2D(curve, axis, angle):
    split = False
    if curve.shape[1]>2: 
        tmp = curve[:,2:]
        curve = curve[:,0:2]
        split = True
    angle = math.radians(angle)
    rotMatrix = ny.array([[math.cos(angle), -math.sin(angle)],[math.sin(angle), math.cos(angle)]])
    rotCurve = ny.dot((curve-axis),rotMatrix) + axis
    if split:
        return ny.hstack([rotCurve,tmp])
    else:
        return rotCurve

def dist2pts(point1, point2):
    dist = 0.
    dimm = ny.min([point1.shape[0],point2.shape[0]])
    for ii in range(dimm):
        dist += (point1[ii] - point2[ii])**2
    return math.sqrt(dist)

def join(curves, reverse):
    N = curves.shape[0]
    if reverse[0]: 
        newCurve = curves[0][-1,:]
    else:
        newCurve = curves[0][0,:]
    for ii in range(N):
        if reverse[ii]: 
            addCurve = ny.flipud(curves[ii])
        else:
            addCurve = curves[ii]
        addCurve = ny.delete(addCurve,0,0)
        newCurve = ny.vstack([newCurve,addCurve])
    return newCurve

def testCurve():
    node = ny.array([[0,0],[2.,2.],[4,1]])
    testVect = ny.array([[0,1.],[1.5,0],[0.5,-1]])
    
    testCurve, Curve1 = pwBezier(node,testVect,25)
    
    curve3 = curveSplit2(testCurve,0.2,1.7)
    
    print testCurve
    plt.plot(testCurve[:,0],testCurve[:,1],'-')
    plt.hold(True)
    plt.plot(curve3[:,0],curve3[:,1],'r-o')
    plt.grid(True)
    plt.show()
    
#testCurve()
    