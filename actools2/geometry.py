# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 11:10:25 2012
geometry 2D toolbox:
    point
    line
    curve
    split
    trim
    rotate
    translate
    scale
@author: Maxim
"""

import matplotlib.pyplot as plt
import numpy as ny
import scipy.interpolate as interpolate
import math
from math import radians, pi,sin,cos
from numpy import flipud, vstack, arange, array, zeros,linspace,hstack,transpose

class xyCurve:
    def __init__(self,pts):
        self.curve = interpolate.interp1d(pts[:,0],pts[:,1],'cubic')
        self.pts = pts
        self.xU = ny.min(pts[:,0])
        self.xL = ny.max(pts[:,0])
    def tanAngle(self,x,dx=1e-3):
        dx,dy = self.tanDirection(x,dx)
        return math.atan(dy/dx)
    def tanDirection(self,x,dx=1e-3):
        dy = self.curve(x+dx) - self.curve(x-dx)
        return [2*dx, dy]
    def tan(self,x,dx=1e-4):
        dx,dy = self.tanDirection(x,dx)
        return dy/dx
    def __call__(self,x):
        return self.curve(x)

def sort_airfoil_coordinates(coordinates,start=0.1,end=0.95):
    crdNew = ny.zeros(coordinates.shape)
    i = 0
    for pt in coordinates:
        if start<=pt[0]<=end:
            crdNew[i] = pt
            i += 1
    return crdNew[:i]


def get_sine_distribution(nPts):
    """
    Returns sine distribution between 0 and 1 with given number of points
    """
    xx = ny.linspace(0.0,1.0,nPts)
    cos_curve = lambda x: 1.0-ny.cos(x*ny.pi/2.0)
    return ny.array([cos_curve(_x) for _x in xx])


def get_cosine_distribution(nPts):
    """
    Returns cosine distribution between 0 and 1 with given number of points
    """
    xx = ny.linspace(0.0,1.0,nPts)
    cos_curve = lambda x: (ny.cos(x*ny.pi)+1.0)/2.0
    return ny.flipud(ny.array([cos_curve(_x) for _x in xx]))


def cross_section_cst(kUp,kLo,width,heightUp,heightLo,numPts=25):
    x = linspace(0,pi,numPts,True)
    xPts = [(cos(xx)+1)/2.0 for xx in x]

    upPts = cs_cst(kUp,xPts)
    loPts = cs_cst(kLo,xPts)
    
    hUp = cs_cst(kUp,[0.5])[0,1]
    hLo = cs_cst(kLo,[0.5])[0,1]
    
    upPts[:,0] = (upPts[:,0]-0.5)*width
    loPts[:,0] = (loPts[:,0]-0.5)*width
    upPts[:,1] = upPts[:,1]*heightUp/hUp
    loPts[:,1] = -loPts[:,1]*heightLo/hLo
    
    plt.figure(1)
    plt.plot(upPts[:,0],upPts[:,1],'*-')
    plt.hold(True)
    plt.plot(loPts[:,0],loPts[:,1],'*r-')
    plt.axis('equal')
    plt.grid(True)
    plt.show()
    
def cs_cst(k,xPts):
    y = zeros(len(xPts))
    i = 0
    for x in xPts:
        y[i] = (x*(1.0-x))**k
        i +=1
    coord = vstack([xPts,y])
    return transpose(coord)

def runTest2():
    cross_section_cst(0.6,0.6,2.5,1.3,1.4)

def ellipse2D(pt1,pt2,center,angleStart=0,angleEnd=360,nPts=10):
    """
    creates 2D ellipse with given edge points
    
    :ivar 2Dpoint pt1: side point
    :ivar 2Dpoint pt2: top point
    :ivar float angleStart: first angle of ellipse
    :ivar float angleEnd: end angle of ellipse
    """
    angleStart = radians(angleStart)
    angleEnd   = radians(angleEnd)
    a = abs(pt1[0] - pt2[0])
    b = abs(pt1[1] - pt2[1])
    step    = (angleEnd-angleStart)/(nPts-1)
    angles  = arange(angleStart,angleEnd+step,step)
    ellipse = zeros([nPts,2])
    for i,angle in enumerate(angles):
        ellipse[i,0] = center[0] + a*cos(angle)
        ellipse[i,1] = center[1] + b*sin(angle)
    return ellipse

def ellipse_circumference(a,b):
    a,b = float(a), float(b)
    C = 4.0*(pi*a*b + (a-b)**2)/(a+b)
    return C

def get_distance(point1,point2=[0.0,0.0]):
    """
    calculates the distance between two points
    
    :param array point1: coordinates of first point
    :param array point2: coordinates of second point
    
    Returns:
        distance between two points
    """
    distance = 0.0
    for i,crd in enumerate(point1):
        distance += (crd-point2[i])**2
    return distance**0.5

def join_coordinates(upPts, loPts):
    """
    joins two arrays of airfoil points into one.
    
    :param array upPts: upper curve points starting from leading edge
    :param array loPts: lower curve points starting from leading edge
    
    Returns:
        airfoil coordinates array starting from upper curve trailing 
        edge point
    """
    if upPts[0,0]<upPts[-1,0]:
        upPts = flipud(upPts)
    if loPts[0,0]>loPts[-1,0]:
        loPts = flipud(loPts)
    d = get_distance(upPts[-1,:],loPts[0,:])
    if d<1e-3:
        loPts = loPts[1:]
    coord = vstack([upPts,loPts])
    return coord

def separate_coordinates(coord):
    """
    separates airfoil coordinates array into array of upper and 
    lower points. Separation point is found as closest point to [0;0]
    
    :param array coord: airfoil coordinates array starting from upper 
        curve trailing edge point
    
    Returns:
        upper and lower curve points starting from leading edge
    """
    minDist = 10.0
    for i,point in enumerate(coord):
        distance = get_distance(point)
        if distance<minDist:
            sepPt = i
            minDist = distance
    upPts = vstack([coord[:sepPt],coord[sepPt]])
    upPts = flipud(upPts)
    loPts = coord[sepPt:]
    return upPts, loPts

def curve_pt_dist_normalized(curvePts):
    length = zeros(len(curvePts))
    dist = 0.0
    for i,pt in enumerate(curvePts[1:]):
        dist += get_distance(pt,curvePts[i])
        length[i+1] = dist
    length = length / dist
    return length

def curve_length(curvePts):
    """
    calculates the length of the curve as sum of distances of all curve points
    """
    length = 0.0
    for i,pt in enumerate(curvePts[1:]):
        pt2 = curvePts[i]
        length += get_distance(pt,pt2)
    return length

def get_pts_in_range(curvePts,lb=[],ub=[]):
    """
    takes curve point in format [x;y] and returns all the points 
    with lb <= x <= ub
    
    :param array curvePts: points of the curve
    :param float lb: lower bound x value
    :param float ub: upper bound x value
    """
    if lb==[]:
        lb = ny.min(curvePts[:,0])
    if ub==[]:
        ub = ny.max(curvePts[:,0])
    newCurvePts = ny.zeros([1,2])
    for pt in curvePts:
        if pt[0]>=lb and pt[0]<=ub:
            newCurvePts = ny.vstack([newCurvePts,pt])
    return newCurvePts[1:]

class LinePtPt():
    #limitation - vertical line is not allowed
    def __init__(self,startPt,endPt):
        self.startPt = startPt
        self.endPt   = endPt
        self.dx = endPt[0] - startPt[0]
        self.dy = endPt[1] - startPt[1]
        self.slope = self.dy / self.dx
        self.b = (endPt[0]*startPt[1] + startPt[0]*endPt[1])
        self.length = (self.dx**2 + self.dy**2)**0.5
    def __call__(self,x):
        y = self.slope*x + self.b
        return y

class SplineCurve():
    def __init__(self,nodes):
        self.type = 'xy'
        x = nodes.x
        y = nodes.y
        self.curve = interpolate.interp1d(x,y,'cubic')
    def __call__(self,x):
        if hasattr(x, '__iter__'):
            x = [float(xx) for xx in x]
        else:
            x = float(x)
        return self.curve(x)

def BPOcoef(order):
            K = ny.ones([order+1],int)
            for ii in range(1,order):
                K[ii] = math.factorial(order)/(math.factorial(ii)*math.factorial(order-ii))
            return K

class ClassFcn():
    def __init__(self,N=[0.5,1.0]):
        self.N1 = N[0]
        self.N2 = N[1]
    def __call__(self,x):
        return x**self.N1 * (1-x)**self.N2

class ShapeFcn():
    def __init__(self,A):
        self.order = len(A)-1
        K = BPOcoef(self.order)
        self.K = A*K
    def __call__(self,x):
        S = 0
        xx = 1-x
        for i,k in enumerate(self.K):
            S+= x**i * xx**(self.order-i) * k
        return S

class CstCurve():
    def __init__(self,A,N=[0.5,1.0]):
        self.type = 'xy'
        self.classCurve = ClassFcn(N)
        self.shapeCurve = ShapeFcn(A)
    def __call__(self,x):
        return self.classCurve(x) * self.shapeCurve(x)
    def get_coordinates(self,x):
        coord = ny.zeros([len(x),2])
        for i,xx in enumerate(x):
            coord[i,0] = xx
            coord[i,1] = self.__call__(xx)
        return coord
    def update_shape(self,Anew):
        self.shapeCurve = ShapeFcn(Anew)

class BezierCurve():
    def __init__(self,nodes):
        self.order = nodes.nPts - 1
        self.nPts = nodes.nPts
        self.dimm = nodes.dimm
        self.K = BPOcoef(self.order)
        self.pts = nodes.coord
        
    def _curve(self,t):
        pt = ny.zeros(self.dimm+1)
        pt[-1] = t
        for j in range(self.dimm):
            tmp = 0.0
            for k in range(self.nPts):
                 tmp += self.K[k]*self.pts[k,j]* (1-t)**(self.order-k) * t**k
            pt[j] = tmp
        return pt
    def __call__(self,t):
        if hasattr(t, '__iter__'):
            nPts = len(t)
            pts = ny.zeros([nPts,self.dimm+1,])
            for i,tt in enumerate(t):
                tt = float(tt)
                pts[i,:] = self._curve(tt)
            return pts
        else:
            t = float(t)
            return self._curve(t)

class PtsArray():
    def __init__(self,pts=[]):
        if pts!=[]:
            self.coord = ny.array(pts)
            self.x = self.coord[:,0]
            self.y = self.coord[:,1]
            self.nPts,self.dimm = self.coord.shape
    def xy(self,x,y):
        coord = ny.transpose(ny.vstack([x,y]))
        self.__init__(coord)

def SplitCurve():
    pass

def JoinCurve():
    pass

def runTest():
    pt1 = [0.0,0.0]
    pt2 = [3. ,2.0]
    pt3 = [5. ,0.0]
    pt4 = [7.,-2]
    
    xx = ny.array([0.,3,5,7])
    yy = ny.array([0.,2,0,-2])
    nodes = PtsArray()
    nodes.xy(xx,yy)
    nodes = PtsArray([pt1,pt2,pt3,pt4])
    curve1 = BezierCurve(nodes)
    curve2 = SplineCurve(nodes)
    
    t1 = ny.linspace(0,1,50)
    x1 = ny.linspace(0,7,50)
    pts = curve1(t1)
    y1 = curve2(x1)
    plt.plot(pts[:,0],pts[:,1],'--')
    plt.hold(True)
    plt.plot(x1,y1,'r-')
    plt.plot(nodes.x,nodes.y,'ro')
    plt.show()
    
def sine_dist_test():
    pts1 = get_cosine_distribution(25)
    pts2 = get_sine_distribution(25)
    plt.hold(True)
    plt.plot(pts1,ny.zeros(pts1.shape)+.2,'ko-')
    plt.plot(pts2,ny.zeros(pts1.shape)+.1,'ro-')
    plt.plot(pts1[0],0.14,'r*')
    plt.plot(pts2[0],0.24,'r*')
    plt.axis([0,1,0,.3])
    plt.show()

if __name__=="__main__":
    runTest()
