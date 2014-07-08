# -*- coding: utf-8 -*-
"""
Created on Tue May 14 19:45:14 2013

@author: Maxim
"""
from numpy import array, zeros, flipud, ones, vstack, linspace, transpose, sin, cos, pi, arctan
from math import factorial
#from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
#from matplotlib.pyplot import figure, plot, axis, grid, show, title,xlim

#class CurveXyt:
#    """
#    parametric curve in format x = x(t), y = y(t)
#    """
#    def __init__(self,x,y,t):
#        self.x = x
#        self.y = y
#        self.t = t
#        self.xtCurve = interp1d(self.t,self.x,'cubic')
#        self.ytCurve = interp1d(self.t,self.y,'cubic')
#    def __call__(self,tnew):
#        xnew = self.xtCurve(tnew)
#        ynew = self.ytCurve(tnew)
#        return array([xnew,ynew])

class CstCurve():
    """
    Class used to create 2D airfoil curves using CST method. 
    The curve represents as:
        :math:`C(x)S(x)`

    where :math:`C(x)=x^{N_0}(1-x)^{N_1}` - class function

    and :math:`S(x)=\sum_{i=0}^nA_ix^(n-i)(1-x)^i` - shape function
    """
    class ClassFcn():
        def __init__(self,N=[0.5,1.0]):
            self.N1 = N[0]
            self.N2 = N[1]
        def __call__(self,x):
            return x**self.N1 * (1-x)**self.N2
    class ShapeFcn():
        def __init__(self,A):
            self.order = len(A)-1
            K = self.BPOcoef(self.order)
            self.K = A*K
        def __call__(self,x):
            S = 0
            xx = 1-x
            for i,k in enumerate(self.K):
                S+= x**i * xx**(self.order-i) * k
            return S
        def BPOcoef(self,order):
            K = ones([order+1],int)
            for ii in range(1,order):
                K[ii] = factorial(order)/(factorial(ii)*factorial(order-ii))
            return K
    def __init__(self,A,N=[0.5,1.0]):
        self.classCurve = self.ClassFcn(N)
        self.shapeCurve = self.ShapeFcn(A)
    def __call__(self,x):
        return self.classCurve(x) * self.shapeCurve(x)


class Airfoil2D():
    def __init__(self,name='sample airfoil'):
        self.upPts = list()
        self.loPts = list()
        self.coord = list()
        self.teGap = 0.0
        self.leRad = 0.0
        self.leCenter = zeros(2)
        self.camberLeSlope = 0.0
        self.name = str(name)
        self.upCurve = 0.0
        self.loCurve = 0.0
        self.coordCurve = 0.0
        self.nPtsUp = 0
        self.nPtsLo = 0

    def read_txt(self,filepath,aftype=1):
        if aftype==1:
            self._read_txt_1(filepath)
        else:
            self._read_txt_2(filepath)
        self.leCenter = zeros(2)
        self.leRad = 0.0
    
    def convert_af(self,af):
        #self.upPts = af.ptsUp
        #self.loPts = af.ptsLo
        self.coord = af.pts
        self.name = af.name
        self._separate_coord()
        self._analyze_geometry()

    def _read_txt_1(self,filepath):
        """
        reads airfoil coordinates from text file stored in xfoil and javafoil 
        like format.
        
        Text format as follows:  
            1. airfoil name
            2. coordinates starting from upper curve trailing edge (xx yy)
        
        Parameters
        ----------
        
        filePath : path
            airfoil coordinates file path
        
        Returns
        -------
        
        coord : 2d array
            airfoil coordinates
        name : string
            airfoil name
        """
        afFile = open(filepath,'rt')
        name = afFile.readline()[:-1]
        coord = array([0.0,0.0])
        for line in afFile:
            if line.strip()!='':
                segLine = line.split()
                point = [float(segLine[0]), float(segLine[1])]
                coord = vstack([coord,point])
        afFile.close()
        self.coord = coord[1:]
        self.name  = name
        self._separate_coord()
        self._analyze_geometry()

    def _read_txt_2(self,filepath):
        """
        reads airfoil coordinates from text file stored in format similar to 
        uiuc db file format.
        
        Text format as follows:
            1. airfoil name
            2. number of upper and lower curve points (n n)
            3. upper curve points starting from leading edge (xx yy)
            4. lower curve points starting from leading edge (xx yy)
        
        Parameters
        ----------
        
        filePath : path
            airfoil coordinates file path

        Returns
        -------
        
        coord : 2d array
            airfoil coordinates
        name : string
            airfoil name
        """
        afFile = open(filepath,'rt')
        name = afFile.readline()[:-1]
        nPts = [int(float(n)) for n in afFile.readline().split()]
        upPts = self.read_lines(afFile,nPts[0])
        loPts = self.read_lines(afFile,nPts[1])
        afFile.close()
        self.upPts = upPts
        self.loPts = loPts
        self.name = name
        self._join_coord()
        self._analyze_geometry()

    def write_txt(self,filepath,tab=True,writeName=True):
        """
        writes airfoil coordinates to text file.
        
        Parameters
        ----------
        
        filePath : path
            airfoil text file path to write coordinates
        tab : bool
            white space parameter. If True then x,y coordinates are separated by tabbing otherwise 
            separated by 4 spacings as required for X-foil input.
        writeName : bool
            if True then airfoil name will be written to the text file header otherwise only 
            coordinates will be stored
        """
        afFile = open(filepath,'wt')
        if writeName:
            afFile.write('%s\n'%self.name)
        if tab:
            whiteSpace = '\t'
        else:
            whiteSpace = '    '
        for point in self.coord:
            afFile.write('%.6f%s%.6f\n'%(point[0],whiteSpace,point[1]))
        afFile.close()

    def _analyze_geometry(self):
        self._calc_gap()
        self.nPtsUp = len(self.upPts)
        self.nPtsLo = len(self.loPts)
    
    def _calc_gap(self):
        self.teGap = self.upPts[-1,1] - self.loPts[-1,1]
        
    def set_trailing_edge_gap(self,gap=0.0):
        dz = (gap - self.teGap)/2.0
        self.upPts[:,1] += self.upPts[:,0]*dz
        self.loPts[:,1] += -self.loPts[:,0]*dz
        self._calc_gap()
        self._join_coord()

    def _separate_coord(self):
        coord = self.coord
        minDist = 10.0
        for i,point in enumerate(coord):
            distance = self._get_distance(point)
            if distance<minDist:
                sepPt = i
                minDist = distance
        upPts = vstack([coord[:sepPt],coord[sepPt]])
        self.upPts = flipud(upPts)
        self.loPts = coord[sepPt:]

    def _join_coord(self):
        upPts = self.upPts
        loPts = self.loPts
        if upPts[0,0]<upPts[-1,0]:
            upPts = flipud(upPts)
        if loPts[0,0]>loPts[-1,0]:
            loPts = flipud(loPts)
        d = self._get_distance(upPts[-1,:],loPts[0,:])
        if d<1e-3:
            loPts = loPts[1:]
        self.coord = vstack([upPts,loPts])

    def _get_distance(self,point1,point2=array([0.0,0.0])):
        distance = 0.0
        for i,crd in enumerate(point1):
            distance += (crd-point2[i])**2
        return distance**0.5

#    def _create_splines(self):
#        lenUp = self._curve_pt_dist_normalized(self.upPts)
#        lenLo = self._curve_pt_dist_normalized(self.loPts)
#        lenCr = self._curve_pt_dist_normalized(self.coord)
#        self.upCurve = CurveXyt(self.upPts[:,0],self.upPts[:,1],lenUp)
#        self.loCurve = CurveXyt(self.loPts[:,0],self.loPts[:,1],lenLo)
#        self.coordCurve = CurveXyt(self.coord[:,0],self.coord[:,1],lenCr)

    def _get_length(self):
        curvePts = self.coord
        length = 0.0
        for i,pt in enumerate(curvePts[1:]):
            pt2 = curvePts[i]
            length += self._get_distance(pt,pt2)
        return length

    def _curve_pt_dist_normalized(self,curvePts):
        length = zeros(len(curvePts))
        dist = 0.0
        for i,pt in enumerate(curvePts[1:]):
            dist += self._get_distance(pt,curvePts[i])
            length[i+1] = dist
        length = length / dist
        return length

    def redim(self,numPts=121):
        def cos_curve(x):
            phi = x*pi
            return (cos(phi)+1.0)/2.0
        if (numPts % 2)==1:
            numPts = int(numPts)/2+1
            tUp = linspace(0.0,1.0,numPts)
            tLo = tUp
        else:
            numPts = int(numPts)/2+1
            tUp = linspace(0.0,1.0,numPts)
            tLo = linspace(0.0,1.0,numPts-1)
        tUp = array([cos_curve(t) for t in tUp])
        tLo = array([cos_curve(t) for t in tLo])
        self._create_splines()
        xUp, yUp = self.upCurve(tUp)
        xLo, yLo = self.loCurve(tLo)
        self.upPts = transpose(vstack([xUp,yUp]))
        self.loPts = transpose(vstack([xLo,yLo]))
        self._join_coord()

    def display(self,linetype='ko-'):
        plt.figure()
        plt.grid(True)
        plt.plot(self.coord[:,0],self.coord[:,1],linetype)
        plt.axis('equal')
        plt.xlim((0.0,1.0))
        plt.title(self.name)
        plt.show()
    
    def _get_cos_distribution(self,npts=30):
        x = linspace(0,pi/2.0,npts)
        return -cos(x) + 1.0

    def naca4(self,thickness=12.0,camber=0.0,camberLoc=0.0):
        t = thickness / 100.0
        m = camber / 100.0
        p = camberLoc / 100.0
        x = self._get_cos_distribution()
        y = t/0.2*(0.2969*x**0.5 - 0.1281*x - 0.3516*x*x + 0.2843*x**3.0 - 0.1015*x**4.0)
        y[-1] = 0.0 #FIXME: y[-1] always led to 1e-17, this is fix
        self.leRad = 1.1019*t*t
        self.name = 'NACA%.0f%.0f%.0f'%(camber,camberLoc/10.0,thickness)
        if not m*p==0.0:
            yc = zeros(len(x))
            tanTheta = zeros(len(x))
            for i,xpt in enumerate(x):
                if xpt<=p:
                    yc[i] = m*xpt/(p*p)*(2.0*p - xpt)
                    tanTheta[i] = 2.0*m/(p*p)*(p-xpt)
                else:
                    yc[i] = m*(1.0-xpt)/((1.0-p)**2.0)*(1.0+xpt-2.0*p)
                    tanTheta[i] = -2.0*m/((1.0-p)**2.0)*(xpt-p)
            theta = arctan(tanTheta)
            xu = x - y*sin(theta)
            xl = x + y*sin(theta)
            yu = yc + y*cos(theta)
            yl = yc - y*cos(theta)
            self.upPts = transpose(vstack([xu,yu]))
            self.loPts = transpose(vstack([xl,yl]))
            self.leSlope = theta[0]
            
        else:
            self.upPts = transpose(vstack([x,y]))
            self.loPts = transpose(vstack([x,-y]))
            self.leSlope = 0.0
        self.leCenter = array([cos(self.leSlope),sin(self.leSlope)]) * self.leRad + self.upPts[0]
        self._join_coord()
        self._calc_gap()
    
    def create_CST(self,Au,Al):
        self.upCurve = CstCurve(Au)
        self.loCurve = CstCurve(Al)
        x = linspace(0,pi/2.0,25)
        x = array([-cos(xx) + 1.0 for xx in x])
        self.upPts = self.upCurve.get_coordinates(x)
        self.loPts = self.loCurve.get_coordinates(x)
        self.join_coordinates()
        self.analyze_geometry()
        self.name = 'CST airfoil'

class Airfoil3D():
    def __init__(self,af2d,plane='xz',chord=1.0e3):
        self.af = af2d
        self.set_plane(plane)
        self.upPts    = None
        self.loPts    = None
        self.leCenter = 0.0
        self.leRad    = 0.0
        self.teGap = af2d.teGap
        self.upPts3D = self._process_2d_pts(self.af.upPts)
        self.loPts3D = self._process_2d_pts(self.af.loPts)
        self.leCenterNondim = self._process_2d_pt(self.af.leCenter)
        self.leRadNondim = af2d.leRad
        self.set_chord(chord)

    def set_plane(self,plane):
        self.plane = plane
        tmp = self._assign_plane(plane[0]) + self._assign_plane(plane[1])
        self.rotAxisDir = zeros(3)
        self.rotAxisDir[3-tmp] = 1.0
        self.planeOrientation = array([self._assign_plane(plane[0]),self._assign_plane(plane[1])])

    def _get_plane_idx(self,plane):
        x = self._assign_plane(plane[0])
        y = self._assign_plane(plane[1])
        return x,y
        
    def _assign_plane(self,axis):
        if axis=='x':
            return 0
        elif axis=='y':
            return 1
        else:
            return 2

    def _process_2d_pts(self,points):
        pts3d = zeros([len(points),3])
        pts3d[:,self.planeOrientation[0]] = points[:,0]
        pts3d[:,self.planeOrientation[1]] = points[:,1]
        return pts3d
    
    def _process_2d_pt(self,point):
        point3d = zeros(3)
        point3d[self.planeOrientation[0]] = point[0]
        point3d[self.planeOrientation[1]] = point[1]
        return point3d

    def _strip_3d_pts(self,points3d):
        xAxis,yAxis = self._get_plane_idx(self.plane)
        points = zeros([len(points3d),2])
        points[:,0] = points3d[:,xAxis]
        points[:,1] = points3d[:,yAxis]
        return points

    def _process_3d_pts(self,points3d,planeNew):
        points = self._strip_3d_pts(points3d)
        return self._process_2d_pts(points,planeNew)

    def set_chord(self,chord):
        self.chord = chord
        self.upPts = self.upPts3D * chord
        self.loPts = self.loPts3D * chord
        self.leCenter = self.leCenterNondim * chord
        self.leRad = self.leRadNondim * chord

    def set_ref_plane(self,plane):
        self.upPts3D = self._process_3d_pts(self.upPts3D,plane)
        self.loPts3D = self._process_3d_pts(self.loPts3D,plane)
        self.plane = plane

    def set_trailing_edge_gap_mm(self,gap):
        pass


def run_test1():
    af = Airfoil2D()
    #af.naca4(12,3,40)
    af.read_txt('exampleAF.txt')
    af3d = Airfoil3D(af,'xz',1e3)
    print af3d.leCenter
    print af3d.leRad
    print af.teGap
    print af3d.teGap
    af.display()

if __name__=="__main__":
    run_test1()