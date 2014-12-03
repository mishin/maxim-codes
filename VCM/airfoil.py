# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 18:14:28 2013
@author: Maxim

Module contains methods for storing and processing arifoils and aerodynamic data. 

"""

from numpy import array, zeros, ones, vstack, transpose, argsort, append, arange, linspace, radians, hstack,sin,cos
from scipy.interpolate import interp1d, RectBivariateSpline
from subprocess import Popen, PIPE
import os
from scipy.optimize import fminbound,bisect
from matplotlib.pyplot import figure, plot, axis, hold, grid, legend, title, show,xlim,ylim
from math import factorial, pi, radians
import shlex
import matplotlib.pyplot as plt
import paths
import FlightConditions as fc
import dbTools
import geometry as geom
import win32com.client

def load(airfoilName,dbPath=''):
    """
    loads airfoil from .xls database
    
    Parameters
    ----------
    
    airfoilName : string
        name of the airfoil to read
    dbPath : path
        file path of the airfoil database. If path is not specified then 
        default db path will be used
    
    Examples
    --------
    
    load airfoil and store in *af* variable
    
    >>> import airfoil
    >>> af = airfoil.load('GA37A315')
    >>> af.plot()
    """
    af = Airfoil()
    af.read_xls(airfoilName,dbPath)
    return af

def cst(Au,Al,nPts=30,dist='sin'):
    af = Airfoil()
    af.create_CST(Au,Al,nPts,dist)
    af.radius_LE = Au[0]**2*50.
    af.camber_slope_LE = 0.0
    return af


def cst_x(A,nPts=25):
    af = Airfoil()
    n = len(A)
    if n%2==0:
        Au = A[:n/2]
        Al = A[n/2:]
    else:
        Au = array(A[:int(n/2)+1])
        Al = hstack([-A[0],A[int(n/2)+1:]])
    af.create_CST(Au,Al,nPts)
    return af

class AirfoilPolar:
    """
    Stores aerodynamic data of airfoil.
    
    Arguments
    ---------
    source : string
        source of aerodynamic data (xfoil, javafoil, CFD, etc.)
    """
    def __init__(self):
        self.source = ''
        self.analysisFlag = False
        self.Mach  = list()
        self.Re    = list()
        self.alpha = list()
        self.cl    = list()
        self.cd    = list()
        self.cd_p  = list()
        self.cm    = list()
        self.TU    = list()
        self.TL    = list()
        self.SU    = list()
        self.SL    = list()
        self.LD    = list()
        self.AC    = list()
        self.CP    = list()
        self.clAlpha = None
    def calc_clmax(self):
        """
        calculates maximum lift coefficient for all Mach numbers stored.
        """
        if not hasattr(self.Mach,'__iter__'):
            clCurve = interp1d(self.alpha,-self.cl,'cubic')
            ub,lb = self.alpha[0],self.alpha[-1]
            self.alphaClmax = fminbound(clCurve,ub,lb)
            self.clmax = -clCurve(self.alphaClmax)
        else:
            n = len(self.Mach)
            self.clmax = zeros(n)
            self.alphaClmax = zeros(n)
            for i in range(n):
                clCurve = interp1d(self.alpha,-self.cl[i])
                ub,lb = self.alpha[0],self.alpha[-1]
                self.alphaClmax[i] = fminbound(clCurve,ub,lb)
                self.clmax[i] = -clCurve(self.alphaClmax[i])
            
        
    def create_splines(self,mach=False):
        """
        creates 1D or 2D spline of lift, drag and moment coefficient using Interp2D
        
        Parameters
        ----------
        
        mach : bool
            if mach=True then generates 2D spline of aero coefficients vs. 
            angle of attack and Mach number, otherwise generates 1D spline 
            for angle of attack at Mach number stored in self.Mach[0]
        """
        if mach:
#            self.clAlpha = Interp2D(self.Mach,self.alpha,self.cl)
#            self.cdAlpha = Interp2D(self.Mach,self.alpha,self.cd)
#            self.cmAlpha = Interp2D(self.Mach,self.alpha,self.cm)
            self.clAlpha = FlatPlateInterp(self.Mach,self.alpha,self.cl,'lift')
            self.cdAlpha = FlatPlateInterp(self.Mach,self.alpha,self.cd,'drag')
        else:
            self.clAlpha = interp1d(self.alpha,self.cl,'cubic')
            self.cdAlpha = interp1d(self.alpha,self.cd,'cubic')
            self.cmAlpha = interp1d(self.alpha,self.cm,'cubic')
    def get_cd_at_cl(self,cl=0.4):
        if self.clAlpha==None:
            self.create_splines()
        def obj(x):
            return self.clAlpha(x)-cl
        alpha = bisect(obj,a=self.alpha[0],b=self.alpha[-1],xtol=1e-3,maxiter=30)
        return self.cdAlpha(alpha)
    def display(self):
        """
        Plots lift and drag coefficient vs. angle of attack for all available 
        Mach numbers
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid(True)
        ax.hold(True)
        ax.plot(self.alpha,self.cl)
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.grid(True)
        ax2.plot(self.alpha,self.cd)
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        ax3.grid(True)
        ax3.plot(self.cd,self.cl,'bo-')
        plt.show()

class Interp2D:
    """
    2D interpolation using rectangular bivariate spline
    scipy.interpolate.RectBivariateSpline
    
    Returns interpolated value if requested value is between bounds otherwise 
    returns value at the boundaryrec
    """
    def __init__(self,x,y,z):
        self.model = RectBivariateSpline(x,y,z)
        self.x = x
        self.xmax = max(x)
        self.xmin = min(x)
        self.y = y
        self.ymax = max(y)
        self.ymin = min(y)
    def __call__(self,x,y):
        if x>self.xmax:
            x = self.xmax
        if x<self.xmin:
            x = self.xmin
        if y>self.ymax:
            y = self.ymax
        if y<self.ymin:
            y = self.ymin
        return self.model(x,y)

class FlatPlateInterp:
    def __init__(self,Mach,alpha,coef,data='lift'):
        self.alphaLim = array([-90.0,90.0])
        self.n = 20
        self.alpha = alpha
        self.Mach = Mach
        self.coef = coef
        self.datatype = data
        self._process()
    def _process(self):
        alphaL = linspace(self.alphaLim[0],self.alpha[0]-5.0,self.n)
        alphaU = linspace(self.alpha[-1]+5.0,self.alphaLim[1],self.n)
        if self.datatype=='lift':
            coefL = self._lift(alphaL)
            coefU = self._lift(alphaU)
        elif self.datatype=='drag':
            coefL = self._drag(alphaL)
            coefU = self._drag(alphaU)
        coef = hstack([coefL,self.coef[0],coefU])
        for i in range(len(self.Mach)-1):
            coef_ = hstack([coefL,self.coef[i+1],coefU])
            coef = vstack([coef,coef_])
        self.alpha = hstack([alphaL,self.alpha,alphaU])
        self.model = RectBivariateSpline(self.Mach,self.alpha,coef)
    def __call__(self,Mach,alpha):
        return self.model(Mach,alpha)[0]
    def _lift(self,alpha):
        return sin(2.0*radians(alpha))
    def _drag(self,alpha):
        return -cos(2.0*radians(alpha))+1.0


class CurveXyt:
    """
    parametric curve in format x = x(t), y = y(t)
    """
    def __init__(self,x,y,t):
        self.x = x
        self.y = y
        self.t = t
        self.xtCurve = interp1d(self.t,self.x,'cubic')
        self.ytCurve = interp1d(self.t,self.y,'cubic')
    def __call__(self,tnew):
        xnew = self.xtCurve(tnew)
        ynew = self.ytCurve(tnew)
        return array([xnew,ynew])

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

class Xfoil:
    """
    class to communicate with Xfoil
    """
    def __init__(self,graphic=False):
        self.path = paths.MyPaths()
        args=shlex.split(self.path.Xfoil,False,os.name=='posix')
        self.ps=Popen(args,stdin=PIPE,stderr=PIPE,stdout=PIPE)
        if graphic==False:
            self.cmd('PLOP\nG\n')
    def cmd(self,command,echo=False):
        self.ps.stdin.write(command+'\n')
        if echo: print command
    def terminate(self):
        self.cmd('\n\n\nQUIT')
        self.ps.stderr.close()
        self.ps.stdout.read()

class AirfoilList:
    #TODO: store only single airfoil object or multiple for design purposes (option)
    #TODO: check if same airfoil exists in list before storing (may cause excel error)
    def __init__(self):
        self.airfoils = list()
        self.uniqueNames = list()
        self.contents = list()
        self.numOfAirfoils = 0
        self.numOfUnique = 0
    def __getitem__(self,k):
        idx = self.contents[k]
        return self.airfoils[idx]
    def add_airfoil_xls(self,name,overwrite=False):
        self.numOfAirfoils += 1
        if overwrite:
            exist = self.check_existing(name)
            if not exist:
                newAirfoil = Airfoil()
                newAirfoil.read_xls(name)
                self.airfoils.append(newAirfoil)
        else:
            newAirfoil = Airfoil()
            newAirfoil.read_xls(name)
            self.airfoils.append(newAirfoil)
            self.contents.append(self.numOfAirfoils-1)
    def save_airfoils(self):
        for af in self.airfoils:
            af.write_xls()
    def check_existing(self,name):
        try:
            idx = self.uniqueNames.index(name)
            self.contents.append(idx)
            return True
        except ValueError:
            self.uniqueNames.append(name)
            self.numOfUnique += 1
            self.contents.append(self.numOfUnique-1)
            return False

class Airfoil:
    """
    main class used for work with airfoil sections. Contains methods to 
    read/write text and xls files, analyze using xfoil and javafoil and 
    some tools to manipulate geometry.
    
    Attributes
    ----------
    
    coord : float array
        coordinates of airfoil starting from upper curve trailing edge
    upPts : float array
        coordinates of upper curve starting from leading edge
    loPts : float array
        coordinates of lower curve starting from leading edge
    polar : Airfoil polar
        aerodynamic data of the airfoil. This data becomes available if 
        it was stored in airfoil db or calculated using *build_aero_table*, 
        *get_X_polar* or *get_J_polar*.

    Examples
    --------
    
    read airfoil text file and plot it.
    
    >>> import airfoil
    >>> af = airfoil.Airfoil()
    >>> af.read_txt('GA37A315.txt')
    >>> af.plot()
    
    read airfoil from xls database
    
    >>> import airfoil
    >>> af = airfoil.load('GA37A315')
    >>> af.plot('rs--')
    
    create NACA0012 airfoil and build aerodynamic coefficients table
    
    >>> import airfoil
    >>> af = airfoil.Airfoil()
    >>> af.naca4(12.0,0,0)
    >>> af.build_aero_table()
    >>> af.polar.display()
    """
    def __init__(self):
        self.dbPath = paths.Database().airfoil
        self.polar  = AirfoilPolar()
        self.path   = paths.MyPaths()
        self.camber_slope_LE = 0.0
    def read_lines(self,fileID,numLines):
        """
        reads lines from open file and returns 2D array
        
        Parameters
        ----------

        fileID : open file id
            file to read
        numLines : integer
            number of lines to read. Empty lines are skipped.
        """
        data = zeros([numLines,2])
        for i in range(numLines):
            line = fileID.readline()
            if line.strip()!='':
                segLine = line.split()
                data[i,0] = float(segLine[0])
                data[i,1] = float(segLine[1])
        return data
    def read_txt(self,filePath,overwrite=True,afType=1):
        """
        reads in airfoil coordinates file in text format.
        
        Parameters
        ----------
        
        filePath : path
            absolute path of text file to read
        overwrite : bool
            if True then overwrites coordinates in self.coord otherwise returns coordinates as 
            numpy array.
        afType : integer
            type of airfoil file. 1 or 2
        """
        if afType==1:
            out = self.read_txt_1(filePath)
            if overwrite:
                self.coord = out[0]
                self.upPts, self.loPts = geom.separate_coordinates(self.coord)
                self.name = out[1]
            else:
                return out[0]
        elif afType==2:
            out = self.read_txt_2(filePath)
            if overwrite:
                self.upPts = out[0]
                self.loPts = out[1]
                self.join_coordinates()
                self.name  = out[2]
            else:
                return geom.join_coordinates(out[0],out[1])
    def read_txt_1(self,filePath):
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
        afFile = open(filePath,'rt')
        name = afFile.readline()[:-1]
        coord = array([0.0,0.0])
        for line in afFile:
            if line.strip()!='':
                segLine = line.split()
                point = [float(segLine[0]), float(segLine[1])]
                coord = vstack([coord,point])
        afFile.close()
        return coord[1:],name
    def read_txt_2(self,filePath):
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
        afFile = open(filePath,'rt')
        name = afFile.readline()[:-1]
        nPts = [int(float(n)) for n in afFile.readline().split()]
        upPts = self.read_lines(afFile,nPts[0])
        loPts = self.read_lines(afFile,nPts[1])
        afFile.close()
        return upPts, loPts, name
    def write_txt(self,filePath,tab=True,writeName=True):
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
        afFile = open(filePath,'wt')
        if writeName:
            afFile.write('%s\n'%self.name)
        if tab:
            whiteSpace = '\t'
        else:
            whiteSpace = '    '
        for point in self.coord:
            afFile.write('%.6f%s%.6f\n'%(point[0],whiteSpace,point[1]))
        afFile.close()
    def read_xls(self,sheetName,dbPath=''):
        r"""
        Parameters
        ----------
        
        sheetName : string
            name of xls sheet with airfoil data
        dbPath : path
            file path of xls format airfoil database. If dbPath='' then default database will be 
            used. 

        Default db path is 'database/airfoil.xls'.
        Default paths can be edited in paths.py
        Data read are: coordinates, optional: aerodynamic data
        """
        if dbPath=='':
            dbPath=self.dbPath
        dbPath = paths.fixPaths(dbPath) #TODO: update using io.path
        db = dbTools.read(dbPath,sheetName)
        self.name = sheetName
        xCoord    = db.readRow(0,1)
        yCoord    = db.readRow(1,1)
        self.coord = transpose([xCoord,yCoord])
        self.separate_coordinates()
        i = db.findHeader('GEOMETRY')
        if i==-1:
            self.analyze_geometry()
        else:
            self.thickness = float(db.readRow(i+1,1))
            self.camber    = db.readRow(i+2,1)
        i = db.findHeader('ANALYSIS')
        if i==-1:
            self.build_aero_table()
        else:
            self.polar.source = db.readRow(i+1,1)
            i = db.findHeader('LIFT COEFFICIENT')
            nAlpha = db.findHeader('DRAG COEFFICIENT') - i-2
            self.polar.Mach  = array(db.readRow(i+1,1))
            self.polar.alpha = array(db.readCol(i+2,0,nAlpha))
            self.polar.cl    = db.readRange(i+2,1,nAlpha)
            i = db.findHeader('DRAG COEFFICIENT')
            self.polar.cd    = db.readRange(i+2,1,nAlpha)
            i = db.findHeader('MOMENT COEFFICIENT')
            self.polar.cm    = db.readRange(i+2,1,nAlpha)
            self.polar.cl   = transpose(self.polar.cl)
            self.polar.cd   = transpose(self.polar.cd)
            self.polar.cm   = transpose(self.polar.cm)
            self.polar.create_splines(mach=True)
        self.polar.calc_clmax()
        self.create_splines()
    def write_xls(self,dbPath='',includePolars=True):
        r"""
        Writes airfoil data (geometry, aerodynamics) to database.

        Parameters
        ----------
        
        dbPath : path
            file path of airfoil database in xls format. If dbPath is not 
            specified then default database will be used
        """
        if dbPath=='':
            dbPath = self.dbPath
        dbPath = paths.fixPaths(dbPath)
        afDB     = dbTools.loadDB(dbPath,mode='w')
        newSheet = afDB.add_sheet(self.name)
        sheet       = dbTools.writeDB(newSheet)
        sheet.writeRow('X',self.coord[:,0])
        sheet.writeRow('Y',self.coord[:,1])
        sheet.writeRow('GEOMETRY')
        sheet.writeRow('thickness',self.thickness)
        sheet.writeRow('camber'   ,self.camber)
        if includePolars:
            sheet.writeRow('ANALYSIS')
            sheet.writeRow('method', self.polar.source)
            sheet.writeRow('LIFT COEFFICIENT')
            sheet.writeRow('Mach list',self.polar.Mach)
            i = sheet._prevRow + 1
            sheet.writeCol(self.polar.alpha)
            sheet.writeRange(transpose(self.polar.cl),i,1)
            sheet.writeRow('DRAG COEFFICIENT')
            sheet.writeRow('Mach list',self.polar.Mach)
            i = sheet._prevRow + 1
            sheet.writeCol(self.polar.alpha)
            sheet.writeRange(transpose(self.polar.cd),i,1)
            sheet.writeRow('MOMENT COEFFICIENT')
            sheet.writeRow('Mach list',self.polar.Mach)
            i = sheet._prevRow + 1
            sheet.writeCol(self.polar.alpha)
            sheet.writeRange(transpose(self.polar.cm),i,1)
        afDB.save_db()
    def create_CST(self,Au,Al,nPts=30,dist='cos'):
        self.name = 'CST airfoil'
        self.upCurve = geom.CstCurve(Au)
        self.loCurve = geom.CstCurve(Al)
        if dist=='cos':
            x = linspace(0,pi/2.0,nPts)
            x = array([-cos(xx) + 1.0 for xx in x])
        else:
            x = linspace(-pi/2.0,pi/2.0,nPts)
            x = array([0.5*(sin(xx)+1) for xx in x])
        self.upPts = self.upCurve.get_coordinates(x)
        self.loPts = self.loCurve.get_coordinates(x)
        self.join_coordinates()
        self.analyze_geometry()

    def naca4(self,thickness=12.0,camber=0.0,camberLoc=0.0,nPts=120,
              closedTE=True):
        """
        creates NACA 4-series airfoil using javafoil.
        
        Parameters
        ----------
        
        thickness : float
            maximum thickness of airfoil in percent of chord
        camber : float
            maximum camber of airfoil in percent of chord
        camberLoc : float
            location of maximum airfoil camber in percent of chord
        nPts : integer
            number of points
        closedTE : bool
            If True - trailing edge thickness is zero, otherwise finite thickness as defined by NACA 
            equation will be used.
        """
        # TODO: replace javafoil with function
        path = paths.MyPaths()
        self.thickness = float(thickness)/100
        self.camber = float(camber)/100
        self.camberLoc = float(camberLoc)/10
        self.name = 'NACA %d%d%d'%(camber,camberLoc/10,thickness)
        path.setRandPrefix()
        tmpJouFile = path.getTmpFile('jfscript')
        tmpAfFile = path.getTmpFile('dat')
        jouFile = open(tmpJouFile,'wt')
        jouFile.write('Options.Country(0)\n')
        if os.name=='nt':
            jouFile.write('Geometry.CreateAirfoil(0;%d;'%nPts)
            jouFile.write('%.2f;30;%.2f;%.2f;0;0;%d)\n'%(thickness,camber,camberLoc,int(closedTE)))
        else:
            jouFile.write('Geometry.CreateAirfoil(0:%d:'%nPts)   
            jouFile.write('%.2f:30:%.2f:%.2f:0:0:%d)\n'%(thickness,camber,camberLoc,int(closedTE)))
        jouFile.write('Geometry.Save(%s)\nExit()'%tmpAfFile)
        jouFile.close()
        cmd = ('%s -cp %s -jar %s Script=%s'%(path.java,path.mhclasses,path.javafoil,tmpJouFile))
        args=shlex.split(cmd,False,os.name=='posix')
        ps=Popen(args,stdin=PIPE,stderr=PIPE,stdout=PIPE)
        output, errors = ps.communicate()
        ps.stderr.close() 
        self.read_txt(tmpAfFile)
        os.remove(tmpAfFile)
        os.remove(tmpJouFile)
    
    def separate_coordinates(self):
        self.upPts, self.loPts = geom.separate_coordinates(self.coord)
    def join_coordinates(self):
        self.coord = geom.join_coordinates(self.upPts, self.loPts)
    def analyze_geometry(self):
        """
        analyzing airfoil geometry by upper and lower curve points.
        Search maximum thickness and maximum camber using cubic spline 
        interpolation and gradient based optimization. To avoid interpolation 
        errors that can occur at leading edge of several airfoil types 
        (mostly NACA cambered airfols) it is assumed that maximum camber and 
        thickness are located between 5 and 95% of airfoil length.
        
        Result is stored in self.thicknessLoc, self.thickness, self.camber, 
        self.camberLoc
        """
        lb = 0.1
        ub = 0.9
        up = geom.get_pts_in_range(self.upPts,lb,ub)
        lo = geom.get_pts_in_range(self.loPts,lb,ub)
        upCurve = interp1d(up[:,0],up[:,1],'cubic')
        loCurve = interp1d(lo[:,0],lo[:,1],'cubic')
        lb = up[0,0]
        ub = up[-1,0]
        def tc(x):
            return loCurve(x) - upCurve(x)
        def camber(x):
            return - (upCurve(x) + loCurve(x)) / 2.0
        self.thicknessLoc = float(fminbound(tc,lb,ub,xtol=0.001))
        self.camberLoc    = float(fminbound(camber,lb,ub,xtol=0.001))
        self.thickness    = -float(tc(self.thicknessLoc))
        self.camber       = -float(camber(self.camberLoc))
    def create_splines(self):
        """
        create cubic splines for x and y coordinate with respect to curve 
        length parameter
        """
        lenUp = geom.curve_pt_dist_normalized(self.upPts)
        lenLo = geom.curve_pt_dist_normalized(self.loPts)
        self.upCurve = CurveXyt(self.upPts[:,0],self.upPts[:,1],lenUp)
        self.loCurve = CurveXyt(self.loPts[:,0],self.loPts[:,1],lenLo)

    def get_length(self):
        """
        calculates the length of the airfoil curve to be used to calculate 
        wing wetted area or similar purposes.
        """
        return geom.curve_length(self.coord)
    def create_by_interp_txt(self,afPath1,afPath2,afName='interpolated',
                             interpFraction=0.5,graphic=False):
        """
        Interpolates or "blending" airfoil in various proportions using Xfoil.
        
        Parameters
        ----------
        
        afPath1 : path
            file path of first airfoil text fle
        afPath2 : path
            file path of second airfoil text file
        afName : string
            name of result airfoil
        interpFraction : float [0,1]
            defines in which proportion airfoils will be merged. 0 - will 
            return 1st airfoil, 1 - will return 2nd airfoil.
        graphic : bool
            enable/disable Xfoil gui appearance
        """
        xfoil = Xfoil()
        self.path.setRandPrefix()
        tmpAfFile = self.path.getTmpFile('dat')
        xfoil.cmd('INTE')
        xfoil.cmd('F')
        xfoil.cmd('%s'%afPath1)
        xfoil.cmd('F')
        xfoil.cmd('%s'%afPath2)
        xfoil.cmd('%.4f'%interpFraction)
        xfoil.cmd('%s'%afName)
        xfoil.cmd('PCOP')
        xfoil.cmd('SAVE')
        xfoil.cmd('%s'%tmpAfFile)
        xfoil.terminate()
        self.read_txt(tmpAfFile)
        os.remove(tmpAfFile)
    def create_by_interp_airfoils(self,af1,af2,afName='interpolated',
                                  interpFraction=0.5,graphic=False):
        """
        Interpolates two Airfoil objects using Xfoil
        
        Parameters
        ----------
        
        af1 : Airfoil
            first airfoil
        af2 : Airfoil
            second airfoil
        afName : string
            name of result airfoil
        interpFraction : float [0,1]
            defines in which proportion airfoils will be merged. 0 - will 
            return 1st airfoil, 1 - will return 2nd airfoil.
        graphic : bool
            enable/disable Xfoil gui appearance
        """
        self.path.setRandPrefix()
        afFile1 = self.path.getTmpFile('dat','1')
        afFile2 = self.path.getTmpFile('dat','2')
        af1.write_txt(afFile1)
        af2.write_txt(afFile2)
        self.create_by_interp_txt(afFile1,afFile2,afName,interpFraction,graphic)
        os.remove(afFile1)
        os.remove(afFile2)
    def redim(self,numPts=121,overwrite=False,updCurves=True):
        """
        Redimension and redistribution of airfoil points using cosine function. 
        More points are located at leading and trailing edge.
        
        Parameters
        ----------
        
        numPts : integer
            number of points for target airfoil. If number of points is same 
            as source airfoil, points will be redistributed to make smooth 
            surface
        overwrite : bool
            If True then self.coord, self.upPts and self.loPts will be 
            overwritten, otherwise will return array of redimensioned points 
            in format of self.coord
        """
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
        if updCurves:
            self.create_splines()
        xUp, yUp = self.upCurve(tUp)
        xLo, yLo = self.loCurve(tLo)
        upPts = transpose(vstack([xUp,yUp]))
        loPts = transpose(vstack([xLo,yLo]))
        
        coord = geom.join_coordinates(upPts,loPts)
        if overwrite:
            self.coord = coord
            self.upPts = upPts
            self.loPts = loPts
        else:
            return coord
    def scale_thickness(self,newThickness,analysis=True):
        """
        Scales airfoil in Y-direction to achieve required thickness. Useful 
        to be used with propeller analysis when same airfoil has different 
        thickness at different *x* stations.
        
        Parameters
        ----------
        
        newThickness: float
            thickness/chord of result airfoil
        analysis : bool
            If True then aerodynamic table will be generated using 
            self.build_aero_table, otherwise only geometry will be changed.
        """
        scale = newThickness / self.thickness
        self.name = self.name + '_tc%.2f'%(newThickness*100.0)
        self.coord[:,1] = self.coord[:,1]*scale
        self.separate_coordinates()
        self.redim(50,overwrite=True)
        if analysis:
            self.build_aero_table()

    def build_aero_table(self,MachSeq=[0.1,0.9,0.1],alphaSeq=[-10,25,1.0],
                         mode='javafoil'):
        """
        builds full table of aerodynamic coefficients using *fast* solvers
        Coefficients to be caclulated are as follows:
            
            - lift coefficient vs. alpha, Mach
            - drag coefficient vs. alpha, Mach
            - moment coefficient vs. alpha, Mach
        
        Parameters
        ----------
        
        MachSeq : array float
            array of Mach number to be analyzed in format [start, end, step]
        alphaSeq : array float
            array of angle of attack sequence to be analyzed in format 
            [start, end, step]
        mode : string
            if mode='javafoil' then coefficients will be generated using 
            javafoil, if 'xfoil' then using xfoil
        """
        Mach = arange(MachSeq[0], MachSeq[1], MachSeq[2])
        Re = []
        for i,M in enumerate(Mach):
            Re.append(fc.FlightConditions(speed=M).Re)
            if mode=='javafoil':
                tmpPolar = self.get_J_polar(M, Re[i], alphaSeq)
            elif mode=='xfoil':
                tmpPolar = self.get_X_polar(M, Re[i], alphaSeq, 200)
            if i==0:
                self.polar.cl = tmpPolar.cl
                self.polar.cd = tmpPolar.cd
                self.polar.cm = tmpPolar.cm
            else:
                self.polar.cl = vstack([self.polar.cl, tmpPolar.cl])
                self.polar.cd = vstack([self.polar.cd, tmpPolar.cd])
                self.polar.cm = vstack([self.polar.cm, tmpPolar.cm])
        self.polar.Re     = Re
        self.polar.source = 'javafoil'
        self.polar.Mach   = Mach
        self.polar.alpha  = tmpPolar.alpha
        self.polar.create_splines(mach=True)
    def get_X_polar(self,Mach,Re,alphaSeq=[-10,25,1.0],nIter=100,
                    graphic=False,smooth=False):
        """
        Calculates aerodynamic coefficients at given flight conditions using Xfoil.
        
        Parameters
        ----------
        
        Mach : float
            Mach number
        Re : float
            Reynolds number
        alphaSeq : array float
            array of angle of attack sequence to be analyzed in format 
            [start, end, step]
        nIter : integer
            number of xfoil iterations. Default value is that is enough for 
            conventional airfoil configurations. For complex airfoil shapes 
            this number should be increased (for example while optimization 
            using genetic algorithm)
        
        Returns
        -------
        
        polar : AirfoilPolar
            airfoil polar with all aerodynamic data at specified flight 
            conditions
        """
        self.path.setRandPrefix()
        tmpAfFile = self.path.getTmpFile('dat')
        tmpPolar  = self.path.getTmpFile('pol')
        self.write_txt(tmpAfFile,False)
        xfoil = Xfoil()
        xfoil.cmd
        xfoil.cmd('LOAD\n%s'%tmpAfFile)
        if smooth:
            xfoil.cmd('GDES\nCADD\n\n\n\n\nPANEL')
        xfoil.cmd('OPER\nVISC\n%.0f\nMACH\n%.4f'%(Re,Mach))
        if nIter>10:
            xfoil.cmd('ITER\n%d' % nIter)
        xfoil.cmd('PACC')
        xfoil.cmd(' ')
        xfoil.cmd(' ')
        if alphaSeq[2] == 0 or alphaSeq[0]==alphaSeq[1]:
            xfoil.cmd('ALFA\n%.2f'%alphaSeq[0])
        elif alphaSeq[0]*alphaSeq[1]<0:
            xfoil.cmd('ASEQ\n%.2f\n%.2f\n%.2f'%(alphaSeq[2], alphaSeq[1], alphaSeq[2]))
            xfoil.cmd('ASEQ\n%.2f\n%.2f\n%.2f'%(0, alphaSeq[0], -alphaSeq[2]))
        elif alphaSeq[0]*alphaSeq[1]>=0 and alphaSeq[0]>=0:
            xfoil.cmd('ASEQ\n%.2f\n%.2f\n%.2f'%(alphaSeq[0], alphaSeq[1], alphaSeq[2]))
        else:
            xfoil.cmd('ASEQ\n%.2f\n%.2f\n%.2f'%(alphaSeq[1], alphaSeq[0], -alphaSeq[2]))
        xfoil.cmd('PWRT')
        xfoil.cmd('%s'%(tmpPolar))
        xfoil.terminate()
        polar = self._read_X_polar(tmpPolar)
        polar.Mach = Mach
        polar.Re = Re
        os.remove(tmpAfFile)
        os.remove(tmpPolar)
        return polar
    def get_J_polar(self,Mach,Re,alphaSeq=[-10,25,1.0]):
        """
        Calculates aerodynamic coefficients at given flight conditions using Xfoil.
        
        Parameters
        ----------
        
        Mach : float
            Mach number
        Re : float
            Reynolds number
        alphaSeq : array float
            array of angle of attack sequence to be analyzed in format 
            [start, end, step]
        nIter : integer
            number of xfoil iterations. Default value is that is enough for 
            conventional airfoil configurations. For complex airfoil shapes 
            this number should be increased (for example while optimization 
            using genetic algorithm)
        
        Returns
        -------
        
        polar : AirfoilPolar
            airfoil polar with all aerodynamic data at specified flight 
            conditions
        
        Note
        ----
        
        Javafoil requires path to java installed in system in 
        :file:`javapath.txt`. For linux OS file should contain *java* 
        keyword only
        """
        self.path.setRandPrefix()
        tmpJournal = self.path.getTmpFile('jfscript')
        tmpAfFile = self.path.getTmpFile('dat')
        tmpPolar  = self.path.getTmpFile('pol')
        self.write_txt(tmpAfFile)
        jouFile = open(tmpJournal,'wt')
        jouFile.write('Options.Country(0)\nGeometry.Clear()\n')
        jouFile.write('Geometry.Open(/%s/)\n' %tmpAfFile)
        jouFile.write('Options.MachNumber(%.4f)\n'%Mach)
        jouFile.write('Options.StallModel(0)\n')
        jouFile.write('Options.TransitionModel(1)\n')
        jouFile.write('Options.GroundEffect(0)\n')
        jouFile.write('Options.HeightOverSpan(0.5)\n')
        jouFile.write('Options.AspectRatio(0)\n')
        if os.name=='nt':
            jouFile.write('Polar.Analyze(%.0f;%.0f;%.0f;%.0f;%.0f;%.2f;100;100;0;0)\n'%(Re,Re,0,alphaSeq[0], alphaSeq[1], alphaSeq[2]))    
        else:                
            jouFile.write('Polar.Analyze(%.0f:%.0f:%.0f:%.0f:%.0f:%.2f:100:100:0:0)\n'%(Re,Re,0,alphaSeq[0], alphaSeq[1], alphaSeq[2]))
        jouFile.write('Polar.Save(/%s/)\n'%tmpPolar)
        jouFile.write('Exit()')
        jouFile.close()
        cmd = ('%s -cp %s -jar %s Script=%s'%(self.path.java,self.path.mhclasses,self.path.javafoil,tmpJournal))
        args=shlex.split(cmd,False,os.name=='posix')
        ps=Popen(args,stdin=PIPE,stderr=PIPE,stdout=PIPE)
        output, errors = ps.communicate()
        ps.stderr.close() 
        polar = self._read_J_polar(tmpPolar)
        polar.Mach  = Mach
        polar.Re = Re
        os.remove(tmpJournal)
        os.remove(tmpAfFile)
        os.remove(tmpPolar)
        return polar
    def _read_X_polar(self,polarPath):
        """
        reads polar file of Xfoil and returns AirfoilPolar object
        """
        polar = AirfoilPolar()
        polarFile = open(polarPath,'rt')
        lines = polarFile.readlines()
        polarFile.close()
        del lines[0:12]
        newLines = [float(s) for s in lines[0].split()]
        for line in lines[1::]:
            addLine = [float(s) for s in line.split()]
            newLines = vstack([newLines,addLine])
        newIdx      = argsort(newLines[:,0])
        polar.alpha = array(newLines[newIdx,0])
        polar.cl    = array(newLines[newIdx,1])
        polar.cd    = array(newLines[newIdx,2])
        polar.cd_p  = array(newLines[newIdx,3])
        polar.cm    = array(newLines[newIdx,4])
        polar.TU    = array(newLines[newIdx,5])
        polar.TL    = array(newLines[newIdx,6])
        return polar
    def _read_J_polar(self,polarPath):
        """
        reads polar file of javafoil and returns AirfoilPolar object
        """
        polar = AirfoilPolar()
        polarFile = open(polarPath,'rt')
        lines = polarFile.readlines()
        polarFile.close()
        del lines[0:5]
        for line in lines:
            if line.strip()!='':
                segLine = line.split()
                try:
                    polar.alpha = append(polar.alpha, float(segLine[0]))
                    polar.cl    = append(polar.cl, float(segLine[1]))
                    polar.cd    = append(polar.cd, float(segLine[2]))
                    polar.cm    = append(polar.cm, float(segLine[3]))
                    polar.TU    = append(polar.TU, float(segLine[4]))
                    polar.TL    = append(polar.TL, float(segLine[5]))
                    polar.SU    = append(polar.SU, float(segLine[6]))
                    polar.SL    = append(polar.SL, float(segLine[7]))
                    polar.LD    = append(polar.LD, float(segLine[8]))
                    polar.AC    = append(polar.AC, float(segLine[9]))
                    polar.CP    = append(polar.CP, float(segLine[10]))
                except ValueError:
                    print 'Polar processing failed'
        return polar
    def create_af_CAT(self,chord = 1.0,batch = True, save = [],filetype = 'igs'):
        CATIA = win32com.client.Dispatch('catia.application')
        partDocument1 = CATIA.Documents.Add("Part")
        partDocument1 = CATIA.ActiveDocument
        Part1 = partDocument1.Part
        hSF = Part1.HybridShapeFactory
        hybridBodies1 = Part1.HybridBodies
        hybridBody1 = hybridBodies1.Item("Geometrical Set.1")
        
        def curve2D(mode):
            le_camber_angle = radians(self.camber_slope_LE)
            if mode == 'up':
                coordinates = self.upPts * chord
                direction = -1
            elif mode == 'lo':
                coordinates = self.loPts * chord
                direction = 1

            Spline1 = hSF.AddNewSpline()
            if self.radius_LE == []:
                start = 0
            else:
                start = 1
                radius_LE = self.radius_LE * chord/100
                Xctr = coordinates[0,0] + radius_LE*cos(le_camber_angle)
                Yctr = coordinates[0,1] + radius_LE*sin(le_camber_angle)
                ctrPC = hSF.AddNewPointCoord(Xctr,Yctr,0)
                LE_pt = hSF.AddNewPointCoord(coordinates[0,0],coordinates[0,1],0)
                originElements = Part1.OriginElements
                plane = originElements.PlaneXY
                ref1 = Part1.CreateReferenceFromObject(ctrPC)
                ref2 = Part1.CreateReferenceFromObject(LE_pt)
                ref3 = Part1.CreateReferenceFromObject(plane)
                LE_circle = hSF.AddNewCircleCtrPt(ref1,ref2,ref3,False)
                LE_circle.SetLimitation(1)
                ref4 = Part1.CreateReferenceFromObject(LE_circle)
                Spline1.AddPointWithConstraintFromCurve(ref2,ref4,1,direction,2)
                
            for pt in coordinates[start:]:
                HSPC = hSF.AddNewPointCoord(pt[0],pt[1],0)
                reference = Part1.CreateReferenceFromObject(HSPC)
                Spline1.AddPointWithConstraintExplicit(reference, None, -1, 1, None, 0)
            hybridBody1.AppendHybridShape(Spline1)
            Part1.update()
            
        curve2D('up')
        curve2D('lo')
        if save != []: partDocument1.ExportData(save,'igs')
        if batch: 
            partDocument1.Close()
            CATIA.quit()

    def write_polar_txt(self, path):
        if self.polar.alpha !=[]:
            fpolar = open(path,'wt')
            fpolar.write('Airfoil: %s\n'%self.name)
            fpolar.write('Re = %.0f\n'%self.polar.Re)
            fpolar.write('Mach = %.4f\n'%self.polar.Mach)
            fpolar.write('alpha\tCL\tCD\tCM\n')
            for ii in range(self.polar.alpha.shape[0]):
                fpolar.write('%.2f\t'%self.polar.alpha[ii])
                fpolar.write('%.6f\t'%self.polar.cl[ii])
                fpolar.write('%.6f\t'%self.polar.cd[ii])
                fpolar.write('%.6f\n'%self.polar.cm[ii])
            fpolar.close()

    def plot(self,linetype='bo-'):
        """
        plots an airfoil using matplotlib pytplot
        
        Parameters
        ----------
        
        linetype : string
            linetype in matplotlb format
        """
        figure()
        grid(True)
        plot(self.coord[:,0],self.coord[:,1],linetype)
        xlim((0.0,1.0))
        axis('equal')
        title(self.name)
        show()
    def report(self):
        #TODO: print out information about airfoil: name, number of points etc.
        pass
    
    def set_trailing_edge(self,zTE=0.0):
        zTEcurrent = self.upPts[-1,1] - self.loPts[-1,1]
        zTEnew = zTE - zTEcurrent
        ptsUpNew = self.upPts[:,1] + self.upPts[:,0]*zTEnew/2.
        ptsLoNew = self.loPts[:,1] - self.loPts[:,0]*zTEnew/2.
        self.upPts[:,1] = ptsUpNew
        self.loPts[:,1] = ptsLoNew
        self.join_coordinates()
        self.analyze_geometry()

def batch_txt_xls():
    dir1 = os.getcwd() + '\\database\\prop_airfoils'
    files = os.listdir(dir1)
    for afFile in files:
        afPath = dir1 + '\\' + afFile
        af = Airfoil()
        af.read_txt(afPath)
        af.build_aero_table()
        af.write_xls(dbPath=r'D:\codes\actools\pyAC\actools\database\prop_airfoil.xls')
        print '%s\t completed'%af.name

def test_function1():
    # create NACA2312 airfoil
    thickness = 12
    camber    = 2
    camberLocation = 30
    af = Airfoil()
    af.naca4(thickness,camber,camberLocation)
    # analyze using X-foil and Javafoil
    Mach = 0.16
    Re   = 2.2e6
    alphaStart = -20.0
    alphaEnd   = 20
    alphaStep  = 1
    alphaSequence = [alphaStart, alphaEnd, alphaStep]
    xFoilPolar    = af.get_X_polar(Mach,Re,alphaSequence, smooth=True)
    javaFoilPolar = af.get_J_polar(Mach,Re,alphaSequence)
    # plot the results of analysis
    xFoilPolar.display()
    javaFoilPolar.display()

def test_function2():
    # load GA37A315 airfoil from Excel airfoil database
    af = load('GA37A315')
    # plot it
    af.plot()

def test_function3():
    # read airfoil from text file
    #airfoilPath = 'GA37A315mod_reworked.txt'
    af = Airfoil()
    #af.read_txt(airfoilPath)
    af.naca4(15,2,30)
    # create aerodynamic characteristics table
    #af.build_aero_table()
    # save to excel database
    af.polar = af.get_X_polar(0.16,3e6,[-20.,20,1.0])
    figure(13)
    hold(True)
    plot(af.polar.alpha, af.polar.TU,'b-')
    plot(af.polar.alpha, af.polar.TL,'r-')
    plot(af.polar.alpha, af.polar.TU+af.polar.TL,'g--')
    show()

def debug1():
    # read airfoil from text file
    airfoilPath = 'GA37A315mod_reworked.txt'
    af = Airfoil()
    af.read_txt(airfoilPath)
    # create aerodynamic characteristics table
    af.build_aero_table()
    Mach = af.polar.Mach[-1]
    alpha = linspace(-30,30)
    cl = zeros(len(alpha))
    for i,a in enumerate(alpha):
        cl[i] = af.polar.clAlpha(Mach,a)
    print af.polar.TU
    figure(13)
    hold(True)
    plot(alpha,cl,'ro')
    plot(af.polar.alpha, af.polar.TU[-1],'b-')
    plot(af.polar.alpha, af.polar.TL[-1],'b-')
    show()

def debug2():
    afList = AirfoilList()
    names = ['GA37A315','GA37A315','GA37A315','NACA0012']
    for name in names:
        afList.add_airfoil_xls(name,False)
    afList[0].scale_thickness(0.25)
    afList[1].scale_thickness(0.20)
#    afList[0].plot()
#    afList[1].plot()
    afList[-1].polar.display()
    print afList.contents
    
def test_04():
    from miscTools import Timer
    timer = Timer()
    af = load('NACA0012')
    Mach = linspace(0.2,0.8,20)
    alpha = linspace(-10,10,20)
    timer.start()
    for i in range(10000):
        for a in alpha:
            for M in Mach:
                cl = af.polar.clAlpha(M,a)
                cd = af.polar.cdAlpha(M,a)
                cm = af.polar.cmAlpha(M,a)
    timer.finish()

def test_05():
    af = Airfoil()
    af.naca4()
    af.set_trailing_edge(0.02)
    pol = af.get_X_polar(0.15,3e6,[-10,20,1.0])
    pol.create_splines()
    print pol.get_cd_at_cl(0.5)
    pol.calc_clmax()
    print pol.clmax
    plt.figure(1)
    plt.plot(pol.alpha,pol.cd,'ro-')
    plt.show()
    af.plot()

def test_06():
    af = load('NACA0012')
    lift = FlatPlateInterp(af.polar.Mach,af.polar.alpha,af.polar.cl,'lift')
    drag = FlatPlateInterp(af.polar.Mach,af.polar.alpha,af.polar.cd,'drag')
    anew = linspace(-90,90,200)
    clnew = lift(0.4,anew)
    plt.figure(1)
    plt.plot(anew,clnew,'ro-')
    plt.show()

def test_cst():
    Au = array([0.2337,0.7288,0.7400,0.0154])
    Al = array([-0.2337,-0.2271,-0.4580,-1.0059])
    af = cst(Au,Al)
    af.plot()
    af.set_trailing_edge(0.02)
    af.plot()
    #af.create_af_CAT(batch=False)

if __name__=="__main__":
    test_05()