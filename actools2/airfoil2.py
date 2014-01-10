# -*- coding: utf-8 -*-
"""
Created on Wed Jan 08 14:02:55 2014

@author: Maxim
"""
from paths import MyPaths
import numpy as np
import matplotlib.pyplot as plt
import db_tools
from scipy.interpolate import interp1d
from scipy.optimize import fminbound
from misc_tools import Timer
import geometry as geom
from airfoil_polar import AirfoilPolar
import xfoil_analysis as xf
import javafoil_analysis as jf

pth = MyPaths()

def load(airfoilName,dbPath=None):
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
    af.read_db(airfoilName,dbPath)
    return af

def cst(Au,Al):
    pass

def naca4(thickness,camber,camberLocation):
    pass


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
        return np.array([xnew,ynew])
        

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
            K = np.ones([order+1],int)
            for ii in range(1,order):
                K[ii] = np.factorial(order)/(np.factorial(ii)*np.factorial(order-ii))
            return K
    def __init__(self,A,N=[0.5,1.0]):
        self.classCurve = self.ClassFcn(N)
        self.shapeCurve = self.ShapeFcn(A)
    def __call__(self,x):
        return self.classCurve(x) * self.shapeCurve(x)


class Airfoil:
    def __init__(self):
        self._db = pth.db.airfoil
        self.name             = None
        self.coord            = None
        self.ptsUp            = None
        self.ptsLo            = None
        self.thickness        = None
        self.camber           = None
        self.camberLocation   = None
        self._curveUp         = None
        self._curveLo         = None
        self._curveParametric = None
        self.polar = AirfoilPolar()

    def read_db(self,name,xlsPath=None):
        """
        Reads airfoil from xls database. If additional geometry information 
        and aerodynamic tables exist then existing tables are used 
        otherwise information will be created using internal tools of module.
        """
        if xlsPath==None:
            xlsPath = pth.db.airfoil
        db = db_tools.ReadDatabase(xlsPath,name)
        self.name = name
        xCoord    = db.read_row(0,1)
        yCoord    = db.read_row(1,1)
        self.coord = np.transpose([xCoord,yCoord])
        self._separate_coordinates()
        # geometry
        i = db.find_header('GEOMETRY')
        if i==-1:
            self._analyze_geometry()
        else:
            self.thickness = db.read_row(i+1,1)
            self.camber    = db.read_row(i+2,1)
            self._analyze_geometry()
        # analysis
        i = db.find_header('ANALYSIS')
        if i==-1:
            self.build_aero_table()
        else:
            self.polar.source = db.read_row(i+1,1)
            i = db.find_header('LIFT COEFFICIENT')
            nAlpha = db.find_header('DRAG COEFFICIENT') - i-2
            self.polar.Mach  = np.array(db.read_row(i+1,1))
            self.polar.alpha = np.array(db.read_column(0,i+2,nAlpha))
            self.polar.cl    = db.read_row_range(i+2,1,nAlpha)
            i = db.find_header('DRAG COEFFICIENT')
            self.polar.cd    = db.read_row_range(i+2,1,nAlpha)
            i = db.find_header('MOMENT COEFFICIENT')
            self.polar.cm    = db.read_row_range(i+2,1,nAlpha)
            self.polar.cl    = np.transpose(self.polar.cl)
            self.polar.cd    = np.transpose(self.polar.cd)
            self.polar.cm    = np.transpose(self.polar.cm)
            self.polar._create_splines()
        self.polar._calc_clmax()

    def write_db(self,xlsPath=None,includePolars=True):
        """
        Writes airfoil to xls format database.
        """
        if xlsPath==None:
            xlsPath=pth.db.airfoil
        db = db_tools.WriteDatabase(xlsPath,self.name)
        db.write_row('X',self.coord[:,0])
        db.write_row('Y',self.coord[:,1])
        db.write_row('GEOMETRY')
        db.write_row('thickness',self.thickness)
        db.write_row('camber'   ,self.camber)
        if includePolars:
            db.write_row('ANALYSIS')
            db.write_row('method', self.polar.source)
            db.write_row('LIFT COEFFICIENT')
            db.write_row('Mach list',self.polar.Mach)
            i = db._irowPrev
            db.write_column(self.polar.alpha)
            db._irowPrev = i
            db.write_row_range(np.transpose(self.polar.cl),-1,1)
            db.write_row('DRAG COEFFICIENT')
            db.write_row('Mach list',self.polar.Mach)
            i = db._irowPrev
            db.write_column(self.polar.alpha)
            db._irowPrev = i
            db.write_row_range(np.transpose(self.polar.cd),-1,1)
            db.write_row('MOMENT COEFFICIENT')
            db.write_row('Mach list',self.polar.Mach)
            i = db._irowPrev
            db.write_column(self.polar.alpha)
            db._irowPrev = i
            db.write_row_range(np.transpose(self.polar.cm),-1,1)
        db.save()

    def read_txt(self,airfoilPath):
        """
        Reads airfoil coordinates from text file. Two formats of coordinates 
        are supported.
        """
        fid = open(airfoilPath,'rt')
        self.name = str(fid.readline()[:-1])
        lines = fid.readlines()
        fid.close()
        seg = lines[0].split()
        if float(seg[0])>3 and float(seg[1])>3:
            iup = int(round(float(seg[0])))
            ilo = int(round(float(seg[1])))
            self._read_txt_type2(lines[1:],iup,ilo)
        else:
            self._read_txt_type1(lines)
        self._analyze_geometry()
        
    def write_txt(self,filePath,tab=True,writeName=True):
        """
        writes airfoil coordinates to text file.
        
        Parameters
        ----------
        
        filePath : path
            airfoil text file path to write coordinates
        tab : bool
            white space parameter. If True then x,y coordinates are separated 
            by tabbing otherwise 
            separated by 4 spacings as required for X-foil input.
        writeName : bool
            if True then airfoil name will be written to the text file header 
            otherwise only coordinates will be stored
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

    def _line_coord_to_float(self,lines):
        output = list()
        for line in lines:
            if line.strip()!='':
                seg = line.split()
                outline = np.zeros(2)
                outline[0] = float(seg[0])
                outline[1] = float(seg[1])
                output.append(outline)
        return np.array(output)
    
    def _separate_coordinates(self):
        self.ptsUp, self.ptsLo = geom.separate_coordinates(self.coord)

    def _join_coordinates(self):
        self.coord = geom.join_coordinates(self.ptsUp, self.ptsLo)

    def _read_txt_type1(self,lines):
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
        self.coord = self._line_coord_to_float(lines)
        self._separate_coordinates()

    def _read_txt_type2(self,lines,nPtsUp,nPtsLo):
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
        coordRaw = self._line_coord_to_float(lines)
        self.ptsUp = coordRaw[:nPtsUp]
        self.ptsLo = coordRaw[nPtsUp:nPtsUp+nPtsLo]
        self._join_coordinates()

    def display(self,linetype='k-'):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid(True)
        ax.axis('equal')
        ax.set_title(self.name)
        ax.plot(self.coord[:,0],self.coord[:,1],linetype)
        plt.show()
    
    def _create_splines(self):
        """
        create cubic splines for x and y coordinate with respect to curve 
        length parameter
        """
        # parametric splines x = x(t), y = y(t)
        lenUp = geom.curve_pt_dist_normalized(self.ptsUp)
        lenLo = geom.curve_pt_dist_normalized(self.ptsLo)
        self.upCurve = CurveXyt(self.ptsUp[:,0],self.ptsUp[:,1],lenUp)
        self.loCurve = CurveXyt(self.ptsLo[:,0],self.ptsLo[:,1],lenLo)
        # non parametric splines y = y(x)
        ptsUpSorted = geom.sort_airfoil_coordinates(self.ptsUp)
        ptsLoSorted = geom.sort_airfoil_coordinates(self.ptsLo)
        self.upCurve2 = interp1d(ptsUpSorted[:,0],ptsUpSorted[:,1],'cubic')
        self.loCurve2 = interp1d(ptsLoSorted[:,0],ptsLoSorted[:,1],'cubic')

    def _analyze_geometry(self):
        """
        calculates airfoil geometry parameters
        """
        self._create_splines()
        _tc  = lambda x: -(self.upCurve2(x) - self.loCurve2(x))
        _cam = lambda x: -(self.upCurve2(x) + self.loCurve2(x))
        xtc = fminbound(_tc,0.1,0.9)
        xcam = fminbound(_cam,0.1,0.9)
        self.thickness = -_tc(xtc)
        self.camber = -_cam(xcam)
        self.camberLocation = xcam

    def redim(self,nPts,overwrite=False,updCurves=False):
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
            phi = x*np.pi
            return (np.cos(phi)+1.0)/2.0
        if (nPts % 2)==1:
            nPts = int(nPts)/2+1
            tUp = np.linspace(0.0,1.0,nPts)
            tLo = tUp
        else:
            nPts = int(nPts)/2+1
            tUp = np.linspace(0.0,1.0,nPts)
            tLo = np.linspace(0.0,1.0,nPts-1)
        tUp = np.array([cos_curve(t) for t in tUp])
        tLo = np.array([cos_curve(t) for t in tLo])
        if updCurves:
            self.create_splines()
        xUp, yUp = self.upCurve(tUp)
        xLo, yLo = self.loCurve(tLo)
        upPts = np.transpose(np.vstack([xUp,yUp]))
        loPts = np.transpose(np.vstack([xLo,yLo]))
        coord = geom.join_coordinates(upPts,loPts)
        if overwrite:
            self.coord = coord
            self.upPts = upPts
            self.loPts = loPts
        else:
            return coord

    def scale_thickness(self,thicknessNew,analysis=False):
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
        scale = thicknessNew / self.thickness
        self.name = self.name + '_tc%.2f'%(thicknessNew*100.0)
        self.coord[:,1] = self.coord[:,1]*scale
        self._separate_coordinates()
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
        Mach = np.arange(MachSeq[0], MachSeq[1], MachSeq[2])
        Re = list()
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
                self.polar.cl = np.vstack([self.polar.cl, tmpPolar.cl])
                self.polar.cd = np.vstack([self.polar.cd, tmpPolar.cd])
                self.polar.cm = np.vstack([self.polar.cm, tmpPolar.cm])
        self.polar.Re     = Re
        self.polar.source = 'javafoil'
        self.polar.Mach   = Mach
        self.polar.alpha  = tmpPolar.alpha
        self.polar.create_splines(mach=True)

    def get_xfoil_polar(self,Mach,Re,alphaSeq=[-15,15,1],nIter=10,graphic=False,smooth=False):
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
        return xf.get_xfoil_analysis(self,Mach,Re,alphaSeq,nIter,graphic,smooth)

    def get_jfoil_polar(self,Mach,Re,alphaSeq=[-15,15,1]):
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
        return jf.get_javafoil_analysis(self,Mach,Re,alphaSeq)


# --- debug section ---
def run_test_geometry():
    af = Airfoil()
    #af.read_db('NACA0012')
    af.read_txt('clarky.txt')
    print af.thickness
    print af.camber
    print af.camberLocation
    af.display()

def run_test_aero_analysis():
    af = Airfoil()
    timer = Timer()
    af.read_db('NACA0012')
    timer.lap('read db')    
    af.get_jfoil_polar(0.2,3e6)
    timer.lap('javafoil')
    af.get_xfoil_polar(0.2,3e6)
    timer.stop('xfoil')

if __name__=="__main__":
    run_test_geometry()