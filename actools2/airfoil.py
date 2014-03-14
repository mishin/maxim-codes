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
import flight_conditions as fc

pth = MyPaths()

# --- wrapper functions ---

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

def cst(Au,Al,nPts=30,zTE=0.0):
    af = Airfoil()
    af.create_cst(Au,Al,nPts,zTE)
    return af

def cst_x(A,nPts=25):
    af = Airfoil()
    n = len(A)
    if n%2==0:
        Au = A[:n/2]
        Al = A[n/2:]
    else:
        Au = np.array(A[:int(n/2)+1])
        Al = np.hstack([-A[0],A[int(n/2)+1:]])
    af.create_cst(Au,Al,nPts)
    return af

def naca4(thickness,camber,camberLocation):
    af = Airfoil()
    af.create_naca4(thickness,camber,camberLocation)
    return af

def read_txt(afPath):
    af = Airfoil()
    af.read_txt(afPath)
    return af

# --- ---

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
        self.pts              = None
        self.ptsUp            = None
        self.ptsLo            = None
        self.thickness        = None
        self.camber           = None
        self.camberLocation   = None
        self.zTrailingEdge    = None
        self._curveUp         = None
        self._curveLo         = None
        self.polar            = None
        self._distribution    = 'cos'

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
        self.pts = np.transpose([xCoord,yCoord])
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
        db.write_row('X',self.pts[:,0])
        db.write_row('Y',self.pts[:,1])
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
        for point in self.pts:
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
        self.ptsUp, self.ptsLo = geom.separate_coordinates(self.pts)

    def _join_coordinates(self):
        self.pts = geom.join_coordinates(self.ptsUp, self.ptsLo)

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
        self.pts = self._line_coord_to_float(lines)
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

    def display(self,linetype='ko-'):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.grid(True)
        ax.axis('equal')
        ax.set_title(self.name)
        ax.plot(self.pts[:,0],self.pts[:,1],linetype)
        plt.show()
    
    def _create_splines(self):
        """
        create cubic splines for x and y coordinate with respect to curve 
        length parameter
        """
        # parametric splines x = x(t), y = y(t)
        lenUp = geom.curve_pt_dist_normalized(self.ptsUp)
        lenLo = geom.curve_pt_dist_normalized(self.ptsLo)
        self._curveUp = CurveXyt(self.ptsUp[:,0],self.ptsUp[:,1],lenUp)
        self._curveLo = CurveXyt(self.ptsLo[:,0],self.ptsLo[:,1],lenLo)

    def _analyze_geometry(self):
        """
        calculates airfoil geometry parameters
        """
        self._create_splines()
        _ymax = lambda t: -self._curveUp(t)[1]
        _ymin = lambda t:  self._curveLo(t)[1]
        t1 = fminbound(_ymax,0,1)
        t2 = fminbound(_ymin,0,1)
        self.thickness = self._curveUp(t1)[1] - self._curveLo(t2)[1]
    
    def _get_point_distribution(self,nPts=30,distribution=None):
        if distribution==None:
            distribution = self._distribution
        if distribution=='sin':
            return geom.get_sine_distribution(nPts)
        elif distribution=='cos':
            return geom.get_cosine_distribution(nPts)
    def redim(self,nPts,overwrite=True,distribution='sin'):
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
            If True then self.pts, self.upPts and self.loPts will be 
            overwritten, otherwise will return array of redimensioned points 
            in format of self.pts
        """
        nPts *= 0.5
        t = self._get_point_distribution(nPts,distribution)
        self._create_splines()
        xUp, yUp = self._curveUp(t)
        xLo, yLo = self._curveLo(t)
        upPts = np.transpose(np.vstack([xUp,yUp]))
        loPts = np.transpose(np.vstack([xLo,yLo]))
        coord = geom.join_coordinates(upPts,loPts)
        if overwrite:
            self.pts = coord
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
        self.pts[:,1] = self.pts[:,1]*scale
        self._separate_coordinates()
        self.redim(50,overwrite=True)
        if analysis:
            self.build_aero_table()
    
    def build_aero_table(self,MachSeq=[0.1,0.9,0.1],alphaSeq=[-20,20,1.0],
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
            Re.append(fc.FlightConditions(M,0.0).Re)
            if mode=='javafoil':
                tmpPolar = self.get_jfoil_polar(M, Re[i], alphaSeq)
            elif mode=='xfoil':
                tmpPolar = self.get_xfoil_polar(M, Re[i], alphaSeq, 200)
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
        self.polar._create_splines()

    def get_xfoil_polar(self,Mach,Re,alphaSeq=[-15,15,1],nIter=10,graphic=False,smooth=False,
                        flapChordRatio=0.3,flapDefl=None):
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
        return xf.get_xfoil_analysis(self,Mach,Re,alphaSeq,nIter,graphic,smooth,flapChordRatio,flapDefl)

    def get_jfoil_polar(self,Mach,Re,alphaSeq=[-15,15,1], 
                        stall='eppler',transition='drelaAfter1991',surface='NACAstandard',
                        flapDefl=None,flapChordRatio=30.0):
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
        stall : string
            stall model. Available options: calcfoil, eppler
        transition : string
            transition model. Availabel options: epplerStandard, epplerExtended, 
            michel1, michel2, H12Re, granville, drelaBefore1991, drelaAfter1991, 
            arnal
        surface : string
            surface type. Available options: smooth, paintedFabrid, 
            NACAstandard, bugsDirt
        
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
        return jf.get_javafoil_analysis(self,Mach,Re,alphaSeq,stall,transition,
                                        surface,flapDefl,flapChordRatio)
    
    def create_cst(self,Au,Al,nPts=25,zTE=None):
        self.upCurve = geom.CstCurve(Au)
        self.loCurve = geom.CstCurve(Al)
        x = self._get_point_distribution(nPts)
        self.ptsUp = self.upCurve.get_coordinates(x)
        self.ptsLo = self.loCurve.get_coordinates(x)
        self._join_coordinates()
        if not zTE==None:
            self.set_trailing_edge(zTE)
        self._analyze_geometry()
        self.name = 'CST airfoil'

    def create_naca4(self,thickness=12.0,camber=0.0,camberLoc=0.0,nPts=25):
        t = thickness / 100.0
        m = camber / 100.0
        p = camberLoc / 100.0
        x = self._get_point_distribution(nPts)
        y = t/0.2*(0.2969*x**0.5 - 0.1281*x - 0.3516*x*x + 0.2843*x**3.0 - 0.1015*x**4.0)
        self.thickness = t
        self.camber = m
        self.camberLocation = p
        self.name = 'NACA%.0f%.0f%.0f'%(camber,camberLoc/10.0,thickness)
        if not m*p==0.0:
            yc = np.zeros(len(x))
            tanTheta = np.zeros(len(x))
            for i,xpt in enumerate(x):
                if xpt<=p:
                    yc[i] = m*xpt/(p*p)*(2.0*p - xpt)
                    tanTheta[i] = 2.0*m/(p*p)*(p-xpt)
                else:
                    yc[i] = m*(1.0-xpt)/((1.0-p)**2.0)*(1.0+xpt-2.0*p)
                    tanTheta[i] = -2.0*m/((1.0-p)**2.0)*(xpt-p)
            theta = np.arctan(tanTheta)
            xu = x - y*np.sin(theta)
            xl = x + y*np.sin(theta)
            yu = yc + y*np.cos(theta)
            yl = yc - y*np.cos(theta)
            self.ptsUp = np.transpose(np.vstack([xu,yu]))
            self.ptsLo = np.transpose(np.vstack([xl,yl]))            
        else:
            self.ptsUp = np.transpose(np.vstack([x,y]))
            self.ptsLo = np.transpose(np.vstack([x,-y]))
        self._join_coordinates()

    def set_trailing_edge(self,zTEnew=0.0):
        zTEcurrent = self.ptsUp[-1,1] - self.ptsLo[-1,1]
        dzTE = zTEnew - zTEcurrent
        ptsUpNew = self.ptsUp[:,1] + self.ptsUp[:,0]*dzTE/2.
        ptsLoNew = self.ptsLo[:,1] - self.ptsLo[:,0]*dzTE/2.
        self.ptsUp[:,1] = ptsUpNew
        self.ptsLo[:,1] = ptsLoNew
        self._join_coordinates()
        self._analyze_geometry()

# --- debug section ---
def run_test_geometry():
    af = Airfoil()
    af = naca4(12,2,30)
    af.display('ko-')
    af.redim(40,True)
    print af.thickness
    af.display('ko-')

def run_test_aero_analysis():
    timer = Timer()
    af = cst_x(np.array([0.18723832, 0.2479892, 0.26252777, 0.31606257, 0.0819584, -0.11217863, -0.14363534, -0.06480575, -0.27817776, -0.08874038]))
    af._analyze_geometry()
    af.display()
    timer.lap('read db')    
    pol1 = af.get_jfoil_polar(0.2,3e6)
    timer.lap('javafoil')
    pol2 = af.get_xfoil_polar(0.2,3e6)
    timer.stop('xfoil')
    plt.figure(1)
    plt.hold(True)
    plt.grid(True)
    plt.plot(pol1.alpha,pol1.cd)
    plt.plot(pol2.alpha,pol2.cd)
    plt.show()

if __name__=="__main__":
    run_test_aero_analysis()