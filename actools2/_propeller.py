# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 17:13:09 2012
@author: Maxim

This module is part of propulsion module that contains propeller input/output 
functions, analysis using BEMT.

TODO: save propeller performance charts (Cp,Ct,effy vs. J), define effy(J=0)=0
create surrogate model via Rbf. Determine first point as [0,0] second at
efficiency approx 20%, to guarantee linear chart region of efficiency at low speed.
"""
import dbTools
import paths
import airfoil
import os
from numpy import array,hstack, zeros, copy, linspace
from scipy.interpolate import interp1d
from numpy import sin, cos, tan, radians, pi, exp, arccos
import FlightConditions
import matplotlib.pyplot as plt
from scipy.optimize import newton, fsolve
from scipy.integrate import simps

def load(name,dbPath='',afDbPath=''):
    prop = propeller()
    prop.read_xls(name,dbPath,afDbPath)
    return prop

def create_from_curves(name, hubDiameter, diameter, chord, pitch, thickness, baseAirfoil,xNew=[]):
    """
    Creates propeller from discrete set of parameters. Usually this set of parameters is obtained 
    by digitizing analog data. Chord, pitch and thickness values are obtained at *x* by cubic 
    spline interpolation. Base airfoil is scaled to fit desired thickness.
    
    Parameters
    ----------

    name : string
        propeller name
    hubDiameter : float, meters
        propeller hub diameter
    diameter : float, meters
        propeller diameter
    chord : float 2-n array, [nondim, meters]
        array of propeller chords and propeller stations in format 
        [[X(1),X(2)...X(n)],[C(1),C(2)...C(n)]]
    pitch : float 2-n array, [nondim, degree]
        array of propeller pitch angle and propeller stations
    thickness : float 2-n array, [nondim, nondim]
        array of propeller section thickness and stations
        [[X(1),X(2)...X(n)],[tc(1),tc(2)...tc(n)]]
    baseAirfoil : airfoil
        base airfoil object.
    x : float array, nondim
        propeller stations for new propeller initialization
    """
    
    if xNew==[]:
        xNew = [0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    
    chordCurve     = interp1d(chord[:,0],chord[:,1],'spline')
    pitchCurve     = interp1d(pitch[:,0],pitch[:,1],'spline')
    thicknessCurve = interp1d(thickness[:,0],thickness[:,1],'spline')
    
    chordNew     = chordCurve(xNew)
    pitchNew     = pitchCurve(xNew)
    thicknessNew = thicknessCurve(xNew)
    
    airfoils = list()
    for tc in thicknessNew:
        afNew = baseAirfoil
        afNew.scale_thickness(thicknessNew)
        airfoils.append(afNew)
    prop = propeller()
    prop.name = name
    prop.diameter = diameter
    prop.hubDiameter = hubDiameter
    prop.chord = chordNew
    prop.beta = pitchNew
    prop.airfoil = afNew
    return prop

class Propeller_new:
    def __init__(self):
        self.name = 'Propeller'
        self.diameter = 0.0
        self.numBlades = 0
        self.radius = 0.0
        self.hubDiam = 0.0
        self.chord = list()
        self.thickness = list()
        self.beta = list()
        self.x = list()
        self.r = list()
        self.airfoil = list()
        self.betaSet = 0.0
        self.defaultPath = paths.Database()

    def create_by_params(self,name,diam,numBlade,hubDiam,chord,beta,thickness,x,baseAirfoil,betaSet):
        self.name = name
        self.diameter  = diam
        self.numBlades = numBlade
        self.hubDiam   = hubDiam
        self.chord     = chord
        self.beta      = beta
        self.thickness = thickness
        self.x         = x
        self.betaSet   = betaSet
        for tc in self.thickness:
            af = baseAirfoil
            af.scale_thickness(tc)
            self.airfoil.append(af)
        self._get_prop_data()

    def create_by_params_unstruct(self,name,diam,numBlade,hubDiam,chord,beta,thickness,xNew,baseAirfoil,betaSet):
        chord = self._interp1d(chord,xNew)
        beta = self._interp1d(beta,xNew)
        thickness = self._interp1d(thickness,xNew)
        self.create_by_params(name,diam,numBlade,hubDiam,chord,beta,thickness,xNew,baseAirfoil,betaSet)

    def create_with_airfoil(self,name,diam,numBlade,hubDiam,chord,beta,airfoil,x,betaSet):
        self.name = name
        self.diameter = diam
        self.numBlades = numBlade
        self.hubDiam = hubDiam
        self.chord = chord
        self.beta = beta
        self.x = x
        self.airfoil = airfoil
        self.betaSet = betaSet
        self._get_prop_data()

    def create_with_airfoil_name(self,name,diam,numBlade,hubDiam,chord,beta,airfoilName,x,betaSet):
        afList = list()
        for afName in airfoilName:
            af = airfoil.load(afName)
            afList.append(af)
        self.create_with_airfoil(self,name,diam,numBlade,hubDiam,chord,beta,afList,x,betaSet)

    def _interp1d(self,data,xNew):
        spline = interp1d(data[:,0],data[:,1],'spline')
        yNew = spline(xNew)
        return yNew

    def _get_prop_data(self):
        self.radius = self.diameter / 2.0
        self.hubRad = self.hubDiam / 2.0
        self.r = self.x * self.radius
        bladeArea = simps(self.chord,self.r)
        self.solidity = bladeArea * self.numBlades / (pi*self.radius**2)
        self.thickness = zeros(len(self.airfoil))
        for i in range(len(self.airfoil)):
            self.thickness[i] = self.airfoil[i].thickness

    def read_xls(self,name,dbPath='',afDbPath=''):
        if dbPath=='':
            dbPath = self.defaultPath.prop
        if afDbPath=='':
            afDbPath = self.defaultPath.propAirfoil
        dbPath = paths.fixPaths(dbPath)
        propDb = dbTools.loadDB(dbPath)
        sheet = propDb.selectByName(name)
        db = dbTools.readDB(sheet)
        self.name = name
        self.diameter = db.readRow(1,1)
        self.hubDiam  = db.readRow(2,1)
        self.numBlades= db.readRow(3,1)
        self.betaRange= db.readRow(4,1)
        self.r        = db.readRow(5,1)
        self.x = 2.0*self.r / self.diameter
        self.chord = db.readRow(7,1)
        self.beta = db.readRow(8,1)
        afname = db.readRow(9,1)
        self.airfoil = list()
        for afn in afname:
            af = airfoil.load(afn)
            self.airfoil.append(af)
        self._get_prop_data()

class propeller:
    """
    Attributes
    ----------
    
    name : string
    diameter : float, meters
    hubDiameter : float, meters
    numBlades : integer
    r : array, meters
        propeler station (distance from propeller center to section)
    x : array
        relative propeller station x=r/R
    chord : array, meters
        array of propeller section chords
    beta : array, degrees
        propeller section pitch angle
    airfoil : list of airfoils
        array of airfoil objects representing given section
    """
    def __init__(self):
        self.name        = 'propeller'
        self.diameter    = 0.0
        self.hubDiameter = 0.0
        self.numBlades   = 0
        self.numStations = 0
        self.r           = list()
        self.x           = list()
        self.chord       = list()
        self.beta        = list()
        self.airfoilName = list()
        self.airfoil     = list()
        self.defaultPath = paths.Database()
    
    def read_xls(self,propName, dbPath='', afDbPath=''):
        """
        reads in excel propeller database. Two database files are needed: propeller, airfoil
        
        Parameters
        ----------
        
        propName : string
            propeller name
        dbPath : path
            path of the propeller database. Default db will be used if not 
            specified.
        afDbPath : path
            path of the airfoil database. Default db will be used if not 
            specified.
        """
        if dbPath=='':
            dbPath = self.defaultPath.prop
        if afDbPath=='':
            afDbPath = self.defaultPath.propAirfoil
        dbPath = paths.fixPaths(dbPath)
        propDb = dbTools.loadDB(dbPath)
        sheet = propDb.selectByName(propName)
        db = dbTools.readDB(sheet)
        self.name        = propName
        self.diameter    = db.readRow(1,1)
        self.radius      = self.diameter / 2.0
        self.hubDiameter = db.readRow(2,1)
        self.numBlades   = db.readRow(3,1)
        self.betaRange   = db.readRow(4,1)
        self.r           = db.readRow(5,1)
        self.x           = 2*self.r / self.diameter
        self.chord       = db.readRow(7,1)
        self.beta        = db.readRow(8,1)
        self.airfoilName = db.readRow(9,1)
        self.numStations = len(self.r)
        self.betaCurrent = 0.0
        self.betaXstation = 0.0
        for afName in self.airfoilName:
            af = airfoil.Airfoil()
            af.read_xls(afName,afDbPath)
            self.airfoil.append(af)
        self.analyze_geometry()

    def write_xls(self,dbPath='',afDbPath=''):
        if dbPath=='':
            dbPath = self.defaultPath.prop
        if afDbPath=='':
            afDbPath = self.defaultPath.propAirfoil
        dbPath   = paths.fixPaths(dbPath)
        print dbPath, afDbPath
        db       = dbTools.loadDB(dbPath,mode='w')
        newSheet = db.add_sheet(self.name)
        sheet    = dbTools.writeDB(newSheet)
        sheet.writeRow('name',self.name)
        sheet.writeRow('diameter',self.diameter)
        sheet.writeRow('hub diameter',self.hubDiameter)
        sheet.writeRow('num blades',self.numBlades)
        sheet.writeRow('beta range deg',self.betaRange)
        sheet.writeRow('r',self.r)
        sheet.writeRow('x',self.x)
        sheet.writeRow('chord',self.chord)
        sheet.writeRow('beta',self.beta)
        sheet.writeRow('airfoils',self.airfoilName)
        db.save_db()
        for af in self.airfoil:
            af.write_xls(afDbPath)
    
    def create_splines(self):
        self.betaAtX = interp1d(self.x,self.beta,'cubic')
        self.chordAtX = interp1d(self.x,self.chord,'cubic')
        
    def analyze_geometry(self):
        self.create_splines()
        self.betaX75 = 0.75
        self.beta75 = self.betaAtX(self.betaX75)
        self._get_solidity_ratio()
#        self.thickness = list()
#        for af in self.airfoil:
#            af.analyze_geometry()
#            self.thickness.append(af.thickness)
#    
    def set_beta(self,betaNew,x=0.75):
        #TODO: if x is same set beta without using interpolation betaAtX
        if x==self.betaX75:
            betaIncrement = betaNew - self.beta75
        else:
            betaCurrent = self.betaAtX(x)
            betaIncrement = betaNew - betaCurrent
            self.betaX75 = x
        self.beta += betaIncrement
        self.beta75 = betaNew
        self._get_solidity_ratio()

    def create_prop(self,propName,diameter,hubDiameter,numBlades,x,chord,beta,airfoil,betaRange):
        """
        creates propeller using given input data.
        This method can be used to initialize propeller using airfoil text files
        """
        self.name = propName
        self.diameter = diameter
        self.hubDiameter = hubDiameter
        self.numBlades = numBlades
        self.x = x
        self.r = diameter * x / 2
        self.chord = chord
        self.beta = beta
        self.airfoil = airfoil
        self.betaRange = betaRange
        for af in airfoil:
            self.airfoilName.append(af.name)
        self.analyze_geometry()

    def _get_solidity_ratio(self):
        x = hstack([0.0,self.x,1.0])
        c = hstack([0.0,self.chord,0.0])
        area = simps(c,x)
        self.sigma = 4.0*area*self.numBlades / (self.diameter**2*pi)

    def analyze_prop(self,betaSet,N,V,rho):
        """
        Parameters
        ----------
        
        betaSet : float, degree
            pitch angle setting at 75% of diameter
        N : float, rpm
            propeller rotation speed
        V : float, m/sec
            true airspeed of aircraft (propeller)
        rho : float, kg/m^3
            air density
        
        Detailed theory explanation will be included here
        
        Performs propeller analysis using Blade-element momentum theory.
        Advance ratio is calculated as 
        
        .. math::
            J = V/nd
        
        Master function to calculate advance ration using combined blade element 
        momentum theory is
        
        .. math::
            J = \frac{\pi x (4F sin^2{\phi}-\lambda_T \sigma)}{4F sin{\phi} cos{\phi}+\lambda_T \sigma}
        
        In case of poor convergence (required J cannot be achieved at given rpm and
        airspeed) function returns advance ratio and CP=None, CT=None, effy=None.
        
        Returns
        -------
        J : float
            advance ratio
        CT : float
            thrust coefficient
        CP : float
            power coefficient
        effy : float
            propeller efficiency
        """
        if N*V*rho==0.0:
            return zeros(6)
        if betaSet!=None:
            self.set_beta(betaSet,0.75)
        Jcorr = 0.97
        n  = N/60.0
        d  = self.diameter
        B  = self.numBlades
        dT = list()
        dQ = list()
        altitude   = FlightConditions.get_density_altitude(rho)
        soundSpeed = FlightConditions.ISAtmosphere(altitude).soundSpeed
        Jreq       = Jcorr * V / (n*d)
        def master_func(alpha,beta,x,r,B,soundSpeed,d,af,c,n):
            phi = radians(beta - alpha)
            # compute Prandtl loss factor
            sinPhiL1 = (1.0 + 1.0/(x*tan(phi))**2)**0.5
            f = B * (1.0/x-1.0) * sinPhiL1/2.0
            F = 2.0 * arccos(exp(-f)) / pi
            vLocal = pi*x*n*d / cos(phi)
            mach = vLocal / soundSpeed
            cl = float(af.polar.clAlpha(mach,alpha))
            cd = float(af.polar.cdAlpha(mach,alpha))
            sinPhi = sin(phi)
            cosPhi = cos(phi)
            lambdaT = cl*cosPhi - cd*sinPhi
            lambdaP = cl*sinPhi + cd*cosPhi
            sigma   = c*B / (2.0*pi*r)
            J = pi*x*(4.0*F*sinPhi**2 - lambdaT*sigma) / (4.0*F*sinPhi*cosPhi + lambdaP*sigma)
            return J,vLocal, lambdaT,lambdaP
        def get_J(alpha,Jreq,beta,x,r,B,soundSpeed,d,af,c,n):
            rslt = master_func(alpha,beta,x,r,B,soundSpeed,d,af,c,n)
            return rslt[0] - Jreq
        for i in range(self.numStations):
            x      = self.x[i]
            c      = self.chord[i]
            beta   = self.beta[i]
            r      = self.r[i]
            alpha0List = [-20,-10,-5,0,5,10,20,30,40,50]
            af     = self.airfoil[i]
            for alpha0 in alpha0List:
                try:
                    alpha  = newton(get_J,alpha0,args=(Jreq,beta,x,r,B,soundSpeed,d,af,c,n),tol=0.01)
                    break
                except RuntimeError:
                    return 0.0,0.0,Jreq, 0.0, 0.0, 0.0
            J,vLocal,lambdaT,lambdaP = master_func(alpha,beta,x,r,B,soundSpeed,d,af,c,n)
            dT.append(  rho*vLocal**2*lambdaT*c*B / 2.0  )
            dQ.append(  rho*vLocal**2*lambdaP*c*B*r / 2.0  )
        dT.insert(0,0.0)
        dQ.insert(0,0.0)
        dT.append(0.0)
        dQ.append(0.0)
        r_ = copy(self.r)
        r_ = hstack([0.0,r_,1.0])
        T  = simps(dT,r_)
        Q  = simps(dQ,r_)
        CT   = T / (rho*n**2*d**4)
        CQ   = Q / (rho*n**2*d**5)
        effy = CT/CQ*Jreq/(2.0*pi)
        effy = effy
        CP   = Jreq*CT/effy
        P    = CP*rho*n**3*d**5
        if CT<=0.0 or CT>1.0:
            CT=0.0
            effy = 0.0
        if CQ<=0.0 or CQ>1.0:
            CQ=0.0
            effy = 0.0
        if CP<=0.0:
            CP=0.0
        if not effy==None:
            effy = effy*0.95
        return T,P,Jreq/Jcorr,CT,CP,effy

    def analyze_at_velocity_range(self,betaSet,V,N,rho):
        T = list()
        P = list()
        J    = list()
        CT   = list()
        CP   = list()
        effy = list()
        for v in V:
            data = self.analyze_prop(betaSet,N,v,rho)
            T.append(data[0])
            P.append(data[1])
            J.append(data[2])
            CT.append(data[3])
            CP.append(data[4])
            effy.append(data[5])
        return T,P,J,CT,CP,effy

    def test_prop(self):
        N = 800.0
        V = linspace(5.0,51.0,50)
        rho = 1.2255
        beta = 45.0
        J    = list()
        CT   = list()
        CP   = list()
        effy = list()
        for v in V:
            data = self.analyze_prop(beta,N,v,rho)
            J.append(data[2])
            CT.append(data[3])
            CP.append(data[4])
            effy.append(data[5])
        plt.figure(1)
        plt.grid(True)
        plt.hold(True)
        plt.plot(J,CT,'o-')
        plt.xlabel('J')
        plt.plot(J,CP,'ro-')
        plt.axis([0,2.5,0,0.15])
        plt.figure(2)
        plt.plot(J,effy,'gs-')
        plt.hold(True)
        plt.grid(True)
        plt.axis([0,2.5,0,1.0])
        plt.xlabel('advance ratio')
        plt.ylabel('efficiency')
        plt.show()
    
    def plot_propeller(self):
        cD = self.chord / self.diameter
        beta = self.beta
        x = self.x
        
        
        fig = plt.figure()
        ax1 = fig.add_subplot(311)
        ax1.grid(True)
        ax1.plot(self.x,self.chord,'bo-')
        ax1.set_ylabel('chord')
        ax2 = fig.add_subplot(312)
        ax2.grid(True)
        ax2.plot(self.x,self.thickness,'ro-')
        ax2.set_ylabel('thickness')
        ax3 = fig.add_subplot(313)
        ax3.grid(True)
        ax3.plot(self.x,self.beta,'go-')
        ax3.set_ylabel('pitch angle')
        plt.show()
#        fig = plt.figure(11)
#        ax1 = fig.add_subplot(111)
#        ax1.hold(True)
#        ax1.grid(True)
#        lns1 = ax1.plot(x,cD,'ks-',label='chord/diameter')
#        lns2 = ax1.plot(x,tc,'bo-',label='thickness')
#        ax1.set_xlabel('x = r/R')
#        ax2 = ax1.twinx()
#        lns3 = ax2.plot(x,beta,'dr-',label='pitch angle')
#        lns = lns1 + lns2 + lns3
#        labs = [l.get_label() for l in lns]
#        ax1.legend(lns, labs, loc=0)
#        ax2.set_ylabel('pitch angle, deg')
#        plt.title('%s propeller geometry'%self.name)
#        ax1.axis([0,1.,0,0.3])
#        ax2.axis([0,1.,0,50])
#        plt.show()

def full_test():
    from miscTools import Timer
    timer = Timer()
    name = 'neuform'
    rho = 1.2255
#    beta = ny.array([15.0,20,25,30,35,40,45])
    #beta = array([15.,30,45])
    beta = linspace(15,45,10)
    N = array([1200.,1200,1000,1000,800,900,800,800,800,800])
    V = linspace(5.0,75.0,100)

    J = zeros([len(beta),len(V)])
    effy = zeros([len(beta),len(V)])
    prop = propeller()
    prop.read_xls(name)
    
    legend = list()
    color = ['b','r','m','y','k','c','b','b','r','m','y','k','c','b']
    for i,b in enumerate(beta):
        legend.append('%.0f deg'%b)
        for j,v in enumerate(V):
            data = prop.analyze_prop(b,N[i],v,rho)
            J[i,j] = data[2]
            effy[i,j] = data[5]

    Jexp, effyExp = read_exp_data()

    plt.figure(13)
    plt.xlabel('J/nD - advance ratio')
    plt.ylabel('efficiency')
    plt.axis([0,2.8,0,1.0])
    plt.title('%s'%name)
    plt.hold(True)
    plt.grid(True)
    for i in range(len(beta)):
        plt.plot(J[i,:],effy[i,:],color[i]+'-')
        #plt.plot(Jexp[i],effyExp[i],color[i]+'o-')
        
    plt.show()

def read_exp_data():
    wbName = r'D:\propeller_results.xls'
    wb = dbTools.loadDB(wbName)
    sh = wb.selectByName('Sheet1')
    sh = dbTools.readDB(sh)
    
    J = list()
    J.append(sh.readCol(3,0,13))
#    J.append(sh.readCol(3,2,17))
#    J.append(sh.readCol(3,4,18))
    J.append(sh.readCol(3,6,21))
#    J.append(sh.readCol(3,8,17))
#    J.append(sh.readCol(3,10,16))
    J.append(sh.readCol(3,12,17))
    
    effy = list()
    effy.append(sh.readCol(3,1,13))
#    effy.append(sh.readCol(3,3,17))
#    effy.append(sh.readCol(3,5,18))
    effy.append(sh.readCol(3,7,21))
#    effy.append(sh.readCol(3,9,17))
#    effy.append(sh.readCol(3,11,16))
    effy.append(sh.readCol(3,13,17))
    
    return J, effy

def create_prop_from_curves():
    from math import atan, pi, degrees
    from numpy import zeros
    from scipy.interpolate import interp1d
    
    x = array([0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95])
    folder = r'D:\light aircraft\V05\propeller analysis'
    chordsFile = folder + '\\chords.txt'
    pitchFile = folder + '\\pitch.txt'
    thicknessFile = folder + '\\thickness.txt'
    
    diameter = 1.73
    chords = dbTools.read_txt_data(chordsFile)
    pitch = dbTools.read_txt_data(pitchFile)
    thickness  = dbTools.read_txt_data(thicknessFile)
    pD = pitch[:,1]
    rR = pitch[:,0]
    phi = zeros(len(pD))
    for i,pd in enumerate(pD):
        phi[i] = atan(pd/(rR[i]*pi))
        phi[i] = degrees(phi[i])

    pitchCurve = interp1d(rR,phi,'cubic')
    chordsCurve = interp1d(chords[:,0],chords[:,1],'cubic')
    thicknessCurve = interp1d(thickness[:,0],thickness[:,1],'cubic')
    
    phiNew = pitchCurve(x)
    tcNew = thicknessCurve(x)
    tcNew[0] = 1.4*tcNew[4]
    tcNew[1] = 1.3*tcNew[4]
    tcNew[2] = 1.2*tcNew[4]
    tcNew[3] = 1.1*tcNew[4]
    chordNew = chordsCurve(x) * diameter

    afName = 'Clark-Y'
    propName = 'NACA5868'
    airfoils = list()
    hubDiameter = 0.18
    numBlades = 2
    for i in range(len(x)):
        af = airfoil.Airfoil()
        af.read_xls(afName)
        af.scale_thickness(tcNew[i])
        af.build_aero_table(alphaSeq=[-90,90,2.0])
        airfoils.append(af)
    
    prop = propeller()
    prop.create_prop(propName,diameter,hubDiameter,numBlades,x,chordNew,phiNew,airfoils,[15.,45.])
    prop.write_xls()

def run_test2():
    propName = 'neuform'
    diameter = 1.736
    hubDiameter = 0.5
    numBlades = 3
    x = array([0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85])
    chord = array([0.06751,0.08984,0.08925,0.08584,0.08017,0.07055,0.05802,0.02807])
    beta = array([36.0,32.29,28.87,25.341,22.4347,19.7793,17.8562,16.2372])
    n = len(x)
    airfoils = list()
    for i in range(n):
        afName = os.getcwd() + '\\prop_airfoils\\prop_%d.txt'%(i+1)
        af = airfoil.Airfoil()
        af.read_txt(afName)
        af.build_aero_table()
        airfoils.append(af)
    prop = propeller()
    prop.create_prop(propName,diameter,hubDiameter,numBlades,x,chord,beta,airfoils)
    prop.write_xls()

def run_test1():
    prop = propeller()
    prop.read_xls('neuform')
    print prop.name
    print prop.diameter
    print prop.r
    print prop.x
    print prop.airfoilName
    prop.full_analysis()

def run_test3():
    prop = propeller()
    prop.read_xls('NACA5868')
    rpm = 1200.0
    V = linspace(1,30.,50)
    rho = 1.2255
    data = prop.analyze_at_velocity_range(15.0,V,rpm,rho)
    T = data[1]
    J = data[2]
    plt.figure(1)
    plt.plot(J,T,'bo-')
    plt.xlim([0.0,1.5])
    plt.show()

def run_test04():
    prop = propeller()
    prop2 = propeller()
    prop2.read_xls('MDO_optimum')
    prop.read_xls('MDO_initial')
    prop.analyze_geometry()
    prop2.analyze_geometry()

    
    V = linspace(1.0,100.0,100)
    N = 1500.0
    rho = 1.2255
    #beta = array([10.0,15.0,20.0,25.0])
    beta = array([10.0,20.0,30.0])
    
    plt.figure()
    plt.grid(True)
    plt.hold(True)
    for b in beta:
        data = prop.analyze_at_velocity_range(b,V,N,rho)
        data2 = prop2.analyze_at_velocity_range(b,V,N,rho)
        plt.plot(data[2],data[-1],'b')
        plt.plot(data2[2],data2[-1],'r')
        #plt.plot(data3[2],data3[-1],'g')
    plt.axis([0.0,1.5,0.0,1.0])
    plt.xlabel('advance ratio')
    plt.ylabel('efficiency')
    plt.legend(['baseline','optimum'],'upper left')
#return T,P,J,CT,CP,effy
    prop.set_beta(19.0)
    prop2.set_beta(28.5)

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.grid(True)
    ax1.hold(True)
    ax1.plot(prop.x,prop.chord,'b')
    ax1.plot(prop2.x,prop2.chord,'r-')
    ax1.axis('equal')
    ax1.set_ylabel('chord')
    ax2 = fig.add_subplot(312)
    ax2.grid(True)
    ax2.hold(True)
    ax2.plot(prop.x,prop.beta,'b')
    ax2.plot(prop2.x,prop2.beta,'r')
    ax2.legend(['baseline','optimum'])
    ax2.set_ylabel('pitch angle,deg')
    ax3 = fig.add_subplot(313)
    ax3.grid(True)
    ax3.hold(True)
#    ax3.plot(prop.x,prop.thickness,'b')
#    ax3.plot(prop2.x,prop2.thickness,'r')
    ax3.set_xlabel('x=r/R')
    plt.show()
    
#T,P,Jreq/Jcorr,CT,CP,effy
if __name__=="__main__":
    run_test04()