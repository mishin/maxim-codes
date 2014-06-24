# -*- coding: utf-8 -*-
"""
Created on Fri Apr 11 12:02:17 2014

@author: Maxim
"""
import numpy as np
from math import asin, log, log10
from scipy.optimize import root

from flight_conditions import FlightConditions
from paths import MyPaths
from db_tools import ReadDatabase
from drag_tools import DragItem, DragList, AircraftDrag

path = MyPaths()

# --- wrapper functions ---
def get_friction_drag_FW(aircraft,velocity,altitude):
    refArea = aircraft.wing.area
    frictionDrag = Friction(refArea)
    frictionDrag.set_flight_conditions(velocity,altitude)
    items = DragList()
    items.name = 'Airframe'
    items.add_item( frictionDrag.analyze_wing(aircraft.wing,'main wing') )
    return items
    
# --- drag data structures ---

class AnalysisFW:
    """
    Wrapper class for drag analysis containing friction drag and components 
    drag using table lookup. Aircraft main wing area is used as reference 
    area. 

    Parameters
    ----------
    
    aircraft : aircraft
        aircraft object
    velocity : float, m/sec
        velocity at which parasite drag will be calculated
    altitude : float, meter
        altitude at which parasite drag will be calculated
    updComponents : bool
        Should be true if components list was updated or at first call. 
        Calculation of components drag takes some time due to xls db 
        interface, so False option should be selected in most of the cases.


    Notes
    -----
    
    Only conventional configuration provided by aircraft class can be analyzed. 
    For unconventional configurations use **Friction** and **PartsDrag** 
    separately.


    Examples
    --------
    
    >>> import aircraft
    >>> ac = aircraft.load('sampleInput2')
    >>> V = 50.0 #m/sec
    >>> h = 2000. #m
    >>> drag = Analysis(ac,V,h,False)
    >>> aircraftDrag = drag.get_aircraft_drag()
    >>> aircraftDrag.display()
    ========================================
    Total drag components breakdown
    ========================================
    Name                      | Drag coef. |
    ----------------------------------------
    main wing                   1.0188e-02
    horizontal stab.            1.8388e-03
    ...
    antenna                     1.0526e-03
    ----------------------------------------
    TOTAL                       2.4408e-02
    ----------------------------------------
    """
    def __init__(self,aircraft,velocity,altitude,updComponents):
        self.ac = aircraft
        self.refArea = self.ac.wing.area
        self.aircraftDrag = AircraftDrag()
        self.get_friction_drag(velocity,altitude)
        if updComponents:
            self.get_components_drag()
        else:
            self.aircraftDrag.components = aircraft.drag.components
        self.aircraftDrag.update_total()
    def get_friction_drag(self,velocity,altitude):
        """
        Calculates friction drag of an aircraft using Friction class 
        at given velocity and altitude.
        
        Returns
        -------
        
        items : DragList
            list of the drag components with drag coefficients used for 
            friction drag calculation (Body, wing, empennage)
        """
        frictionDrag = Friction(self.refArea)
        frictionDrag.set_flight_conditions(velocity,altitude)
        items = DragList()
        items.add_item( frictionDrag.analyze_wing(self.ac.wing,'main wing') )
        self.aircraftDrag.friction = items
        return items
    def get_components_drag(self):
        """
        Calculates components drag using PartsDrag class methods using table 
        lookup methods.
        
        Returns
        -------
        
        items : DragList
            list of the drag components with drag coefficients obtained by 
            PartsDrag (langing gear, antenna etc.)
        """
        drag = PartsDrag(self.refArea)
        items = drag.get_items_drag(self.ac.drag.components)
        self.aircraftDrag.components = items
        return items
    def get_aircraft_drag(self):
        """
        Returns list of full drag components.
        """
        self.aircraftDrag.update_total()
        return self.aircraftDrag





# --- drag calculation section ---


class Friction(object):
    """
    calculates friction and form drag of the full aircraft
     or separate components. The code is based on Mason's friction fortran code.
    """
    def __init__(self,refArea):
        self.refArea = float(refArea)
        self.gamma   = 1.4
        self.Pr      = 0.72
        self.Te      = 390.0
        self.K       = 200.0
        self.TwTaw   = 1.0 # adiabatic wall condition
        self.wetArea      = list()
        self.refLength    = list()
        self.tc           = list()
        self.bodyType     = list()
        self.transitionPt = list()
    def set_reference_area(self,refArea=None):
        """
        reinitialize with new reference area
        """
        if refArea==None:
            self.__init__(self.refArea)
        else:
            self.__init__(refArea)
    def set_flight_conditions(self,velocity,altitude):
        """
        setting flight conditions to be analyzed. Flight conditions are calculated
         using FlightConditions module. Reynolds number is calculated using ref 
        length 1.0 that is stored in self.Re1
        """
        self.fc = FlightConditions(velocity,altitude)
        self.Re1 = self.fc.Re # Re for L=1

    def analyze_wing(self,wing,name):
        """
        analyzing friction and form drag of the wing-type bodies. Creates Component 
        object returns with all fields filled including geometry required for 
        drag analysis and drag value itself. 
        
        Parameters
        ----------
        
        wing : aircraft.wing
            wing, horizontal or vertical stabilizer object
        name : string
            name of the wing. Can be arbitrary and will be used for 
            report generation
        """
        self.set_reference_area()
        wingDrag = DragItem()
        wingDrag.name = name
        self.add_wing(wing)
        wingDrag.CD = self.analyze()
        wingDrag.quantity = 1
        wingDrag.CDtotal = wingDrag.CD
        return wingDrag

    def analyze_body(self,body,name):
        """
        analyzing revolution body object (fuselages with non-revolution body 
            shape are also analyzed)
        
        Parameters
        ----------
        
        body : aircraft.fuselage
            fuselage or nacelle object
        name : string
            name of analyzed body. Can be arbitrary and will be used for 
            report generation
    
        Returns
        -------
        
        bodyDrag : DragItem
            drag of the fuselage with reference area specified while initialization
        """
        self.set_reference_area()
        bodyDrag = DragItem()
        bodyDrag.name = name
        self.add_body(body)
        bodyDrag.CD = self.analyze()
        bodyDrag.quantity = 1
        bodyDrag.CDtotal = bodyDrag.CD
        return bodyDrag
    def add_wing(self,wing,transitionPoint=0.0):
        """
        adds wing type body into stack to be analyzed.
        """
        self.wetArea = np.hstack([self.wetArea,wing.wettedArea])
        for i, span in enumerate(wing.segSpans[1:]):
            self.bodyType.append(0)
            thickness1 = wing.secThickness[i]
            thickness2 = wing.secThickness[i+1]
            self.tc.append((thickness1+thickness2)/2.0) #average between two sections
            self.transitionPt.append(transitionPoint) #TODO update to get transition point from polar
            self.refLength.append(wing.segMAC[i])
    def add_body(self,fuselage):
        self.wetArea = np.hstack([self.wetArea,fuselage.wettedArea])
        self.bodyType.append(1)
        self.tc.append(fuselage.diameter / fuselage.length)
        self.transitionPt.append(0.0)
        self.refLength.append(fuselage.diameter)
    def analyze(self):
        n = len(self.wetArea)
        cfsw   = np.zeros(n)
        cfswff = np.zeros(n)
        for i,Swet in enumerate(self.wetArea):
            ff = self._get_form_factor(i)
            cf = self._get_transition_cf(i)
            cfsw[i]   = cf * self.wetArea[i]
            cfswff[i] = cfsw[i] * ff
        sum1 = np.sum(cfsw)
        sum2 = np.sum(cfswff)
        cdFriction = sum1 / self.refArea
        cdForm = (sum2 - sum1) / self.refArea
        cdTotal = cdFriction + cdForm
        return cdTotal

    def _get_transition_cf(self,i):
        Rel = self.Re1 * self.refLength[i]
        cfTurbL = self._get_turbulent_cf(Rel)
        if self.transitionPt[i]==0:
            cf = cfTurbL
        else:
            Rec = Rel * self.transitionPt[i]
            cfTurbC = self._get_turbulent_cf(Rec)
            cfLamC = self._get_laminar_cf(Rec)
            cf = cfTurbL - self.transitionPt[i]*(cfTurbC - cfLamC)
        return cf
    def _get_laminar_cf(self,Rex):
        r = self.Pr**0.5
        Mach = self.fc.Mach
        TwTe = self.TwTaw*(1.0+r*(self.gamma-1.0)*Mach**2/2)
        TstarTe = 0.5 + 0.039*Mach**2+0.5*TwTe
        KTe = self.K / self.Te
        Cstar = TstarTe**0.5*(1.0+KTe)/(TstarTe+KTe)
        cf = 2.0*0.664*Cstar**0.5/Rex**0.5
        return cf
    def _get_turbulent_cf(self,Rex):
        tol = 1.0e-8
        r = 0.88
        Te = 222.0
        Mach = self.fc.Mach
        m = (self.gamma - 1.0)*Mach*Mach/2
        TawTe = 1.0 + r*m
        F = self.TwTaw*TawTe
        Tw = F*Te
        A = (r*m/F)**0.5
        B = (1.0+r*m-F)/F
        denom = (4*A*A + B**2)**0.5
        alpha = (2.0*A*A-B)/denom
        beta = B / denom
        if Mach>0.1:
            Fc = r*m/((asin(alpha)+asin(beta))**2)
        else:
            Fc = ((1.0+F**0.5)/2)**2
        xden = 1.0 + 122.0**(-5/Tw)/Te
        xnum = 1.0 + 122.0**(-5/Tw)/Tw
        Ftheta = xnum/xden * (1/F)**0.5
        Fx = Ftheta / Fc
        RexBar = Fx * Rex
        Cfb = 0.074/(RexBar**0.2)
        itr = 0
        err = tol+1.0
        while err>tol and itr<100:
            Cf0 = Cfb
            xnum = 0.242 - Cfb**0.5*log10(RexBar*Cfb)
            xden = 0.121 + Cfb**0.5 / log(10.0)
            Cfb = Cfb*(1.0 + xnum/xden)
            err = abs(Cf0-Cfb)
            itr += 1
        return Cfb/Fc
    def _get_form_factor(self,i):
        if self.bodyType[i]==0:
            tc = self.tc[i]
            ff = 1.0 + 2.7*tc + 100*tc**4
        else:
            dl = self.tc[i]
            ff = 1.0 + 1.5*dl**1.5 + 50*dl**3 #TODO: change 7 to 50 to follow ref
        return ff


class PartsDrag:
    """
    drag of specific components of aircraft
    
    ref: "Light aircraft design" - Badyagin, Mukhamedov
    :ivar relCD: relative drag coefficient
    """
    def __init__(self,refArea):
        self.refArea = float(refArea)   # wing area
        self.dbPath = path.db.drag
        self.name  = list()
        self.desc  = list()
        self.relCD = list()
        self.mode  = list()
        self._read_db()
    
    def _read_db(self,sheetName='Sheet1'):
        db = ReadDatabase(self.dbPath,sheetName)
        nComp = db.nrows - 1
        self.name = db.read_column(1,0,nComp,True)
        self.desc = db.read_column(1,1,nComp,True)
        self.relCD = db.read_column(1,2,nComp,True)
        mode = db.read_column(1,3,nComp,True)
        for m in mode:
            if m=='*':
                self.mode.append(True)
            else:
                self.mode.append(False)        

    def get_item_drag(self,item):
        idx = self.name.index(item.name)
        if self.mode[idx]:
            item.CDtimesArea = self.relCD[idx]*item.frontArea
        else:
            item.CDtimesArea = self.relCD[idx]
        item.CD = item.CDtimesArea/self.refArea
        item.CDtotal = item.CD*item.quantity
        return item
    def get_items_drag(self,items):
        itemsNew = DragList()
        for item in items:
            item = self.get_item_drag(item)
            itemsNew.add_item(item)
        return itemsNew


def _calc_ARe(AR,Kp):
    v = lambda x: 2.0/(1.0/x+x**3.0) + 3.33/(2.0/(x**3.0)+1.0)
    f = lambda ARe: v(ARe) - Kp*v(AR)
    sol = root(f,AR)
    return sol.x[0]


class WaveDrag(object):
    def __init__(self,refArea):
        self.refArea = float(refArea)
    
    def analyze_wing(self,wing):
        #if self.fc.Mach>1.0:
        tc  = wing.equivThickness
        xtc = wing.equivThicknessLoc
        r0t = wing.equivLEradius/tc
        ht  = wing.equivCamber/tc
        TR = wing.taperRatio
        cosLE = np.cos(wing.equivSweepLErad)
        cosTE = np.cos(wing.equivSweepTErad)
        tanLE = np.tan(wing.equivSweepLErad)
        tanTE = np.tan(wing.equivSweepTErad)
        # airfoil factors
        Kb = 1.069 # for curved sections; Kb=1.0 for wedge type
        Kw = 1.0 # for curved sections; Kw=1.2 for wedge type
        # thickness factor
        r0t2 = r0t**0.5
        Kt = 1.0+4.0*(0.5-xtc*(1.0+0.5*r0t2))**2.0 - 0.25*r0t2*(1.0-xtc)**2.0
        Kc = 1.0+2.5*(ht)**2. # airfoil camber factor
        Kp = (cosLE + 0.5/((1.0+TR)**2.0)*(tanLE**2.0 - tanTE**2.0))
        Kp = Kp/(1.0 + 1.0/((1.0+TR*TR)**2.0)*(tanLE + tanTE)**2.0)
        Fb = 0.3+0.7*Kp**(1.0+2.0*TR**(1./3.))
        
        betaLim = abs(tanLE)
        M = self.fc.Mach*tanLE
        beta = (M**2.0-1.0)**0.5
        betaBar = beta/(tc**(1.0/3.0))
        if beta>betaLim:
            z = cosLE + cosTE
            m = 0.5*(1.0+TR**2.0*(2.0-TR)**3.0)*(betaLim/beta)**z
        else:
            m = 0.5*(1.0+TR**2.0*(2.0-TR)**3.0)
        ARe = _calc_ARe(wing.aspectRatio,Kp)
        # main equation
        Fbm = Fb**m
        var1 = Kt*Kw*Kc*Kb
        var2 = (1.0/ARe)/(1.0+(1.0+TR)*Fb*betaBar**(1.0+Kw))
        var3 = ARe**3.0/(1.0+1.0/3.0*ARe**3.0*betaBar**4.0)
        CDbar = 2.0*var1/(betaBar*Kb*Kw*Fbm + (var2+var3))
        var4 = betaBar*Kb*Kw**3.8*Fbm
        var5 = (2.0/ARe**3.0)/(1.0+(2.0/3.0+TR)*Fb*betaBar**(1.0+Kw**3.8))
        var6 = 1.0/(1.0+3.0*ARe*betaBar**4.0)
        CDbar += 3.33*var1/(var4+var5+var6)
        CDwave = CDbar*tc**(5./3)
        return CDwave
#        else:
#            return 0.0
    
    def set_flight_conditions(self,velocity,altitude):
        """
        setting flight conditions to be analyzed. Flight conditions are calculated
         using FlightConditions module. Reynolds number is calculated using ref 
        length 1.0 that is stored in self.Re1
        """
        self.fc = FlightConditions(velocity,altitude)
        self.Re1 = self.fc.Re # Re for L=1

if __name__=="__main__":
    import aircraft_FW
    import matplotlib.pyplot as plt
    ac = aircraft_FW.load('X45C')
    wd = WaveDrag(ac.wing.area)
    fd = Friction(ac.wing.area)
    M = np.arange(0.1,2.1,0.01)
    cd = np.zeros(len(M))
    for i,mach in enumerate(M):
        wd.set_flight_conditions(mach,10000)
        fd.set_flight_conditions(mach,10000)
        item = fd.analyze_wing(ac.wing,'main wing')
        cdf = item.CDtotal
        cdw = wd.analyze_wing(ac.wing)
        cd[i] = cdw + cdf
        print '%.4f\t%.6f\t%.6f'%(mach,cdw,cdf)
    plt.figure(1)
    plt.grid(True)
    plt.plot(M,cd,'b-')
    plt.xlabel('Mach')
    plt.ylabel('Friction+Form+Wave drag')
    plt.hold(True)
#    plt.axis([1,2,0,0.08])
    plt.show()