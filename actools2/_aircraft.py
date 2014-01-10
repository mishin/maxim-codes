# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 15:06:38 2012
Modules for loading and processing aircraft data files and storing all 
data for a given aircraft configuration
@author: Daniel
"""
import numpy
import dbTools
import engine as eng
import propeller as prop
import airfoil as Airfoil_Class
import paths
import miscTools
import drag
import mass as weight
import aerodynamics
import AircraftAnalysis
from warnings import warn
from copy import copy
from math import tan, radians, sin, cos
import FlightConditions as FC
from matplotlib.mlab import find
import matplotlib.pyplot as plt
#import AircraftAnalysis
import geometry as geom
from scipy import interp
from scipy.integrate import simps, trapz

def load(name):
    """
    function to load aircraft from xls database.

    Parameters
    ----------

    name : string
        name of the aircraft configuration (sheet name in db)

    Examples
    --------

    load aircraft configuration from db and store it into *ac* variable

    >>> import aircraft
    >>> ac = aircraft.load('sampleInput1')
    >>> ac.display()
    >>> print ac.get_mass_empty() # print empty mass
    >>> print ac.get_mass_total() # print total mass
    >>> print ac.get_drag() # print total drag
    """
    pth=paths.Database()
    acDB =dbTools.loadDB(pth.aircraft)
    return aircraft(acDB.selectByName(name))


class aircraft(object):
    """
    Main class describing aircraft configuration.
    """
    def __init__(self,aircraftInputSheet):
        self._inputSheet    =aircraftInputSheet
        self.entryName      =str(self._inputSheet.name)
        self.name           =self.entryName#str(self._inputSheet.cell(0,1).value)
        self._readData()
        self._get_volume_ratios()
        # components mass and drag are calculated here
        self.mass = weight.get_general_aviation_mass(self)
        self.drag = drag.get_total_drag(self,self.designGoals.designSpeed,
                                        self.designGoals.designAltitude,True)
        # TODO: separate analysis into different methods
        self.analysis       =AircraftAnalysis.analysis(self)
    def display(self):
        """
        Displays 3D aircraft with all mass comonents.
        """
        #TODO: mesh to be displayed should be in main class (wing, fus)
        import plot_aircraft
        plot_aircraft.plot_aircraft2(self)
    def _readData(self):
        db=dbTools.readDB(self._inputSheet)
        # Reading in Design Quantities Section
        i=db.findHeader("SECTION: DESIGN QUANTITIES")
        self.type                  =db.readRow(i+1,1)
        fuelMass                   =db.readRow(i+2,1)
        designGrossMass            =db.readRow(i+3,1)
        designSpeed                =db.readRow(i+4,1)
        designAltitude             =db.readRow(i+5,1)
        designLoadFactor           =db.readRow(i+6,1)
        designLandingLoadFactor    =db.readRow(i+7,1)
        occupantNumber             =db.readRow(i+8,1)
        payloadMass                =db.readRow(i+9,1)
        self.designGoals=designGoals(fuelMass,designGrossMass,designSpeed,designAltitude,designLoadFactor,designLandingLoadFactor,occupantNumber,payloadMass)
        # Reading in Main Wing Section        
        i=db.findHeader("SECTION: MAIN WING")
        segmentSpans               =db.readRow(i+1,1,iterable=True)
        chords                     =db.readRow(i+2,1,iterable=True)
        airfoilName                =db.readRow(i+3,1,iterable=True)
        segmentOffsets             =db.readRow(i+4,1,iterable=True)
        segmentDihedrals           =db.readRow(i+5,1,iterable=True)
        twistAngles                =numpy.concatenate((numpy.zeros(1),db.readRow(i+6,1,iterable=True)))
        incidence                  =db.readRow(i+7,1)
        aapex                      =db.readRow(i+8,1)
        aileronLocation            =db.readRow(i+9,1)
        flapLocation               =db.readRow(i+10,1)
        flapType                   =db.readRow(i+11,1)
        fuelWingCG                 =db.readRow(i+12,1)
        material                   =db.readRow(i+13,1)
        self.wing=wing(segmentSpans,chords,airfoilName,segmentOffsets,
                       segmentDihedrals,twistAngles,incidence,aapex,material,1)
        self.wing.setFuelTank(fuelMass/2.0,1,fuelWingCG[0],fuelWingCG[1])
        self.wing.addAileron(aileronLocation)
        self.wing.addFlap(flapLocation,flapType)
        # Reading in Horizontal Stab Section
        i=db.findHeader("SECTION: HORIZONTAL STABILIZER")
        segmentLengths             =db.readRow(i+1,1,iterable=True)
        chords                     =db.readRow(i+2,1,iterable=True)
        airfoilName                =db.readRow(i+3,1,iterable=True)
        segmentOffsets             =db.readRow(i+4,1,iterable=True)
        segmentDihedrals           =db.readRow(i+5,1,iterable=True)
        twistAngles                =numpy.zeros(len(chords))
        incidence                  =db.readRow(i+6,1)
        aapex                      =db.readRow(i+7,1)
        elevatorLocation           =db.readRow(i+8,1)
        material                   =db.readRow(i+8,1)
        #self.hStab=hStab(segmentLengths,chords,airfoilName,segmentOffsets,segmentDihedrals,incidence,aapex,elevatorLocation,material)
        self.hStab=wing(segmentLengths,chords,airfoilName,segmentOffsets,segmentDihedrals,twistAngles,incidence,aapex,material,1,False)
        self.hStab.addElevator(elevatorLocation)
        # Reading in Vertical Stab Section
        i=db.findHeader("SECTION: VERTICAL STABILIZER")
        number                     =db.readRow(i+1,1)
        lateralDistance            =db.readRow(i+2,1)
        segmentLengths             =db.readRow(i+3,1,iterable=True)
        chords                     =db.readRow(i+4,1,iterable=True)
        airfoilName                =db.readRow(i+5,1,iterable=True)
        segmentOffsets             =db.readRow(i+6,1,iterable=True)
        segmentDihedrals           =db.readRow(i+7,1,iterable=True)
        twistAngles                =numpy.zeros(len(chords))
        aapex                      =db.readRow(i+8,1)
        rudderLocation             =db.readRow(i+9,1)
        material                   =db.readRow(i+10,1)
        self.vStab=wing(segmentLengths,chords,airfoilName,segmentOffsets,segmentDihedrals,twistAngles,incidence,aapex,material,0,True)
        self.vStab.addRudder(rudderLocation)
        self.vStabNum              =number
        self.vStabSeparation       =lateralDistance 
        # Reading in Fuselage Data
        i=db.findHeader("SECTION: FUSELAGE")
        length                     =db.readRow(i+1,1)
        diameter                   =db.readRow(i+2,1)
        wettedArea                 =db.readRow(i+3,1)
        sideProfile_X              =db.readRow(i+4,1)
        sideProfile_Y              =db.readRow(i+5,1)
        topProfile_X               =db.readRow(i+6,1)
        topProfile_Y               =db.readRow(i+7,1)
        material                   =db.readRow(i+8,1)        
        self.fuselage=fuselage(length,diameter,wettedArea,sideProfile_X,sideProfile_Y,topProfile_X,topProfile_Y,material)
        # Reading in Engine Data
        i=db.findHeader("SECTION: PROPULSION")
        name                       =db.readRow(i+1,1)
        CG_X                       =db.readRowAsArray(i+2,1)
        CG_Y                       =db.readRowAsArray(i+3,1)
        CG_Z                       =db.readRowAsArray(i+4,1)
        propName                   =db.readRow(i+5,1)
        propHub_X                  =db.readRowAsArray(i+6,1)
        propHub_Y                  =db.readRowAsArray(i+7,1)
        propHub_Z                  =db.readRowAsArray(i+8,1)
        minBeta                    =db.readRowAsArray(i+9,1)
        maxBeta                    =db.readRowAsArray(i+10,1)
        self.engine=engine(name,CG_X,CG_Y,CG_Z)
        self.propeller=propeller(propName,propHub_X,propHub_Y,propHub_Z,minBeta,maxBeta)
        # Reading in Landing Gear Data
        i=db.findHeader("SECTION: LANDING GEAR")
        gearType                   =db.readRow(i+1,1)
        tireWidth                  =db.readRow(i+2,1)
        tireDiameter               =db.readRow(i+3,1)
        strutLength                =db.readRow(i+4,1)
        strutType                  =db.readRow(i+5,1)
        fairingType                =db.readRow(i+6,1)
        groundContact_X            =db.readRow(i+7,1)
        groundContact_Y            =db.readRow(i+8,1)
        groundContact_Z            =db.readRow(i+9,1)
        self.landingGear=landingGear(gearType,tireWidth,tireDiameter,strutLength,strutType,fairingType,groundContact_X,groundContact_Y,groundContact_Z)
        # Reading in VLM Operating Parameters        
        i=db.findHeader("SECTION: OPTIONAL VLM PARAMS")
        defaultFlag                =int(db.readRow(i+1,1))
        if defaultFlag==1:
            wing_spanwise          =db.readRow(i+2,1)
            wing_chordwise         =db.readRow(i+3,1)
            hStab_spanwise         =db.readRow(i+4,1)
            hStab_chordwise        =db.readRow(i+5,1)
            vStab_spanwise         =db.readRow(i+6,1)
            vStab_chordwise        =db.readRow(i+7,1)
            cosDistribution        =db.readRow(i+8,1)
            sinDistribution        =db.readRow(i+9,1)
            self.VLMsetup=VLMsetup(wing_spanwise,wing_chordwise,hStab_spanwise,hStab_chordwise,vStab_spanwise,vStab_chordwise,cosDistribution,sinDistribution)
        else:
            self.VLMsetup=VLMsetup()
        # Reading in Mass List
        massSec=db.read_section('MASS SOURCES',0)
        self.mass = weight.AircraftMass(self.name)
        for line in massSec:
            name=str(line[0])
            values  =numpy.array(line[1:])
            mass    =values[0]
            CG      =values[1:4]
            MOM     =values[5:8]
            self.mass.payload.add_item(name,mass,CG,MOM)
        self.mass.fuel.add_item('Fuel tank Right',fuelMass/2.0)
        self.mass.fuel.add_item('Fuel tank Left',fuelMass/2.0)
        self.mass.update_total()
        # Reading in Drag List
        dragSec=db.read_section('DRAG SOURCES',0)
        self.drag = drag.AircraftDrag()
#        self.drag=drag.DragList()
        for line in dragSec:
            drg=drag.DragItem()
            drg.name=str(line[0])
            drg.quantity=int(line[1])
            drg.frontArea=float(line[2])
            self.drag.components.add_item(drg)
    
    def save(self,aircraftName):
        name = aircraftName

    def _get_volume_ratios(self):
        pass
    def update_mass(self):
        """
        Update airframe mass of the aircraft
        """
        self.mass = weight.get_general_aviation_mass(self)
    def get_mass_empty(self, updMass=True):
        """
        Returns empty mass of current configuration
        """
        if updMass:
            self.update_mass()
        emptyMass = self.mass.airframe.get_total_mass()
        return emptyMass
    def get_mass_total(self, updMass=False):
        """
        Returns total mass of current configuration
        """
        if updMass:
            self.update_mass()
        totalMass = self.mass.total.get_total_mass()
        return totalMass
    def get_CG(self, updMass=False):
        """
        Returns CG location of current configuration
        """
        if updMass:
            self.update_mass()
        CG = self.mass.total.get_CG()
        if CG[0]>self.landingGear.groundContact_X[1]:
            warn('Aircraft CG should be in front of main landing gear')
        return CG
    def get_inertia(self, updMass=False):
        """
        Returns inertia of current configuration
        """
        if updMass:
            self.update_mass()
        return self.mass.inertia
    def update_drag(self,velocity=0.0,altitude=0.0):
        """
        Updates parasite drag of current configuration
        
        Parameters
        ----------
        
        velocity : float, m/sec
        altitude : float, m
        """
        if velocity<=0.0:
            velocity=self.designGoals.designSpeed
        if altitude<=0.0:
            altitude=self.designGoals.designAltitude
        self.drag = drag.get_total_drag(self,velocity,altitude)
    def get_drag(self,velocity=0.0,altitude=0.0):
        """
        Returns parasite drag of current configuration
        
        Parameters
        ----------
        
        velocity : float, m/sec
        altitude : float, m
        """
        self.update_drag(velocity,altitude)
        return self.drag.total.get_total_drag()
    def update_aero_trim(self,velocity=0.0,altitude=0.0,loadFactor=1.0,CmTrim=0.0,flapDefl=0.0):
        """
        Updates aerodynamic information stored in self.aeroResults at given 
        flight conditions
        
        Parameters
        ----------
        
        velocity : float, m/sec
            aircraft true airspeed
        altitude : float, m
            altitude at which calculation is performed
        loadFactor : float
            load factor
        CmTrim : float
            trim pitch moment coefficient
        """
        if velocity<=0.0:
            velocity=self.designGoals.designSpeed
        aero = aerodynamics.aerodynamics(self)
        atmosphere = FC.ISAtmosphere_simple(altitude)
        rho = atmosphere.density
        aero.update(velocity,rho,loadFactor,CmTrim,flapDefl)
        self.aeroResults = aero.results
    def analyze_aero_trim(self,velocity=0.0,altitude=0.0,loadFactor=1.0,CmTrim=0.0,flapDefl=0.0):
        """
        Updates aerodynamic information using self.update_aero_trim and 
        returns result
        """
        self.update_aero_trim(velocity,altitude,loadFactor,CmTrim,flapDefl)
        return self.aeroResults

class designGoals():
    def __init__(self,fuelMass,designGrossMass,designSpeed,designAltitude,designLoadFactor,designLandingLoadFactor,occupantNumber,payloadMass):
        self.fuelMass                           =fuelMass            
        self.designGrossMass                    =designGrossMass
        self.designSpeed                        =designSpeed
        self.designAltitude                     =designAltitude
        self.designLoadFactor                   =designLoadFactor
        self.designLandingLoadFactor            =designLandingLoadFactor
        self.occupantNumber                     =occupantNumber
        self.payloadMass                        =payloadMass


class wing():
    """
    Class contains wing configuration.
    """
    def __init__(self,segmentLengths,chords,airfoilName,segmentOffsets,
                 segmentDihedrals,twistAngles,incidence,aapex,material,symFlag,
                 vertical=False):
        pth=paths.Database()
        self._foilDBpath=pth.airfoil        
        self.setWing(segmentLengths,chords,airfoilName,segmentOffsets,
                     segmentDihedrals,twistAngles,incidence,aapex,material,symFlag,vertical)                
    def __str__(self):
        self._getWingData()
        return "Wing Span = %4.3f m"%(self.span)
    def setWing(self,segmentSpans,chords,airfoilName,segmentOffsets,
                segmentDihedrals,twistAngles,incidence,aapex,material,
                symFlag,vertical):
        foils=list()
        for name in airfoilName:
            airfoil=Airfoil_Class.Airfoil()
            airfoil.read_xls(str(airfoilName[0]),self._foilDBpath)
            airfoil.polar.calc_clmax()
            foils.append(airfoil)
        self.numSegments                        =len(segmentSpans)
        self.segmentSpans                       =segmentSpans
        self.chords                             =chords
        self.airfoilName                        =airfoilName
        self.airfoil                            =foils
        self.segmentOffsets                     =segmentOffsets
        self.segmentDihedrals                   =segmentDihedrals
        self.segmentDihedrals_rad               =segmentDihedrals*numpy.pi/180.
        self.twistAngles                        =twistAngles
        self.incidence                          =incidence
        self.aapex                              =aapex
        self.material                           =material
        self.symmetryFlag                       =symFlag  
        self.vertical                           =vertical
        self.controlSurfaces                    =list()
        self._get_segment_aapex()
        self._getWingData()
    def create_mesh(self):
        if self.vertical:
            self._create_mesh_vertical('b',0.5,20)
        else:
            self._create_mesh_horizontal('b',0.5,20)
    def _create_mesh_horizontal(self,color='y',opacity=0.8,numPts=20):
        self.meshX = list()
        self.meshY = list()
        self.meshZ = list()
        nSeg = len(self.segmentSpans)
        refPt = copy(self.aapex)
        refPt[1] += self.segmentSpans[0]
        for i in range(nSeg):
            airfoil = self.airfoil[i].redim(numPts,overwrite=False)
            section = airfoil*self.chords[i]
            axis1 = [self.chords[i]*0.25,0]
            if i==0:
                angle = self.incidence
                section = self.rotate2D(section,axis1,angle)
                wingX = section[:,0] + refPt[0]
                wingY = numpy.zeros(numPts)
                wingZ = section[:,1] + refPt[2]
            else:
                angle = self.incidence + self.twistAngles[i-1]
                section = self.rotate2D(section,axis1,angle)
                refPt[0] += self.segmentOffsets[i-1]
                refPt[1] += self.segmentSpans[i]
                refPt[2] += self.segmentSpans[i]*tan(radians(self.segmentDihedrals[i-1]))
                secX = section[:,0] + refPt[0]
                secY = numpy.zeros(numPts)  + refPt[1]
                secZ = section[:,1] + refPt[2]
                wingX = numpy.vstack([secX,wingX,secX])
                wingY = numpy.vstack([-secY,wingY,secY])
                wingZ = numpy.vstack([secZ,wingZ,secZ])
        self.meshX = wingX
        self.meshY = wingY
        self.meshZ = wingZ
        self.color = color
        self.opacity = opacity
    def _create_mesh_vertical(self,color='y',opacity=0.8,numPts=20):
        nSeg = len(self.segmentSpans)
        refPt = copy(self.aapex)
        refPt[2] += self.segmentSpans[0]
        for i in range(nSeg):
            airfoil = self.airfoil[i].redim(numPts,overwrite=False)
            section = airfoil*self.chords[i]
            axis1 = [self.chords[i]*0.25,0.0]
            if i==0:
                angle = self.incidence
                section = self.rotate2D(section,axis1,angle)
                vstabX = section[:,0] + refPt[0]
                vstabY = section[:,1] + refPt[1]
                vstabZ = numpy.zeros(numPts) + refPt[2]
            else:
                angle = self.incidence + self.twistAngles[i-1]
                section = self.rotate2D(section,axis1,angle)
                refPt[0] += self.segmentOffsets[i-1]
                refPt[1] += self.segmentSpans[i]*tan(radians(self.segmentDihedrals[i-1]))
                refPt[2] += self.segmentSpans[i]
                secX = section[:,0] + refPt[0]
                secY = section[:,1] + refPt[1]
                secZ = numpy.zeros(numPts)  + refPt[2]
                vstabX = numpy.vstack([vstabX,secX])
                vstabY = numpy.vstack([vstabY,secY])
                vstabZ = numpy.vstack([vstabZ,secZ])
        self.meshX = vstabX
        self.meshY = vstabY
        self.meshZ = vstabZ
        self.color = color
        self.opacity = opacity
    def rotate2D(self,curve, axis, angle):
        """
        rotates 2D curve around 2D axis point to a given angle
        
        :ivar array curve: array of curve points
        :ivar array axis: rotation axis in format array([x,y])
        :ivar float angle: rotation anlge in degrees
        
        :math:`x_{new} = cos(\theta)x - sin(\theta)y`
        
        :math:`y_{new} = sin(\theta)x + cos(\theta)y`
        """
        angle = radians(angle)
        rotMatrix = numpy.array([[cos(angle), -sin(angle)],[sin(angle), cos(angle)]])
        return numpy.dot( (curve-axis), rotMatrix) + axis
    def set_aapex(self,aapex):
        #TODO: horizontal tail aapex location moves to the right
        self.aapex = aapex
        fuelMass = self.fuelTank.fuelMass
        self.setFuelTank(fuelMass,1,self.fuelTankCGspanRatio,self.fuelTankCGchordRatio)
        self._getWingData()
    def _get_segment_aapex(self):
        """
        calculates location of each segment's leading edge point.
        Twist and incidence angles are not considered
        """
        x = 0.0
        y = self.segmentSpans[0]
        z = 0.0
        n = len(self.segmentOffsets)
        segmentAapex = numpy.array([x,y,z])
        for i in range(n):
            x  = self.segmentOffsets[i]
            y += self.segmentSpans[i+1]
            z += self.segmentSpans[i+1]*numpy.tan(self.segmentDihedrals_rad[i])
            segmentAapex = numpy.vstack([segmentAapex,[x,y,z]])
        self.segmentAapex = segmentAapex
    def airfoil_clmax(self,airfoilNum,Mach):
        clmax_list=self.airfoil[airfoilNum].polar.clmax
        Mach_list =self.airfoil[airfoilNum].polar.Mach
        return interp(Mach,Mach_list,clmax_list)
    def setFuelTank(self,fuelMass,mirrorFlag,CG_spanRatio,CG_chordRatio):
        CG=self.locateOnWing(CG_spanRatio,CG_chordRatio)
        self.fuelTankCGspanRatio = CG_spanRatio
        self.fuelTankCGchordRatio = CG_chordRatio
        self.fuelTank=_wingFuelTank(fuelMass,mirrorFlag,CG,CG_spanRatio,CG_chordRatio)
    def addAileron(self,location):
        self.aileron=_controlSurface(self,'plain',location,0)
    def addElevator(self,location):
        self.elevator=_controlSurface(self,'plain',location,1)
    def addRudder(self,location):
        self.rudder=_controlSurface(self,'plain',location,0)
    def addFlap(self,location,configuration):
        self.flap=_controlSurface(self,configuration,location,1)
    def _getWingData(self):
        self._get_segment_aapex()
        if self.vertical:
            self.span=numpy.sum(self.segmentSpans)
        else:
            self.span=numpy.sum(self.segmentSpans) * 2.0
        self.MAC = self._getMAC()
        self.area = self._getWingArea()
        self.segCoords=self._getSegCoords()
        self.aspectRatio=self.span**2/self.area
        self._calc_wetted_area()
        self._getBasicGeometryParams()
        self._getCentroid()
        self.create_mesh()
    def _getCentroid(self):
        xLE=self.segCoords[0]
        yLE=self.segCoords[1]
        xTE=xLE+self.chords[::-1]
        yTE=yLE
        x=numpy.concatenate((xLE,xTE[::-1]))
        y=numpy.concatenate((yLE,yTE[::-1]))
        Cx,Cy,A=miscTools.polyCentroid(x,y)
        b=numpy.cumsum(self.segmentSpans)
        Cz=0.0
        i=find(b>Cy)[0]-1
        R=(Cy-b[i])/(b[i+1]-b[i])
        Cz=R*self.segCoords[2][i]+(1.-R)*self.segCoords[2][i+1]
        if self.vertical:
            self.centroid = numpy.array([Cx,0.0,Cy+self.aapex[2]])
        else:
            self.centroid=numpy.array([Cx,Cy,Cz])
    def locateOnWing(self,spanRatio,chordRatio):
        x = self.aapex[0]
        z = self.aapex[2]
        spanLocation = spanRatio*self.span/2.0
        if spanRatio==0.0:
            i=0
        else:
            for idx, bb in enumerate(self.segmentAapex[:,1]):
                if spanLocation<=bb:
                    i=idx-1
                    break
        locSpan = spanLocation - self.segmentAapex[i,1]
        locSpanRatio = locSpan/self.segmentSpans[i+1]
        tmp1 = self.chords[i] - self.chords[i+1]
        tmp2 = self.segmentAapex[i+1,1]-spanLocation
        locChord = tmp1*tmp2 / self.segmentSpans[i+1] + self.chords[i+1]
        dx = locSpanRatio*self.segmentOffsets[i] + locChord*chordRatio
        dz = locSpan * numpy.tan(self.segmentDihedrals_rad[i])
        x += self.segmentAapex[i,0] + dx
        z += self.segmentAapex[i,2] + dz
        return numpy.array([x,spanLocation,z])
    def _getBasicGeometryParams(self):
        b=self.segmentSpans
        o=numpy.concatenate((numpy.array([0.0]),self.segmentOffsets))
        c4line=o+self.chords/4.        
        d =numpy.zeros(len(c4line)-1)
        t =numpy.zeros(len(c4line)-1)
        tc=numpy.zeros(len(c4line)-1)
        for i in range(1,len(c4line),1):         
            d[i-1]=c4line[i]-c4line[i-1]
            t[i-1]=self.chords[i]/self.chords[i-1]
            tc[i-1]=self.airfoil[i].thickness
        theta=numpy.arctan(d/b[1:])
        self.sweep=(numpy.sum(theta*b[1:])/numpy.sum(b))
        self.taper=(numpy.sum(t*b[1:])/numpy.sum(b))
        self.thicknessRatio=(numpy.sum(tc*b[1:])/numpy.sum(b))
    def _getSegCoords(self):
        N=len(self.segmentSpans)
        xOffset=numpy.concatenate((numpy.zeros(1),self.segmentOffsets))
        zOffset=numpy.zeros(N)
        for i in range(1,N,1):
            #zOffset[i]=zOffset[i-1]+self.segmentSpans[i]*numpy.tan(self.segmentDihedrals_rad[i-1])
            zOffset[i]=self.segmentSpans[i]*numpy.tan(self.segmentDihedrals_rad[i-1])
        x=self.aapex[0]+xOffset
        y=self.aapex[1]+self.segmentSpans.cumsum()
        z=self.aapex[2]+zOffset.cumsum()
        return x,y,z
    def _getWingArea(self):
        self.segAreas=numpy.zeros(len(self.segmentSpans)-1)
        for i in range(len(self.segmentSpans)-1):
            self.segAreas[i]=self.segmentSpans[i+1]*(self.chords[i]+self.chords[i+1])/2.0
        area=numpy.sum(self.segAreas)*2.0
        if self.vertical:
            area = area / 2.0
        return area
    def _getMAC(self):
        b2=numpy.sum(self.segmentSpans)
        S2=self._getWingArea()/2
        MAC = 1.0/S2*miscTools.integrateRomberg(self._cATy_sq,0,b2,10.0e-6)
        if self.vertical:
            MAC *= 0.5 #TODO check MAC calculation - smth strange
        return MAC
    def _calc_wetted_area(self):
        self.wettedArea = list()
        for i,span in enumerate(self.segmentSpans[1:]):
            a = self.airfoil[i].get_length() * self.chords[i]
            b = self.airfoil[i+1].get_length() * self.chords[i+1]
            area = (a+b)*span/2
            if self.symmetryFlag:
                area = 2*area
            self.wettedArea.append(area)
        self.wettedAreaFull = numpy.sum(self.wettedArea)
    def _cATy_sq(self,y):
        b2=numpy.cumsum(self.segmentSpans)
        i=numpy.where(b2>=y)[0][0]
        c1=self.chords[i-1]
        c2=self.chords[i]     
        r=(y-b2[i-1])/(b2[i]-b2[i-1])
        c=r*c2+(1-r)*c1
        return c**2.0


class _wingFuelTank():
    def __init__(self,fuelMass,mirrorFlag,CG,chordRatio,spanRatio):
        self.fuelMass   = fuelMass
        self.mirrorFlag = mirrorFlag
        self.CG         = CG
        self.chordRatio = chordRatio
        self.spanRatio  = spanRatio


class _controlSurface():
    def __init__(self,wing,configuration,location,symDeflectionFlag):
        # configuration = plain, slotted, fowler
        allow   =numpy.array(['plain','slotted','split'])
        if any(allow==configuration)==False:
            print 'Control type must be either plain, slotted, or split'
        assert any(allow==configuration)==True
        startSec=numpy.argmax(location<1)
        stopSec =len(location)-numpy.argmax(location[::-1]<1)-1
        self.flappedArea                        =numpy.sum(wing.segAreas[startSec:stopSec])*2.
        self.flappedChord                       =numpy.average(wing.chords[startSec:stopSec])
        self.flapChord                          =self.flappedChord*numpy.average(location[startSec:stopSec])
        self.avgChordRatio                      =1.-numpy.average(location[startSec:stopSec])
        self.avgChord                           =wing.MAC*self.avgChordRatio
        self.configuration                      =configuration
        self.location                           =location
        self.symDeflectionFlag                  =symDeflectionFlag
        chords = wing.chords*(1.0-location)
        for i,c in enumerate(chords):
            if not c==0:
                idx=i
                break
            else:
                idx=0
        self.area = 0.5*(chords[idx]+chords[idx+1])*wing.segmentSpans[idx+1]


class fuselage():
    def __init__(self,length,diameter,wettedArea,sideProfile_X,sideProfile_Y,topProfile_X,topProfile_Y,material):
        self.setFuselage(length,diameter,wettedArea,sideProfile_X,sideProfile_Y,topProfile_X,topProfile_Y,material)
    def setFuselage(self,length,diameter,wettedArea,sideProfile_X,sideProfile_Y,topProfile_X,topProfile_Y,material):
        self.length                             =length
        self.diameter                           =diameter
        self.sideProfile_X                      =sideProfile_X
        self.sideProfile_Y                      =sideProfile_Y
        self.coords                             =numpy.array([sideProfile_X,sideProfile_Y]).transpose()
        self.topProfile_X                       =topProfile_X
        self.topProfile_Y                       =topProfile_Y
        self._fit_profile(length,diameter)
        self.material                           =material
        self.crossSecRatio = 1.00    # 1.0 if cs is ellipse, >1.0 if close to square
        self.separate_points()
        tmpArea = self.calc_wetted_area()
        if wettedArea<=0.5*tmpArea or wettedArea>=1.5*tmpArea:
            self.wettedArea = tmpArea
        else:
            self.wettedArea = wettedArea
        self._getCentroid()
    def write_txt(self,filePath):
        with open(filePath, "w") as coordData:
            for coord in self.coords:
                coordData.write('%.6f\t%.6f\n'%(coord[0],coord[1]))
    def separate_points(self):
        """
        from fuselage.sideProfile and fuselage.topProfile creates
        upper, lower and side curve to be used with geometry modules
        """
        pts = numpy.vstack([self.sideProfile_X,self.sideProfile_Y])
        pts = numpy.transpose(pts)
        nosePtIdx = numpy.argmin(pts[:,0])
        self.profileUp = numpy.flipud(pts[:nosePtIdx+1])
        self.profileLo = pts[nosePtIdx:]
        self.profileSide = numpy.vstack([self.topProfile_X,self.topProfile_Y])
        self.profileSide = numpy.transpose(self.profileSide)
    def _fit_profile(self,length,diameter):
        if length>=0.0:
            lengthExact = max(self.sideProfile_X) - min(self.sideProfile_X)
            lengthMult = length / lengthExact
            self.sideProfile_X *= lengthMult
            self.topProfile_X *= lengthMult
        if diameter>=0.0:
            diamExact   = max(self.sideProfile_Y) - min(self.sideProfile_Y)
            diamMult   = diameter / diamExact
            self.sideProfile_Y *= diamMult
            self.topProfile_Y *= diamMult
    def _getCentroid(self):
        Cy=0.
        Cx1,Cz,A1=miscTools.polyCentroid(self.sideProfile_X,self.sideProfile_Y)
        Cx2,Cy2,A2=miscTools.polyCentroid(self.topProfile_X,self.topProfile_Y)
        A1=numpy.abs(A1)
        A2=numpy.abs(2.*A2)
        Cx=(A1*Cx1+A2*Cx2)/(A1+A2)
        self.centroid=numpy.array([Cx,Cy,Cz])
    def calc_wetted_area(self):
        """
        calculates wetted area of the fuselage assuming elliptic cross sections.
        """
        length1 = list()
        offset = list()
        for i in range(len(self.profileUp)-1):
            h = self.profileUp[i+1,0] - self.profileUp[i,0]
            a1 = self.profileUp[i,1] - self.profileLo[i,1]
            a2 = self.profileUp[i+1,1] - self.profileLo[i+1,1]
            a = (a1+a2)/4
            b = (self.profileSide[i,1] + self.profileSide[i+1,1])/2
            offset.append(h)
            length1.append(geom.ellipse_circumference(a,b))
        area = simps(length1,numpy.cumsum(offset))
        return area * self.crossSecRatio

class engine():
    """
    Describes engine configuration.
    """
    def __init__(self,nameList,CG_X,CG_Y,CG_Z):
        """
        Parameters
        ----------
        
        nameList : list, string
            list of engine names
        CG_X : float, meters
            x coordinate of engine CG
        CG_Y : float, meters
            y coordinate of engine CG
        CG_Z : float, meters
            z coordinate of engine CG
        """
        pth=paths.Database()                
        self.CG_X                               =CG_X
        self.CG_Y                               =CG_Y
        self.CG_Z                               =CG_Z
        self.engineList                         =list()
        database                                =dbTools.loadDB(pth.engine)
        if type(nameList)==type(list()):
            self.nameList=nameList
        else:
            self.nameList=list()
            self.nameList.append(nameList)
        for name in self.nameList:
            dataSheet=database.selectByName(name)
            self.engineList.append(eng.pistonEngine(dataSheet))


class landingGear():
    def __init__(self,gearType,tireWidth,tireDiameter,strutLength,strutType,fairingType,groundContact_X,groundContact_Y,groundContact_Z):
        self.type                               =gearType
        self.tireWidth                          =tireWidth
        self.tireDiameter                       =tireDiameter
        self.strutLength                        =strutLength
        self.strutType                          =strutType
        self.fairingType                        =fairingType
        self.groundContact_X                    =groundContact_X
        self.groundContact_Y                    =groundContact_Y
        self.groundContact_Z                    =groundContact_Z


class propeller(aircraft):
    def __init__(self,nameList,propHub_X,propHub_Y,propHub_Z,minBeta,maxBeta):
        pth=paths.Database()    
        self.hub_X                              =propHub_X
        self.hub_Y                              =propHub_Y
        self.hub_Z                              =propHub_Z
        self.betaRange                          =[minBeta,maxBeta]
        self._database                          =dbTools.loadDB(pth.prop)
        if type(nameList)==type(list()):
            self.nameList=nameList
        else:
            self.nameList=list()
            self.nameList.append(nameList)
        self.propellerList                      =list()
        for i,name in enumerate(self.nameList):
            P=prop.propeller()
            P.read_xls(name)
            betaLo=self.betaRange[0][i]
            betaHi=self.betaRange[1][i]
            P.betaRange=[betaLo,betaHi]
            self.propellerList.append(P)


class VLMsetup():
    """
    Data structure contains data for running aerodynamic analysis using VLM
    
    Parameters
    ----------
    
    wing_spanwise : integer
        number of wing spanwise panels
    wing_chordwise : integer
        number of wing chordwise panels
    hStab_spanwise : integer
        number of horizontal tail spanwise panels
    hStab_chordwise : integer
        number of horizontal tail chordwise panels
    vStab_spanwise : integer
        number of vertical tail spanwise panels
    vStab_chordwise : integer
        number of vertical tail chordwise panels
    cosDistribution : float
        parameter of AVL panels cosine distribution
    sinDistribution : float
        parameter of AVL panels sine distribution
    """
    def __init__(self,wing_spanwise=20,wing_chordwise=6,hStab_spanwise=10,hStab_chordwise=6,vStab_spanwise=7,vStab_chordwise=6,cosDistribution=0,sinDistribution=0):
        self.wing_spanwise                      =wing_spanwise
        self.wing_chordwise                     =wing_chordwise
        self.hStab_spanwise                     =hStab_spanwise
        self.hStab_chordwise                    =hStab_chordwise
        self.vStab_spanwise                     =vStab_spanwise
        self.vStab_chordwise                    =vStab_chordwise
        self.cosDistribution                    =cosDistribution
        self.sinDistribution                    =sinDistribution

# --- test section ---
def runTest1():
    pth=paths.Database()
    acDB =dbTools.loadDB(pth.aircraft)
    #acDB=dbTools.loadDB("//database//aircraft.xls")
    workbook = acDB.workbook
    worksheets = workbook.sheet_names()
    for worksheet_name in worksheets:
        worksheet = workbook.sheet_by_name(worksheet_name)
        print worksheet
    sheetList,sheetNames=acDB.loadData()
    print sheetList
    print sheetNames
    acData=acDB.selectByName("sampleInput1")
    print acData
    ac=aircraft(acData)
    ac=load("sampleInput2")
    print ac.wing.airfoilName
    print ac.wing.aapex
    print ac.wing.incidence
    print ac.engine.CG_X
    print ac.landingGear.groundContact_X
    print ac.VLMsetup.wing_spanwise  
    print ac.mass[0].name
    print ac.mass.nameList()
    print 'HT'
    print ac.hStab.segmentDihedrals
    print ac.hStab.segmentOffsets
    print ac.vStab.segmentSpans

    ac.mass.removeItem("passenger")
    print ac.wing.airfoilName
    print ac.fuselage.sideProfile_X
    print ac.mass.nameList()
    print ac.engine.nameList
    print ac.propeller.nameList
    print ac.engine.engineList[0].RPM
def runTest2():
    ac=load('V0510')
    print '{0:12} = {1:.4f}'.format('Wing Span',ac.wing.span)
    print '{0:12} = {1:.4f}'.format('Wing Area',ac.wing.area)
    print '{0:12} = {1:.4f}'.format('Wing MAC',ac.wing.MAC)
    print '{0:12} = {1:.4f}'.format('Wing Sweep',ac.wing.sweep)
    print '{0:12} = {1:.4f}'.format('Wing AR',ac.wing.aspectRatio)
    print '{0:12} = {1:.4f}'.format('hStab Span',ac.hStab.span)
    print '{0:12} = {1:.4f}'.format('hStab Area',ac.hStab.area)
    print '{0:12} = {1:.4f}'.format('hStab MAC',ac.hStab.MAC)
    print '{0:12} = {1:.4f}'.format('hStab Sweep',ac.hStab.sweep)
    print '{0:12} = {1:.4f}'.format('hStab AR',ac.hStab.aspectRatio)
    print '{0:12} = {1:.4f}'.format('vStab Span',ac.vStab.span)
    print '{0:12} = {1:.4f}'.format('vStab Area',ac.vStab.area)
    print '{0:12} = {1:.4f}'.format('vStab MAC',ac.vStab.MAC)
    print '{0:12} = {1:.4f}'.format('vStab Sweep',ac.vStab.sweep)
    print '{0:12} = {1:.4f}'.format('vStab AR',ac.vStab.aspectRatio)
    print ac.wing.locateOnWing(0.8,0.3)
    print ac.fuselage.centroid
    print ac.wing.airfoil_clmax(0,0.1)
    ac.display()    

def debug1():
    ac = load('V0510')
    ac.mass.airframe.display()
    ac.display()
if __name__=="__main__":
#    ac = load('F01')
#    print ac.get_CG()
#    ac.display()
    runTest2()
