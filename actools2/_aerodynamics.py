# -*- coding: utf-8 -*-
r"""
This module contains classes for performing 3d aerodynamcis analysis on
aircraft geometry models.

Module provides the following functions:
    - input/output: sets flight conditions to be analyzed
f
    - aerodynamic analysis: performs analysis using AVL

    - stability analysis: interprets stability derivatives to assess 
    dynamic motions
"""
import paths
import subprocess as sp
import os
import numpy
import shlex
import FlightConditions as fc
import re
from StringIO import StringIO

class aerodynamics():
    """
    Simplified wrapper for the aerodynamics library
    
    Parameters
    ----------
    
    aircraft : aircraft
        aircraft configuration
    
    Examples
    --------
    >>> import aircraft
    >>> import FlightConditions as FC
    >>> ac = aircraft.load('V05')
    >>> aero = aerodynamics(ac)
    >>> V = 50.0
    >>> alt = 1500.0
    >>> dT = 0.0
    >>> n = 1.0
    >>> cruise = FC.FlightConditions(V,alt,dT,ac.wing.MAC)
    >>> rho = cruise.atmosphere.density
    >>> aero.update(V,rho,n)
    >>> aero.results.display()
    ... 'displays full list of aerodynamic coefficients'
    """
    def __init__(self,aircraft):
        self.aircraft         =aircraft
        self.flightConditions =flightConditions()
        self.results          =None
    def update(self,velocity,density,loadFactor=1.0,CmTrim=0.0,flapDefl=0.0,
               mass=None,CG=None,inertia=None,name='default'):
        """
        Updates aerodynamics results in self.results
        
        Parameters
        ----------
        
        velocity : float, m/sec
            aircraft TAS
        density : float, kg/m**3
            air density
        loadFactor : float
            load factor
        CmTrim : float
            trim pitch moment coefficient
        flapDefl : float, deg
            plain flap deflection
        """
        self.flightConditions =flightConditions()
        densityAltitude       =fc.get_density_altitude(density)
        if mass==None:
            mass                  =self.aircraft.get_mass_total()
        if CG==None:
            CG                    =self.aircraft.get_CG()
        if inertia==None:
            inertia               =self.aircraft.get_inertia()
        Cd0                   =self.aircraft.get_drag(velocity,densityAltitude)
        self.flightConditions.addTrimmedFlightCondition(name,mass,CG,inertia,velocity,density,loadFactor,Cd0,CmTrim,flapDefl)
        vlm                   =aero3d_VLM(self.aircraft)
        out                   =vlm.runVLM(self.flightConditions)
        self.results          =out.results[0]
        self.results.name = name
        NP = self.results.xNP
        self.results.SM = (NP-CG[0])/self.aircraft.wing.MAC
        self.results.mass = mass
        self.results.CG = CG
        self.results.inertia = inertia
        self.results.velocity = velocity
        self.results.density = density
    def update_notrim(self,velocity,density,loadFactor=1.0,flap=0.0,aileron=0.0,
                      elevator=0.0,rudder=0.0,alpha=0.0,beta=0.0,mass=None,
                      CG=None,inertia=None):
        """
        Performs aerodynamic analysis at specified flight conditions.
        
        Parameters
        ----------
        
        velocity : float, m/s
        density : float, kg/m**3
        loadFactor : float
        flap : float,deg
            flap deflection angle. Note that plain flap is analyzed, slotted flap 
            analysis is not available
        aileron : +float, deg
            aileron deflection angle. Deflection is mirrored.
        elevator : float, deg
            elevator deflection angle
        rudder : float, deg
            rudder deflection angle
        alpha : float, deg
            angle of attack. Note that solver may fail on high angles of attack. 
            Use values [-5;5].
        beta : float, deg
            sideslip angle.
        mass : float, kg
            if mass is not defined, than it will be obtained from current 
            aircraft configuration
        CG : float vector, m
            aircraft center of gravity in format [x,y,z]. 
            If CG is not defined it will be obtained from current aircraft configuration
        inertia : float vector, kg*m**2
            aircraft moment of inertia in format [Ixx,Iyy,Izz].
        """
        if mass==None:
            mass = self.aircraft.get_mass_total()
        if CG==None:
            CG = self.aircraft.get_CG()
        if inertia==None:
            inertia = self.aircraft.get_inertia()
        altitude = fc.get_density_altitude(density)
        CD0 = self.aircraft.get_drag(velocity,altitude)
        flightCond = flightConditions()
        flightCond.add_flight_condition('untrimmed',mass,CG,inertia,velocity,density,
                                loadFactor,CD0,alpha,beta,flap,elevator,
                                aileron,rudder)
        vlm = aero3d_VLM(self.aircraft)
        out = vlm.runVLM(flightCond)
        return out.results[0]
    def Clmax(self,Mach,flapSetting_deg):
        if self.results==None:
            self.runDefaultCase()
        CLmax                 =self.CLmax_clean(Mach)
        dCLmax                =self._dCLmax_flap(flapSetting_deg)#*Sf/S      
        return CLmax+dCLmax
    def CLmax_clean(self,Mach):
        clmax2d               =self.aircraft.wing.airfoil_clmax(0,Mach)
        avgTaper              =self.aircraft.wing.taper
        k                     =self._kFcn(avgTaper)
        return k*clmax2d
    def _kFcn(self,tr):
        """
        Returns k, the max CL reduction factor in a finite tapered wing 
        Polynomial fit of a table from Roskam
        """
        return (-0.2917*tr**2+0.2917*tr+0.88)
    def _dCLmax_flap(self,flapSetting_deg):
        if flapSetting_deg>5.0:
            Sf                 =self.aircraft.wing.flap.flappedArea    
            S                  =self.aircraft.wing.area
            a                  =self.results.a
            cr                 =self.aircraft.wing.flap.avgChordRatio
            thetaf             =numpy.arccos(2*cr-1)
            tau                =1+numpy.sin(thetaf)/numpy.pi-thetaf/numpy.pi    
            N                  =self._getEffy(flapSetting_deg)
            dCL                =a*numpy.pi/180.0*tau*N*flapSetting_deg*Sf/S
        else: dCL=0.0
        return dCL
    def _getEffy(self,df):
        cfg                    =str(self.aircraft.wing.flap.configuration)
        if cfg=='plain':      N=-89.9015/df**2+15.8051/df+0.1191
        if cfg=='split':      N=0.4228/df-0.0031*df+0.539
        if cfg=='slotted':    N=(df/35.27)**5-(df/24.21)**4+(df/18.03)**3-(df/14.53)**2+(df/17.54)+0.570
        else:print 'Specified flap type is not recognized.  Should be either plain, split, or slotted'
        return N
    def runDefaultCase(self):
        V                      =self.aircraft.designGoals.designSpeed        
        atm                    =fc.ISAtmosphere(self.aircraft.designGoals.designAltitude)
        rho                    =atm.density
        self.update(V,rho)    


class flightConditions:
    """
    Sets and stores a list of flight conditions to be analyzed by the
    aerodynamics solvers.
    """
    def __init__(self):
        self.flightConditionList=list()
    def __getitem__(self,k):
        return self.flightConditionList[k]
    def __len__(self):
        return len(self.flightConditionList)
    def addTrimmedFlightCondition(self,name,mass,CG,inertia,velocity,density,
                                  loadFactor=1.0,Cd0=0.0,CmTrim=0.0,flapDefl=0.0,
                                  aileron=0.0,rudder=0.0):
        """
        Sets and stores a list of flight conditions to be analyzed by the
        aerodynamics solvers.
        
        Parameters
        ----------
        
        name : string
            name of the flight condition case
        mass : float, kg
            total mass of the aircraft
        CG : float array, m
            array of CG coordinates [x,y,z]
        inertia : float array
            array of moments of inertia [Ixx, Iyy, Izz]
        velocity : float, m/sec
            aircraft TAS
        density : float, kg/m**3
            atmospheric density
        loadFactor : float
            load factor for steady pull-out flight conditions
        Cd0 : float
            zero lift drag coefficient of aircraft
        CmTrim : float
            trim pitch moment coefficient. Set CmTrim=None if non-trimmed flight 
            condition is required
        flapDefl : float, deg
            plain flap deflection angle
        """
        newFC                      =self._flightCondition()
        newFC.pitchTrim            =True
        newFC.CG                   =CG
        newFC.name                 =name
        newFC.mass                 =float(mass)
        newFC.inertia              =inertia   
        newFC.velocity             =float(velocity)
        newFC.density              =float(density)
        newFC.densityAltitude      =fc.get_density_altitude(newFC.density)
        atm                        =fc.ISAtmosphere(newFC.densityAltitude,0.0)
        newFC.soundSpeed           =atm.soundSpeed
        newFC.pressure             =atm.pressure
        newFC.temperature          =atm.temperature
        newFC.Mach                 =newFC.velocity/newFC.soundSpeed
        newFC.loadFactor           =loadFactor
        newFC.Cd0                  =Cd0
        newFC.CmTrim               =CmTrim
        newFC.flap                 =flapDefl
        newFC.aileron              =aileron
        newFC.rudder               =rudder
        self.flightConditionList.append(newFC)

    def add_flight_condition(self,name,mass,CG,inertia,velocity,density,
                             loadFactor=1.0,CD0=0.0,alpha=0.0,beta=0.0,flap=0.0
                             ,elevator=0.0,aileron=0.0,rudder=0.0):
        """
        Non trimmed flight conditions
        """
        newFC                      =self._flightCondition()
        newFC.pitchTrim            =False
        newFC.CG                   =CG
        newFC.name                 =name
        newFC.mass                 =float(mass)
        newFC.inertia              =inertia   
        newFC.velocity             =float(velocity)
        newFC.density              =float(density)
        newFC.densityAltitude      =fc.get_density_altitude(newFC.density)
        atm                        =fc.ISAtmosphere(newFC.densityAltitude,0.0)
        newFC.soundSpeed           =atm.soundSpeed
        newFC.pressure             =atm.pressure
        newFC.temperature          =atm.temperature
        newFC.Mach                 =newFC.velocity/newFC.soundSpeed
        newFC.loadFactor           =loadFactor
        newFC.Cd0                  =CD0
        newFC.alpha                =alpha
        newFC.beta                 =beta
        newFC.flap                 =flap
        newFC.elevator             =elevator
        newFC.aileron              =aileron
        newFC.rudder               =rudder
        self.flightConditionList.append(newFC)
        
    def removeItem(self,itemName):
        try:
            k=self.nameList().index(itemName)
            self.flightConditionList.remove(self.flightConditionList[k])
        except:
            print "Specified flight condition not found"
    class _flightCondition:
        def __init__(self):
            atm0=fc.ISAtmosphere(0.0)
            self.name                           =str()                
            self.mass                           =0.0
            self.CG                             =numpy.array(3)
            self.inertia                        =numpy.array(3)
            self.density                        =atm0.density
            self.velocity                       =0.0
            self.aileron                        =0.0
            self.elevator                       =0.0
            self.flap                           =0.0
            self.rudder                         =0.0
            self.pitchTrim                      =True
            self.denistyAltitude                =0.0
            self.g                              =atm0.g
            self.loadFactor                     =1.0
            self.Cd0                            =0.0
            self.CmTrim                         =0.0
        def display(self):
            for attr, value in self.__dict__.iteritems():
                if type(value)==type(0.0):
                    print '{:<16} = {:<+12.5f}  '.format(attr,value)
        def setRudder(self,angle):
            """
            Sets the rudder angle for a given flight condition
            
            Parameters
            ----------
            
            angle : float, deg
                rudder angle in degrees
            """
            self.rudder=angle
        def setFlap(self,angle):
            """
            Sets the flap angle for a given flight condition
            
            Parameters
            ----------
            
            angle : float, deg
                flap angle in degrees
            """
            self.flap=angle
        def setAileron(self,angle):
            """
            Sets the aileron angle for a given flight condition
            
            Parameters
            ----------
            
            angle : float, deg
                aileron angle in degrees
            """
            self.aileron=angle
        def setElevator(self,angle):
            """
            Sets the elevator angle for a given flight condition
            
            Parameters
            ----------
            
            angle : float, deg
                elevator angle in degrees
            """
            self.elevator=angle
    def nameList(self):
        """
        returns a list of the flight condition names
        """
        outList=list()
        for k in range(self.flightConditionList.__len__()):
            outList.append(self.flightConditionList[k].name)
        return outList


class aero3d_VLM():
    """
    Main class for performing VLM analysis.  Writes out the required input
    files, runs the VLM code and a parasite drag calculation.  
    Performs cleanup, and returns results.
    
    Parameters
    ----------
    
    aircraft : aircraft
        aircraft configuration model
    flightConditions : flightCoditions
        list of flight conditions to be analyzed
    """
    def __init__(self,aircraft):
        pth=paths.MyPaths()
        self._tempFile=pth.getTmpFile()
        self.aircraft=aircraft
    def runVLM(self,flightConditions,runway=False):
        self.writeAVLinputFile(runway)
        results=runAVL(self._tmp_avl,flightConditions)
        i = 0
        for rslt,fc in zip(results.results,flightConditions):
            hm = self._calc_hinge_moments(rslt.hingeMoments,fc)
            results.results[i].hingeMoments.Mflap = hm[0]
            results.results[i].hingeMoments.Maileron = hm[1]
            results.results[i].hingeMoments.Melevator = hm[2]
            results.results[i].hingeMoments.Mrudder = hm[3]
            i += 1
        self.cleanup()
        return results
    def _calc_hinge_moments(self,hingeCoef,flightCondition):
        def _get_data(control,sym=False):
            S = control.area
            c = control.avgChord
            if sym: 
                S *= 2.0
            return S, c
        velocity = flightCondition.velocity
        q = flightCondition.density*velocity*velocity/2.0
        Se, ce = _get_data(self.aircraft.hStab.elevator,True)
        Sa, ca = _get_data(self.aircraft.wing.aileron)
        Sr, cr = _get_data(self.aircraft.vStab.rudder)
        Sf, cf = _get_data(self.aircraft.wing.flap)
#        S = self.aircraft.wing.area
#        c = self.aircraft.wing.MAC
#        Se, Sa, Sr, Sf = S,S,S,S
#        ce, ca, cr, cf = c,c,c,c
        flap = hingeCoef.flap* q*Sf*cf
        ail  = hingeCoef.aileron* q*Sa*ca
        elev = hingeCoef.elevator* q*Se*ce
        rud  = hingeCoef.rudder* q*Sr*cr
        return flap, ail, elev, rud
    def writeAVLinputFile(self,runway):
        """
        Writes a geometry input file for AVL
        """
        ac=self.aircraft  
        vlm=ac.VLMsetup
        self._tmpFileList=list()        
        self._tmp_avl=self._tempFile+'.avl'        
        self._tmp_fus=self._tempFile+'_fus.dat'
        
        self._tmpFileList.append(self._tmp_avl)
        self._tmpFileList.append(self._tmp_fus)   
        ac.fuselage.write_txt(self._tmp_fus)
        wingFiles =self.writeAirfoilFiles(ac.wing,'W')
        hStabFiles=self.writeAirfoilFiles(ac.hStab,'H')
        vStabFiles=self.writeAirfoilFiles(ac.vStab,'V')
        with open(self._tmp_avl, "w") as avlin:
            avlin.write('AVL INPUT FILE - GENERATED BY AERODYNAMICS.PY\n')
            avlin.write('0\n')
            if runway:
                groundLocation = ac.landingGear.groundContact_Z[1]
                avlin.write('0 1 %.4f\n'%(groundLocation))
            else:
                avlin.write('0 0 0\n')
            avlin.write('%f %f %f\n'%(ac.wing.area,ac.wing.MAC,ac.wing.span))
            avlin.write('0 0 0\n')
            avlin.write('0\n')
            avlin.write('#====================BODY======================#\n')
            avlin.write('BODY\n')
            avlin.write('FUSELAGE\n')
            avlin.write('28 2\n')
            avlin.write('BFIL\n')
            avlin.write(self._tmp_fus+'\n')
            avlin.write('#==================MAIN WING===================#\n')
            avlin.write('SURFACE\n')
            avlin.write('WING\n')
            avlin.write('%d %d %d %d\n'%(vlm.wing_chordwise,vlm.cosDistribution,vlm.wing_spanwise,vlm.cosDistribution))
            avlin.write('YDUPLICATE\n')
            avlin.write('0\n')
            avlin.write('ANGLE\n')
            avlin.write('%f\n'%(ac.wing.incidence))
            self._writeWingSections(ac.wing,wingFiles,'WING',avlin)
            avlin.write('#====================HSTAB=====================#\n')
            avlin.write('SURFACE\n')
            avlin.write('HSTAB\n')
            avlin.write('%d %d %d %d\n'%(vlm.hStab_chordwise,vlm.cosDistribution,vlm.hStab_spanwise,vlm.cosDistribution))
            avlin.write('YDUPLICATE\n')
            avlin.write('0\n')
            avlin.write('ANGLE\n')
            avlin.write('%f\n'%(ac.hStab.incidence))
            self._writeWingSections(ac.hStab,hStabFiles,'HSTAB',avlin)
            avlin.write('#====================VSTAB=====================#\n')
            avlin.write('SURFACE\n')
            avlin.write('VSTAB\n')
            avlin.write('%d %d %d %d\n'%(vlm.vStab_chordwise,vlm.cosDistribution,vlm.vStab_spanwise,vlm.cosDistribution))
            self._writeWingSections(ac.vStab,vStabFiles,'VSTAB',avlin,True)
    def _writeWingSections(self,wing,wingFiles,name,fileObj,vStabFlag=False):
        secNum=1;
        for i in range(len(wing.segmentSpans)):
            x,y,z=wing.segCoords
            fileObj.write('#------------------%s SEC %d-----------------#\n'%(name,secNum))
            fileObj.write('SECTION\n')
            if vStabFlag:
                fileObj.write('%f %f %f %f %f\n'%(x[i],z[i],y[i],wing.chords[i],wing.twistAngles[i]))
            else:                   
                fileObj.write('%f %f %f %f %f\n'%(x[i],y[i],z[i],wing.chords[i],wing.twistAngles[i]))
            fileObj.write('AFIL\n')
            fileObj.write(wingFiles[i]+'\n')
            if hasattr(wing,'aileron'):
                if wing.aileron.location[i]<>1:
                    fileObj.write('CONTROL\n')
                    fileObj.write('aileron 1 %f 0 0 0 -1\n'%(wing.aileron.location[i]))
            if hasattr(wing,'flap'):
                if wing.flap.location[i]<>1:
                    fileObj.write('CONTROL\n')
                    fileObj.write('flap 1 %f 0 0 0 +1\n'%(wing.flap.location[i]))            
            if hasattr(wing,'elevator'):
                if wing.elevator.location[i]<>1:
                    fileObj.write('CONTROL\n')
                    fileObj.write('elevator 1 %f 0 0 0 +1\n'%(wing.elevator.location[i]))   
            if hasattr(wing,'rudder'):
                if wing.rudder.location[i]<>1:
                    fileObj.write('CONTROL\n')
                    fileObj.write('rudder 1 %f 0 0 0 +1\n'%(wing.rudder.location[i]))
            secNum+=1  
    def writeAirfoilFiles(self,wing,name):
        fileList=list()
        for i,foil in enumerate(wing.airfoil):            
            pth=self._tempFile+'_'+name+str(i)+'.dat'
            foil.write_txt(pth)
            self._tmpFileList.append(pth)
            fileList.append(pth)
        return fileList
    def cleanup(self):
        for fname in self._tmpFileList:
            os.remove(fname)


class runAVL():
    """
    Executes AVL VLM code, reads and returns the results
    
    Parameters
    ----------
    
    inputFilePath : path
        full path of the AVL geometry input file
    flightConditions : list
        list of flight conditions to be analyzed
    Cd0 : float
        parasite drag constant
    """
    def __init__(self,inputFilePath,flightConditions):
        self.results=list()
        for i in range(len(flightConditions)):
            self.results.append(_results())
        self.inputFilePath=paths.fixPaths(inputFilePath)
        self.flightConditions=flightConditions
        pth=paths.MyPaths()
        cmd=pth.avl+' '+inputFilePath
        args=shlex.split(cmd,False,os.name=='posix')
        self._ps=sp.Popen(args,stdin=sp.PIPE,stderr=sp.PIPE,stdout=sp.PIPE)
        self._ps.stderr.close()                
        self._runFlightConditions(self.flightConditions)        
        self.rawOutput=str(self._ps.stdout.read()).replace('\r','')
        self._processOutput()
        self._setStabResults()
    def _issueCmd(self,cmd):
        self._ps.stdin.write(str(cmd)+'\n')
    def _runFlightConditions(self,FCs):
        """
        run AVL for all flight conditions specified in FCs
        """
        # General Options
        self._issueCmd(' ')
        self._issueCmd('oper')
        self._issueCmd('O')
        self._issueCmd('P')
        self._issueCmd('T,T,T,T')
        self._issueCmd('H')
        self._issueCmd('D') 
        self._issueCmd(' ')        
        # Set each run case
        for i in range(len(FCs)):
            FC=FCs[i]
            if i>0:
                self._issueCmd('+')
            self._issueCmd('N')
            self._issueCmd(FC.name)
            self._issueCmd('M')
            self._issueCmd('M')
            self._issueCmd(FC.mass*FC.loadFactor)
            self._issueCmd('MN')
            self._issueCmd(FC.Mach)
            self._issueCmd('V')
            self._issueCmd(FC.velocity)
            self._issueCmd('D')
            self._issueCmd(FC.density)
            self._issueCmd('G')
            self._issueCmd(FC.g)
            self._issueCmd('CD')
            self._issueCmd(FC.Cd0)
            self._issueCmd('IX')
            self._issueCmd(FC.inertia[0])
            self._issueCmd('IY')
            self._issueCmd(FC.inertia[1])
            self._issueCmd('IZ')
            self._issueCmd(FC.inertia[2])
            self._issueCmd('X')
            self._issueCmd(FC.CG[0])
            self._issueCmd('Y')
            self._issueCmd(FC.CG[1])
            self._issueCmd('Z')
            self._issueCmd(FC.CG[2])
            self._issueCmd(' ')
            self._issueCmd('C1')
            self._issueCmd('V')
            self._issueCmd(FC.velocity)
            self._issueCmd('M')
            self._issueCmd(FC.mass*FC.loadFactor)
            self._issueCmd('D')
            self._issueCmd(FC.density)
            self._issueCmd('G')
            self._issueCmd(FC.g)
            self._issueCmd('X')
            self._issueCmd(FC.CG[0])
            self._issueCmd('Y')
            self._issueCmd(FC.CG[1])
            self._issueCmd('Z')
            self._issueCmd(FC.CG[2]) 
            self._issueCmd(' ')
            self._issueCmd('D1')
            self._issueCmd('D1')
            self._issueCmd(FC.flap)
            self._issueCmd('D2')
            self._issueCmd('D2')
            self._issueCmd(FC.aileron)
            self._issueCmd('D4')
            self._issueCmd('D4')
            self._issueCmd(FC.rudder)
            if FC.pitchTrim:
                self._issueCmd('D3')
                self._issueCmd('PM')
                self._issueCmd(FC.CmTrim)
            else:
                self._issueCmd('A')
                self._issueCmd('A')
                self._issueCmd(FC.alpha)
                self._issueCmd('B')
                self._issueCmd('B')
                self._issueCmd(FC.beta)
                self._issueCmd('D3')
                self._issueCmd('D3')
                self._issueCmd(FC.elevator)
        # Run each case 
        for i in range(1,len(FCs)+1,1):
            self._issueCmd(i)
            self._issueCmd('x')
            self._issueCmd('ST')        
            self._issueCmd(' ')
        # Quit AVL
        self._issueCmd(' ')
        self._issueCmd('mode')
        for i in range(1,len(FCs)+1,1):
            self._issueCmd(i)
            self._issueCmd('M')
            self._issueCmd('M')
            self._issueCmd(FC.mass*FC.loadFactor)
            self._issueCmd('MN')
            self._issueCmd(FC.Mach)
            self._issueCmd('V')
            self._issueCmd(FC.velocity)
            self._issueCmd('D')
            self._issueCmd(FC.density)
            self._issueCmd('G')
            self._issueCmd(FC.g)
            self._issueCmd('CD')
            self._issueCmd(FC.Cd0)
            self._issueCmd('IX')
            self._issueCmd(FC.inertia[0])
            self._issueCmd('IY')
            self._issueCmd(FC.inertia[1])
            self._issueCmd('IZ')
            self._issueCmd(FC.inertia[2])
            self._issueCmd('X')
            self._issueCmd(FC.CG[0])
            self._issueCmd('Y')
            self._issueCmd(FC.CG[1])
            self._issueCmd('Z')
            self._issueCmd(FC.CG[2]) 
            self._issueCmd(' ')
            self._issueCmd('S')
            self._issueCmd(' ')
        self._issueCmd(' ')
        self._issueCmd('quit')
    def _processOutput(self):
        coefficients_raw=self._splitResults(self.rawOutput,'Enter filename, or <return> for screen output   s>','(  > 1 if spirally stable )')
        forces_raw =self._splitResults(self.rawOutput,'(referred to Sref,Cref,Bref about Xref,Yref,Zref)','Control Hinge Moments')        
        matrix_raw =self._splitResults(self.rawOutput,'u         w         q        the',': \n  Run-case parameters for eigenmode analyses ...')        
        hingeMoments_raw=self._splitResults(self.rawOutput,'Control Hinge Moments',' ---------------------------------------------------------------\n\n Operation of run case')               
        self._setCoefficients(coefficients_raw)
        self._setDerivs(coefficients_raw)
        self._setMatrix(matrix_raw)
        self._setForces(forces_raw)
        self._setHingeMoments(hingeMoments_raw)
    def _setStabResults(self):
        for i,result in enumerate(self.results):
            self.results[i].dynamicStability=dynamicStability(result.sysMatrix)
    def _setHingeMoments(self,hingeMoments_raw):
        for i in range(len(hingeMoments_raw)):
            self.results[i].hingeMoments.flap=self._selNum(hingeMoments_raw[i],'flap     ',20)  
            self.results[i].hingeMoments.aileron=self._selNum(hingeMoments_raw[i],'aileron  ',20)  
            self.results[i].hingeMoments.elevator=self._selNum(hingeMoments_raw[i],'elevator ',20)  
            self.results[i].hingeMoments.rudder=self._selNum(hingeMoments_raw[i],'rudder   ',20)  
    def _setForces(self,forces_raw):
        for i in range(len(forces_raw)):
            n1=forces_raw[i].find('Strip Forces referred to Strip Area, Chord')+42
            n2=forces_raw[i].find('  Surface # 2')
            s=StringIO(forces_raw[i][n1:n2])
            self.results[i].force_wing=numpy.genfromtxt(s,skip_header=2,skip_footer=0)
    def _setMatrix(self,matrix_raw):
        for i in range(len(matrix_raw)):
            s=StringIO(matrix_raw[i])
            self.results[i].sysMatrix=numpy.genfromtxt(s,skip_header=1,skip_footer=1)
            self.results[i].eigenvals=numpy.linalg.eig(self.results[i].sysMatrix[:,0:12])      
    def _setCoefficients(self,coeff_raw):
        for i in range(len(coeff_raw)):
            self.results[i].Cd0  =self.flightConditions[i].Cd0
            self.results[i].alpha=self._selTxt(coeff_raw[i],'Alpha =','pb/2V')
            self.results[i].beta =self._selTxt(coeff_raw[i],'Beta  =','qc/2V')
            self.results[i].Mach =self._selTxt(coeff_raw[i],'Mach  =','rb/2V')
            self.results[i].Sref =self._selTxt(coeff_raw[i],'Sref =','Cref =')
            self.results[i].Cref =self._selTxt(coeff_raw[i],'Cref =','Bref =')
            self.results[i].Bref =self._selTxt(coeff_raw[i],'Bref =','    \n  Xref')
            self.results[i].ARref=self.results[i].Bref**2/self.results[i].Sref
            self.results[i].coefficients.CX=self._selTxt(coeff_raw[i],'CXtot =','Cltot =')
            self.results[i].coefficients.CY=self._selTxt(coeff_raw[i],'CYtot =','Cmtot =')
            self.results[i].coefficients.CZ=self._selTxt(coeff_raw[i],'CZtot =','Cntot =')
            self.results[i].coefficients.Cl=self._selTxt(coeff_raw[i],'Cltot =','Cl\'tot =')
            self.results[i].coefficients.Cm=self._selTxt(coeff_raw[i],'Cmtot =','\n  CZtot =')
            self.results[i].coefficients.Cn=self._selTxt(coeff_raw[i],'Cntot =','Cn\'tot =')
            self.results[i].coefficients.CL=self._selTxt(coeff_raw[i],'CLtot =','\n  CDtot')
            self.results[i].coefficients.CD=self._selTxt(coeff_raw[i],'CDtot =','\n  CDvis')
            self.results[i].coefficients.CDind=self._selTxt(coeff_raw[i],'CDff  =','\| Trefftz\n  CYff ')
            self.results[i].k=self.results[i].coefficients.CDind/(self.results[i].coefficients.CL**2)         
            self.results[i].e=1.0/(self.results[i].k*numpy.pi*self.results[i].ARref)            
            self.results[i].a=self._selTxt(coeff_raw[i],'CLa =','CLb =')            
            self.results[i].xNP=self._selTxt(coeff_raw[i],'Neutral point  Xnp =','\n\n Clb Cnr')
            self.results[i].CL0=self.results[i].coefficients.CL-self.results[i].a*numpy.pi/180*self.results[i].alpha
            self.results[i].elevator=self._selNum(coeff_raw[i],'elevator        =',10)
    def _setDerivs(self,coeff_raw):
        for i in range(len(coeff_raw)):
            self.results[i].derivs.CLa=self._selTxt(coeff_raw[i],'CLa =','CLb =')            
            self.results[i].derivs.CYa=self._selTxt(coeff_raw[i],'CYa =','CYb =')
            self.results[i].derivs.Cla=self._selTxt(coeff_raw[i],'Cla =','Clb =')
            self.results[i].derivs.Cma=self._selTxt(coeff_raw[i],'Cma =','Cmb =')
            self.results[i].derivs.Cna=self._selTxt(coeff_raw[i],'Cna =','Cnb =')
            
            self.results[i].derivs.CLb=self._selNum(coeff_raw[i],'CLb =',11)
            self.results[i].derivs.CYb=self._selNum(coeff_raw[i],'CYb =',11)            
            self.results[i].derivs.Clb=self._selNum(coeff_raw[i],'Clb =',11)
            self.results[i].derivs.Cmb=self._selNum(coeff_raw[i],'Cmb =',11)
            self.results[i].derivs.Cnb=self._selNum(coeff_raw[i],'Cnb =',11)            
            
            self.results[i].derivs.CLp=self._selTxt(coeff_raw[i],'CLp =','CLq =')      
            self.results[i].derivs.CYp=self._selTxt(coeff_raw[i],'CYp =','CYq =')      
            self.results[i].derivs.Clp=self._selTxt(coeff_raw[i],'Clp =','Clq =')      
            self.results[i].derivs.Cmp=self._selTxt(coeff_raw[i],'Cmp =','Cmq =')      
            self.results[i].derivs.Cnp=self._selTxt(coeff_raw[i],'Cnp =','Cnq =')                  
            
            self.results[i].derivs.CLp=self._selTxt(coeff_raw[i],'CLq =','CLr =')  
            self.results[i].derivs.CLq=self._selTxt(coeff_raw[i],'CYq =','CYr =')  
            self.results[i].derivs.CYq=self._selTxt(coeff_raw[i],'Clq =','Clr =')  
            self.results[i].derivs.Cmq=self._selTxt(coeff_raw[i],'Cmq =','Cmr =')  
            self.results[i].derivs.Cnq=self._selTxt(coeff_raw[i],'Cnq =','Cnr =')  

            self.results[i].derivs.CLr=self._selNum(coeff_raw[i],'CLr =',11)
            self.results[i].derivs.CYr=self._selNum(coeff_raw[i],'CYr =',11)
            self.results[i].derivs.Clr=self._selNum(coeff_raw[i],'Clr =',11)
            self.results[i].derivs.Cmr=self._selNum(coeff_raw[i],'Cmr =',11)
            self.results[i].derivs.Cnr=self._selNum(coeff_raw[i],'Cnr =',11)

            self.results[i].derivs.CLdf=self._selTxt(coeff_raw[i],'CLd1 =','CLd2 =')
            self.results[i].derivs.CYdf=self._selTxt(coeff_raw[i],'CYd1 =','CYd2 =')
            self.results[i].derivs.Cldf=self._selTxt(coeff_raw[i],'Cld1 =','Cld2 =')
            self.results[i].derivs.Cmdf=self._selTxt(coeff_raw[i],'Cmd1 =','Cmd2 =')
            self.results[i].derivs.Cndf=self._selTxt(coeff_raw[i],'Cnd1 =','Cnd2 =')
            self.results[i].derivs.CDff=self._selTxt(coeff_raw[i],'CDffd1 =','CDffd2 =')
            self.results[i].derivs.edf =self._selTxt(coeff_raw[i],'ed1 =','ed2 =')
            
            self.results[i].derivs.CLda=self._selTxt(coeff_raw[i],'CLd2 =','CLd3 =')
            self.results[i].derivs.CYda=self._selTxt(coeff_raw[i],'CYd2 =','CYd3 =')
            self.results[i].derivs.Clda=self._selTxt(coeff_raw[i],'Cld2 =','Cld3 =')
            self.results[i].derivs.Cmda=self._selTxt(coeff_raw[i],'Cmd2 =','Cmd3 =')
            self.results[i].derivs.Cnda=self._selTxt(coeff_raw[i],'Cnd2 =','Cnd3 =')
            self.results[i].derivs.CDfa=self._selTxt(coeff_raw[i],'CDffd2 =','CDffd3 =')
            self.results[i].derivs.eda =self._selTxt(coeff_raw[i],'ed2 =','ed3 =')
           
            self.results[i].derivs.CLde=self._selTxt(coeff_raw[i],'CLd3 =','CLd4 =')
            self.results[i].derivs.CYde=self._selTxt(coeff_raw[i],'CYd3 =','CYd4 =')
            self.results[i].derivs.Clde=self._selTxt(coeff_raw[i],'Cld3 =','Cld4 =')
            self.results[i].derivs.Cmde=self._selTxt(coeff_raw[i],'Cmd3 =','Cmd4 =')
            self.results[i].derivs.Cnde=self._selTxt(coeff_raw[i],'Cnd3 =','Cnd4 =')
            self.results[i].derivs.CDfe=self._selTxt(coeff_raw[i],'CDffd3 =','CDffd4 =')
            self.results[i].derivs.ede =self._selTxt(coeff_raw[i],'ed3 =','ed4 =')
            
            self.results[i].derivs.CLdr=self._selNum(coeff_raw[i],'CLd4 =',11)
            self.results[i].derivs.CYdr=self._selNum(coeff_raw[i],'CYd4 =',11)
            self.results[i].derivs.Cldr=self._selNum(coeff_raw[i],'Cld4 =',11)
            self.results[i].derivs.Cmdr=self._selNum(coeff_raw[i],'Cmd4 =',11)
            self.results[i].derivs.Cndr=self._selNum(coeff_raw[i],'Cnd4 =',11)
            self.results[i].derivs.CDfr=self._selNum(coeff_raw[i],'CDffd4 =',11)
            self.results[i].derivs.edr =self._selNum(coeff_raw[i],'ed4 =',11)
            
            self.results[i].derivs.CNa = self._get_normal_force_coef(i,self.flightConditions[i].Cd0)
    def _selNum(self,source,fromStr,num):
        a=self._findall(source,fromStr)
        val=source[a[0][1]:a[0][1]+num]
        return float(val)
    def _selTxt(self,source,fromStr1,toStr2):
        a=self._findall(source,fromStr1)
        b=self._findall(source,toStr2)
        assert len(a)==1
        assert len(b)==1
        num=source[a[0][1]:b[0][0]]
        return float(num)
    def _splitResults(self,source,fromStr1,toStr2):
        a=self._findall(source,fromStr1)
        b=self._findall(source,toStr2)
        assert len(a)==len(b)
        out=list()
        for i in range(len(a)):
            out.append(source[a[i][1]:b[i][0]])
        return out
    def _findall(self,source,string):
        return [(a.start(), a.end()) for a in list(re.finditer(string, source))]
    def _get_normal_force_coef(self,i,CD0):
        alpha = numpy.radians(self.results[i].alpha)
        cosA = numpy.cos(alpha)
        sinA = numpy.sin(alpha)
        k = self.results[i].k
        a = self.results[i].a
        CLa = self.results[i].derivs.CLa
        CL = self.results[i].CL0 + alpha*a
        CNa = a*cosA - CL*sinA + 2.*k*CL*CLa*sinA + (CD0+k*CL*CL)*cosA
        return CNa

def print_header(header):
    print '\n'+header
    print '='*len(header)

class parasiteDrag():
    """
    Class for estimating aircraft parasite drag
    
    Parameters
    ----------
    
    aircraft : aircraft
        aircraft configuration object
    """
    def __init(self,aircraft):
        self.aircraft=aircraft
    def getCd0(self,V,rho):
        """
        Executes AeroCD0 skin friction / form factor drag analysis code
        Returns parasite drag (Cd0)
        
        Parameters
        ----------
        
        V : float, m/sec
            velocity
        rho : float, kg/m**3
            density
        """
        alt = fc.get_density_altitude(rho)
        Cd0= self.aircraft.get_drag(V,alt)
        return Cd0


class _results():
    """
    Private class for holding the results of the aerodynamics analysis  
    Holds a structured set of all results including basic aerodynamic
    coefficients, stability derivatives, stability analysis results, 
    hinge moments, and force distributions
    """
    def __init__(self):
        self.alpha=0.0
        self.beta=0.0
        self.Mach=0.0
        self.Sref=0.0
        self.Cref=0.0
        self.Bref=0.0
        self.ARref=0.0
        self.k=0.0
        self.e=0.0
        self.a=0.0
        self.xNP=0.0
        self.coefficients=_coefficients()
        self.derivs=_derivs()
        self.hingeMoments=_hingeMoments()
    def display(self):

        print_header('BASIC TERMS')
        for attr, value in self.__dict__.iteritems():
            if type(value)==type(0.0):
                print '{:<8} = {:<+12.5f}  '.format(attr,value)
        print_header("COEFFICIENTS")
        self.coefficients.display()
        print_header("DERIVATIVES")
        self.derivs.display()
        print_header('HINGE MOMENTS')
        self.hingeMoments.display()
        


class _hingeMoments():
    def __init__(self):
        self.flap    =0.0
        self.aileron =0.0
        self.elevator=0.0
        self.rudder  =0.0
        self.Mflap = 0.0
        self.Maileron = 0.0
        self.Melevator = 0.0
        self.Mrudder = 0.0
    def display(self):        
        for attr, value in self.__dict__.iteritems():
            print '{:<8} = {:<+12.5f}  '.format(attr,value)


class _coefficients():
    def __init__(self):
        self.CX=0.0
        self.CY=0.0
        self.CZ=0.0
        self.Cl=0.0
        self.Cm=0.0
        self.Cn=0.0
        self.CL=0.0
        self.CD=0.0
        self.CDind=0.0
    def display(self):        
        for attr, value in self.__dict__.iteritems():
            print '{:<8} = {:<+12.5f}  '.format(attr,value)


class _derivs():
    def __init__(self):
        self.CLa=0.0
        self.CYa=0.0
        self.Cla=0.0
        self.Cma=0.0
        self.Cna=0.0
        
        self.CLb=0.0
        self.CYb=0.0
        self.Clb=0.0
        self.Cmb=0.0
        self.Cnb=0.0
        
        self.CLp=0.0 
        self.CYp=0.0
        self.Clp=0.0
        self.Cmp=0.0
        self.Cnp=0.0
        
        self.CLp=0.0
        self.Clq=0.0
        self.CYq=0.0
        self.Cmq=0.0
        self.Cnq=0.0

        self.CLr=0.0
        self.CYr=0.0
        self.Clr=0.0
        self.Cmr=0.0
        self.Cnr=0.0

        self.CLdf=0.0
        self.CYdf=0.0
        self.Cldf=0.0
        self.Cmdf=0.0
        self.Cndf=0.0
        self.CDdf=0.0
        self.edf =0.0
        
        self.CLda=0.0
        self.CYda=0.0
        self.Clda=0.0
        self.Cmda=0.0
        self.Cnda=0.0
        self.CDfa=0.0
        self.eda =0.0
       
        self.CLde=0.0
        self.CYde=0.0
        self.Clde=0.0
        self.Cmde=0.0
        self.Cnde=0.0
        self.CDde=0.0
        self.ede =0.0
        
        self.CLdr=0.0
        self.CYdr=0.0
        self.Cldr=0.0
        self.Cmdr=0.0
        self.Cndr=0.0
        self.CDdr=0.0
        self.edr =0.0
    def display(self):        
        for attr, value in self.__dict__.iteritems():
            print '{:<8} = {:<+12.5f}  '.format(attr,value)
            #print '\\texttt{%s} & $C_'%attr


class dynamicStability():
    """
    Interprets the aircraft system matrix computed by AVL.  Calculates the 
    eigenvalues and the associated frequency and damping associated with each
    considered dynamic mode.  Catagorizes the results according to londitudinal
    and lateral modes
    
    Parameters
    ----------
    
    sysMatrix : matrix
        full aircraft matrix including control vectors
    aircraftMatrix:  matrix
        aircraft derivative matrix without control vectors
    longMatrix : matrix
        Decoupled longitudinal matrix
    latMatrix : matrix
        Decoupled lateral matrix
    """
    def __init__(self,sysMatrix):
        self.sysMatrix      =sysMatrix
        self.aircraftMatrix =sysMatrix[:,0:11]
        self.controlMatrix  =sysMatrix[:,12:]
        self.longMatrix     =sysMatrix[0:4,0:4]
        self.latMatrix      =sysMatrix[4:8,4:8]
        #print self.longMatrix
        longParams          =self._getStabParams(self.longMatrix)
        latParams           =self._getStabParams(self.latMatrix)
        # Longitudinal Params
        self.phugoid     =longParams[longParams.w.argmin()]
        self.shortPeriod =longParams[longParams.w.argmax()]
        # Lateral Params
        self.dutchRoll   =latParams[numpy.array(latParams.w<>0).argmax()]
        spiralRollParams =latParams[latParams.w==0]
        self.spiral      =spiralRollParams[spiralRollParams.n.argmin()]
        self.roll        =spiralRollParams[spiralRollParams.n.argmax()]
    def _getStabParams(self,M):
        params=_stabilityParams()
        c=numpy.log(2)/numpy.pi/2
        eig,v=numpy.linalg.eig(M)
        eig=eig[eig<>0]
        eig=eig[eig.imag>=0]                 
        n=eig.real       
        w =numpy.abs(eig.imag)
        ii=w<>0
        jj=w==0    
        Wn=(w**2+n**2)**.5
        Z =-n/Wn
        T=numpy.zeros([len(n)])
        T[ii] =2*numpy.pi/w[ii]
        T[jj] =numpy.inf
        N =numpy.abs(c*w/n)
        T2=numpy.log(2)/numpy.abs(n)
        params.T=T
        params.Z=Z
        params.N=N
        params.w=w
        params.Wn=Wn
        params.n=n
        params.T2=T2
        return params

class _stabilityParams():
    def __init__(self):
        self.T =numpy.array([])
        self.Z =numpy.array([])
        self.N =numpy.array([])
        self.w =numpy.array([])
        self.Wn=numpy.array([])
        self.n =numpy.array([])
        self.T2=numpy.array([])
    def __getitem__(self,key):
        out=_stabilityParams()
        out.T  =self.T[key]
        out.Z  =self.Z[key]
        out.N  =self.N[key]
        out.w  =self.w[key]
        out.Wn =self.Wn[key]
        out.n  =self.n[key]
        out.T2 =self.T2[key]
        return out
    def display(self):        
        for attr, value in self.__dict__.iteritems():
            print '{:<8} = {:<12s}  '.format(attr,str(value))

# --- debug section ---
def example1():
    # run trim analysis
    import aircraft
    ac = aircraft.load('V0510')
    fc = flightConditions()
    CG = numpy.array([1.7,0,0])
    I = numpy.array([500., 1000, 1000])
    M = 600.
    V = 50.
    density = 1.2255
    Cd0 = ac.get_drag(V,0.0)
    fc.addTrimmedFlightCondition('cruise1',M,CG,I,V,density,CmTrim=0.0,Cd0=Cd0)
    a3d = aero3d_VLM(ac)
    rslt = a3d.runVLM(fc,False)
    print 'GENERAL'
    rslt.results[0].display()
    print 'PHUGOID'
    rslt.results[0].dynamicStability.phugoid.display()
    print 'SHORT PERIOD'
    rslt.results[0].dynamicStability.shortPeriod.display()
    print 'DUTCH ROLL'
    rslt.results[0].dynamicStability.dutchRoll.display()
    print 'SPIRAL'
    rslt.results[0].dynamicStability.spiral.display()
    print 'ROLL'
    rslt.results[0].dynamicStability.roll.display()

def example2():
    # run trim using modified function
    import aircraft
    ac = aircraft.load('V0510')
    aero = aerodynamics(ac)
    aero.update(65.,0.8,3.0)
    aero.results.display()
    ac.display()

def runTest2():
    import aircraft
    AC =aircraft.load("Xtail_sample")
    FC =flightConditions()
    CG =numpy.array([0.65,0,0])
    I  =numpy.array([500.0,1000.0,1000.0])
    M  =75.0
    V  =50.0
    rho=1.225
    CD0=AC.drag.total.get_total_drag()
    FC.addTrimmedFlightCondition('cruise1',M,CG,I,V,rho,CmTrim=0.0,Cd0=CD0)
    FC.flightConditionList[0].flap = 0.0
    #FC.addTrimmedFlightCondition('cruise2',M/2,CG,I,V/2,rho,Cd0=CD0)
    a3d=aero3d_VLM(AC)
    R=a3d.runVLM(FC,False)
    R.results[0].display()
    print_header('PHUGOID')
    R.results[0].dynamicStability.phugoid.display()
    print_header('SHORT PERIOD')
    R.results[0].dynamicStability.shortPeriod.display()
    print_header('DUTCH ROLL')
    R.results[0].dynamicStability.dutchRoll.display()
    print_header('SPIRAL')
    R.results[0].dynamicStability.spiral.display()
    print_header('ROLL')
    R.results[0].dynamicStability.roll.display()
    AC.display()

def run_test():
    import aircraft
    ac = aircraft.load('V0510')
    fc = flightConditions()
    CG = ac.get_CG(False)
    I = ac.get_inertia(False)
    M = ac.get_mass_total(False)
    V = 1.0
    rho = 1.2255
    Cd0 = ac.get_drag(V,0.0)
    fc.addTrimmedFlightCondition('trim',M, CG, I, V, rho)
    fc.flightConditionList[0].setRudder(0.)
    aero3D = aero3d_VLM(ac)
    rslt = aero3D.runVLM(fc)
    rslt.results[0].display()

def run_test2():
    import aircraft
    ac = aircraft.load('V0510')
    aero = aerodynamics(ac)
    aero.update(65.,0.8,3.0)
    aero.results.display()
    ac.display()

def run_test3():
    import aircraft
    import matplotlib.pyplot as plt
    ac = aircraft.load('V0510')
    fc = flightConditions()
    CG = ac.get_CG()
    I = ac.get_inertia()
    m = ac.get_mass_total()
    V = 50.0
    rho = 1.2255
    CD0 = ac.get_drag(V,0.0)
    ail = [-10,0,2,5,7,10]
    ailDisp = list()
    r = 15
    for a in ail:
        fc.add_flight_condition('untrimmed',m,CG,I,V,rho,CD0=CD0,aileron=a,rudder=r)
        ailDisp.append(a)
    aero3d = aero3d_VLM(ac)
    rslt = aero3d.runVLM(fc)
    cl = [rslt.results[i].coefficients.Cl for i in range(len(ail))]
    plt.plot(ailDisp,cl,'ro-')
    plt.show()

def calc_hinges():
    import aircraft
    ac = aircraft.load('V0510')
    fc = flightConditions()
    CG = [1.824,0.0,-0.125]
    I = ac.get_inertia()
    m = 620.0
    V = 75.0
    rho = 1.2255
    CD0 = ac.get_drag(V,0.0)
    n = 1.0
    fc.add_flight_condition('hinge',m,CG,I,V,rho,n,CD0,alpha=20.0,elevator=30.0)
    aero3d = aero3d_VLM(ac)
    rslt = aero3d.runVLM(fc)
    rslt.results[0].hingeMoments.display()

if __name__=="__main__":
    example1()