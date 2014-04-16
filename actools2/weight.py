# -*- coding: utf-8 -*-
"""
Set of classes and functions for easy weight and balance analysis. 
Contains equations for light aircraft mass estimation. Can be easy extended 
for different types of aircraft. Analysis uses constants.xls.
Mass calculation routines are available from aircraft.py.
This module can be used for explicit weight and balance calculation. For example 
stability and control calculation with different payload and fuel mass. 
"""
import numpy as np
import constants
import quantities as qu
#import FlightConditions as fc
from time import ctime

def get_general_aviation_mass(aircraft):
    """
    Function returns aircraft mass object that contains lists of the mass 
    components. Mass is obtained using textbook methods by D.Raymer.
    """
    mass = GeneralAviationMass(aircraft)
    return mass.aircraftMass

class MassComponent:
    """
    Describing mass component by parameters described below
    
    Parameters
    ----------
    
    name : string
        name of the mass component
    mass : float, kg
        mass of the component
    coords : 3d array float, m
        CG coordinates of the mass component
    inertia : 3d array float, kg*m^2
        moment of inertia of the mass component
    
    Examples
    --------
    
    >>> pilot   = MassComponent('pilot',   86.0, [1.86,-0.25,-0.075])
    >>> baggage = MassComponent('baggage', 70.0, [2.85,    0,-0.125])
    >>> baggage.display(header=True)
    Name              |Mass,kg   |X,m       |Y,m       |Z,m       |
    baggage               70.00     2.8500     0.0000    -0.1250
    """
    numberOfItems = 0
    def __init__(self,name='',mass=0.0,CG=np.zeros(3),inertia=np.zeros(3)):
        MassComponent.numberOfItems += 1
        if name=='':
            self.name = 'component %d'%MassComponent.numberOfItems
        else:
            self.name = str(name)
        self.mass = mass
        self.coords  = np.array(CG)
        self.inertia = np.array(inertia)
    def display(self,header=True,disp=False):
        """
        prints out data about given item
        
        Parameters
        ----------

        header : bool, optional
            if True then header of the table is displayed
        """
        report = ''
        if header:
            line = '{0:18}|{1:10}|{2:10}|{3:10}|{4:10}|'
            report += line.format('Name','Mass,kg','X,m','Y,m','Z,m')+'\n'
        line = '{0:18} {1:8.2f} {2:10.4f} {3:10.4f} {4:10.4f}'
        x, y, z = self.coords[0], self.coords[1], self.coords[2]
        report += line.format(self.name,self.mass,x,y,z)+'\n'
        if disp:
            print report
        else:
            return report


class MassList:
    """
    List of mass components.
    
    Parameters
    ----------
    
    name : string
        name of the mass list. For example *payload* or *airframe*.


    Examples
    --------
    
    >>> payload = MassList('payload')
    >>> payload.add_item('passenger',86,[1.86,0.25,-0.075],[0,0,0])
    >>> baggage    = MassComponent('baggage',  70.0,[2.85,0,-0.125])
    >>> payload.add_component(baggage)
    >>> print payload.get_total_mass()
    156.0
    >>> print payload.get_CG()
    [ 2.30423077  0.13782051 -0.0974359 ]
    """
    n = 0
    def __init__(self,name=''):
        MassList.n += 1
        if name=='':
            self.name = 'mass list %d'%MassList.n
        else:
            self.name = name
        self.items = list()
        self.totalMass = 0.0
        self.CG = np.zeros(3)
        self.names = list()

    def __getitem__(self,k):
        return self.items[k]

    def __call__(self):
        self.update_totals()
        return self.totalMass

    def add_item(self,name,mass,CG=np.zeros(3),inertia=np.zeros(3)):
        """
        Add new item by set of parameters to the list. After the item is added 
        total CG and mass are updated.
        
        Parameters
        ----------
        name : string
            name of the mass item
        mass : float
            mass of new item
        CG : array
            center of gravity of new item in format [x,y,z] in meters
        inertia : array
            inertial of new mass item in format [Ixx, Iyy, Izz]
        
        Examples
        --------
        
        >>> payload = MassList('payload')
        >>> payload.add_item('passenger',86,[1.86,0.25,-0.075],[0,0,0])
        """
        newItem = MassComponent()
        newItem.name = name
        self.names.append(name)
        newItem.mass = mass
        newItem.coords = CG
        newItem.inertia = inertia
        self.items.append(newItem)
        self.update_totals()

    def add_component(self,component):
        """
        Add new mass component. After the item is added total CG and mass are 
        updated.
        
        Paramters
        ---------
        component : MassComponent
        """
        self.names.append(component.name)
        self.items.append(component)
        self.update_totals()

    def add_mass_list(self,newMassList):
        """
        Add all items from newMassList to current one. After the item is added 
        total CG and mass are updated.

        Parameters
        ----------
        
        newMassList : MassList
            mass list to be added
        """
        for newItem in newMassList:
            self.items.append(newItem)
        self.update_totals()

    def get_item_by_name(self,name,remove=False):
        """
        returns mass component by name
        """
        try:
            k=self.names.index(name)
            item = self[k]
            if remove:
                self.remove_item(name)
            return item
        except:
            print 'Specified mass item not found'

    def get_item_mass_by_name(self,name):
        """ returns mass of the component by name
        """
        item = self.get_item_by_name(name)
        return item.mass

    def get_item_cg_by_name(self,name):
        """ returns CG of the component by name
        """
        item = self.get_item_by_name(name)
        return item.coords

    def get_item_inertia_by_name(self,name):
        """ returns inertia of the component by name
        """
        item = self.get_item_by_name(name)
        return item.inertia

    def check_item(self,name):
        """
        Returns True if item with name specified exists in mass list
        
        Parameters
        ----------
        name : string
            name of the item to be checked
        """
        try:
            self.names.index(name)
            exist = True
        except:
            exist = False
        return exist
        
    def remove_item(self,name):
        """ remove mass component with specified name from mass list
        """
        try:
           k=self.names.index(name)
           del self.items[k]
        except:
            print 'Specified mass item not found: %s'%name
    def update_item(self,name,mass=None,cg=None):
        """ update mass and CG of the component with given name
        
        Parameters
        ----------
        
        name : string
            name of the component to update
        mass : float, kg
            new mass of the component
        cg : array, m
            new center of gravity of the component
        """
        if not mass==None:
            self.update_item_mass(name,mass)
        if not cg==None:
            self.update_item_cg(name,cg)
    def update_item_mass(self,name,mass):
        try:
            k = self.names.index(name)
            self.items[k].mass = float(mass)
            self.update_totals()
        except:
            print 'Specified mass item not found'
    def update_item_cg(self,name,cgX,cgY,cgZ):
        try:
            k = self.names.index(name)
            self.items[k].coords = np.array([cgX,cgY,cgZ])
            self.update_totals()
        except:
            print 'Specified mass item not found: %s'%name
    def __len__(self):
        return len(self.items)
    def update_totals(self):
        totalMass = 0.0
        CG = np.zeros(3)
        for item in self.items:
            totalMass += item.mass
            CG += item.mass*item.coords
        self.totalMass = totalMass
        if totalMass>0.0:
            self.CG = CG / totalMass
        else:
            self.CG = np.zeros(3)
    def get_total_mass(self):
        """ returns total mass of the list
        """
        self.update_totals()
        return self.totalMass
    def get_CG(self,upd=True):
        """ returns total CG of the list
        """
        if upd:
            self.update_totals()
        return self.CG
    def get_name_list(self):
        outList = list()
        for item in self.items:
            outList.append(item.name)
        return outList
    def _report(self):
        """
        Displays mass components breakdown report: name, mass, CG
        """
        limiter1 = '-'*63
        limiter2 = '='*63
        report = limiter2+'\n'
        report += '%s mass components breakdown\n'%self.name
        report += limiter2+'\n'
        lineFormat = '{0:20}|{1:8}|{2:10}|{3:10}|{4:10}|\n'
        report += lineFormat.format('Name','Mass,kg','CG:X,m','CG:Y,m','CG:Z,m')
        report += limiter1+'\n'
        for item in self.items:
            report += item.display(False)
        report += limiter1+'\n'
        lineFormat = '{0:18} {1:8.2f} {2:10.4f} {3:10.4f} {4:10.4f}\n'
        report += lineFormat.format('TOTAL',self.totalMass,
                                self.CG[0],self.CG[1],self.CG[2])
        report += limiter1 +'\n'
        return report
    def display(self):
        """ prints out tabulated information about mass list
        """
        print self._report()
    def save_txt(self,filePath='',acname=''):
        """
        saves mass list in tabulated form to text file
        
        Parameters
        ----------
        
        filePath : string
            output text file path
        acname : string
            name of the aircraft that will be written in file header. If not 
            specified then empty field will be written.
        """
        if filePath=='':
            filePath = '%s_mass_report.txt'%self.name
        if acname!='':
            filePath = '%s_'%acname + filePath
        lineFormat = '{0:20}:{1:20}\n'
        fid = open(filePath,'wt')
        fid.write(lineFormat.format('Generated on',ctime()))
        fid.write(self._report())
        fid.close()


class Fuel(MassList):
    """ Class describing fuel tank as point mass with given CG location and
    fuel mass. Provides methods for setting fuel mass in different tanks or total.
    """
    def _init_fuel(self):
        self.density = 750.0
        self.maxMass = 0.0
        self.volume = 0.0
        self._dict = {'right':'Fuel tank Right',
                      'left':'Fuel tank Left',
                      'all':'all',
                      'both':'all',
                      'total':'total'}
    def set_fuel_mass(self,fuelMass,tankName='total'):
        n = len(self.items)
        tankName = self._dict[tankName]
        if tankName=='total':
            fuelPerTank = fuelMass/float(n)
            for i in range(n):
                self.items[i].mass = fuelPerTank
        elif tankName=='all':
            for i in range(n):
                self.items[i].mass = fuelMass
        else:
            self.update_item_mass(tankName,fuelMass)
    def set_fuel_mass_burned(self,burnedFuelMass,tankName='total'):
        mass = self.get_fuel_mass_total()
        newMass = mass - burnedFuelMass
        self.set_fuel_mass(newMass,tankName)
    def get_fuel_mass_total(self):
        mass = 0.0
        for item in self.items:
            mass += item.mass
        return mass
    def set_fuel_fraction(self,fuelFraction,tankName):
        n = len(self.items)
        if tankName=='all':
            for i in range(n):
                self.items[i].mass *= fuelFraction
        else:
            m = self.get_item_by_name(tankName).mass
            self.update_item_mass(tankName,m*fuelFraction)


class AircraftMass:
    """
    Class describing set of aircraft mass components. It is used in aircraft.mass
    Contains several major mass lists: airframe, payload and fuel.
    Total mass is stored in separate mass list and should be updated manually 
    after any mass list was changed.
    
    AircraftMass object is created automatically when loading aircraft
    
    >>> import aircraft
    >>> ac = aircraft.load('V0510')
    >>> ac.mass.display() # AircraftMass object is aircraft.mass
    ===============================================================
    Total mass components breakdown
    ===============================================================
    Name                |Mass,kg |CG:X,m    |CG:Y,m    |CG:Z,m    |
    ---------------------------------------------------------------
    Left Wing             37.23     1.9300    -2.3750    -0.3103
    ...
    baggage               70.00     2.8500     0.0000    -0.1250
    Fuel tank Left        47.50     1.8100    -1.4250    -0.3366
    ---------------------------------------------------------------
    TOTAL                581.74     1.8705    -0.0370    -0.1535
    ---------------------------------------------------------------
    
    Update baggage mass
    
    >>> ac.mass.update_item_mass('Payload','baggage',20.0)
    >>> ac.mass.display()
    ===============================================================
    Total mass components breakdown
    ===============================================================
    Name                |Mass,kg |CG:X,m    |CG:Y,m    |CG:Z,m    |
    ---------------------------------------------------------------
    Left Wing             37.23     1.9300    -2.3750    -0.3103
    ...
    baggage               20.00     2.8500     0.0000    -0.1250
    Fuel tank Left        47.50     1.8100    -1.4250    -0.3366
    ---------------------------------------------------------------
    TOTAL                531.74     1.8705    -0.0370    -0.1535
    ---------------------------------------------------------------    
    
    Update fuel mass in right wing fuel tank
    
    >>> ac.mass.update_item_mass('Fuel','Fuel Tank Right',2.5)
    
    Parameters
    ----------
    
    name : string
        aircraft name. Name is used in text reports only.
    """
    def __init__(self,name='noname aircraft'):
        self.aircraftName = name
        self.airframe = MassList('Airframe')
        self.payload = MassList('Payload')
        self.fuel = Fuel('Fuel')
        self.fuel._init_fuel()
        self.total = MassList('Total')
        self.inertia = np.zeros(3)
        self.totalMass = 0.0
        self.totalCG = np.zeros(3)
    def update_total(self):
        """
        Updates total mass list.
        """
        self.total = MassList('Total')
        self.total.add_mass_list(self.airframe)
        self.total.add_mass_list(self.payload)
        self.total.add_mass_list(self.fuel)
        self.total.update_totals()
        self.totalMass = self.total.totalMass
        self.totalCG = self.total.CG
    def display(self):
        """ prints out tabulated information about all mass components
        """
        self.total.display()

    def save_txt(self,filePath=''):
        """ saves total mass list in tabulated form to text file
        """
        self.total.save_txt(filePath,self.aircraftName)

    def __call__(self):
        massTotal = self.total.get_total_mass()
        return massTotal

    def update_item_mass(self,listName,itemName,mass):
        """
        Update item mass in a list.
        
        Paramteres
        ----------
        listName : string
            item will be changed in listName list if exist. If list does not 
            exits nothing will be changed
        itemName : string
        """
        if listName=='Payload':
            self.payload.update_item_mass(itemName,mass)
        elif listName=='Airframe':
            self.airframe.update_item_mass(itemName,mass)
        self.update_total()
        
    def update_item_cg(self,listName,itemName,cgX,cgY,cgZ):
        """ updates CG of the item in specified list
        
        Parameters
        ----------
        
        listName : string
            name of the list (example: payload)
        itemName : string
            name of the item to be modified
        cgX : float
            X coordinate of CG
        cgY : float 
            Y coordinate of CG
        cgZ : float
            Z coordinate of CG
        """
        if listName=='Payload':
            self.payload.update_item_cg(itemName,cgX,cgY,cgZ)
        elif listName=='Airframe':
            self.airframe.update_item_cg(itemName,cgX,cgY,cgZ)
        self.update_total()
    def set_fuel_mass(self,fuelMass,tankName):
        """
        Sets fuel mass in specified fuel tank if several of them exists
        
        Parameters
        ----------
        
        fuelMass : float, kg
            new mass of the fuel
        tankName : string
            fuel tank name
        """
        self.fuel.set_fuel_mass(fuelMass,tankName='total')
        self.update_total()



class BlendedWingBodyMass(object):
    def __init__(self,aircraft):
        self.ac = aircraft
    def analyze(self):
        pass



class GeneralAviationMass:
    def __init__(self,aircraft,overrideColumn=0):
        self.aircraft = aircraft
        self.constants = constants.load('mass',overrideColumn)
        
        self.aircraftMass = aircraft.mass
        self.aircraftMass.airframe = MassList('Airframe')
        #atm                    =fc.ISAtmosphere(aircraft.designGoals.designAltitude)
        self._designQ          =aircraft.designGoals.fc.dynamicPressure
        self._designM          =aircraft.designGoals.fc.Mach
        self._Wdg              =float((aircraft.designGoals.designGrossMass*qu.kg).rescale('pound'))    
        self._Nz               =aircraft.designGoals.designLoadFactor*1.5
        self._q                =float((self._designQ*qu.Pa).rescale('pound_force/ft**2'))
        self._Nen              =float(len(aircraft.engine.engineList))

        self._wing()
        self._hStab()
        self._vStab()
        self._fuselage()
        self._landingGear()
        self._engine()
        self._fuel()
        self._misc()
#        self.emptyMass    =0.0
#        self.payloadMass  =0.0
#        self.totalMass    =0.0
        self.CG=np.zeros(3)
        self.get_inertia()
        self.aircraftMass.update_total()
    def _wing(self):
        #TODO: include equations here
        ac        =self.aircraft
        Sw        =float((ac.wing.area*qu.m**2).rescale('ft**2'))
        Wfw       =float((ac.designGoals.fuelMass*qu.kg).rescale('pound'))
        A         =ac.wing.aspectRatio                                      
        W_swept   =np.radians(ac.wing.sweep)
        taper     =ac.wing.taper
        t_c       =ac.wing.thicknessRatio
        Ww        =0.036*(Sw**0.758)*(Wfw**0.0035)*((A/((np.cos(W_swept))**2))**0.6)*(self._q**0.006)*((taper)**0.04)*(((100*(t_c ))/np.cos(W_swept))**(-0.3))*((self._Nz*self._Wdg)**0.49) 
        Mw        =float((Ww*qu.pound).rescale('kg'))
        CGR       =ac.wing.centroid
        CGR[0] = ac.wing.aapex[0]+0.4*ac.wing.MAC
        CGL       =np.array([1.,-1.,1.])*ac.wing.centroid
        inertia   =np.zeros(3)  #TODO - calculate inertia based on geometry
        if ac.wing.material=='composite':
            Mw *= 0.85
        self.aircraftMass.airframe.add_item('Left Wing',Mw/2.,CGL,inertia)
        self.aircraftMass.airframe.add_item('Right Wing',Mw/2.,CGR,inertia)

    def _hStab(self):
        ac        =self.aircraft
        A         =ac.hStab.aspectRatio
        Sht       =float((ac.hStab.area*qu.m**2).rescale('ft**2'))
        W_swept   =float(ac.hStab.sweep*qu.rad)
        taper_ht  =ac.hStab.taper
        t_c       =ac.hStab.thicknessRatio
        Wht       =0.016*((self._Nz*self._Wdg)**0.414)*(self._q**0.0168)*(Sht**0.896)*(((100*(t_c ))/np.cos(W_swept))**(-0.12))*((A/((np.cos(W_swept))**2))**0.043)*(taper_ht**(-0.02));
        Mht       =float((Wht*qu.pound).rescale('kg'))
        CGR       =ac.hStab.centroid
        CGL       =np.array([1.,-1.,1.])*ac.hStab.centroid
        inertia   =np.zeros(3)  #TODO - calculate inertia based on geometry
        if ac.hStab.material=='composite':
            Mht *= 0.83
        self.aircraftMass.airframe.add_item('Left Hstab',Mht/2.,CGL,inertia)
        self.aircraftMass.airframe.add_item('Right Hstab',Mht/2.,CGR,inertia)
    def _vStab(self):
        ac        =self.aircraft
        A         =ac.vStab.aspectRatio
        Svt       =float((ac.vStab.area*qu.m**2).rescale('ft**2'))
        W_swept   =float(ac.vStab.sweep*qu.rad)
        taper_vt  =ac.vStab.taper
        t_c       =ac.vStab.thicknessRatio
        Ht        =ac.hStab.aapex[2]-ac.vStab.aapex[2]
        Hv        =ac.vStab.aapex[2]+ac.vStab.span
        Ht = abs(Ht)
        Wvt       =0.073*(1.+(0.2*(Ht/Hv)))*((self._Nz*self._Wdg)**0.376)*(self._q**0.122)*(Svt**0.3873)*(((100*(t_c))/np.cos(W_swept))**(-0.49))*((A/((np.cos(W_swept))**2))**0.375)*(taper_vt**(0.039))
        Mvt       =float((Wvt*qu.pound).rescale('kg'))
        CG        =ac.vStab.centroid
        inertia   =np.zeros(3)  #TODO - calculate inertia based on geometry
        if ac.vStab.material=='composite':
            Mvt *= 0.83
        #Mvt *= 2.0 #TODO: find reliable mass calculation equation than this one
        Mvt *= 2.0
        self.aircraftMass.airframe.add_item('vStab',Mvt,CG,inertia)
    def _fuselage(self):
        #FIXME: fuselage mass prediction difference is due to wrong wetted area calculation
        ac        =self.aircraft
        Sf        =float((ac.fuselage.wettedArea*qu.m**2).rescale('ft**2'))
        xMW       =ac.wing.centroid[0]+ac.wing.MAC/4.
        xHT       =ac.hStab.centroid[0]+ac.hStab.MAC/4.
        Lt_si     =np.abs(xHT-xMW)
        Lt        =float((Lt_si*qu.m).rescale('ft'))
        L         =float((ac.fuselage.length*qu.m).rescale('ft'))
        D         =float((ac.fuselage.diameter*qu.m).rescale('ft'))
        Wpress    =0.
        Wfuselage =0.052*(Sf**1.086)*((self._Nz*self._Wdg)**0.177)*(Lt**(-0.051))*((L/D)**(-0.072))*(self._q**0.241)+Wpress
        Mfuselage =float((Wfuselage*qu.pound).rescale('kg'))
        if ac.fuselage.material=='composite':
            Wfuselage *= 0.9
        #CG        =ac.fuselage.centroid
        #TODO: check fuselage center of gravity prediction
        CG = self.constants.getValue('fuseCGratio') * self.aircraft.fuselage.length
        CG = np.array([CG,0.0,0.0])
        inertia   =np.zeros(3)
        self.aircraftMass.airframe.add_item('fuselage',Mfuselage,CG,inertia)
    def _landingGear(self):
        ac        =self.aircraft
        iNose     =np.argmin(ac.landingGear.groundContact_X)
        Nl        =ac.designGoals.designLandingLoadFactor*1.5        
        Wl        =self._Wdg        
        ln        =float((ac.landingGear.strutLength[iNose]*qu.m).rescale('ft'))
        WnoseGear =0.125*((Nl*Wl)**0.566)*((ln/12.)**0.845)
        MnoseGear =float((WnoseGear*qu.pound).rescale('kg'))
        CGnoseGear=np.array([ac.landingGear.groundContact_X[iNose],ac.landingGear.groundContact_Y[iNose],ac.landingGear.groundContact_Z[iNose]+ac.landingGear.tireDiameter[iNose]/2.])
        InoseGear =np.zeros(3)
        self.aircraftMass.airframe.add_item('nose gear',MnoseGear,CGnoseGear,InoseGear)
        nMG       =float(len(ac.landingGear.groundContact_X)-1.)
        for i in range(len(ac.landingGear.groundContact_X)):
            if i<>iNose:
                lm      =ac.landingGear.strutLength[i]
                Wmlg    =0.095*((Nl*Wl)**0.768)*((lm/12.)**0.409)
                Mmlg    =float((Wmlg*qu.pound).rescale('kg'))/nMG
                CG      =np.array([ac.landingGear.groundContact_X[i],ac.landingGear.groundContact_Y[i],ac.landingGear.groundContact_Z[i]+ac.landingGear.tireDiameter[i]/2.])
                inertia =np.zeros(3)
                name    ='landing gear '+str(i)
                self.aircraftMass.airframe.add_item(name,Mmlg,CG,inertia)
    def _engine(self):
        ac        =self.aircraft
        for i,engine in enumerate(ac.engine.engineList):
            Men  =engine.mass
            Wen  =float((Men*qu.kg).rescale('lb'))
            Weng =2.3*(Wen**0.922)
            Meng =float((Weng*qu.pound).rescale('kg'))
            CG   =np.array([ac.engine.CG_X[i],ac.engine.CG_Y[i],ac.engine.CG_Z[i]])
            I    =np.zeros(3)
            name ='engine '+str(i+1)
            self.aircraftMass.airframe.add_item(name,Meng,CG,I)
    def _fuel_(self,fuelFraction=1.0):
        # will be replaced, kept here for reference. Do not use!!!
        ac        =self.aircraft
        if ac.wing.fuelTank.mirrorFlag<>1:M=ac.wing.fuelTank.fuelMass*fuelFraction
        else:M=2.*ac.wing.fuelTank.fuelMass*fuelFraction
        CGR       =ac.wing.fuelTank.CG
        CGL       =CGR*np.array([1.,-1.,1])
        I         =np.zeros(3)
        self.aircraftMass.fuel.add_item('Fuel tank Right',M/2.,CGR,I)
        self.aircraftMass.fuel.add_item('Fuel tank Left' ,M/2.,CGL,I)
    def _fuel(self):
        sr = self.aircraft.wing.fuelTank.spanRatio
        cr = self.aircraft.wing.fuelTank.chordRatio
        CG = self.aircraft.wing.locateOnWing(sr,cr)
        self.aircraftMass.fuel.update_item_cg('Fuel tank Right',CG[0], CG[1],CG[2])
        self.aircraftMass.fuel.update_item_cg('Fuel tank Left', CG[0],-CG[1],CG[2])
    def _fuel_system(self):
        ac        =self.aircraft
        fuseCG    =(self.aircraftMass.airframe.get_item_by_name('fuselage')).coords
        I         =np.zeros(3)
        # Fuel System
        Vt        =float((ac.fuselage.diameter*qu.m).rescale('ft'))
        Vi        =0.0
        Nt        =2.0
        Wfs       =2.49*(Vt**0.726)*((1./(1.+(Vi/Vt)))**0.363)*(Nt**0.242)*(self._Nen**0.157)
        Mfs       =float((Wfs*qu.pound).rescale('kg'))
        name      ='fuel system'
        self.aircraftMass.airframe.add_item(name,Mfs,fuseCG,I)
    def _flight_controls(self):
        ac        =self.aircraft
        fuseCG    =(self.aircraftMass.airframe.get_item_by_name('fuselage')).coords
        I         =np.zeros(3)
        L         =float((ac.fuselage.length*qu.m).rescale('ft'))
        Bw        =float((ac.wing.span*qu.m).rescale('ft'))
        Wfc       =0.053*(L**0.1536)*(Bw**0.371)*((self._Nz*self._Wdg*(10.**(-4)))**0.8)
        Mfc       =float((Wfc*qu.pound).rescale('kg'))
        name      ='flight controls'
        self.aircraftMass.airframe.add_item(name,Mfc,fuseCG,I)
    def _hydraulics(self):
        ac        =self.aircraft
        fuseCG    =(self.aircraftMass.airframe.get_item_by_name('fuselage')).coords
        I         =np.zeros(3)
        Kh        =0.013
        W         =float((ac.fuselage.diameter*qu.m).rescale('ft'))
        M         =self._designM
        Whyd      =Kh*(W**0.8)*(M**0.5)
        Mhyd      =float((Whyd*qu.pound).rescale('kg'))
        name      ='hydraulics'
        self.aircraftMass.airframe.add_item(name,Mhyd,fuseCG,I)
    def _misc(self):
        self._fuel_system()
        self._flight_controls()
        self._hydraulics()
    def get_inertia(self):
        #Temporary statistical method until more detailed analysis becomes available
        #TODO: add detailed intertia analysis
        ac=self.aircraft 
        self.aircraftMass.airframe.update_totals()
        Mtot=self.aircraftMass.airframe.get_total_mass()
        R=[.25,.38,.39]
        Ixx=ac.wing.span**2*Mtot*R[0]**2/4.
        Iyy=ac.fuselage.length**2*Mtot*R[1]**2/4.
        Izz=(ac.wing.span/2+ac.fuselage.length/2)**2*Mtot*R[2]**2/4.
        self.inertia=np.array([Ixx,Iyy,Izz])
        self.aircraftMass.inertia = self.inertia



def debug1():
    import aircraft
    ac = aircraft.load('V0510')
    ac.mass.fuel.set_fuel_mass(30,'total')
    ac.mass.display()
    ac.mass.fuel.set_fuel_mass_burned(1.5,'total')
    ac.mass.display()
    
    baggage = MassComponent('baggage', 70.0, [2.85,0,-0.125])
    baggage.display(True,True)

def kla100_mass():
    import aircraft
    ac = aircraft.load('V0510')
    print 'Aft CG\n======'
    ac.mass.update_item_mass('Payload','baggage',70.)
    ac.mass.display()
    ac.mass.update_item_mass('Fuel','Fuel Tank Right',47.5/2.)
    ac.mass.update_item_mass('Fuel','Fuel Tank Left',47.5/2.)
    ac.mass.total.display()
    ac.mass.update_item_mass('Fuel','Fuel Tank Right',2.5)
    ac.mass.update_item_mass('Fuel','Fuel Tank Left',2.5)
    ac.mass.total.display()
    print 'front CG\n========'
    ac.mass.update_item_mass('Payload','passenger',86.)
    ac.mass.update_item_mass('Payload','pilot',86.)
    ac.mass.update_item_mass('Payload','baggage',0.0)

    ac.mass.update_item_mass('Fuel','Fuel Tank Right',47.5)
    ac.mass.update_item_mass('Fuel','Fuel Tank Left',47.5)
    ac.mass.total.display()

    ac.mass.update_item_mass('Fuel','Fuel Tank Right',47.5/2)
    ac.mass.update_item_mass('Fuel','Fuel Tank Left',47.5/2)
    ac.mass.total.display()

    ac.mass.update_item_mass('Fuel','Fuel Tank Right',2.5)
    ac.mass.update_item_mass('Fuel','Fuel Tank Left',2.5)
    ac.mass.total.display()
    ac.display()

if __name__=="__main__":
    debug1()