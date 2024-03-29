# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 21:27:51 2014

@author: Maxim
"""
import numpy as np
from time import ctime


def get_total_cg(mass,CGx,CGy,CGz):
    mTotal = 0.0
    M  = np.zeros(3)
    for m,x,y,z in zip(mass,CGx,CGy,CGz):
        M = m*np.array([x,y,z])
        mTotal += m
    CG = M/mTotal
    return mTotal, CG

class AircraftMass:
    def __init__(self,name='Aircraft',MAC=None,xMAC=None):
        self.name = str(name)
        self.empty      = MassList(self.name + ' empty',MAC,xMAC)
        self.payload    = MassList(self.name + ' payload',MAC,xMAC)
        self.fuel       = Fuel()
        self.total      = MassList(self.name + ' total',MAC,xMAC)
        self.totalMass = 0.0
        self.totalCG   = np.zeros(3)
        self.totalMOI  = np.zeros(3)
        self.MAC = MAC
        self.xMAC = xMAC
        self.update_mac(MAC,xMAC)
    
    def set_fuel_mass(self,mass,CG):
        self.fuel = Fuel(mass,CG)
        self.update_total() 

    def __call__(self):
        self.update_total()
        return self.total.totalMass

    def update_total(self):
        """
        Updates total mass list.
        """
        self.total = MassList(self.name +' total')
        self.total.add_mass_list(self.empty)
        self.total.add_mass_list(self.payload)
        self.total.add_component(self.fuel)
        #self.total.add_mass_list(self.fuel)
        self.total.update_totals()
        self.totalMass = self.total.totalMass
        self.totalCG = self.total.CG
        self.update_mac(self.MAC,self.xMAC)
    
    def update_mac(self,MAC,xMAC):
        if not (MAC==None and xMAC==None):
            MAC  = float(MAC)
            xMAC = float(xMAC)
        self.total.MAC    = MAC
        self.total.xMAC   = xMAC
        self.empty.MAC    = MAC
        self.empty.xMAC   = xMAC
        self.payload.MAC  = MAC
        self.payload.xMAC = xMAC
        self.fuel.MAC     = MAC
        self.fuel.xMAC    = xMAC


    def display(self):
        """ prints out tabulated information about all mass components
        """
        self.update_total()
        self.total.display()


class MassComponent(object):
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
    def __init__(self,name='',mass=0.0,CG=None,inertia=np.zeros(3)):
        MassComponent.numberOfItems += 1
        if name=='':
            self.name = 'component %d'%MassComponent.numberOfItems
        else:
            self.name = str(name)
        self.cgIsKnown = True
        if CG==None:
            self.cgIsKnown = False
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
        
        if self.cgIsKnown:
            line = '{0:18} {1:8.2f} {2:10.4f} {3:10.4f} {4:10.4f}'
            x, y, z = self.coords[0], self.coords[1], self.coords[2]
        else:
            line = '{0:18} {1:>8.2f} {2:>10} {3:>10} {4:>10}'
            x, y, z = 'N/A', 'N/A', 'N/A'
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
    def __init__(self,name='',MAC=None,xMAC=None):
        MassList.n += 1
        if name=='':
            self.name = 'mass list %d'%MassList.n
        else:
            self.name = name
        self.items = list()
        self.totalMass = 0.0
        self.CG = np.zeros(3)
        self.names = list()
        self.MAC  = None
        self.xMAC = None
        if not self.MAC==None:
            self.MAC = float(MAC)
        if not self.xMAC==None:
            self.xMAC = float(xMAC)

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
        newItem = MassComponent(name,mass,CG,inertia)
        self.names.append(name)
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
        return self.item_exists(name)
        
    def remove_item(self,name):
        """ remove mass component with specified name from mass list
        """
        if name in self.names:
           k=self.names.index(name)
           del self.items[k]
           del self.names[k]
           self.update_totals()
        else:
            print 'Specified mass item not found: %s'%name
            print self.names
    
    def item_exists(self,name):
        return name in self.names
            
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
        tmpMass = 0.0
        CG = np.zeros(3)
        for item in self.items:
            totalMass += item.mass
            if item.cgIsKnown:
                tmpMass += item.mass
                CG += item.mass*item.coords
        self.totalMass = totalMass
        if totalMass>0.0:
            self.CG = CG / tmpMass
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
        if self.MAC!=None and self.xMAC!=None:
            cgMAC = (self.CG[0]-self.xMAC)/self.MAC*100.
            report += '{0:29} {1:+8.2f} % of MAC\n'.format('',cgMAC)
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


class Fuel(MassComponent):
    def __init__(self,mass=0.0,CG=np.zeros(3),inertia=np.zeros(3),
                 density=750.,maximumMass=0.0,tankVolume=0.0):
        super(Fuel,self).__init__('Fuel',mass,CG,inertia)
        self.density = float(density)
        self.maximumMass = float(maximumMass)
        self.tankVolume = float(tankVolume)
    
    def set_fuel_mass(self,fuelMass):
        self.mass = float(fuelMass)

    def set_fuel_burned(self,burnedFuelMass):
        self.mass += -float(burnedFuelMass)

    def set_fuel_fraction(self,fuelFraction):
        self.mass = self.maximumMass*fuelFraction
    
    def __repr__(self):
        out = '{0:10} = {1:8.4f}kg'.format('FuelMass',self.mass)
        return out


class _AircraftMass:
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
    def __init__(self,name='noname aircraft',MAC=None,xMAC=None):
        self.MAC = MACcomponents
        self.xMAC = xMAC
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
        self.update_total()
        self.total.display(self.MAC,self.xMAC)

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


def run_test1():
    acmass = AircraftMass('sample aircraft')
    baggage = MassComponent('baggage',25,[1.0,0.0,0.0])
    acmass.payload.add_component(baggage)
    acmass.set_fuel_mass(50,[1,2,3])
    acmass.display()
    acmass.payload.add_item('Weapons',200,[2.,0,0])
    acmass.display()


if __name__=="__main__":
    run_test1()