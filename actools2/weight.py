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

from weight_ga import GeneralAviationMass
from weight_fw import BlendedWingBodyMass


def get_general_aviation_mass(aircraft):
    """
    Function returns aircraft mass object that contains lists of the mass 
    components. Mass is obtained using textbook methods by D.Raymer.
    """
    mass = GeneralAviationMass(aircraft)
    return mass.aircraftMass

def get_flying_wing_mass(aircraft):
    mass = BlendedWingBodyMass(aircraft)
    mass.analyze()



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