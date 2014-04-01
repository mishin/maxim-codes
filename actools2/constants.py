# -*- coding: utf-8 -*-
"""
Created on Wed Jan 08 13:37:27 2014

@author: Maxim
"""
from paths import MyPaths
from db_tools import ReadDatabase

GRAVITY_ACCEL = 9.80665
EARTH_RADIUS = 6371.0e3
DENSITY_GASOLINE = 710.0 #kg/m3
DENSITY_KEROSENE = 800.0

def load(sheetName):
    """
    loads list of constants from xls database
    """
    pth    = MyPaths()
    output = {}
    sh = ReadDatabase(pth.db.constants,sheetName)
    for irow in range(1,sh.nrows):
        name  = sh.read_cell(irow,1)
        value = sh.read_row(irow,2,False)
        output[str(name)] = value
    return output

def get_gravity_acceleration(altitude):
    c = (EARTH_RADIUS/(EARTH_RADIUS+altitude))**2.0
    return GRAVITY_ACCEL*c


mass = load('mass')
fuelDensity = load('fuel_density')

# --- debug ---
def run_test1():
    mass = load('mass')
    print mass['compCorr']
    print mass['fuseCGratio']
    print mass['wingCGratio']

if __name__=="__main__":
    run_test1()