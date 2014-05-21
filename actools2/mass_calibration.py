# -*- coding: utf-8 -*-
"""
Created on Mon May 12 16:01:31 2014

@author: Maxim
"""

import aircraft_FW
import numpy as np
import matplotlib.pyplot as plt

def run():
    mass1 = np.array([1740., 3600, 4900, 6350])
    mass4 = np.array([2678.,8165,7000,20215])
    name = ['X47A','X45C','nEUROn','X47B']

    ac1 = aircraft_FW.load(name[0])
    ac2 = aircraft_FW.load(name[1])
    ac3 = aircraft_FW.load(name[2])
    ac4 = aircraft_FW.load(name[3])
    m1 = ac1.mass.empty.get_total_mass()
    m2 = ac2.mass.empty.get_total_mass()
    m3 = ac3.mass.empty.get_total_mass()
    m4 = ac4.mass.empty.get_total_mass()
    m11 = ac1.mass.total.get_total_mass()
    m12 = ac2.mass.total.get_total_mass()
    m13 = ac3.mass.total.get_total_mass()
    m14 = ac4.mass.total.get_total_mass()

    ac1.mass.display()
    ac2.mass.display()
    ac3.mass.display()
    ac4.mass.display()

#    ac2.display()
#    ac4.display()
    print ac1.designGoals.grossMass
    print ac2.designGoals.grossMass
    print ac3.designGoals.grossMass
    print ac4.designGoals.grossMass

    mass2 = np.array([m1,m2,m3,m4])
    mass3 = np.array([m11,m12,m13,m14])
    loc = np.array([1,2,3,4])
    width = 0.4

    print 'exact\tcalculated\t\terror,%'
    for me,mc in zip(mass1,mass2):
        print '%.2f\t%.2f\t\t%+.4f'%(me,mc,100.*(mc-me)/me)

    plt.figure(1)
    plt.hold(True)
    plt.bar(loc-width,mass1,width=width,color='b')
    plt.bar(loc,mass2,width=width,color='r')
    plt.legend(['data','calculated'],'upper left')
    plt.ylabel('Empty mass, kg')
    plt.xticks(loc,name)
    plt.axis([0.4,4.6,0,7000])
    
    plt.figure(2)
    plt.hold(True)
    plt.bar(loc-width,mass4,width=width,color='b')
    plt.bar(loc,mass3,width=width,color='r')
    plt.legend(['data','calculated'],'upper left')
    plt.ylabel('Total mass, kg')
    plt.xticks(loc,name)
    plt.axis([0.4,4.6,0,22000])
    plt.show()

if __name__=="__main__":
    run()