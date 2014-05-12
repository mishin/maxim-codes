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
    name = ['X47A','X45C','nEUROn','X47B']
    
    ac1 = aircraft_FW.load(name[0])
    ac2 = aircraft_FW.load(name[1])
    ac3 = aircraft_FW.load(name[2])
    ac4 = aircraft_FW.load(name[3])
    m1 = ac1.mass.empty.get_total_mass()
    m2 = ac2.mass.empty.get_total_mass()
    m3 = ac3.mass.empty.get_total_mass()
    m4 = ac4.mass.empty.get_total_mass()
    
    ac1.mass.empty.display()

    mass2 = np.array([m1,m2,m3,m4])
    loc = np.array([1,2,3,4])
    width = 0.4

    print 'exact\tcalculated\t\terror,\%'
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
    plt.show()

if __name__=="__main__":
    run()