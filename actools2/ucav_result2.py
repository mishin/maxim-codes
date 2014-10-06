# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 10:30:40 2014

@author: Maxim
"""

from ucav_design_1 import DesignFormulation
import numpy as np
import matplotlib.pyplot as plt

ac = DesignFormulation()
ac.load_xls('Baseline1 (2)')
#cg1 = ac.get_cg()
#cg1[0] = cg1[0]+0.15
#print cg1
#ac.get_aero_single_point(0.7,1e4).display()
#ac.get_aero_trim(70,0,0.0,cg=cg1).display()
#ac.get_aero_trim(0.7,1e4,cg=cg1).display()
#print '--->'
#raw_input()

ac.setup()
ac.aero.display()


x0 = ac.wing.x2d
y0 = ac.wing.y2d
cg0 = ac.get_cg()
print ac.maxSectionLength
print ac.analysisData

# aero, performance, mission, stability
xopt1 = np.array([ 0.20520747,-0.39226611,-0.47326701,-1.,0.14291925,0.98650447,  0.37346922,  0.37756372,  0.65654722])
ac.set_x(xopt1)
ac.aero.display()

x1 = ac.wing.x2d
y1 = ac.wing.y2d
cg1 = ac.get_cg()

# added Cmde constraint
xopt2 = np.array([0.95866629, -0.07903226, -0.30488744, -0.2990727, -0.9974649, -0.59189978,
                  -0.99881686, 0.86857917, 0.97005337])
ac.set_x(xopt2)
x2 = ac.wing.x2d
y2 = ac.wing.y2d
cg2 = ac.get_cg()

print ac.analysisData

fig = plt.figure(1)
ax1 = fig.add_subplot(111)
ax1.hold(True)
ax1.axis('equal')
ax1.plot(x0,y0,linestyle='-',linewidth=2,color='b')
ax1.plot(x1,y1,linestyle='--',linewidth=2,color='k')
#plt.plot(x2,y2,'k-')
ax1.plot(0,-cg0[0],'bo')
ax1.plot(0,-cg1[0],'ks')
ax1.legend(['Baseline','Low-fi optimum'])
#plt.plot(0,-cg2[0],'ko')
plt.show()