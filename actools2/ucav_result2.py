# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 10:30:40 2014

@author: Maxim
"""

from ucav_design_1 import DesignFormulation
import numpy as np

ac = DesignFormulation()
ac.load_xls('Baseline1')
ac.setup()


xopt = np.array([-0.6510831,0.10411176,0.18537688,-0.99810656,-0.99857516,
                 0.73880065,0.44773939,0.99990615,0.22948609])
ac.set_x(xopt)

ac.mass.display()
ac.display()
print ac.get_cg()
print ac.get_mass()
print ac.wing.area
print ac.analysisData
print ac.g(xopt)
print ac.propulsion.engine.length
print ac.wing.chords[0]
print ac.wing.MAC
print ac.get_cg()
ac.mass.empty.update_item_cg('engine',5,0,0)
ac.set_engine_cg(5,0,0)
ac.mass.display()
print ac.get_cg()