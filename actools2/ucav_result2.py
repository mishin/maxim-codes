# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 10:30:40 2014

@author: Maxim
"""

from ucav_design_1 import DesignFormulation
import numpy as np
import matplotlib.pyplot as plt


def plot_results():
    ac = DesignFormulation()
    ac.load_xls('Baseline1')

    # baseline
    ac.setup()
    ac.gvfmAero = True
    x0 = ac.wing.x2d
    y0 = ac.wing.y2d
    cg0 = ac.get_cg()
    
    xbase = ac.x0norm
    ac.set_x(xbase)
    ac.show_results()
    print ac.get_drag()
    # low fidelity optimum
    xopt1 = np.array([ 0.20520747,-0.39226611,-0.47326701,-1.,0.14291925,0.98650447, 
                      0.37346922,  0.37756372,  0.65654722])

    print 'Low fidelity optimum\n==========='
    ac.set_x(xopt1)
    ac.show_results()
    print ac.get_drag()
  
    x1 = ac.wing.x2d
    y1 = ac.wing.y2d
    cg1 = ac.get_cg()
    
    xgvfm1 = np.array([-0.10127195,-0.10127187,-0.60597138,-0.99999966,-0.71421434,0.97300896,  1.00000803, 0.57455958,0.3130853])
    xgvfm2 = np.array([-0.04819471,-0.04819471,-0.59713618,-0.99999966,-1.,0.94601792,1.00000803,0.81889676,-0.3738294])
    xgvfm4 = np.array([-0.10676065, -0.10676065, -0.61007048, -1.00003353, -0.57133293,  0.9857291,  1.00000803,  0.56576939,  0.08997523])
    print np.linalg.norm(xgvfm1-xgvfm2)
    ac.set_x(xgvfm4)
    print ac.get_cg()
    print ac.wing.area
    print ac.wing.MAC
    print ac.wing.span
    ac.show_results()
    print ac.get_drag()
  
    x3 = ac.wing.x2d
    y3 = ac.wing.y2d
    cg3 = ac.get_cg()
    

    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.hold(True)
    ax1.axis('equal')
    ax1.plot(x0,y0,'ko-',linewidth=2)
    ax1.plot(x1,y1,'b^:',linewidth=2)
    ax1.plot(x3,y3,'rs--',linewidth=2)
    #plt.plot(x2,y2,'k-')
    ax1.plot(0,-cg0[0],'ko')
    ax1.plot(0,-cg1[0],'b^')
    ax1.plot(0,-cg3[0],'rs')
    ax1.legend(['Baseline','Low-fi optimum','GVFM optimum'])
    #plt.plot(0,-cg2[0],'ko')
    plt.show()


if __name__=="__main__":
    plot_results()