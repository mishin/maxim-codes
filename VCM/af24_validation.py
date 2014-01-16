# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 16:00:01 2014

@author: Maxim
"""
from numpy import array
from scipy.interpolate import interp1d
#from matplotlib.pyplot import figure, plot, hold, show, axis, title, set_set_xlabel, set_ylabel, legend
import matplotlib.pyplot as plt

class Polar:
    def __init__(self):
        self.alpha = None
        self.cl = None
        self.cd = None
        self.LD = None
    def create_curves(self):
        self.clalpha = interp1d(self.alpha,self.cl,'linear')
        self.cdalpha = interp1d(self.alpha,self.cd,'linear')
        

def airfoil_validation():
    labels = ['experiment','fluent(k-epsilon)','fluent(SA)','Xfoil','javafoil']
    linetype = ['ro','b-','b--','k.-','k-.']
    exp = Polar()
    cfd1 = Polar()
    cfd2 = Polar()
    xfoil = Polar()
    javafoil = Polar()

    exp.alpha = array([-5.82, -4.86, -3.85, -2.79, -1.73, -1.28, -0.74, -0.24, 0.26, 0.73, 1.24, 1.8, 2.34, 2.86, 3.38, 4.36, 5.38, 6.43, 7.45, 8.4, 9.46, 10.45, 11.49])
    exp.cl = array([-0.376, -0.286, -0.165, -0.033, 0.082, 0.131, 0.184, 0.23, 0.28, 0.33, 0.382, 0.441, 0.492, 0.548, 0.601, 0.703, 0.798, 0.901, 0.987, 1.062, 1.121, 1.169, 1.179])
    exp.cd = array([0.0179, 0.0128, 0.0091, 0.0075, 0.0069, 0.0065, 0.0063, 0.0061, 0.006, 0.0061, 0.0066, 0.0068, 0.0072, 0.0075, 0.008, 0.0091, 0.0105, 0.0124, 0.0149, 0.0186, 0.024, 0.0313, 0.0543])
    
    cfd1.alpha = array([-10, -5, 0, 5, 10, 12, 14, 16, 18])
    cfd1.cl = array([-0.800332, -0.29402, 0.239798, 0.763995, 1.21425, 1.3484, 1.4016, 1.09993, 1.04271])
    cfd1.cd = array([0.116982, 0.0367235, 0.022792, 0.0386915, 0.0879292, 0.123326, 0.17428, 0.210824, 0.276682])
    
    cfd2.alpha = array([ -10, -5, 0, 5, 10, 12, 14, 16, 18])
    cfd2.cl = array([-0.809755, -0.293005, 0.227697, 0.745432, 1.19817, 1.34888, 1.19873, 1.01502, 0.936905])
    cfd2.cd = array([0.109056, 0.0332904, 0.0191247, 0.0344199, 0.0827258, 0.118184, 0.182808, 0.243119, 0.27792])
    
    xfoil.alpha = array([-10, -8, -4, -2, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13])
    xfoil.cl = array([-0.3999, -0.4038, -0.2062, 0.087, 0.2913, 0.4014, 0.511, 0.6198, 0.7275, 0.8335, 0.9366, 1.0341, 1.1209, 1.1953, 1.2485, 1.204, 1.1271])
    xfoil.cd = array([0.10403, 0.06903, 0.012, 0.00734, 0.006, 0.00633, 0.00688, 0.00767, 0.00869, 0.00998, 0.01171, 0.01426, 0.01852, 0.02395, 0.04571, 0.0606, 0.08885])

    
    javafoil.alpha = array([-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    javafoil.cl = array([-0.412, -0.388, -0.285, -0.153, 0.077, 0.311, 0.545, 0.775, 0.976, 0.881, 0.868, 0.76, 0.612, 0.472, 0.36, 0.276])
    javafoil.cd = array([0.08063, 0.06094, 0.04475, 0.01091, 0.00894, 0.00804, 0.00891, 0.01097, 0.01502, 0.05345, 0.07205, 0.09717, 0.12667, 0.16308, 0.19865, 0.24155])
    
    aNew = exp.alpha
    exp.create_curves()
    cfd1.create_curves()
    cfd2.create_curves()
    xfoil.create_curves()
    javafoil.create_curves()
    
    clerr2 = cfd1.clalpha(aNew) - exp.cl
    cderr2 = cfd1.cdalpha(aNew) - exp.cd
    clerr3 = cfd2.clalpha(aNew) - exp.cl
    cderr3 = cfd2.cdalpha(aNew) - exp.cd
    clerr4 = xfoil.clalpha(aNew) - exp.cl
    cderr4 = xfoil.cdalpha(aNew) - exp.cd
    clerr5 = javafoil.clalpha(aNew) - exp.cl
    cderr5 = javafoil.cdalpha(aNew) - exp.cd

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.hold(True)
    ax1.set_xlabel('AoA')
    ax1.set_ylabel('Cl')
    ax1.plot(exp.alpha, exp.cl, linetype[0])
    ax1.plot(cfd1.alpha, cfd1.cl, linetype[1])
    ax1.plot(cfd2.alpha, cfd2.cl, linetype[2])
    ax1.plot(xfoil.alpha, xfoil.cl, linetype[3])
    ax1.plot(javafoil.alpha, javafoil.cl, linetype[4])
    ax1.legend(labels,'upper left')
    
    
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111)
    ax2.hold(True)
    ax2.set_xlabel('AoA')
    ax2.set_ylabel('Cd')
    ax2.plot(exp.alpha, exp.cd, linetype[0])
    ax2.plot(cfd1.alpha, cfd1.cd, linetype[1])
    ax2.plot(cfd2.alpha, cfd2.cd, linetype[2])
    ax2.plot(xfoil.alpha, xfoil.cd, linetype[3])
    ax2.plot(javafoil.alpha, javafoil.cd, linetype[4])
    ax2.legend(labels,'upper left')
    
    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot(111)
    ax3.hold(True)
    ax3.set_xlabel('AoA')
    ax3.set_ylabel('L/D')
    ax3.plot(exp.alpha, exp.cl/exp.cd, linetype[0])
    ax3.plot(cfd1.alpha, cfd1.cl/cfd1.cd, linetype[1])
    ax3.plot(cfd2.alpha, cfd2.cl/cfd2.cd, linetype[2])
    ax3.plot(xfoil.alpha, xfoil.cl/xfoil.cd, linetype[3])
    ax3.plot(javafoil.alpha, javafoil.cl/javafoil.cd, linetype[4])
    ax3.legend(labels,'upper left')
    
    fig4 = plt.figure(4)
    ax4 = fig4.add_subplot(111)
    ax4.set_xlabel('Angle of attack, deg')
    ax4.set_ylabel('Lift coefficient error')
    ax4.hold(True)
    ax4.grid(True)
    #ax4.set_yscale('log')
    ax4.plot(aNew, clerr2, linetype[1], lw=1.5)
    ax4.plot(aNew, clerr3, linetype[2], lw=1.5)
    ax4.plot(aNew, clerr4, linetype[3], lw=1.5)
    ax4.plot(aNew, clerr5, linetype[4], lw=1.5)
    ax4.legend(labels[1:],'lower left')
    
    
    fig5 = plt.figure(5)
    ax5 = fig5.add_subplot(111)
    ax5.set_xlabel('Angle of attack, deg')
    ax5.set_ylabel('Drag coefficient error')
    ax5.hold(True)
    ax5.grid(True)
    #ax5.set_yscale('log')
    ax5.plot(aNew, cderr2, linetype[1], lw=1.5)
    ax5.plot(aNew, cderr3, linetype[2], lw=1.5)
    ax5.plot(aNew, cderr4, linetype[3], lw=1.5)
    ax5.plot(aNew, cderr5, linetype[4], lw=1.5)
    ax5.legend(labels[1:],'upper left')
    
    
    
    plt.show()


if __name__=="__main__":
    airfoil_validation()