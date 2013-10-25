# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 12:35:21 2013

@author: Maxim
"""
import airfoil
import matplotlib.pyplot as plt

def test_plot_airfoil():
    af1 = airfoil.Airfoil()
    af1.read_txt('GA37A315.txt')
    Au = [0.1823,0.3380,0.2650,0.2016]
    Al = [-0.1823,-0.0708,-0.1885,-0.0384]
    af2 = airfoil.cst(Au,Al)
    
    rslt1 = af1.get_X_polar(0.18,4.4e6,[-10,20,1.0])
    rslt2 = af2.get_X_polar(0.18,4.4e6,[-10,20,1.0])
    
    alpha1 = [0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19,20]
    cl1 = [0.243296,0.453613,0.558244,0.661459,0.763304,0.859225,0.951732,1.04177,1.12103,1.19799,1.26788,1.32786,1.38623,1.41998,1.42465,1.23664,1.18999,1.17667]
    
    alpha2 = [0,1,2,5,7,8,12,13,14,15,16,17,18,19,20]
    cl2 = [0.264347,0.372483,0.481921,0.801596,0.998552,1.09719,1.41826,1.46123,1.52317,1.54066,1.4781,1.35913,1.31233,1.27969,1.26525]    
    print af2.thickness
    af1.analyze_geometry()
    print af1.thickness
    
    plt.figure(1)
    plt.grid(True)
    plt.hold(True)
    plt.plot(af1.coord[:,0],af1.coord[:,1],'b-')
    plt.plot(af2.coord[:,0],af2.coord[:,1],'r-')
    plt.legend(['baseline','optimum'])
    plt.figure(2)
    plt.grid(True)
    plt.hold(True)
    plt.plot(rslt1.cd,rslt1.cl,'b-')
    plt.plot(rslt2.cd,rslt2.cl,'r-')
    plt.axis([0.002,0.012,-0.5,1.0])
    plt.legend(['baseline','optimum'])
    plt.xlabel('drag coefficient')
    plt.ylabel('lift coefficient')
    
    plt.figure(3)
    plt.grid(True)
    plt.hold(True)
    plt.plot(rslt1.alpha,rslt1.cl,'b-')
    plt.plot(alpha1,cl1,'ko-')
    plt.plot(rslt2.alpha,rslt2.cl,'r-')
    plt.plot(alpha2,cl2,'go-')
    plt.axis([0,20,0,1.8])
    plt.legend(['baseline-low','baseline-high','optimum-low','optimum-high'],'lower right')
    plt.show()


def test_plot_convergence():
    nIter = [1,2,3]
    f = [0.004461,0.004679,0.004703]
    ghi = [1.4492,1.4983,1.5000]
    glow = [0.9016,1.1613,1.1615]
    gsc = [1.5000,1.5000,1.5004]
    rho = [0.3603,0.9670,0.8105]
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.grid(True)
    ax1.plot(nIter,f,'bo-')
    ax1.set_ylabel('drag coefficient')
    ax2 = fig.add_subplot(312)
    ax2.grid(True)
    ax2.plot(nIter,ghi,'b-')
    ax2.plot(nIter,glow,'k*-')
    ax2.plot(nIter,gsc,'ro')
    ax2.set_ylabel('c_lmax')
    ax3 = fig.add_subplot(313)
    ax3.grid(True)
    ax3.plot(nIter,rho,'ko-')
    ax3.set_xlabel('iterations')
    ax3.set_ylabel('rho')
    plt.show()
if __name__=="__main__":
    test_plot_airfoil()
    #test_plot_convergence()