# -*- coding: utf-8 -*-
"""
Created on Sun May 04 23:08:08 2014

@author: Maxim
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from misc_tools import RbfMod, read_tabulated_data_without_header
from ga_new import get_fitness, selection, cross_over


#func = lambda x: np.sin(2.*x[0]) + np.sin(5.*x[1]) #[0;1] [0;1]
def func(x):
    x = 10*np.asarray(x)-5.0
    return (1.0-x[0])**2.0 + 100.*(x[1]-x[0]**2.0)**2.0

class Function:
    def __init__(self,func):
        self.func = func
    def __call__(self,x):
        x = np.asarray(x)
        n = x.shape[0]
        err = np.zeros(n)
        errSum = 0.0
        fval = np.zeros(n)
        for i,xval in enumerate(x):
            fval[i] = self.func(xval)

        for i in range(n):
            xx = np.copy(x)
            ff = np.copy(fval)
            xx = np.delete(xx,i,0)
            ff = np.delete(ff,i,0)
            model = RbfMod(xx,ff)
            err[i] = model(x[i]) - fval[i]
            errSum += err[i]**2.0
        return -err, (errSum**0.5)/n
            

def test_plot():
    lb = np.array([0,0.])
    ub = np.array([1.,1])
    n = 20
    x = y = np.linspace(lb[0],ub[0],n)
    X, Y = np.meshgrid(x,y)
    zs1 = np.array([func([x,y]) for x,y in zip(np.ravel(X),np.ravel(Y))])
    Z1 = zs1.reshape(X.shape)
    
    fig2 = plt.figure(1)
    ax2 = Axes3D(fig2)
    ax2.hold(True)
    ax2.plot_wireframe(X,Y,Z1,color='r')
    ax2.legend(['High fidelity function','Low fidelity function'])
    ax2.set_xlabel('x1')
    ax2.set_ylabel('x2')
    ax2.set_zlabel('f (x1,x2)')
    plt.show()

def run_test1():
    xdoe = np.array([[0.,0],[0,0.5],[0,1],[0.5,0],[0.5,0.5],[0.5,1],[1.,0],[1.,0.5],[1,1]])
    fdoe = [func(xx) for xx in xdoe]
    newfunc = Function(func)
    err1,errSq1 = newfunc(xdoe)
    fitness = get_fitness(err1)
    xnew = np.zeros(xdoe.shape)
    for i in range(4):
        idx1 = selection(fitness)
        idx2 = selection(fitness,idx1)
        xnew[2*i],xnew[2*i+1] = cross_over(xdoe[idx1],xdoe[idx2])
    idx1 = selection(fitness)
    idx2 = selection(fitness,idx1)
    xnew[-1],_t = cross_over(xdoe[idx1],xdoe[idx2])
    fnew = [func(xx) for xx in xnew]
    xnew1 = np.vstack([xdoe,xnew])
    fnew1 = [func(xx) for xx in xnew1]
    err,errSq2 = newfunc(xnew1)
    print errSq1, errSq2
    xrand = read_tabulated_data_without_header('LHC18_2.txt')
    xrand = (xrand+1.0)/2.
    #xrand = np.random.random([18,2])
    frand = [func(xx) for xx in xrand]
    print newfunc(xrand)[1]
    
    model1 = RbfMod(xnew1,fnew1)
    model2 = RbfMod(xrand,frand)
    lb = np.array([0,0.])
    ub = np.array([1.,1])
    n = 40
    xs = ys = np.linspace(lb[0],ub[0],n)
    X, Y = np.meshgrid(xs,ys)
    zs1 = np.array([func([xs,ys]) for xs,ys in zip(np.ravel(X),np.ravel(Y))])
    Z1 = zs1.reshape(X.shape)
    zs2 = np.array([model1([xs,ys]) for xs,ys in zip(np.ravel(X),np.ravel(Y))])
    Z2 = zs2.reshape(X.shape)
    zs3 = np.array([model2([xs,ys]) for xs,ys in zip(np.ravel(X),np.ravel(Y))])
    Z3 = zs3.reshape(X.shape)
    fig2 = plt.figure(1)
    ax2 = Axes3D(fig2)
    ax2.hold(True)
    ax2.plot_wireframe(X,Y,Z1,color='b')
    ax2.plot_wireframe(X,Y,Z2,color='r')
    ax2.plot_wireframe(X,Y,Z3,color='y')
    ax2.plot(xdoe[:,0],xdoe[:,1],fdoe,'ro')
    ax2.plot(xnew[:,0],xnew[:,1],fnew,'go')
    #ax2.plot(xdoe[:,0],xdoe[:,1],err1,'ko')
    ax2.plot(xrand[:,0],xrand[:,1],frand,'bo')
    ax2.legend(['High fidelity function','Low fidelity function'])
    ax2.set_xlabel('x1')
    ax2.set_ylabel('x2')
    ax2.set_zlabel('f (x1,x2)')
    plt.show()


if __name__=="__main__":
    run_test1()