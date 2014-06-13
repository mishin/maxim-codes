# -*- coding: utf-8 -*-
"""
Created on Tue Jun 03 21:15:26 2014

@author: Maxim
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from mpl_toolkits.mplot3d import axes3d
from misc_tools import SaveTextData

out = SaveTextData('output1.py')

tol = 0.01
err = tol+1.

def f(x):
    f=x[0]+x[1]
    return f
def g1(x):
    C1=x[0]**2*x[1]/20.0-1.
    return C1

def g2(x):
    C2=(x[0]+x[1]-5.0)**2/30.+(x[0]-x[1]-12.0)**2/120.0-1.
    return C2
    
def g3(x):
    C3=80.0/(x[0]**2+8.0*x[1]+5.0)-1
    return C3

def g2l(x):
    return 0.05*x[0]+0.8*x[1]-1.0



dg = 0.0
x0 = np.array([5.,5])
i = 0
print 'iter\tx1\tx2\tdg\tg2high\tg2low\tg2adapt\tfval'

x1 = list()
x2 = list()
fval = list()
g1val = list()
g2hval = list()
g2lval = list()
g2aval = list()
g3val = list()
errorVal = list()
offset = list()
itr = list()

neval = 0


cnstr = ({'type':'ineq','fun':g1},
             {'type':'ineq','fun':g2},
             {'type':'ineq','fun':g3},)
rslt = minimize(f,x0,method='SLSQP',constraints=cnstr,bounds=((0,10),(0,10)))
xexact = rslt.x
fexact = rslt.fun
while err>tol:
    g2new = lambda x: g2l(x) + dg
    itr.append(i)
    x1.append(x0[0])
    x2.append(x0[1])
    fval.append(f(x0))
    g1val.append(g1(x0))
    g2hval.append(g2(x0))
    g2lval.append(g2l(x0))
    g2aval.append(g2new(x0))
    g3val.append(g3(x0))
    offset.append(dg)
    errorVal.append(err)
    i += 1
    
   
    cnstr = ({'type':'ineq','fun':g1},
             {'type':'ineq','fun':g2new},
             {'type':'ineq','fun':g3},)
    rslt = minimize(f,x0,method='SLSQP',constraints=cnstr,bounds=((0,10),(0,10)))
    xnew = rslt.x
    neval += rslt.nfev
    err = np.linalg.norm(x0-xnew)
    
    x0 = xnew
    dg += g2(rslt.x) - g2new(rslt.x)

    #print '%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f'%(i,xopt[-1][0],xopt[-1][1],dg,g2hval[-1],g2lval[-1],g2scval[-1],fopt[-1])

print np.linalg.norm(xnew-xexact)
print fval[-1]-fexact
xopt = np.transpose(np.vstack([x1,x2]))
print neval
print i

out.write_array(itr,'itr')
out.write_array(x1,'x1')
out.write_array(x2,'x2')
out.write_array(fval,'f')
out.write_array(errorVal,'err')
out.write_array(offset,'dg')
out.write_array(g1val,'g1')
out.write_array(g2hval,'g2h')
out.write_array(g2lval,'g2l')
out.write_array(g2aval,'g2new')
out.write_array(g3val,'g3')
out.close()










delta = 0.1
x = np.arange(0, 10+delta, delta)
y = np.arange(0, 10+delta, delta)
X, Y = np.meshgrid(x, y)

Z1 = np.array([g1([xx,yy]) for xx,yy in zip(X,Y)])
Z2 = np.array([g2([xx,yy]) for xx,yy in zip(X,Y)])
Z3 = np.array([g3([xx,yy]) for xx,yy in zip(X,Y)])
Z4 = np.array([g2l([xx,yy]) for xx,yy in zip(X,Y)])
Z1 = np.reshape(Z1,X.shape)
Z2 = np.reshape(Z2,X.shape)
Z3 = np.reshape(Z3,X.shape)
Z4 = np.reshape(Z4,X.shape)
Z5 = np.array([g2new([xx,yy]) for xx,yy in zip(X,Y)])
Z5 = np.reshape(Z5,X.shape)

plt.figure(1)
plt.hold(True)
cs1 = plt.contour(X,Y,Z1,levels=[0])
cs2 = plt.contour(X,Y,Z2,levels=[0],colors=['r'],label='g2high')
cs3 = plt.contour(X,Y,Z3,levels=[0])
cs4 = plt.contour(X,Y,Z4,levels=[0],colors=['y'],label='g2low')
cs5 = plt.contour(X,Y,Z5,levels=[0],colors=['g'],label='g2low+dg')
plt.plot(xopt[:,0],xopt[:,1],'ro-')

pth1 = cs1.collections[0].get_paths()[0].vertices
pth2 = cs2.collections[0].get_paths()[0].vertices
pth3 = cs3.collections[0].get_paths()[0].vertices
pth4 = cs4.collections[0].get_paths()[0].vertices
pth5 = cs5.collections[0].get_paths()[0].vertices

plt.figure(2)
plt.hold(True)
plt.plot(pth1[:,0],pth1[:,1],'k-',label='g1')
plt.plot(pth2[:,0],pth2[:,1],'r-',label='g2high',linewidth=2.0)
plt.plot(pth3[:,0],pth3[:,1],'k-',label='g3')
plt.plot(pth4[:,0],pth4[:,1],'b.-',label='g2low')
plt.plot(pth5[:,0],pth5[:,1],'g--',label='g2adapted')
plt.plot(xopt[:,0],xopt[:,1],'mo-',label='X optimum')
plt.xlabel('X1')
plt.ylabel('X2')
plt.legend()

#fig2 = plt.figure(2)
#ax = fig2.add_subplot(111, projection='3d')
#ax.hold(True)
#ax.plot_wireframe(X, Y, Z2,color='r')
#ax.plot_wireframe(X, Y, Z4,color='g')
plt.show()