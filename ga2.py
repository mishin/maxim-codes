# -*- coding: utf-8 -*-
"""
Created on Sun May 27 21:45:36 2012
the bounds are controlled in Mutation only (the only source)
@author: Maxim
"""

import numpy as ny
import matplotlib.pyplot as plt
import sys
import math

class gaOptions:
    def __init__(self,Lbound,Ubound):
        Lbsize, Ubsize = Lbound.shape[0],Ubound.shape[0]
        if  Lbsize == Ubsize:
            self.lb = ny.zeros([Lbsize])
            self.ub = ny.zeros([Lbsize])
            self.Nvar = Lbsize
            for ii in range(Lbsize):
                self.lb[ii] = float(Lbound[ii])
                self.ub[ii] = float(Ubound[ii])
        else:
            sys.exit("Error: Bound dimmensions should match")
            
        self.PopSize = 5*self.Nvar
        self.IterMax = 50
        self.LagIterMax = int(0.5*self.IterMax)
        self.MutRate = 0.2
        self.sigma = 1.0
        self.sigmaEnd = 0.001
        self.EliteRatio = 0.2
        self.rho = 1e10
        self.CorrNan = False
        self.HistFile = []
    def MaxIterations(self,iterMax,lagIterMax=[]):
        self.IterMax = iterMax
        if lagIterMax == []:
            self.LagIterMax = int(0.4*self.IterMax)
        else:
            self.LagIterMax = lagIterMax
    def getSigma(self,CurIter):
        CurIter = float(CurIter)
        IterMax = float(self.IterMax)
        return (1-CurIter/IterMax)*(self.sigma-self.sigmaEnd)+self.sigmaEnd

def realMax(X):
    maximum = -float('inf')
    for x in X:
        if math.isnan(x)==False and maximum<x: maximum = x
    return maximum

def getCost(X,CostFcn,lb,ub,rho=100):
    X = denormalize(X,lb,ub)
    cost = []
    for x in X:
        f,g,h = CostFcn(x)
        P = 0
        if hasattr(g,"__iter__"):
            for gg in g: P = P + ny.max([0,gg])**2
        else:
            P = P + ny.max([0,g])**2
        if hasattr(h,"__iter__"):
            for hh in h: P = P + hh**2
        else:
            P = P + h**2
        cost = ny.append(cost,f+rho*P)
        
    return cost

def getFitness(cost,CorrNan):
    if CorrNan==True:
        fworst = -float('inf')
        for fval in cost:
            if math.isnan(fval)==False and fworst<fval: fworst = fval
        idx = 0
        for fval in cost:
            if math.isnan(fval)==True: cost[idx] = fworst
            idx +=1
    else: fworst = ny.max(cost)
    diff = fworst - cost
    totalDiff = ny.sum(diff)
    return diff/totalDiff

def selection(fitness,prev=-1.):
    Selected = prev
    CumulSum = ny.cumsum(fitness)
    LagRunMax = 25
    LagRun = 0
    while Selected == prev:
        K = ny.random.rand()
        for ii in range(1,CumulSum.shape[0]):
            if K>CumulSum[ii-1] and K <=CumulSum[ii]: Selected = ii
        if K <=CumulSum[0]: Selected = 0
        LagRun +=1
        if LagRunMax<=LagRun: 
            print 'Fitness:\n',fitness
            sys.exit('Execution stopped!\nIncrease population size')
    return Selected

def mutation(x, sigma):
    outOfBound = True
    while outOfBound == True:
        xMut = ny.random.normal(x, sigma, x.shape)
        if ny.all(xMut>=0) and ny.all(xMut<=1): outOfBound = False
    return xMut

def elite(X, Nelite, fitness):
    idx = ny.argsort(1-fitness)[0:Nelite]
    return X[idx]

def NextGen(X,fitness, MutRate,sigma, EliteRatio):
    PopSize = X.shape[0]
    Nelite = int(EliteRatio*PopSize)
    Nmut = int(MutRate*PopSize)
    Xnew = elite(X,Nelite,fitness)
    for ii in range(Nmut):
        idx1 = selection(fitness)
        Xmut = mutation(X[idx1],sigma)
        Xnew = ny.vstack([Xnew,Xmut])
    while Xnew.shape[0] < PopSize:
        idx1 = selection(fitness)
        idx2 = selection(fitness,idx1)
        C1,C2 = crossOver(X[idx1],X[idx2])
        Xnew = ny.vstack([Xnew,C1,C2])
    if Xnew.shape[0] > PopSize: Xnew = Xnew[0:PopSize]
    return Xnew

def crossOver(P1,P2):
    C1 = ny.zeros(P1.shape)
    C2 = ny.zeros(P1.shape)
    for ii in range(P1.shape[0]):
        K = ny.random.rand()
        C1[ii] = K * P1[ii] + (1-K) * P2[ii]
        C2[ii] = (1-K) * P1[ii] + K * P2[ii]
    return C1, C2

def normalize(X,lb,ub):
    return (X-lb) / (ub-lb)

def denormalize(X,lb,ub):
    return X*(ub-lb)+lb
    
def writeHist(Options,Iter,xBest,fBest):
    if Iter==1:
        Hfile = open(Options.HistFile,'wt')
        Hfile.write('Iter\tFval\t')
        for ii in range(xBest.shape[0]): Hfile.write('x%d\t'%ii)
        Hfile.write('\n%d\t%.6f\t' %(Iter,fBest))
        for xx in xBest: Hfile.write('%.6f\t'%xx)
        Hfile.close()
    else:
        Hfile = open(Options.HistFile,'at')
        Hfile.write('\n%d\t%.6f\t' %(Iter,fBest))
        for xx in xBest: Hfile.write('%.6f\t'%xx)
        Hfile.close()
   

def gaMain(CostFcn, Options):
    Nvar    = Options.Nvar
    PopSize = Options.PopSize
    IterMax = Options.IterMax
    lb      = Options.lb
    ub      = Options.ub
    MutRate = Options.MutRate
    LagIterMax = Options.LagIterMax
    EliteRatio = Options.EliteRatio
    rho = Options.rho
    X = ny.random.random([PopSize,Nvar])
    fHist = []
    xHist = ny.zeros([1,Nvar])
    Iter = 1
    LagIter = 0
    fBest = []
    while Iter < IterMax and LagIter < LagIterMax:
        cost = getCost(X,CostFcn,lb,ub,rho)
        fitness = getFitness(cost,Options.CorrNan)
        sigma = Options.getSigma(Iter)
        fBestNew = ny.min(cost)
        if fBestNew == fBest: LagIter+=1 
        else: LagIter = 0
        fBest = fBestNew
        idxBest = ny.argmin(cost)
        xBest = X[idxBest]
        xBest = denormalize(xBest,lb,ub)
        xHist = ny.vstack([xHist,xBest])
        fHist = ny.append(fHist,fBest)
        if Options.HistFile != []: writeHist(Options,Iter,xBest,fBest)
        X = NextGen(X,fitness, MutRate,sigma,EliteRatio)
        Iter +=1
    xHist = ny.delete(xHist,0,0)
    return fBest,xBest,fHist,xHist, Iter-1

#Test section

def testFcn1(x):
    #return x+1/x
    return (1-x[0])**2+100*(x[1]-x[0]**2)**2
    #return x[0]+x[1]

def testRun():
    def CostFcn(x):
        f = testFcn1(x)
        h=[]
        g=[]
        return f,g,h

    lb = ny.array([0,0])
    ub = ny.array([6,6])
    
    opt = gaOptions(lb,ub)
    opt.PopSize = 25
    opt.MaxIterations(100)
    
    fBest,xBest,fHist,xHist = gaMain(CostFcn, opt)
    
    print '\n',fBest, xBest
    plt.figure(1)
    plt.plot(xHist[:,0],xHist[:,1],'o')
    plt.hold(True)
    plt.plot(xBest[0],xBest[1],'ro')
    plt.grid(True)
    plt.axis([0,6,0,6])
    plt.xlabel('generation')
    plt.ylabel('function value')
    plt.title('convergence history')
    plt.show()

#testRun()