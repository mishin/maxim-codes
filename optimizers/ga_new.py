# -*- coding: utf-8 -*-
"""
Created on Sun May 04 13:48:30 2014

@author: Maxim
"""
import numpy as np
import math
from misc_tools import Normalization

class CostFunction(object):
    def __init__(self,costFcn,lb,ub):
        """ cost function should return f,g"""
        self.costFcn = costFcn
        self.nEval = 0
        self.penalty = 1.0e5
        self.norm = Normalization(lb,ub,0,1)

    def evaluate(self,x,penalty=None):
        x = self.norm.denormalize(x)
        if penalty==None:
            penalty = self.penalty
        f,g = self.costFcn(x)
        self.nEval += 1
        cost = f
        for cnstr in g:
            cost += max([0.0,cnstr]) * penalty
        return cost,f,g

    def __call__(self,x,penalty=None):
        x = np.asarray(x)
        if x.ndim==1:
            return self.evaluate(x,penalty)
        elif x.ndim==2:
            n    = x.shape[0]
            cost = np.zeros(n)
            f    = np.zeros(n)
            cost[0],f[0],gtmp = self.evaluate(x[0],penalty)
            g = np.zeros([n,len(gtmp)])
            g[0] = gtmp
            for i in range(1,n):
                cost[i],f[i],g[i] = self.evaluate(x[i],penalty)
            return cost,f,g



class GAresult:
    def __init__(self):
        self.xBest = None
        self.fBest = None
        self.costBest = None
        self.gBest = None
        self.xHistory = None
        self.fHistory = None
        self.favgHistory = None
        self.iterations = 0
        self.lagiterations = 0
        self.message = 0
        self.success = True

def get_sigma(nIter,iterMax):
    ratio = float(nIter)/float(iterMax)
    sigma = 0.2 + 0.8*ratio
    return sigma

def get_fitness(cost):
    nanRatio = 0.5
    nNan = 0
    n = len(cost)
    fworst = np.max(cost)
    for i,fval in enumerate(cost):
        if math.isnan(fval):
            cost[i] = fworst
            nNan += 1
    if nNan>n*nanRatio:
        raise ValueError('Cost function returns too many NaN values')
    diff = fworst - cost
    totalDiff = np.sum(diff)
    return diff/totalDiff
        

def selection(fitness,idxPrevious=None):
    idxSelected = idxPrevious
    cumulSum = np.cumsum(fitness)
    lagrunMax = 25
    lagrun = 0
    n = cumulSum.shape[0]
    while idxSelected==idxPrevious:
        k = np.random.rand()
        for i in range(1,n):
            if k>cumulSum[i-1] and k<=cumulSum[i]:
                idxSelected = i
        if k<=cumulSum[0]:
            idxSelected=0
        lagrun += 1
        if lagrunMax<=lagrun:
            raise ValueError('Increase population size!')
    return idxSelected

def cross_over(p1,p2):
    n = p1.shape[0]
    c1 = np.zeros(n)
    c2 = np.zeros(n)
    for i in range(n):
        k = np.random.rand()
        c1[i] = k*p1[i] + (1.0-k)*p2[i]
        c2[i] = (1.0-k)*p1[i] + k*p2[i]
    return c1, c2

#def cross_over(p1,p2):
#    k1 = np.random.rand()
#    c1 = k1*p1 + (1.0-k1)*p2
#    k2 = np.random.rand()
#    c2 = k2*p1 + (1.0-k2)*p2
#    return c1,c2

def mutation(x, sigma):
    outOfBounds = True
    while outOfBounds:
        xmut = np.random.normal(x,sigma,x.shape)
        if np.all(xmut>=0) and np.all(xmut<=1):
            outOfBounds=False
    return xmut

def new_generation(x,fitness,nMut,nElite,nPop,sigma=0.2):
    idx = np.argsort(fitness)
    idx = np.flipud(idx)
    idxElite = idx[:nElite]
    xElite = x[idxElite]
    #nPop = x.shape[0]
    nVar = x.shape[1]
    nNonElite = nPop - nElite
    
    xNonElite = np.zeros([nPop-nElite,nVar])
    for i in range(nMut):
        idxSel = selection(fitness)
        xmut   = mutation(x[idxSel],sigma)
        xNonElite[i] = xmut
    nCrossOver = int((nNonElite-nMut)/2)+1
    for i in range(nCrossOver):
        idx1 = selection(fitness)
        idx2 = selection(fitness,idx1)
        c1, c2 = cross_over(x[idx1],x[idx2])
        if nNonElite-nMut-2*i>1:
            xNonElite[2*i+nMut] = c1
            xNonElite[2*i+nMut+1] = c2
        elif nNonElite-nMut-2*i==1:
            xNonElite[2*i+nMut] = c1

    return xElite, idxElite, xNonElite
        
def get_summary(x,cost,f):
    idx = np.min(cost)
    xBest = x[idx]
    fBest = f[idx]
    costBest = cost[idx]
    fAverage = f.mean()
    return xBest,fBest,fAverage,costBest

def ga(costFcn,lb,ub,populationSize,eliteRatio=0.1,mutationRate=0.2,
       iterMax=100,iterLag=50,initialPopulation=None,display=False):
    assert len(lb)==len(ub)
    costFcn = CostFunction(costFcn,lb,ub)
    nVar = len(lb)
    nPop = int(populationSize)
    nElite = int(nPop*eliteRatio)
    nMut   = int(nPop*mutationRate)

    iterMax = int(iterMax)
    iterLag = int(iterLag)
    
    if initialPopulation==None:
        initialPopulation = np.random.random([nPop,nVar])
    
    # solve
    nIter = 0
    nIterLag = 0
    
    costBest = float('inf')
    fBest = 0.0
    xBest = np.zeros(nVar)
    fAverage = np.zeros(nVar)
    
    xElite    = initialPopulation[:nElite]
    xNonElite = initialPopulation[nElite:]
    costElite,fElite,gElite = costFcn(xElite)
    
    xHist = np.zeros([iterMax,nVar])
    fHist = np.zeros(iterMax)
    fAvgHist = np.zeros(iterMax)
    while nIter<iterMax and nIterLag<iterLag:
        costNonElite,fNonElite,gNonElite = costFcn(xNonElite)
        xTotal    = np.vstack([xElite,xNonElite])
        fTotal    = np.hstack([fElite,fNonElite])
        costTotal = np.hstack([costElite,costNonElite])
        
        fitness = get_fitness(costTotal)
        sigma = get_sigma(nIter,iterMax)
        xElite, idxElite, xNonElite = new_generation(xTotal,fitness,nMut,nElite,nPop,sigma)
        
        costElite = costTotal[idxElite]
        fElite    = fTotal[idxElite]
        
        xBest,fBest,fAverage,costBestNew = get_summary(xTotal,costTotal,fTotal)
        xHist[nIter] = costFcn.norm.denormalize(xBest)
        fHist[nIter] = fBest
        fAvgHist[nIter] = fAverage
        print '%.4f\t%.4f\t%.4f'%(xHist[nIter,0],xHist[nIter,1],fBest)
        nIter += 1
        if costBestNew<costBest:
            nIterLag = 0
            costBest = costBestNew
        else:
            nIterLag += 1
    xBest = costFcn.norm.denormalize(xBest)
    return xBest,fBest
        
        


def rosenbrock(x):
    f = (1.0-x[0])**2.0 + 100.*(x[1]-x[0]**2.0)**2.0
    return f,[0]
    
def run_test1():
    xInit = np.random.random([100,2])
    print ga(rosenbrock, [-2,-2],[2,2], 20, eliteRatio=0.05, mutationRate=0.3,initialPopulation=xInit)

if __name__=="__main__":
    run_test1()