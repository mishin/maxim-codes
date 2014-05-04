# -*- coding: utf-8 -*-
"""
Created on Fri Apr 04 22:18:19 2014

@author: Maxim
"""

#from numpy import array, zeros, sin, cos, random, vstack, hstack, cumsum
import numpy as np
from misc_tools import Normalization
import math
import sys

def solve(costFcn, lb, ub, populationSize, initialPopulation=None, 
                 eliteRatio=0.1, mutationRate=0.2, iterMax=100, iterLag=20, 
                 display=False):
    ga = GeneticAlgorithm(costFcn, lb, ub, populationSize, initialPopulation, 
                 eliteRatio, mutationRate, iterMax, iterLag, display)
    result = ga.solve()
    return result

class GeneticAlgorithm(object):
    def __init__(self,costFcn, lb, ub, populationSize, initialPopulation=None, 
                 eliteRatio=0.1, mutationRate=0.2, iterMax=100, iterLag=20, 
                 display=False):
        self.costFcn = costFcn
        assert len(lb)==len(ub)
        self.lb            = lb
        self.ub            = ub
        self.norm          = Normalization(lb,ub,0,1)
        self.nVar          = len(lb)
        self.nPop          = int(populationSize)
        self.nElite        = int(self.nPop*eliteRatio)
        self.nMutation     = int(self.nPop*mutationRate)
        self.nCrossover    = self.nPop-self.nElite-self.nMutation
        self.penaltyFactor = 1.0e5
        self.iterMax       = int(iterMax)
        self.iterLag       = int(iterLag)
        self._lagRunSelection = 20
        if initialPopulation==None:
            self.initialPopulation = np.random.random([self.nPop,self.nVar])
        else:
            self.initialPopulation = initialPopulation

    def solve(self):
        nIter   = 0
        lagIter = 0
        
        xElite    = self.initialPopulation[:self.nElite]
        xNonElite = self.initialPopulation[self.nElite:]
        fElite, gElite, costElite = evaluate_cost(costFcn, xElite)
        
        while nIter<self.iterMax and lagIter<self.iterLag:
            fNonElite, gNonElite, costNonElite = self._evaluate_cost(xNonElite)
            xTotal    = np.vstack([xElite,xNonElite])
            fTotal    = vstack([fElite,fNonElite])
            gTotal    = vstack([gElite,gNonElite])
            costTotal = vstack([costElite,costNonElite])
            fitness = self.get_fitness(costTotal)
            xElite, idxElite, xNonElite = self._new_generation(fitness,xTotal)
    
    def _evaluate_cost(self, xArray):
        n = xArray.shape[0]
        f1,g1 = self.costFcn(xArray[0])
        f = zeros(n)
        cost = zeros(n)
        g = zeros([n,len(g1)])
        f[0] = f1
        g[0] = array(g1)
        cost[0] = self._penalty(f[0],g[0])
        for i in range(n-1):
            fnew,gnew = self.costFcn(xArray[i+1])
            f[i+1] = fnew
            g[i+1] = array(gnew)
            cost[i+1] = self._penalty(f[i+1], g[i+1])
        return f,g,cost

    def _penalty(self,f,g):
        cost = f
        for gval in g:
            cost += self.penaltyFactor*max([gval,0.0])
        return cost

    def _get_fitness(self, cost, correctNan=True):
        bigNumber = sys.float_info.max / (1e3*len(cost))
        if correctNan:
            for i,val in enumerate(cost):
                if math.isnan(val):
                    cost[i] = bigNumber
        fWorst = max(cost)
        diff = fWorst - cost
        diffTotal = sum(diff)
        fitness   = diff/diffTotal
        return fitness

    def _mutate(self,xinput,sigma=0.2):
        xnew = 2.0
        while xnew>1.0 or xnew<0.0:
            xnew = np.random.normal(xinput,sigma)
        return xnew

    def _crossover(self):
        pass

    def _selection(self,cumFitness,idxPrevious=None):
        idxSelected = idxPrevious
        lagRun = 0
        while idxSelected==idxPrevious:
            k = np.random.rand()
            for i in range(cumFitness.shape[0]):
                if k>cumFitness[i] and k<=cumFitness[i+1]:
                    idxSelected = i+1
            if k<=cumFitness[0]:
                idxSelected = 0
            lagRun += 1
            if lagRun>=self._lagRunSelection:
                raise ValueError('Increase population size')
        return idxSelected

    def _new_generation(self, fitness, xTotal):
        idx = fitness.argsort()
        xTotal = xTotal[idx]
        idxElite = idx[:self.nElite]
        xElite = xTotal[idxElite]
        xMutate    = zeros([self.nMutation,self.nVar])
        xCrossover = zeros([self.nCrossover,self.nVar])
        idxPrevious = None
        cumFitness = np.cumsum(fitness)
        for i in range(self.nMutation):
            idxm = self._select(cumFitness,idxPrevious)
            xMutate[i] = self._mutate(xTotal[idxm])
        for i in range(self.n):
            pass
        return xElite, idxElite


def ga(costFcn, lb, ub, populationSize, initialPopulation=None, eliteRatio=0.1,
       mutationRate=0.2, iterMax=100, iterLag=20, xtol=1e-6, ftol=1e-6, display=False):
           
    random.seed(1394) #NOTE: this one for debugging only
    rho = 1e5
    xbest = list()
    fbest = list()
    gbest = list()
    
    assert len(lb)==len(ub)
    
    normalizer = Normalization(lb,ub,0.0,1.0)
    
    xerr = xtol+1.
    ferr = ftol+1.
    nIter = 0
    
    nVariables = len(lb)
    nPopulation = int(populationSize)
    nElite = int(populationSize*eliteRatio)
    nMutation = int(populationSize*mutationRate)
    
    if initialPopulation==None:
        initialPopulation=random.random([nPopulation,nVariables])
    xElite = initialPopulation[:nElite]
    xNonElite = initialPopulation[nElite:]
    fElite, gElite,costElite = evaluate_cost(costFcn, xElite, rho)
    
    while nIter<iterMax or xerr>xtol or ferr>ftol:
        fNonElite, gNonElite, costNonElite = evaluate_cost(costFcn, xNonElite, rho)
        xTotal = vstack([xElite,xNonElite])
        fTotal = vstack([fElite,fNonElite])
        gTotal = vstack([gElite,gNonElite])
        costTotal = vstack([costElite,costNonElite])
        fitness = get_fitness(costTotal)
        xElite, xNonElite = new_generation(fitness,xTotal)


def get_fitness(cost, correctNan=True):
    bigNumber = sys.float_info.max / (1e3*len(cost))
    if correctNan:
        for i,val in enumerate(cost):
            if math.isnan(val):
                cost[i] = bigNumber
    fWorst = max(cost)
    diff = fWorst - cost
    diffTotal = sum(diff)
    fitness = diff/diffTotal
    return fitness
 
def evaluate_cost(costFunction,xArray,rho=1.0e5):
    n = xArray.shape[0]
    f1,g1 = costFunction(xArray[0])
    f = zeros(n)
    cost = zeros(n)
    g = zeros([n,len(g1)])
    f[0] = f1
    g[0] = array(g1)
    cost[0] = penalty(f[0],g[0],rho)
    for i in range(n-1):
        fnew,gnew = costFunction(xArray[i+1])
        f[i+1] = fnew
        g[i+1] = array(gnew)
        cost[i+1] = penalty(f[i+1], g[i+1],rho)
    return f,g,cost

def penalty(f,g,rho):
    cost = f
    for gval in g:
        cost += rho*max([gval,0.0])
    return cost

def costFcn(x):
    f = (1.+cos(12.*(x[0]*x[0]+x[1]*x[1])**0.5))/(0.5*(x[0]*x[0]+x[1]*x[1])+2.)
    g = []
    return f,g

def run_test():
    lb = array([-5.0,-5.0])
    ub = array([5.0, 5.0])
    result = solve(costFcn,lb,ub,20)
    print result

if __name__=="__main__":
    run_test()