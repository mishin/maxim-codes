# -*- coding: utf-8 -*-
"""
Created on Mon May 07 13:28:33 2012
.. moduleauthor:: Daniel Neufeld

Multi-processor/core enabled, constrained, real coded genetic algorithm with crossover,
mutation, elitism, and attrition operators.  Run with no input to execute
several example problems for single and multi-core algorithms.

Module provides the following functions:
    - input/output: set GA run params and const function

    - returns best post-optimization design variable and function value
"""

import numpy as np
import pp as pp
#import pylab as pl
import pylab
class ga:
    """
    Sets up optimzation problem and GA run params
    
    :ivar costFcn: cost function - needs to returns objetive function value
        to be minimized and a list of constraints that are defined such that 
        values of less than or equal to zero are feasible.  Must import all 
        necessary modules to run if multi-core optimization is used
    :ivar lb: lower boundary numpy array
    :ivar ub: upper boundary numpy array
    :ivar params: optional parameters that will be passed into the cost
        function
    :ivar nCPU: number of CPUs / CPU cores to be used
    :ivar N: population size
    :ivar muRate: mutation rate
    :ivar sigma1: initial mutation standard deviation
    :ivar sigma2: final mutation standard deviation
    :ivar eliteRatio: ratio of elite values
    :ivar attrition1: initial attrition rate
    :ivar attrition2: final attrition rate
    :ivar maxIter: maximum number of iterations
    :ivar displayFlag: toggle [0 or 1] to display best value at each gen
    :ivar tolerance: constraint tolerance
    """
    def __init__(self,costFcn,lb,ub,params=0.0,nCPU=1,N=30,muRate=0.2,sigma1=0.3,sigma2=0.1,eliteRatio=0.2,attrition1=0.2,attrition2=0.2,maxLagIter=10,maxIter=100,displayFlag=1,tolerance=0.000001):
        self.costFcn=costFcn
        self.lb=np.float64(np.array(lb))
        self.ub=np.float64(np.array(ub))
        self.populationSize=np.int(N)
        self.displayFlag=np.int(displayFlag)
        self.muRate=np.float64(muRate)
        self.Nmut=np.int(np.ceil(N*self.muRate))
        self.sigma1=np.float64(sigma1)
        self.sigma2=np.float64(sigma2)
        self.maxLagIter=np.int(maxLagIter)
        self.maxIter=np.int(maxIter)
        self.N=np.int(N);
        self.nVar=self.lb.shape[0]
        self.params=params;
        self.tolerance=np.float64(tolerance)
        self.eliteRatio=np.float64(eliteRatio)
        self.Nelite=np.int(np.ceil(N*eliteRatio))
        self.attrition1=np.float64(attrition1)
        self.attrition2=np.float64(attrition2)
        self.fHist=list()
        self.gHist=list()
        self.xHist=list()        
        self.fOpt=np.array([])
        self.gOpt=np.array([])
        self.xOpt=np.array([])
        self.funcCount=0
        self.iteration=0;
        f,g=costFcn(lb,params)
        self.numCons=g.shape[0]
        self.setCPU(nCPU)
        self._disp=_display()
    def solve(self):
        xBestList=list()
        fBestList=list()
        gBestList=list()
        sigmaArray=np.linspace(self.sigma1,self.sigma2,self.maxIter)
        attrArray=np.linspace(self.attrition1,self.attrition2,self.maxIter)
        i=0        
        lagIter=0        
        x=np.random.rand(self.N,self.nVar)        
        f,g=self._systemEval(x)
        P=self._penaltyEval(f,g)
        xElite,fElite,gElite,Pelite=self._elite(x,f,g,P) 
        while i<self.maxIter and lagIter<self.maxLagIter:
            xBestList.append(self._denorm(xElite[0,:]))
            fBestList.append(fElite[0])
            gBestList.append(gElite[0,:])
            self._display(fBestList[i],xBestList[i],gBestList[i])
            x,xElite,fElite,gElite,Pelite=self._buildPopulation(x,f,g,xElite,fElite,gElite,Pelite,sigmaArray[i],attrArray[i])
            f,g=self._systemEval(x)            
            P=self._penaltyEval(f,g)
            if i>=1 and fBestList[i]==fBestList[i-1]:
                lagIter+=1
            else:
                lagIter=0
            if self.displayFlag==1:self._disp.update(fBestList[i])
            i+=1
            self.iteration=i
            self.fHist=fBestList
            self.gHist=gBestList
            self.xHist=xBestList              
        self.fOpt=fBestList[fBestList.__len__()-1]
        self.xOpt=xBestList[xBestList.__len__()-1]
        self.gOpt=gBestList[gBestList.__len__()-1]  
    def setCPU(self,N):
        self._nCPU=np.int(N) 
        if N>1:
            self._job_server=pp.Server(self._nCPU)        
        else:
            self._job_server=list()
    def _display(self,f,x,g):
        if self.displayFlag==1:
            s='F= {0:8.5f}, X=['.format(f)
            for xi in x:s+='{0:8.5f}'.format(xi)
            s+=']'            
            print s
    def plot(self):
        pylab.plot(range(self.fHist.__len__()),self.fHist)
        pylab.title("Solution History")
        pylab.xlabel("Generations")
        pylab.ylabel("f(x)")
        pylab.show()
    def _denorm(self,x):
        X=np.zeros(x.shape);
        if np.size(x.shape)>1:
            for i in np.arange(x.shape[0]):
                X[i,:]=x[i,:]*(self.ub-self.lb)+self.lb            
        else:
            X=x*(self.ub-self.lb)+self.lb
        return X
    def _systemEval(self,x):
        X=self._denorm(x)
        if self._nCPU<2:
            f,g=self._systemEval_lp(X)
        else:
            f,g=self._systemEval_pp(X)
        return f,g
    def _systemEval_lp(self,X):
        N=X.shape[0]
        f=np.zeros(N)
        g=np.zeros((N,self.numCons))
        for i in range(X.shape[0]):
            f[i],g[i,:]=self.costFcn(X[i,:],self.params)
        self.funcCount+=X.shape[0]
        return (f,g)          
    def _systemEval_pp(self,X):
        jobs=[]
        N=X.shape[0]
        f=np.zeros(N)
        g=np.zeros((N,self.numCons))        
        for i in range(N):
            jobs.append(self._job_server.submit(self.costFcn, (X[i,:], self.params)))
        i=0
        for job in jobs:
            completedJob=job()
            f[i]=completedJob[0]
            g[i,:]=completedJob[1]                                    
            i+=1
        self.funcCount+=X.shape[0]
        return f,g
    def _penaltyEval(self,f,g):
        P=np.zeros([f.shape[0]])   
        for i in range(f.shape[0]):
            r=f[i]/self.tolerance
            P[i]=f[i]+r*np.sum(g[i,:][g[i,:]>0])
        return P
    def _fitnessEval(self,P):
        realMin=np.finfo(np.double).tiny 
        D=P.max()-P        
        F=D/(D.sum()+realMin)
        return F  
    def _select(self,F,Idisq=-1):
        k=np.float64(np.random.rand())
        I=0
        for i in range(F.shape[0]):
            S1=F[0:i].sum()
            S2=F[0:i+1].sum()
            if k>S1 and k<=S2:I=i
        if I==Idisq:I=self._select(F,Idisq)
        return I
    def _mutate(self,xi,sigma):
        xiMut=xi.copy()
        i=np.random.randint(0,xi.shape[0])
        tmpMut=xi[i]+np.float64(np.random.normal(0,sigma))
        if tmpMut >1:
            xiMut[i]=1
        elif tmpMut<0:
            xiMut[i]=0
        else:
            xiMut[i]=tmpMut
        return xiMut 
    def _crossover(self,P1,P2):
        C1=np.zeros(P1.shape)
        C2=np.zeros(P2.shape)
        for i in range(P1.shape[0]):
            wt=np.float64(np.random.rand())
            C1[i]=wt*P1[i]+(1-wt)*P2[i]
            C2[i]=(1-wt)*P1[i]+wt*P2[i]
        return C1,C2
    def _elite(self,x,f,g,P):
        fitness=self._fitnessEval(P)
        I=fitness.argsort()[::-1]
        xElite=np.zeros([self.Nelite,x.shape[1]])
        fElite=np.zeros([self.Nelite])
        gElite=np.zeros([self.Nelite,g.shape[1]])
        pElite=np.zeros([self.Nelite])
        for i in range(self.Nelite):
            xElite[i,:]=x[I[i],:]
            fElite[i]=f[I[i]]
            gElite[i,:]=g[I[i],:]
            pElite[i]=P[I[i]]
        return xElite,fElite,gElite,pElite
    def _buildPopulation(self,x,f,g,xElite,fElite,gElite,Pelite,sigma,attRate):
        x_=np.concatenate((xElite,x))              
        f_=np.concatenate((fElite,f))
        g_=np.concatenate((gElite,g))              
        P_=self._penaltyEval(f_,g_)                 
        F=self._fitnessEval(P_)
        xMut=np.zeros([self.Nmut,x.shape[1]])
        xCrs=np.zeros([self.N-self.Nmut+1,x.shape[1]])
        for i in range(self.Nmut):
            I=self._select(F)
            xMut[i,:]=self._mutate(x_[I,:],sigma)            
        i=0
        while i < (self.N-self.Nmut):
            I1=self._select(F)
            I2=self._select(F,I1)
            C1,C2=self._crossover(x_[I1,:],x_[I2,:])
            xCrs[i,:]=C1
            xCrs[i+1,:]=C2
            i+=2
        xCrs=xCrs[0:xCrs.shape[0]-1,:]
        xNew=np.concatenate((xMut,xCrs))
        xEliteTmp,fEliteTmp,gEliteTmp,PeliteTmp=self._elite(x_,f_,g_,P_)
        xe=np.concatenate((xElite,xEliteTmp))        
        fe=np.concatenate((fElite,fEliteTmp))
        ge=np.concatenate((gElite,gEliteTmp))
        Pe=np.concatenate((Pelite,PeliteTmp))
        xEliteNew,fEliteNew,gEliteNew,PeliteNew=self._elite(xe,fe,ge,Pe)
        return xNew,xEliteNew,fEliteNew,gEliteNew,PeliteNew
class _display:
    def __init__(self):
        pylab.ion()       # Turn on interactive mode.
        pylab.hold(False)
        self.f=list()
        self.i=list()
        self._i=0
    def update(self,f):
        self._i+=1
        self.f.append(f)
        self.i.append(self._i)
        pylab.plot(self.i,self.f)
        pylab.pause(10**-5)
def runExamples():
    import time
    def costFcn(x,params):
        import numpy as np
        g=np.array([0.0,0.0])
        for i in range(100):
            f=100*(x[1]-x[0]**2)**2+(1-x[0])**2
            g[0]=2.0*x[0]-x[1]
            g[1]=1.0-x[1]
        return (f,g)
    def peaksFcn(x,params):
        import numpy as np
        g=np.array([0])
        f=3*(1-x[0])**2*np.exp(-(x[0]**2)-(x[1]+1)**2)-10*(x[0]/5-x[0]**3-x[1]**5)*np.exp(-x[0]**2-x[1]**2)-1/3*np.exp(-(x[0]+1)**2-x[1]**2)        
        return f,g
    def waveDropFcn(x,params):
        import numpy as np
        g=np.array([0])        
        f=-20*np.sin(((x[0]-4)**2+(x[1]-4)**2+0.1)**0.5)/((x[0]-4)**2+(x[1]-4)**2+0.1)**0.5
        return f,g
    def rosenFcn(x,params):
        import numpy as np
        g=np.array([0])
        f=100*(x[1]-x[0]**2)**2+(1-x[0])**2
        return f,g
        
    lb=np.array([-2.0,-5.0])
    ub=np.array([3.0,5.0])

    GA=ga(costFcn,lb,ub)
    GA.maxIter=200
    GA.maxLagIter=200
    
    
    print "Running test function - constrained Rosenbrock function"
    GA.setCPU(8)  
    GA.displayFlag=1;         
    GA.solve()
    
    print "Running single and multi-core comparions..."
    GA.setCPU(1)  
    GA.displayFlag=0;      
    time1=time.time()      
    GA.solve()
    
    
    GA.setCPU(8)     
    GA.displayFlag=0;
    time2=time.time()      
    GA.solve()
    
    print "Single core time = ",time.time()-time1, "s"  
    print "Multi  core time = ",time.time()-time2, "s"  
        
    print "Solving peaks function      - should be f~6.54  at x=[+0.227,-1.630]"
    GA=ga(peaksFcn,lb,ub)
    GA.setCPU(8)  
    GA.displayFlag=1;         
    GA.solve()
    print 'f=%.4f x1=%.4f x2=%.4f, Neval=%d' %(GA.fOpt,GA.xOpt[0],GA.xOpt[1],GA.funcCount)
    
    print "Solving rosenbrock function - should be f=0.00  at x=[+1.000,+1.000]"
    GA=ga(rosenFcn,lb,ub)
    GA.setCPU(8)  
    GA.displayFlag=1;
    GA.maxIter=200
    GA.maxLagIter=100
    GA.solve()
    print 'f=%.4f x1=%.4f x2=%.4f, Neval=%d' %(GA.fOpt,GA.xOpt[0],GA.xOpt[1],GA.funcCount)
    
    print "Solving wavedrop function   - should be f~16.53 at x=[+3.000,+4.000]"
    GA=ga(waveDropFcn,lb,ub)
    GA.setCPU(8)  
    GA.displayFlag=1;         
    GA.solve()
    print 'f=%.4f x1=%.4f x2=%.4f, Neval=%d' %(GA.fOpt,GA.xOpt[0],GA.xOpt[1],GA.funcCount)
    GA.plot()
if __name__=="__main__":
    runExamples()