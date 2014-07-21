# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 20:58:10 2014

@author: Maxim
"""

import design
from scipy.optimize import minimize


class DesignFormulation(design.Design):
    def setup(self):
        self.lb = None
        self.ub = None
        self.x0 = None
        self.xCurrent = None

    def check_x(self,x):
        """
        To minimize function evaluation, this function re-calculates analysis if 
        new x is given, otherwise precalculated results are used.
        """
        if x==self.xCurrent:
            self.xCurrent = x
            self._upd_analysis(x)


def run_optimization():
    pass



if __name__=="__main__":
    run_optimization()