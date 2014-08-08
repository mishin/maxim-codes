# -*- coding: utf-8 -*-
"""
Created on Wed Aug 06 15:11:19 2014

@author: Maxim
"""

# Optional - turn off bytecode (.pyc files)
import sys
sys.dont_write_bytecode = True

from SALib.sample import saltelli, morris_oat, fast_sampler
from SALib.analyze import sobol, morris, extended_fast
from SALib.test_functions import Sobol_G, Ishigami
from SALib.util import scale_samples, read_param_file
import numpy as np
import random as rd
from ucav_design_1 import DesignFormulation


param_file = 'ucav_sensitivity_input.txt'
pf = read_param_file(param_file)

param_values = saltelli.sample(25, pf['num_vars'], calc_second_order = False)


scale_samples(param_values, pf['bounds'])


ac = DesignFormulation()
ac.load_xls('Baseline1')
ac.setup()

Y = np.zeros([len(param_values),8])
fid = open('SGOutput.txt','wt')
fid.close()
for i,param in enumerate(param_values):
    out = ac.run_full_analysis(param)
    Y[i] = out
    fid = open('SGOutput.txt','at')
    for val in out:
        fid.write('%.10e '%val)
    fid.write('\n')
    fid.close()


for i in range(len(out)):
    Si = sobol.analyze(param_file, 'SGOutput.txt', column = i, calc_second_order = False, conf_level = 0.95)
