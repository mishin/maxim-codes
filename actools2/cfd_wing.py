# -*- coding: utf-8 -*-
"""
Created on Sat Aug 30 11:28:44 2014

@author: Maxim
"""

from ucav_design_1 import DesignFormulation
from catia_wing import create_catia_wing
from pointwise_fw import create_fw_cmesh
import paths
import os

def run_cfd_wing_analysis():
    ac = DesignFormulation()
    ac.load_xls('Baseline1')
    ac.setup()
    
    paths.myPaths.set_file_prefix('fw_debug')
    igsPath = paths.myPaths.get_tmp_file('igs')
    symCasPath = paths.myPaths.get_tmp_file('cas','_sym')
    nonsymCasPath = paths.myPaths.get_tmp_file('cas','_nonsym')

    create_catia_wing(ac,igsPath)
    create_fw_cmesh(ac,igsPath,symCasPath, nonsymCasPath)
    os.remove(igsPath)


if __name__=="__main__":
    run_cfd_wing_analysis()
    