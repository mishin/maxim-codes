# -*- coding: utf-8 -*-
"""
Created on Tue Jul 08 13:19:19 2014

@author: Maxim
"""

from part_designer import CadDesigner
import numpy as np

import aircraft_FW
from airfoil_cad import Airfoil2D
import paths


def create_catia_wing(ac, fileIgs):
    #ac = aircraft_FW.load('Baseline1')
    scale = 1.
    symmetric = True
    show = False
    Rle = 0.1207**2./2
    #paths.myPaths.set_file_prefix('fw2')
    #fileIgs = paths.myPaths.get_tmp_file('igs')

    pd = CadDesigner(show)
    pd.add_new_part()
    zte = 0.0 # trailing edge gap
    xref, yref, zref = 0.0, 0.0, 0.0
    pd.refPoint = np.array([xref,yref,zref])
    chord = ac.wing.chords[0] *scale
    incidence = ac.wing.incidence
    af = Airfoil2D()
    af.convert_af(ac.wing.airfoils[0])
    af.leRad = Rle
    af.leCenter = [Rle, 0]
    af.leSlope = 0.0
    af.set_trailing_edge_gap(zte/chord)
    pd.create_airfoil(af,chord,incidence,0.25,True)
    for i in range(ac.wing.nSec-1):
        chord = ac.wing.chords[i+1]*scale
        incidence += ac.wing.secTwist[i]
        span = ac.wing.segSpans[i]*scale
        sweep = ac.wing.segSweepLEdeg[i]
        dihedral = ac.wing.segDihedral[i]
#        af = Airfoil2D()
#        af.convert_af(ac.wing.airfoils[i+1])
#        af.leRad = Rle
#        af.leCenter = [Rle, 0]
#        af.leSlope = 0.0
        af.set_trailing_edge_gap(zte/chord)
        pd.create_tip_section(af,chord,incidence,0.25,span,sweep,dihedral,False,True)
    pd.rot_sym_cfd()

    pd.save_part_igs(fileIgs)
    pd.close_part()


if __name__=="__main__":
    create_catia_wing()