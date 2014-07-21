# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 15:46:30 2014

@author: Maxim
"""

from aircraft_FW import FlyingWing
import numpy as np

def load(name):
    ac = Design()
    ac.load_xls(name)
    return ac

class Design(FlyingWing):
    def set_chords(self,chords):
        assert len(chords)==self.wing.nSec
        self.wing.chords = chords
        self.wing._process_data()
    
    def set_chord_by_index(self,chord,idx):
        assert idx<=self.wing.nSec
        self.wing.chords[idx] = chord
        self.wing._process_data()
    
    def set_sweep_angles(self,sweep):
        """ NOTE: spans should be specified before sweeps """
        assert len(sweep)==self.wing.nSeg
        offset = np.tan(sweep)*self.wing.segSpans
        self.wing.secOffset = offset
        self.wing._process_data()
    
    def set_sweep_by_index(self,sweep,idx):
        pass

    def set_spans(self,spans):
        spans = np.asarray(spans)
        assert len(spans)==len(self.wing.segSpans)
        self.wing.segSpans = spans
        newOffsets = np.tan(self.wing.segSweepLErad)*spans
        self.wing.secOffset = newOffsets
        self.wing._process_data()
    
    def set_spans_by_index(self,spans,idx):
        pass
    
    def set_taper_ratio_by_index(self,TR,idx):
        """
        Parameters
        ----------
        
        idx: int
            index segment where taper ratio will be applied
        """
        tipChord = self.wing.chords[idx]*TR
        self.set_chord_by_index(tipChord,idx+1)

    def _upd_mass(self):
        self._update_mass()
    
    def _upd_drag(self):
        self._update_parasite_drag()
    
    def _upd_aero(self):
        pass


def run_test1():
    ac = load('Baseline1')
    print ac.get_mass()
    aero = ac.get_aero_trim(250.,1e4)
    print aero.derivs.Cnb, aero.derivs.Clb
    #ac.display()
    ac.set_chord_by_index(2,2.0)
    ac._upd_mass()
    aero = ac.get_aero_trim(250.,1e4)
    print aero.derivs.Cnb, aero.derivs.Clb
    print ac.get_mass()
    #ac.display()
    ac.set_spans([1.5,5.0])
    ac._upd_mass()
    print ac.get_mass()
    #ac.display()
    aero = ac.get_aero_trim(250.,1e4)
    print aero.derivs.Cnb, aero.derivs.Clb
    ac.wing.secTwist[-1] = -2.0
    ac.wing._process_data()
    aero = ac.get_aero_trim(250.,1e4)
    print aero.derivs.Cnb, aero.derivs.Clb

if __name__=="__main__":
    run_test1()
