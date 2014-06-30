# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 15:46:30 2014

@author: Maxim
"""

from aircraft_FW import FlyingWing

class Design(FlyingWing):
    def set_chords(self,chords):
        assert len(chords)==self.wing.nSec
        self.wing.chords = chords
        self.wing._process_data()
    
    def set_chord_by_index(self,chord,idx):
        assert idx<=len(self.wing.chords)
        self.wing.chords[idx] = chord
        self.wing._process_data()
    
    def set_sweep_angles(self,sweep):
        assert len(sweep)==self.wing.nSeg
        
    
    def set_sweep_by_index(self,sweep,idx):
        pass

    def set_spans(self,spans):
        assert len(spans)==len(self.wing.segSpans)
    
    def set_spans_by_index(self,spans,idx):
        pass
    
    def set_taper_ratio_by_index(self,TR,idx):
        pass
    
    def _upd_mass(self):
        pass
    
    def _upd_drag(self):
        pass
    
    def _upd_aero(self):
        pass