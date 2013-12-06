# -*- coding: utf-8 -*-
"""
Created on Fri Dec 06 15:30:56 2013

@author: Maxim
"""

class FileOutput(object):
    def __init__(self,path):
        self.path = path
    
    def write_array(self,array,name=''):
        fid = open(self.path,'at')
        fid.write('%s = array(['%name)
        fid.write('%.6f'%array[0])
        for val in array[1:]:
            fid.write(' ,%.6f\t'%val)
        fid.write('])\n')
        fid.close()
    
    def write_string(self,string):
        fid = open(self.path,'at')
        fid.write(string)
        fid.write('\n')
        fid.close()
        